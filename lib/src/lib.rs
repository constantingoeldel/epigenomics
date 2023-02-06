use crate::{arguments::Args, error::Error};

use files::*;
pub use methylation_site::*;
use rayon::prelude::*;
use setup::set_up_output_dir;
use std::{
    ffi::OsString,
    fs::{self, File},
    io::{self, BufRead},
    path::PathBuf,
    sync::Mutex,
};
use structs::*;
use windows::*;
pub mod arguments;
pub mod error;
pub mod files;
pub mod methylation_site;
pub mod setup;
pub mod structs;
pub mod windows;

pub fn extract(args: Args) -> Result<()> {
    let start = std::time::Instant::now();
    let mut args = args;

    // Adj ust window_step to default value
    if args.window_step == 0 {
        args.window_step = args.window_size;
    }

    let methylome_files = load_methylome(&args.methylome)?;
    let annotation_lines = lines_from_file(&args.genome)?;

    let mut genes: Vec<Gene> = Vec::new();

    // Parse annotation file to extract genes
    for line in annotation_lines {
        let line = line?;
        let gene = Gene::from_annotation_file_line(&line, args.invert);
        if let Some(gene) = gene {
            genes.push(gene)
        }
    }

    // number of different chromosomes assuming they are named from 1 to highest
    let chromosome_count = genes
        .iter()
        .max_by_key(|g| g.chromosome)
        .unwrap()
        .chromosome;

    genes.sort_by_key(|g| g.start); // Sort genes by start bp (propably already the case), needed for binary search

    // Structure genes first by chromosome, then by + and - strand => [Chromosome_1(+ Strand, - Strand), Chromosome_2(+,-), ..]
    let mut structured_genes: Vec<GenesByStrand> = vec![
        GenesByStrand {
            sense: Vec::new(),
            antisense: Vec::new()
        };
        chromosome_count.into()
    ];
    // Put genes into their correct bucket
    let mut gene_length_sum = 0;
    let mut sense_gene_count = 0;
    genes.iter().for_each(|g| {
        let chromosome = &mut structured_genes[(g.chromosome - 1) as usize];
        let strand = match &g.strand {
            Strand::Sense => &mut chromosome.sense,
            Strand::Antisense => &mut chromosome.antisense,
            _ => return,
        };
        if g.strand == Strand::Sense {
            sense_gene_count += 1;
        }
        gene_length_sum += g.end - g.start;
        strand.push(g.to_owned());
    });
    let average_gene_length = gene_length_sum / genes.len() as u32;
    println!(
        "Average gene length: {} bp, {} genes, of which {} are on the sense strand and {} on the antisense strand",
        average_gene_length,
        genes.len(),
        sense_gene_count,
        genes.len() - sense_gene_count
    );

    // Determine the maximum gene length by iterating over all genes
    let mut max_gene_length: u32 = 100; // if not using absolute window sizes, the maximum gene length will be 100%
    if args.absolute {
        for gene in &genes {
            let length = gene.end - gene.start;
            if length > max_gene_length {
                max_gene_length = length
            }
        }
        println!("The maximum gene length is {max_gene_length} bp");
    }

    set_up_output_dir(args.clone(), max_gene_length)?;

    let distributions = Mutex::new(Vec::new());
    let steady_state_methylations = Mutex::new(Vec::new());

    methylome_files.par_iter().try_for_each_with(
        structured_genes,
        |genome, (path, filename)| -> Result<()> {
            let file = open_file(path, filename)?;
            let mut windows =
                extract_windows(file, genome.to_vec(), max_gene_length, args.clone())?;
            if args.invert {
                windows = windows.inverse();
            }
            windows.save(&args, filename)?;
            distributions.lock().unwrap().push(windows.distribution());
            steady_state_methylations
                .lock()
                .unwrap()
                .push(windows.steady_state_methylation());

            Ok(())
        },
    )?;
    let sample_size = steady_state_methylations.lock().unwrap().len();
    let mut average_methylation = vec![0.0; steady_state_methylations.lock().unwrap()[0].len()];
    for source in steady_state_methylations.lock().unwrap().iter() {
        for (i, window) in source.iter().enumerate() {
            average_methylation[i] += window / sample_size as f64
        }
    }

    let distribution_file = format!("{}/distribution.txt", &args.output_dir);
    let methylation_file = format!("{}/steady_state_methylation.txt", &args.output_dir);

    fs::write(
        distribution_file,
        Windows::print_distribution(&distributions.lock().unwrap()[0]),
    )?;
    fs::write(
        methylation_file,
        Windows::print_steady_state_methylation(&average_methylation),
    )?;

    println!("Done in: {:?}", start.elapsed());

    Ok(())
}
