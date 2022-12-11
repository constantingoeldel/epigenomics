use crate::{arguments::Args, error::Error};
use clap::Parser;

use methylation_site::*;
use rayon::prelude::*;
use setup::set_up_output_dir;
use std::{
    ffi::OsString,
    fs::{self, File},
    io::{self, BufRead},
    path::PathBuf,
};
use structs::*;
use windows::*;

mod arguments;
mod error;
mod methylation_site;
mod setup;
mod structs;
mod windows;

fn main() {
    let args = Args::parse();

    let window_step = if args.window_step.is_some() {
        args.window_step.unwrap()
    } else {
        args.window_size
    };

    let methylome_files = load_methylome(args.methylome);

    if let Err(err) = methylome_files {
        println!("{}", err);
        return;
    }

    let annotation_lines = lines_from_file(args.genome);

    if let Err(err) = annotation_lines {
        println!("{}", err);
        return;
    }

    let mut genes: Vec<Gene> = Vec::new();

    for line in annotation_lines.unwrap() {
        let gene = Gene::from_annotation_file_line(&line.unwrap()).unwrap();
        genes.push(gene)
    }

    let chromosome_count = genes // number of different chromosomes assuming they are named from 1 to highest
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
    genes.iter().for_each(|g| {
        let chromosome = &mut structured_genes[(g.chromosome - 1) as usize];
        let strand = match &g.strand {
            Strand::Sense => &mut chromosome.sense,
            Strand::Antisense => &mut chromosome.antisense,
        };
        strand.push(g.to_owned());
    });

    let mut max_gene_length: i32 = 100; // if not using absolute window sizes, the maximum gene length will be 100%
    if args.absolute {
        for gene in &genes {
            let length = gene.end - gene.start;
            if length > max_gene_length {
                max_gene_length = length
            }
        }
        println!("The maximum gene length is {} bp", max_gene_length);
    }

    let result = set_up_output_dir(
        &args.output_dir,
        args.force,
        window_step as usize,
        args.absolute,
        args.cutoff,
        max_gene_length,
    );

    if let Err(err) = result {
        println!("{}", err);
        return;
    }

    let result = methylome_files.unwrap().par_iter().try_for_each_with(
        structured_genes,
        |genome, (path, filename)| -> Result<()> {
            let file = File::open(path).map_err(|_| {
                Error::File(
                    String::from("methylome file"),
                    String::from(filename.to_str().unwrap()),
                )
            })?;

            let result = extract_windows(
                file,
                genome.to_vec(),
                args.window_size,
                window_step,
                max_gene_length as i32,
                args.ignore_strand,
                args.cutoff,
            );
            match result {
                Ok(windows) => windows.save(&args.output_dir, filename, window_step as usize),
                Err(error) => Err(error),
            }
        },
    );
    if let Err(err) = result {
        println!("{}", err)
    } else {
        println!("Done!")
    }
}

fn lines_from_file(filename: String) -> Result<io::Lines<io::BufReader<File>>> {
    let file = File::open(&filename)
        .map_err(|_| Error::File(String::from("Annotation file"), filename))?;
    Ok(io::BufReader::new(file).lines())
}

fn load_methylome(methylome: String) -> Result<Vec<(PathBuf, OsString)>> {
    let methylome_dir = fs::read_dir(&methylome)
        .map_err(|_| Error::File(String::from("Methylome directory"), methylome.clone()))?;
    let methylome_files = methylome_dir
        .map(|f| (f.as_ref().unwrap().path(), f.unwrap().file_name()))
        .collect();
    Ok(methylome_files)
}
