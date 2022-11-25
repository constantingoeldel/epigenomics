use crate::error::Error;
use clap::Parser;
use itertools::Itertools;
use rayon::prelude::*;
use std::{
    ffi::OsString,
    fmt::Display,
    fs::{self, File},
    io::{self, BufRead, Write},
    path::PathBuf,
};

mod error;

type Result<T> = std::result::Result<T, error::Error>;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path of directory containing the methlyome files from which to extract the CG-sites
    #[arg(short, long)]
    methylome: String,

    /// Path of the annotation file containing information about beginning and end of gbM-genes
    #[arg(short, long)]
    annotation: String,

    /// Size of the window in percent of the gbM-gene length
    #[arg(short, long, default_value_t = 5)]
    window_size: u32,

    /// Path of the directory where extracted segments shall be stored
    #[arg(short, long)]
    output_dir: String,

    /// Overwrite current content of the output directory?
    #[arg(short, long, default_value_t = false)]
    force: bool,
}

struct GbmGene {
    chromosome: u8,
    start: u32,
    end: u32,
    name: String,
    strand: bool, // true for + strand, false for - strand
}

struct CgSite {
    chromosome: u8,
    location: u32,
    strand: bool, // true for + strand, false for - strand
    original: String,
}

impl GbmGene {
    fn from_annotation_file_line(s: &str) -> Option<Self> {
        s.split('\t')
            .collect_tuple()
            .map(|(chromosome, start, end, name, _, strand)| GbmGene {
                chromosome: chromosome.parse::<u8>().unwrap(),
                start: start.parse::<u32>().unwrap(),
                end: end.parse::<u32>().unwrap(),
                name: String::from(name),
                strand: strand == "+",
            })
    }
}

impl CgSite {
    fn from_methylome_file_line(s: &str) -> Option<Self> {
        s.split('\t')
            .collect_tuple()
            .filter(|(chromosome, _, _, _, _, _, _, _, _)| chromosome != &"seqnames") // Filter out header row
            .map(|(chromosome, location, strand, _, _, _, _, _, _)| CgSite {
                chromosome: chromosome.parse::<u8>().unwrap(),
                location: location.parse::<u32>().unwrap(),
                strand: strand == "+",
                original: s.to_owned().clone() + "\n",
            })
    }
}

impl Display for GbmGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Gene {} is located on the {} strand of chromosome {} ranging from bp {} to bp {} ",
            self.name,
            if self.strand { "+" } else { "-" },
            self.chromosome,
            self.start,
            self.end
        )
    }
}

impl Display for CgSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CG site is located on the {} strand of chromosome {} at bp {}",
            if self.strand { "+" } else { "-" },
            self.chromosome,
            self.location
        )
    }
}

fn main() {
    let args = Args::parse();
    let methylome_files = load_methylome(args.methylome);

    if let Err(err) = methylome_files {
        println!("{}", err);
        return;
    }

    let annotation_lines = lines_from_file(args.annotation);

    if let Err(err) = annotation_lines {
        println!("{}", err);
        return;
    }

    let result = set_up_output_dir(&args.output_dir, args.window_size as usize, args.force);
    if let Err(err) = result {
        println!("{}", err);
        return;
    }

    let mut genes: Vec<GbmGene> = Vec::new();

    for line in annotation_lines.unwrap() {
        let gene = GbmGene::from_annotation_file_line(&line.unwrap()).unwrap();
        genes.push(gene)
    }

    let result = extract_windows(
        methylome_files.unwrap(),
        genes,
        args.window_size,
        args.output_dir,
    );

    if let Err(e) = result {
        println!("{}", e);
        return;
    }
}

fn lines_from_file(filename: String) -> Result<io::Lines<io::BufReader<File>>> {
    let file = File::open(&filename)
        .map_err(|_| Error::FileError(String::from("Annotation file"), String::from(filename)))?;
    Ok(io::BufReader::new(file).lines())
}

fn load_methylome(methylome: String) -> Result<Vec<(PathBuf, OsString)>> {
    let methylome_dir = fs::read_dir(&methylome)
        .map_err(|_| Error::FileError(String::from("Methylome directory"), methylome.clone()))?;
    let methylome_files = methylome_dir
        .map(|f| (f.as_ref().unwrap().path(), f.unwrap().file_name()))
        .collect();
    Ok(methylome_files)
}

fn extract_windows(
    methylome_files: Vec<(PathBuf, OsString)>,
    genes: Vec<GbmGene>,
    window_size: u32,
    output_dir: String,
) -> Result<()> {
    methylome_files
        .par_iter()
        .try_for_each_with(&genes, |genes, (path, filename)| -> Result<()> {
            let mut last_gene = &genes[0];
            let mut last_gene_index = 0;
            let mut searched_count: i64 = 0;
            println!("Extracting windows for file {:#?} ", filename);
            let file = File::open(&path).map_err(|_| {
                Error::FileError(
                    String::from("methylome file"),
                    String::from(filename.to_str().unwrap()),
                )
            })?;
            let lines = io::BufReader::new(file).lines();

            for (i, line_result) in lines.enumerate() {
                if i % 50000 == 0 {
                    println!("Extracting cg_site {i}");
                }

                if let Ok(line) = line_result {
                    let Some(cg) = CgSite::from_methylome_file_line(&line) else {continue;}; // If cg site could not be extracted, continue with the next line. Happens on header rows, for example.
                    let cg_site_in_gene = |cg: &CgSite, gene: &GbmGene| {
                        cg.strand == gene.strand
                            && cg.chromosome == gene.chromosome
                            && gene.start <= cg.location
                            && cg.location <= gene.end
                    };

                    if cg_site_in_gene(&cg, last_gene) {
                        place_site(cg, last_gene, window_size, &output_dir, &filename)?;
                        continue;
                    }
                    let next_gene = &genes[last_gene_index + 1];
                    if cg_site_in_gene(&cg, next_gene) {
                        place_site(cg, next_gene, window_size, &output_dir, &filename)?;
                        continue;
                    }

                    for (i, gene) in genes.iter().enumerate() {
                        searched_count += 1;
                        if cg_site_in_gene(&cg, gene) {
                            place_site(cg, gene, window_size, &output_dir, &filename)?;
                            last_gene = gene;
                            last_gene_index = i;
                            break;
                        }
                    }
                }
            }
            println!("Searched through a total of {} genes", searched_count);
            Ok(())
        })
}

fn place_site(
    cg: CgSite,
    gene: &GbmGene,
    window_size: u32,
    output_dir: &str,
    filename: &OsString,
) -> Result<()> {
    let mut percentile = (cg.location - gene.start) * 100 / (gene.end - gene.start);
    if percentile == 100 {
        // very unelegant fix but otherwise the last cg-site would land in its own category
        percentile = 99
    }
    let window = (percentile / window_size) * window_size; // integer division rounds down
    let output_file = format!("{}/{}/{}", output_dir, window, filename.to_str().unwrap());
    let mut file = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(&output_file)
        .expect(&format!("Could not open output file, {output_file}"));

    file.write(cg.original.as_bytes())?;
    Ok(())
}

fn set_up_output_dir(output_path: &str, window_size: usize, overwrite: bool) -> Result<()> {
    fs::read_dir(&output_path).map_err(|_| {
        Error::FileError(String::from("Output directory"), String::from(output_path))
    })?; // Throw error if base output dir does not exist

    if overwrite {
        fs::remove_dir_all(&output_path).unwrap();
        fs::create_dir(&output_path).unwrap();
    }

    for i in (0..100).step_by(window_size) {
        let path = format!("{}/{}", &output_path, i);
        let window_dir = fs::read_dir(&path);

        match window_dir {
            Ok(_) => continue,
            Err(_) => fs::create_dir(path).unwrap(),
        };
    }
    Ok(())
}
