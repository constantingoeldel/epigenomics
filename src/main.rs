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

/// simple tool to separate a methylome by position within a gene
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

    /// Use absolute length in basis-pairs for window size instead of percentage of gene length
    #[arg(long, default_value_t = false)]
    absolute: bool,

    /// Match strands? If supplied, only CG sites on the same strand as the gene will be included in the result
    #[arg(short, long, default_value_t = false)]
    ignore_strand: bool,
}
#[derive(Clone, PartialEq)]
enum Strand {
    Sense,
    Antisense,
}

#[derive(Clone)]
struct GbmGene {
    chromosome: u8,
    start: u32,
    end: u32,
    name: String,
    strand: Strand,
}
#[derive(Clone)]
struct CgSite {
    chromosome: u8,
    location: u32,
    strand: Strand,
    original: String,
}
#[derive(Clone)]
struct GenesByStrand {
    sense: Vec<GbmGene>,
    antisense: Vec<GbmGene>,
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
                strand: if strand == "+" {
                    Strand::Sense
                } else {
                    Strand::Antisense
                },
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
                strand: if strand == "+" {
                    Strand::Sense
                } else {
                    Strand::Antisense
                },
                original: s.to_owned().clone(),
            })
    }
}

impl Display for GbmGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Gene {} is located on the {} strand of chromosome {} ranging from bp {} to bp {} ",
            self.name, self.strand, self.chromosome, self.start, self.end
        )
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Sense => write!(f, "+"),
            Strand::Antisense => write!(f, "-"),
        }
    }
}

impl Display for CgSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CG site is located on the {} strand of chromosome {} at bp {}",
            self.strand, self.chromosome, self.location
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

    let mut genes: Vec<GbmGene> = Vec::new();

    for line in annotation_lines.unwrap() {
        let gene = GbmGene::from_annotation_file_line(&line.unwrap()).unwrap();
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

    let mut max_gene_length: u32 = 100; // if not using absolute window sizes, the maximum gene length will be 100%
    if args.absolute {
        for gene in &genes {
            let length = gene.end - gene.start;
            if length > max_gene_length {
                max_gene_length = length
            }
        }
        println!("The maximum gene length is {} bp", max_gene_length);
    }
    let window_count = (max_gene_length / args.window_size) as usize;

    let result = set_up_output_dir(
        &args.output_dir,
        args.force,
        window_count,
        args.window_size as usize,
    );

    if let Err(err) = result {
        println!("{}", err);
        return;
    }

    let result = extract_windows(
        methylome_files.unwrap(),
        structured_genes,
        args.window_size,
        &args.output_dir,
        window_count,
        args.ignore_strand,
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
    genome: Vec<GenesByStrand>,
    window_size: u32,
    output_dir: &str,
    window_count: usize,
    ignore_strand: bool,
) -> Result<()> {
    methylome_files.par_iter().try_for_each_with(
        &genome,
        |genome, (path, filename)| -> Result<()> {
            let mut last_gene = &genome[0].sense[0];
            // let mut last_gene_index = 0;
            let mut searched_count: i64 = 0;
            let mut skipped_count: i64 = 0;
            let mut windows: Vec<Vec<CgSite>> = vec![Vec::new(); window_count];
            println!("Extracting windows for file {:#?}", filename);

            let file = File::open(&path).map_err(|_| {
                Error::FileError(
                    String::from("methylome file"),
                    String::from(filename.to_str().unwrap()),
                )
            })?;
            let lines = io::BufReader::new(file).lines();
            for line_result in lines {
                if let Ok(line) = line_result {
                    let Some(cg) = CgSite::from_methylome_file_line(&line) else {continue;}; // If cg site could not be extracted, continue with the next line. Happens on header rows, for example.
                    let cg_site_in_gene = |cg: &CgSite, gene: &GbmGene| {
                        cg.chromosome == gene.chromosome
                            && gene.start <= cg.location
                            && cg.location <= gene.end
                            && (ignore_strand || cg.strand == gene.strand)
                    };

                    if cg_site_in_gene(&cg, last_gene) {
                        skipped_count += 1;
                        place_site(cg, last_gene, window_size, &mut windows)?;
                        continue;
                    }

                    let mut results: Vec<&GbmGene> = Vec::new(); // Collection of closest genes on all chromosomes and strands
                    for chromosome in genome.iter() {
                        for strand in [&chromosome.sense, &chromosome.antisense] {
                            let result =
                                strand.binary_search_by_key(&cg.location, |gene| gene.start);

                            let gene = match result {
                                Ok(i) => &strand[i],
                                Err(i) => {
                                    if strand.len() > i + 1 {
                                        &strand[i + 1]
                                    } else {
                                        &strand[i - 1] // Ugly hack #2, when cg site outside of chromosome range, jsut return any gene that will be thrown out in the next step
                                    }
                                }
                            };
                            results.push(gene);
                        }
                    }

                    for gene in results.iter() {
                        searched_count += 1;
                        if cg_site_in_gene(&cg, gene) {
                            last_gene = &gene;
                            place_site(cg, gene, window_size, &mut windows)?;
                            break;
                        }
                    }
                }
            }
            println!(
                "Done with {}! Searched through a total of {} genes, skipped {} times",
                filename.to_str().unwrap(),
                searched_count,
                skipped_count
            );
            write_windows(windows, output_dir, &filename, window_size as usize)?;
            Ok(())
        },
    )
}

fn place_site(
    cg: CgSite,
    gene: &GbmGene,
    window_size: u32,
    windows: &mut Vec<Vec<CgSite>>,
) -> Result<()> {
    // Offset from start for + strand, offset from end for - strand
    let offset = match &cg.strand {
        Strand::Sense => cg.location - gene.start,
        Strand::Antisense => gene.end - cg.location,
    };

    let mut percentile = offset * 100 / (gene.end - gene.start);
    // very unelegant fix but otherwise the last cg-site would land in its own category
    if percentile != 0 {
        percentile -= 1
    }
    let window = (percentile / window_size) as usize; // integer division rounds down
    windows[window].push(cg);
    Ok(())
}

fn set_up_output_dir(
    output_path: &str,
    overwrite: bool,
    window_count: usize,
    window_size: usize,
) -> Result<()> {
    fs::read_dir(&output_path).map_err(|_| {
        Error::FileError(String::from("Output directory"), String::from(output_path))
    })?; // Throw error if base output dir does not exist

    if overwrite {
        fs::remove_dir_all(&output_path).unwrap();
        fs::create_dir(&output_path).unwrap();
    }
    let edgelist = String::from(
        "from to
0_0 1_2
1_2 2_2
2_2 3_2
3_2 4_2
4_2 5_2
5_2 6_2
6_2 7_2
7_2 8_2
8_2 9_2
9_2 10_2
10_2 11_2
0_0 1_8
1_8 2_8
2_8 3_8
3_8 4_8
4_8 5_8
5_8 6_8
6_8 7_8
7_8 8_8
8_8 9_8
9_8 10_8
10_8 11_8",
    );
    for i in 0..window_count {
        let window = i * window_size;
        let nodelist = format!(
            "filename,node,gen,meth
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G0_All.txt,0_0,0,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G1_L2_All.txt,1_2,1,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G1_L8_All.txt,1_8,1,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G2_L2_All.txt,2_2,2,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G2_L8_All.txt,2_8,2,Y
-,3_2,3,N
-,3_8,3,N
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G4_L2_All.txt,4_2,4,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G4_L8_All.txt,4_8,4,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G5_L2_All.txt,5_2,5,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G5_L8_All.txt,5_8,5,Y
-,6_2,6,N
-,6_8,6,N
-,7_2,7,N
-,7_8,7,N
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G8_L2_All.txt,8_2,8,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G8_L8_All.txt,8_8,8,Y
-,9_2,9,N
-,9_8,9,N
-,10_2,10,N
-,10_8,10,N
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G11_L2_All.txt,11_2,11,Y
/mnt/extStorage/constantin/windows/{window}/methylome_Col0_G11_L8_All.txt,11_8,11,Y
"
        );

        let path = format!("{}/{}", &output_path, window);
        let window_dir = fs::read_dir(&path);

        match window_dir {
            Ok(_) => continue,
            Err(_) => {
                fs::create_dir(&path).unwrap();
                fs::write(path.to_owned() + "/nodelist.fn", nodelist)
                    .expect("Nodelist not writable at ");
                fs::write(path.to_owned() + "/edgelist.fn", &edgelist).expect("msg");
            }
        };
    }
    Ok(())
}

fn write_windows(
    windows: Vec<Vec<CgSite>>,
    output_dir: &str,
    filename: &OsString,
    window_size: usize,
) -> Result<()> {
    for (window, cg_sites) in windows.iter().enumerate() {
        let output_file = format!(
            "{}/{}/{}",
            output_dir,
            window * window_size,
            filename.to_str().unwrap()
        );
        let mut file = fs::OpenOptions::new()
            .append(true)
            .create(true)
            .open(&output_file)
            .expect(&format!("Could not open output file, {output_file}"));
        let metadata = file.metadata();
        if metadata.unwrap().len() == 0 {
            // On first write to file, create header line
            file.write("seqnames\tstart\tstrand\tcontext\tcounts.methylated\tcounts.total\tposteriorMax\tstatus\trc.meth.lvl\n".as_bytes())?;
        }
        file.write(cg_sites.iter().map(|e| &e.original).join("\n").as_bytes())?;
    }
    Ok(())
}
