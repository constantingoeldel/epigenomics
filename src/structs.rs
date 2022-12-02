use std::{ffi::OsString, fmt::Display, fs::OpenOptions, io::Write};

use itertools::Itertools;

use crate::error;

pub type Result<T> = std::result::Result<T, error::Error>;

#[derive(Clone, PartialEq)]
pub enum Strand {
    Sense,
    Antisense,
}

#[derive(Clone)]
pub struct Gene {
    pub chromosome: u8,
    pub start: u32,
    pub end: u32,
    pub name: String,
    pub strand: Strand,
}
#[derive(Clone, PartialEq)]
pub struct MethylationSite {
    pub chromosome: u8,
    pub location: u32,
    pub strand: Strand,
    pub original: String,
}
#[derive(Clone)]
pub struct GenesByStrand {
    pub sense: Vec<Gene>,
    pub antisense: Vec<Gene>,
}
pub type Window = Vec<MethylationSite>;

pub struct Windows {
    pub upstream: Vec<Window>,
    pub gene: Vec<Window>,
    pub downstream: Vec<Window>,
}

impl Gene {
    pub fn from_annotation_file_line(s: &str) -> Option<Self> {
        s.split('\t')
            .collect_tuple()
            .map(|(chromosome, start, end, name, _, strand)| Gene {
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

impl MethylationSite {
    pub fn from_methylome_file_line(s: &str) -> Option<Self> {
        s.split('\t')
            .collect_tuple()
            .filter(|(chromosome, _, _, _, _, _, _, _, _)| chromosome != &"seqnames") // Filter out header row
            .map(
                |(chromosome, location, strand, _, _, _, _, _, _)| MethylationSite {
                    chromosome: chromosome.parse::<u8>().unwrap(),
                    location: location.parse::<u32>().unwrap(),
                    strand: if strand == "+" {
                        Strand::Sense
                    } else {
                        Strand::Antisense
                    },
                    original: s.to_owned().clone(),
                },
            )
    }
}

impl Windows {
    pub fn new(window_count: usize) -> Self {
        Windows {
            upstream: vec![Vec::new(); window_count],
            gene: vec![Vec::new(); window_count],
            downstream: vec![Vec::new(); window_count],
        }
    }
    pub fn save(&self, output_dir: &str, filename: &OsString, step: usize) -> Result<()> {
        for windows in vec![&self.upstream, &self.gene, &self.downstream].iter() {
            for (window, cg_sites) in windows.iter().enumerate() {
                let output_file = format!(
                    "{}/{}/{}",
                    output_dir,
                    window * step,
                    filename.to_str().unwrap()
                );
                let mut file = OpenOptions::new()
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
        }
        Ok(())
    }
}

impl Display for Gene {
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

impl Display for MethylationSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CG site is located on the {} strand of chromosome {} at bp {}",
            self.strand, self.chromosome, self.location
        )
    }
}
