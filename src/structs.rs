use std::fmt::Display;

use itertools::Itertools;

use crate::error;

pub type Result<T> = std::result::Result<T, error::Error>;

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Strand {
    Sense,
    Antisense,
}

#[derive(Debug)]
pub enum Location {
    Upstream,
    Gene,
    Downstream,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gene {
    pub chromosome: u8,
    pub start: i32,
    pub end: i32,
    pub name: String,
    pub strand: Strand,
}

#[derive(Clone)]
pub struct GenesByStrand {
    pub sense: Vec<Gene>,
    pub antisense: Vec<Gene>,
}

impl GenesByStrand {
    pub fn new() -> Self {
        GenesByStrand {
            sense: Vec::new(),
            antisense: Vec::new(),
        }
    }

    pub fn insert(&mut self, gene: Gene) {
        match gene.strand {
            Strand::Sense => self.sense.push(gene),
            Strand::Antisense => self.antisense.push(gene),
        }
    }

    pub fn sort(&mut self) {
        self.sense.sort_by(|a, b| a.start.cmp(&b.start));
        self.antisense.sort_by(|a, b| a.start.cmp(&b.start));
    }
}

impl Gene {
    pub fn from_annotation_file_line(s: &str) -> Option<Self> {
        s.split('\t')
            .collect_tuple()
            .map(|(chromosome, start, end, name, _, strand)| Gene {
                chromosome: chromosome.parse::<u8>().unwrap(),
                start: start.parse::<i32>().unwrap(),
                end: end.parse::<i32>().unwrap(),
                name: String::from(name),
                strand: if strand == "+" {
                    Strand::Sense
                } else {
                    Strand::Antisense
                },
            })
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
