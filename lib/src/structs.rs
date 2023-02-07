use std::fmt::Display;

use itertools::Itertools;

use crate::error;

pub type Result<T> = std::result::Result<T, error::Error>;

#[derive(Clone, Debug, sqlx::Type)]
#[sqlx(type_name = "strandness")]
pub enum Strand {
    #[sqlx(rename = "+")]
    Sense,
    #[sqlx(rename = "-")]
    Antisense,
    #[sqlx(rename = "*")]
    Unknown,
}

impl PartialEq for Strand {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Sense, Self::Antisense) => false,
            (Self::Antisense, Self::Sense) => false,
            _ => true, // If one of them is unknown, they are equal. WARNING: This leads to unexpected behaviour if you want to check for one specific strand! in that case use `if let Strand::X = self.strand {}
        }
    }
}

impl Eq for Strand {}

#[derive(Debug, Clone, sqlx::Type)]
#[sqlx(type_name = "region")]
pub enum Region {
    Upstream,
    Gene,
    Downstream,
}

impl Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Upstream => write!(f, "upstream"),
            Self::Gene => write!(f, "gene"),
            Self::Downstream => write!(f, "downstream"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gene {
    pub chromosome: u8,
    pub start: u32,
    pub end: u32,
    pub name: String,
    pub annotation: String,
    pub strand: Strand,
}

#[derive(Clone)]
pub struct GenesByStrand {
    pub sense: Vec<Gene>,
    pub antisense: Vec<Gene>,
}

impl Default for GenesByStrand {
    fn default() -> Self {
        GenesByStrand::new()
    }
}

impl GenesByStrand {
    pub fn new() -> Self {
        GenesByStrand {
            sense: Vec::new(),
            antisense: Vec::new(),
        }
    }

    pub fn chain(&self) -> Vec<&Gene> {
        let sense_iter = self.sense.iter();
        let antisense_iter = self.antisense.iter();

        sense_iter.chain(antisense_iter).collect()
    }

    pub fn insert(&mut self, gene: Gene) {
        match gene.strand {
            Strand::Sense => self.sense.push(gene),
            Strand::Antisense => self.antisense.push(gene),
            _ => (),
        }
    }

    pub fn sort(&mut self) {
        self.sense.sort_by(|a, b| a.start.cmp(&b.start));
        self.antisense.sort_by(|a, b| a.start.cmp(&b.start));
    }
}

impl Gene {
    pub fn from_annotation_file_line(s: &str, invert_strand: bool) -> Option<Self> {
        s.split('\t')
            .collect_tuple()
            .map(|(chromosome, start, end, name, annotation, strand)| Gene {
                chromosome: chromosome.parse::<u8>().unwrap(),
                start: start.parse::<u32>().unwrap(),
                end: end.parse::<u32>().unwrap(),
                name: String::from(name),
                annotation: String::from(annotation),
                strand: if (strand == "+") ^ invert_strand {
                    // XOR: if the strand is + and we don't want to invert it, or if the strand is - and we do want to invert it -> Sense
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
            Strand::Unknown => write!(f, "*"),
        }
    }
}
