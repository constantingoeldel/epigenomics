use std::fmt::Display;

use itertools::Itertools;

use crate::*;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct MethylationSite {
    pub chromosome: u8,
    pub location: i32,
    pub strand: Strand,
    pub original: String,
}

impl MethylationSite {
    pub fn from_methylome_file_line(s: &str) -> Option<Self> {
        s.split('\t')
            .collect_tuple()
            .filter(|(_, _, _, context, _, _, _, _, _)| context == &"CG")
            .map(
                |(chromosome, location, strand, _, _, _, _, _, _)| MethylationSite {
                    chromosome: chromosome.parse::<u8>().unwrap(),
                    location: location.parse::<i32>().unwrap(),
                    strand: if strand == "+" {
                        Strand::Sense
                    } else {
                        Strand::Antisense
                    },
                    original: s.to_owned(),
                },
            )
    }
    pub fn is_in_gene(&self, gene: &Gene, ignore_strand: bool, cutoff: i32) -> bool {
        self.chromosome == gene.chromosome
            && gene.start <= self.location + cutoff
            && self.location <= gene.end + cutoff
            && (ignore_strand || self.strand == gene.strand)
    }

    pub fn find_gene<'a>(
        &'a self,
        genome: &'a Vec<GenesByStrand>,
        ignore_strand: bool,
        cutoff: i32,
    ) -> Option<&Gene> {
        let chromosome = &genome[(self.chromosome - 1) as usize];
        let strand = match self.strand {
            Strand::Sense => &chromosome.sense,
            Strand::Antisense => &chromosome.antisense,
        };
        let result = strand.binary_search_by_key(&self.location, |gene| gene.end + cutoff);
        let gene = match result {
            Ok(i) => &strand[i],
            Err(i) => {
                if strand.len() > i + 1 {
                    &strand[i]
                } else {
                    &strand[i - 1] // Ugly hack #2, when cg site outside of chromosome range, jsut return any gene that will be thrown out in the next step
                }
            }
        };
        if self.is_in_gene(gene, ignore_strand, cutoff) {
            return Some(gene);
        }

        None
    }

    pub fn place_in_windows(
        &self,
        gene: &Gene,
        window_size: i32,
        window_step: i32,
        windows: &mut Windows,
    ) -> Result<()> {
        let gene_length = gene.end - gene.start;
        // Offset from start for + strand, offset from end for - strand. Can be negative for upstream sites
        let offset: i32 = match &self.strand {
            Strand::Sense => self.location - gene.start,
            Strand::Antisense => gene.end - self.location,
        };

        let mut normalized_offset = offset.abs() % (gene_length);

        if offset == gene_length {
            normalized_offset = 100
        }

        // TODO update for abolute values
        for window in 0..100 / window_step {
            if normalized_offset >= window * window_step
                && normalized_offset <= window * window_step + window_size
            {
                // Clone necessary as each site can be part of multiple windows if step of sliding window < window size
                if offset < 0 {
                    windows.upstream[(window) as usize].push(self.clone())
                } else if offset > 100 {
                    windows.downstream[(window) as usize].push(self.clone())
                } else {
                    windows.gene[window as usize].push(self.clone())
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {

    use super::MethylationSite;
    use crate::*;

    #[test]
    fn test_extract_gene() {}
    #[test]
    fn test_place_site() {
        const WINDOW_SIZE: i32 = 2;

        let mut windows = Windows::new(100_usize);

        let cg_a = MethylationSite {
            chromosome: 1,
            location: 80,
            strand: Strand::Sense,
            original: String::new(),
        };
        let cg_b = MethylationSite {
            chromosome: 1,
            location: 100,
            strand: Strand::Sense,
            original: String::new(),
        };
        let cg_c = MethylationSite {
            chromosome: 1,
            location: 123,
            strand: Strand::Sense,
            original: String::new(),
        };
        let cg_d = MethylationSite {
            chromosome: 1,
            location: 200,
            strand: Strand::Sense,
            original: String::new(),
        };
        let cg_e = MethylationSite {
            chromosome: 1,
            location: 201,
            strand: Strand::Sense,
            original: String::new(),
        };
        let cg_f = MethylationSite {
            chromosome: 1,
            location: 220,
            strand: Strand::Sense,
            original: String::new(),
        };

        let gene = Gene {
            chromosome: 1,
            start: 100,
            end: 200,
            strand: Strand::Sense,
            name: String::new(),
        };

        cg_a.place_in_windows(&gene, WINDOW_SIZE, WINDOW_SIZE / 2, &mut windows)
            .unwrap();
        cg_b.place_in_windows(&gene, WINDOW_SIZE, WINDOW_SIZE / 2, &mut windows)
            .unwrap();
        cg_c.place_in_windows(&gene, WINDOW_SIZE, WINDOW_SIZE / 2, &mut windows)
            .unwrap();
        cg_d.place_in_windows(&gene, WINDOW_SIZE, WINDOW_SIZE / 2, &mut windows)
            .unwrap();
        cg_e.place_in_windows(&gene, WINDOW_SIZE, WINDOW_SIZE / 2, &mut windows)
            .unwrap();
        cg_f.place_in_windows(&gene, WINDOW_SIZE, WINDOW_SIZE / 2, &mut windows)
            .unwrap();

        println!("{}", windows);
        assert!(windows.upstream[18].contains(&cg_a));
        assert!(windows.upstream[19].contains(&cg_a));
        assert!(windows.upstream[20].contains(&cg_a));
        assert!(windows.gene[0].contains(&cg_b));
        assert!(windows.gene[21].contains(&cg_c));
        assert!(windows.gene[22].contains(&cg_c));
        assert!(windows.gene[23].contains(&cg_c));
        assert!(windows.gene[98].contains(&cg_d));
        assert!(windows.gene[99].contains(&cg_d));
        assert!(windows.downstream[0].contains(&cg_e));
        assert!(windows.downstream[1].contains(&cg_e));
        assert!(windows.downstream[18].contains(&cg_f));
        assert!(windows.downstream[19].contains(&cg_f));
        assert!(windows.downstream[20].contains(&cg_f));
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
