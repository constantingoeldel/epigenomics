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

    pub fn place_in_windows(&self, gene: &Gene, windows: &mut Windows, args: &Args) -> Result<()> {
        let gene_length = gene.end - gene.start;
        // Offset from start for + strand, offset from end for - strand. Can be negative for upstream sites
        let offset: i32 = match &self.strand {
            Strand::Sense => self.location - gene.start,
            Strand::Antisense => gene.end - self.location,
        };

        let max = if args.absolute { gene_length } else { 100 };
        let mut normalized_offset = offset as f32; // In basepairs, absolute distance to 5'
        if normalized_offset < 0.0 {
            normalized_offset *= -1.0;
        } else if normalized_offset > gene_length as f32 {
            normalized_offset -= gene_length as f32;
        }
        if !args.absolute && (offset < 0 || offset > max) {
            // CGs at the end of the 2kb should be placed in the last window of upstream/downstream if in relative mode, even if the gene is longer/shorter than 2kb
            normalized_offset *= gene_length as f32 / args.cutoff as f32;
        }

        if offset == gene_length {
            normalized_offset = max as f32
        }

        for window in 0..max / args.window_step {
            let lower_bound = window * args.window_step;
            let upper_bound = lower_bound + args.window_size;
            if normalized_offset >= lower_bound as f32 && normalized_offset <= upper_bound as f32
            // Check that site is in window
            {
                // Clone necessary as each site can be part of multiple windows if step of sliding window < window size
                if offset < 0 {
                    windows.upstream[(window) as usize].push(self.clone())
                } else if offset > max {
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
            location: 512 + 100 + 100,
            strand: Strand::Sense,
            original: String::new(),
        };
        let cg_g = MethylationSite {
            chromosome: 1,
            location: 1024 + 100 + 100,
            strand: Strand::Sense,
            original: String::new(),
        };
        let cg_h = MethylationSite {
            chromosome: 1,
            location: 2048 + 100 + 100,
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

        let args = Args {
            absolute: false,
            cutoff: 2048,
            force: true,
            genome: String::from("not relevant"),
            ignore_strand: false,
            methylome: String::from("also not relevant"),
            output_dir: String::from("also not relevant"),
            window_size: 2,
            window_step: 1,
        };
        let mut windows = Windows::new(1000, &args);

        cg_a.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_b.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_c.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_d.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_e.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_f.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_g.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_h.place_in_windows(&gene, &mut windows, &args).unwrap();

        println!("{}", windows);
        assert!(windows.upstream[0].contains(&cg_a));
        assert!(windows.gene[0].contains(&cg_b));
        assert!(windows.gene[21].contains(&cg_c));
        assert!(windows.gene[22].contains(&cg_c));
        assert!(windows.gene[23].contains(&cg_c));
        assert!(windows.gene[98].contains(&cg_d));
        assert!(windows.gene[99].contains(&cg_d));
        assert!(windows.downstream[0].contains(&cg_e));
        assert!(windows.downstream[24].contains(&cg_f));
        assert!(windows.downstream[49].contains(&cg_g));
        assert!(windows.downstream[99].contains(&cg_h));
    }
    #[test]
    fn test_place_site_absolute() {
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

        let args = Args {
            absolute: true,
            cutoff: 2048,
            force: true,
            genome: String::from("not relevant"),
            ignore_strand: false,
            methylome: String::from("also not relevant"),
            output_dir: String::from("also not relevant"),
            window_size: 2,
            window_step: 1,
        };
        let mut windows = Windows::new(100, &args);

        cg_a.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_b.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_c.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_d.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_e.place_in_windows(&gene, &mut windows, &args).unwrap();
        cg_f.place_in_windows(&gene, &mut windows, &args).unwrap();

        println!("{}", windows);
        assert!(windows.upstream[2].contains(&cg_a));
        assert!(windows.upstream[3].contains(&cg_a));
        assert!(windows.upstream[4].contains(&cg_a));
        assert!(windows.gene[0].contains(&cg_b));
        assert!(windows.gene[3].contains(&cg_c));
        assert!(windows.gene[4].contains(&cg_c));
        assert!(windows.gene[18].contains(&cg_d));
        assert!(windows.gene[19].contains(&cg_d));
        assert!(windows.downstream[0].contains(&cg_e));
        assert!(windows.downstream[2].contains(&cg_f));
        assert!(windows.downstream[3].contains(&cg_f));
        assert!(windows.downstream[4].contains(&cg_f));
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
