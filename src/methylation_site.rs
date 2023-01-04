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
    pub fn from_methylome_file_line(s: &str) -> Result<Self> {
        s.split('\t')
            .collect_tuple()
            .filter(|(_, _, _, context, _, _, _, _, _)| context == &"CG")
            .map(|(chromosome, location, strand, _, _, _, _, _, _)| {
                Ok(MethylationSite {
                    chromosome: chromosome.parse::<u8>()?,
                    location: location.parse::<i32>()?,
                    strand: if strand == "+" {
                        Strand::Sense
                    } else {
                        Strand::Antisense
                    },
                    original: s.to_owned(),
                })
            })
            .ok_or(Error::CGSite)?
    }
    pub fn is_in_gene(&self, gene: &Gene, cutoff: i32) -> bool {
        self.chromosome == gene.chromosome
            && gene.start <= self.location + cutoff
            && self.location <= gene.end + cutoff
            && self.strand == gene.strand
    }

    pub fn find_gene<'short, 'long>(
        // The lifetime of the genome is longer than the lifetime of the cg site
        // cg sites exist only while a single methylation file is being processed
        // the genome is loaded once and exists for the entire program
        &'short self,
        genome: &'long [GenesByStrand],
        cutoff: i32,
    ) -> Option<&'long Gene> {
        let chromosome = &genome[(self.chromosome - 1) as usize];
        let strand = match self.strand {
            Strand::Sense => &chromosome.sense,
            Strand::Antisense => &chromosome.antisense,
        };
        let first_matching_gene_index = strand
            .binary_search_by_key(&self.location, |gene| gene.end + cutoff)
            .unwrap_or_else(|x| x); // Collapse exact match on gene end and closest previous match into one, as both are valid
        if strand.len() < first_matching_gene_index + 1 {
            return None;
        }

        let gene = &strand[first_matching_gene_index];
        if self.is_in_gene(gene, cutoff) {
            return Some(gene);
        }

        None
    }

    pub fn place_in_windows(
        &self,
        gene: &Gene,
        windows: &mut Windows,
        args: &Args,
    ) -> Vec<(Region, usize)> // Return a vector of (strand, window) tuples for each window the CG site is in
    {
        const E: f32 = 0.1; // Epsilon for floating point comparison
        let location = self.location as f32;
        let cutoff = args.cutoff as f32;
        let step = args.window_step as f32;
        let size = args.window_size as f32;
        let start = gene.start as f32;
        let end = gene.end as f32;
        let length = end - start;

        //
        // // Offset from start for + strand, offset from end for - strand. Can be negative for upstream sites
        let offset = match &self.strand {
            Strand::Sense => location - start,
            Strand::Antisense => end - location,
        };
        let mut windows_in = Vec::new();

        let region = match offset {
            x if x < 0.0 => Region::Upstream,
            x if x >= length => Region::Downstream, // CG site exactly on the end of the gene is considered downstream
            _ => Region::Gene,
        };
        let local_windows = windows.get_mut(&region);

        // let max = if args.absolute { gene_length } else { 100 };
        let mut position = match region {
            // Position within the region of the gene
            Region::Upstream => location - start + cutoff,
            Region::Gene => location - start,
            Region::Downstream => location - end,
        };

        if !args.absolute {
            position = match region {
                Region::Upstream => position / cutoff,
                Region::Gene => position / length,
                Region::Downstream => position / cutoff,
            };
            position *= 100.0; // Normalize to 0-100%
        }

        for (i, window) in local_windows.iter_mut().enumerate() {
            let lower_bound = i as f32 * step - E;
            let upper_bound = lower_bound + size + E;

            if position >= lower_bound && position <= upper_bound {
                window.push(self.clone());
                windows_in.push((region.clone(), i));
            }
        }
        windows_in
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

#[cfg(test)]
mod tests {

    use super::MethylationSite;
    use crate::*;

    const GENE: Gene = Gene {
        chromosome: 1,
        start: 50,
        end: 100,
        strand: Strand::Sense,
        name: String::new(),
    };
    const WITHIN_CG: MethylationSite = MethylationSite {
        chromosome: 1,
        location: 80,
        strand: Strand::Sense,
        original: String::new(),
    };

    const OPPOSITE_STRAND_CG: MethylationSite = MethylationSite {
        chromosome: 1,
        location: 80,
        strand: Strand::Antisense,
        original: String::new(),
    };

    const HIGHER_CG: MethylationSite = MethylationSite {
        chromosome: 1,
        location: 150,
        strand: Strand::Sense,
        original: String::new(),
    };
    const LOWER_CG: MethylationSite = MethylationSite {
        chromosome: 1,
        location: 0,
        strand: Strand::Sense,
        original: String::new(),
    };

    #[test]
    fn test_instantiate_from_methylome_file_line() {
        let line = "1	23151	+	CG	0	8	0.9999	U	0.0025";
        let cg = MethylationSite::from_methylome_file_line(line).unwrap();
        assert_eq!(cg.chromosome, 1);
    }

    #[test]
    fn test_instantiate_from_methylome_file_line_invalid_line() {
        let line = "1	23151	+	CG	0	8	0.9999	";
        let cg = MethylationSite::from_methylome_file_line(line);
        assert!(cg.is_err());
    }

    #[test]
    fn test_instantiate_from_methylome_file_line_invalid_chromosome() {
        let line = "X	23151	+	CG	0	8	0.9999	U	0.0025";
        let cg = MethylationSite::from_methylome_file_line(line);
        assert!(cg.is_err());
    }

    #[test]
    fn test_is_in_gene() {
        assert!(WITHIN_CG.is_in_gene(&GENE, 0));
        assert!(!HIGHER_CG.is_in_gene(&GENE, 0));
        assert!(HIGHER_CG.is_in_gene(&GENE, 50));
        assert!(!LOWER_CG.is_in_gene(&GENE, 0));
        assert!(LOWER_CG.is_in_gene(&GENE, 50));
    }

    #[test]
    fn test_find_gene() {
        let mut genes = GenesByStrand::new();
        for i in 0..100 {
            genes.insert(Gene {
                chromosome: 1,
                start: i,
                end: i + 50,
                strand: Strand::Sense,
                name: String::new(),
            });
        }

        let genome = vec![genes.clone()];
        assert!(OPPOSITE_STRAND_CG.find_gene(&genome, 0).is_none());
        assert_eq!(
            Some(WITHIN_CG.find_gene(&genome, 0)),
            Some(genes.sense.get(30))
        );
        assert_eq!(
            Some(HIGHER_CG.find_gene(&genome, 50)),
            Some(genes.sense.get(50))
        );
        assert_eq!(
            Some(LOWER_CG.find_gene(&genome, 50)),
            Some(genes.sense.get(0))
        );
    }

    #[test]
    fn test_extract_gene() {}

    #[test]
    fn test_place_site_absolute() {
        let args = Args {
            absolute: true,
            cutoff: 1000,
            genome: String::from("not relevant"),
            methylome: String::from("also not relevant"),
            output_dir: String::from("also not relevant"),
            window_size: 2,
            window_step: 1,
        };
        let all_within_gene = Gene {
            chromosome: 1,
            start: 1000,
            end: 2000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            chromosome: 1,
            start: 2000,
            end: 3000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            chromosome: 1,
            start: 0,
            end: 1000,
            strand: Strand::Sense,
            name: String::new(),
        };

        let mut windows = Windows::new(1000, &args);
        for i in 0..1000 {
            let cg = MethylationSite {
                chromosome: 1,
                location: i + 1000,
                strand: Strand::Sense,
                original: String::new(),
            };
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

            println!("Placing {}", i);
            println!("Upstream: {:?}", upstream);
            println!("Gene: {:?}", gene);
            println!("Downstream: {:?}", downstream);
            assert!(windows.upstream[i as usize].contains(&cg));
            assert!(windows.gene[i as usize].contains(&cg));
            assert!(windows.downstream[i as usize].contains(&cg));

            if i > 3 {
                assert!(windows.upstream[i as usize - 1].contains(&cg));
                assert!(windows.gene[i as usize - 1].contains(&cg));
                assert!(windows.downstream[i as usize - 1].contains(&cg));
                assert!(windows.upstream[i as usize - 2].contains(&cg));
                assert!(windows.gene[i as usize - 2].contains(&cg));
                assert!(windows.downstream[i as usize - 2].contains(&cg));
            }
        }
    }
    #[test]
    fn test_place_site_relative_acting_like_absolute() {
        let args = Args {
            absolute: false,
            cutoff: 100,
            genome: String::from("not relevant"),
            methylome: String::from("also not relevant"),
            output_dir: String::from("also not relevant"),
            window_size: 2,
            window_step: 1,
        };
        let all_within_gene = Gene {
            chromosome: 1,
            start: 100,
            end: 200,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            chromosome: 1,
            start: 200,
            end: 300,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            chromosome: 1,
            start: 0,
            end: 100,
            strand: Strand::Sense,
            name: String::new(),
        };

        let mut windows = Windows::new(100, &args);
        for i in 0..100 {
            let cg = MethylationSite {
                chromosome: 1,
                location: i + 100,
                strand: Strand::Sense,
                original: String::new(),
            };
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

            println!("Placing {}", i);
            println!("Upstream: {:?}", upstream);
            println!("Gene: {:?}", gene);
            println!("Downstream: {:?}", downstream);
            assert!(windows.upstream[i as usize].contains(&cg));
            assert!(windows.gene[i as usize].contains(&cg));
            assert!(windows.downstream[i as usize].contains(&cg));
        }
    }
    #[test]
    fn test_place_site_relative() {
        let args = Args {
            absolute: false,
            cutoff: 1000,
            genome: String::from("not relevant"),
            methylome: String::from("also not relevant"),
            output_dir: String::from("also not relevant"),
            window_size: 2,
            window_step: 1,
        };
        let all_within_gene = Gene {
            chromosome: 1,
            start: 1000,
            end: 2000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            chromosome: 1,
            start: 2000,
            end: 3000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            chromosome: 1,
            start: 0,
            end: 1000,
            strand: Strand::Sense,
            name: String::new(),
        };

        let mut windows = Windows::new(1000, &args);
        assert!(windows.upstream.len() == 100);
        for i in 0..1000 {
            let cg = MethylationSite {
                chromosome: 1,
                location: i + 1000,
                strand: Strand::Sense,
                original: String::new(),
            };
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

            println!("Placing {}", i);
            println!("Upstream: {:?}", upstream);
            println!("Gene: {:?}", gene);
            println!("Downstream: {:?}", downstream);
            println!("{}: {}", i / 10, windows.upstream[(i / 10) as usize].len());
            assert!(windows.upstream[(i / 10) as usize].contains(&cg));
            assert!(windows.gene[(i / 10) as usize].contains(&cg));
            assert!(windows.downstream[(i / 10) as usize].contains(&cg));
        }
    }

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
            genome: String::from("not relevant"),
            methylome: String::from("also not relevant"),
            output_dir: String::from("also not relevant"),
            window_size: 2,
            window_step: 1,
        };
        let mut windows = Windows::new(1000, &args);

        cg_a.place_in_windows(&gene, &mut windows, &args);
        cg_b.place_in_windows(&gene, &mut windows, &args);
        cg_c.place_in_windows(&gene, &mut windows, &args);
        cg_d.place_in_windows(&gene, &mut windows, &args);
        cg_e.place_in_windows(&gene, &mut windows, &args);
        cg_f.place_in_windows(&gene, &mut windows, &args);
        cg_g.place_in_windows(&gene, &mut windows, &args);
        cg_h.place_in_windows(&gene, &mut windows, &args);

        println!("{}", windows);
        assert!(windows.upstream[98].contains(&cg_a));
        assert!(windows.upstream[99].contains(&cg_a));
        assert!(windows.gene[0].contains(&cg_b));
        assert!(windows.gene[21].contains(&cg_c));
        assert!(windows.gene[22].contains(&cg_c));
        assert!(windows.gene[23].contains(&cg_c));
        assert!(windows.downstream[0].contains(&cg_d));
        assert!(windows.downstream[0].contains(&cg_e));
        assert!(windows.downstream[24].contains(&cg_f));
        assert!(windows.downstream[49].contains(&cg_g));
        assert!(windows.downstream[99].contains(&cg_h));
    }
    #[test]
    fn test_place_site_absolute_2() {
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
            genome: String::from("not relevant"),
            methylome: String::from("also not relevant"),
            output_dir: String::from("also not relevant"),
            window_size: 2,
            window_step: 1,
        };
        let mut windows = Windows::new(100, &args);

        cg_a.place_in_windows(&gene, &mut windows, &args);
        cg_b.place_in_windows(&gene, &mut windows, &args);
        cg_c.place_in_windows(&gene, &mut windows, &args);
        cg_d.place_in_windows(&gene, &mut windows, &args);
        cg_e.place_in_windows(&gene, &mut windows, &args);
        cg_f.place_in_windows(&gene, &mut windows, &args);

        println!("{}", windows);
        // assert!(windows.upstream[2].contains(&cg_a));
        // assert!(windows.upstream[3].contains(&cg_a));
        // assert!(windows.upstream[4].contains(&cg_a));
        assert!(windows.gene[0].contains(&cg_b));
        assert!(windows.gene[21].contains(&cg_c));
        assert!(windows.gene[22].contains(&cg_c));
        assert!(windows.gene[23].contains(&cg_c));
        assert!(windows.downstream[0].contains(&cg_d));
        assert!(windows.downstream[0].contains(&cg_e));
        assert!(windows.downstream[1].contains(&cg_e));
        assert!(windows.downstream[18].contains(&cg_f));
        assert!(windows.downstream[19].contains(&cg_f));
        assert!(windows.downstream[20].contains(&cg_f));
    }
}
