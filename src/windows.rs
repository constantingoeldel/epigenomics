use std::{
    fs::File,
    io::{self, BufRead},
};

use crate::*;

pub fn extract_windows(
    methylome_file: File,
    genome: Vec<GenesByStrand>,
    window_size: i64,
    window_step: i64,
    window_count: usize,
    ignore_strand: bool,
) -> Result<Windows> {
    let mut last_gene = &genome[0].sense[0];
    // let mut last_gene_index = 0;
    let mut searched_count: i64 = 0;
    let mut skipped_count: i64 = 0;
    let mut windows = Windows::new(window_count);

    let lines = io::BufReader::new(methylome_file).lines();
    for (i, line_result) in lines.enumerate() {
        if let Ok(line) = line_result {
            if i % 1000000 == 0 {
                println!("Done with methylation site {i}");
            }

            let Some(cg) = MethylationSite::from_methylome_file_line(&line) else {continue;}; // If cg site could not be extracted, continue with the next line. Happens on header rows, for example.
            let cg_site_in_gene = |cg: &MethylationSite, gene: &Gene| {
                cg.chromosome == gene.chromosome
                    && gene.start <= cg.location
                    && cg.location <= gene.end
                    && (ignore_strand || cg.strand == gene.strand)
            };

            if cg_site_in_gene(&cg, last_gene) {
                skipped_count += 1;
                place_site(&cg, last_gene, window_size, window_step, &mut windows)?;
                continue;
            }

            let mut results: Vec<&Gene> = Vec::new(); // Collection of closest genes on all chromosomes and strands
            for chromosome in genome.iter() {
                for strand in [&chromosome.sense, &chromosome.antisense] {
                    let result = strand.binary_search_by_key(&cg.location, |gene| gene.start);

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
                    place_site(&cg, gene, window_size, window_step, &mut windows)?;
                    break;
                }
            }
        }
    }
    println!(
        "Done! Searched through a total of {} genes, skipped {} times",
        searched_count, skipped_count
    );
    Ok(windows)
}

fn place_site(
    cg: &MethylationSite,
    gene: &Gene,
    window_size: i64,
    window_step: i64,
    windows: &mut Windows,
) -> Result<()> {
    let gene_length = (gene.end - gene.start) as i64;
    // Offset from start for + strand, offset from end for - strand. Can be negative for upstream sites
    let offset: i64 = match &cg.strand {
        Strand::Sense => (cg.location - gene.start) as i64,
        Strand::Antisense => (gene.end - cg.location) as i64,
    };

    let normalized_offset = offset.abs() % gene_length;

    for window in (0..100).step_by(window_step as usize) {
        if normalized_offset > window * window_step
            && normalized_offset < window * window_step + window_size
        {
            // Clone necessary as each site can be part of multiple windows if step of sliding window < window size
            if offset < 0 {
                windows.upstream[((window - 1) * -1) as usize].push(cg.clone())
            } else if offset / gene_length >= 1 {
                windows.downstream[(window - 100) as usize].push(cg.clone())
            } else {
                windows.gene[window as usize].push(cg.clone())
            }
        }
    }

    Ok(())
}

// #[cfg(test)]
// mod tests {
//     use super::{extract_windows, place_site};
//     use crate::{structs::Windows, *};

//     #[test]
//     fn test_place_site() {
//         const WINDOW_SIZE: u32 = 5;

//         let mut windows = Windows::new(20);

//         let cg_a = MethylationSite {
//             chromosome: 1,
//             location: 80,
//             strand: Strand::Sense,
//             original: String::new(),
//         };
//         let cg_b = MethylationSite {
//             chromosome: 1,
//             location: 100,
//             strand: Strand::Sense,
//             original: String::new(),
//         };
//         let cg_c = MethylationSite {
//             chromosome: 1,
//             location: 123,
//             strand: Strand::Sense,
//             original: String::new(),
//         };
//         let cg_d = MethylationSite {
//             chromosome: 1,
//             location: 200,
//             strand: Strand::Sense,
//             original: String::new(),
//         };
//         let cg_e = MethylationSite {
//             chromosome: 1,
//             location: 220,
//             strand: Strand::Sense,
//             original: String::new(),
//         };

//         let gene = Gene {
//             chromosome: 1,
//             start: 100,
//             end: 200,
//             strand: Strand::Sense,
//             name: String::new(),
//         };

//         assert_eq!(
//             place_site(&cg, &gene, WINDOW_SIZE, &mut windows).unwrap(),
//             ()
//         );
//         assert!(windows[3].contains(&cg))
//         // assert_eq!(calculate_total(10, 10), 20);
//         // assert_eq!(calculate_total(-5, 10), 5);
//     }
// }
