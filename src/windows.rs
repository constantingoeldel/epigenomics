use std::{
    fmt::Display,
    fs::{File, OpenOptions},
    io::{self, BufRead, Write},
};

use itertools::Itertools;

use crate::*;

pub type Window = Vec<MethylationSite>;
#[derive(Debug, Eq, PartialEq)]
pub struct Windows {
    pub upstream: Vec<Window>,
    pub gene: Vec<Window>,
    pub downstream: Vec<Window>,
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
        for windows in vec![
            (&self.upstream, "upstream"),
            (&self.gene, "gene"),
            (&self.downstream, "downstream"),
        ]
        .iter()
        {
            for (window, cg_sites) in windows.0.iter().enumerate() {
                let output_file = format!(
                    "{}/{}/{}/{}",
                    output_dir,
                    windows.1,
                    window * step,
                    filename.to_str().unwrap()
                );
                let mut file = OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(&output_file)?;
                // .map_err(|_| {
                //     Error::File(output_file, output_dir.to_owned());
                // })?;
                let metadata = file.metadata();
                if metadata.unwrap().len() == 0 {
                    // On first write to file, create header line
                    file.write_all("seqnames\tstart\tstrand\tcontext\tcounts.methylated\tcounts.total\tposteriorMax\tstatus\trc.meth.lvl\n".as_bytes())?;
                }
                file.write_all(cg_sites.iter().map(|e| &e.original).join("\n").as_bytes())?;
            }
        }
        Ok(())
    }
}

pub fn extract_windows(
    methylome_file: File,
    genome: Vec<GenesByStrand>,
    window_size: i32,
    window_step: i32,
    max_gene_length: i32,
    ignore_strand: bool,
    cutoff: i32,
    absolute: bool,
) -> Result<Windows> {
    let mut last_gene: Gene = genome[0].sense[0].clone();

    let mut windows = Windows::new((max_gene_length / window_step) as usize);

    let lines = io::BufReader::new(methylome_file).lines();
    for (i, line_result) in lines.enumerate().skip(1) {
        // skip header row
        if let Ok(line) = line_result {
            if i % 100_000 == 0 {
                println!("Done with methylation site {i} ");
            }

            // If cg site could not be extracted, continue with the next line. Happens on header rows, for example.
            let Some(cg) = MethylationSite::from_methylome_file_line(&line) else {continue;};
            if cg.is_in_gene(&last_gene, ignore_strand, cutoff) {
                cg.place_in_windows(&last_gene, window_size, window_step, &mut windows, absolute)?;
                continue;
            }
            let gene = cg.find_gene(&genome, ignore_strand, cutoff);
            if let Some(gene) = gene {
                last_gene = gene.clone();
                cg.place_in_windows(gene, window_size, window_step, &mut windows, absolute)?;
            }
        }
    }
    Ok(windows)
}

impl Display for Windows {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Upstream: {:?}\n\nGene: {:?}\n\nDownstream: {:?}\n\n",
            self.upstream, self.gene, self.downstream
        )
    }
}
