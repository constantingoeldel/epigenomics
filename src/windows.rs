use std::{
    fmt::Display,
    fs::{File, OpenOptions},
    io::{self, BufRead, Write},
};

use itertools::Itertools;

use crate::*;

pub type Window = Vec<MethylationSite>;
#[derive(Debug, PartialEq)]
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

pub fn extract_windows(
    methylome_file: File,
    genome: Vec<GenesByStrand>,
    window_size: i32,
    window_step: i32,
    max_gene_length: i32,
    ignore_strand: bool,
) -> Result<Windows> {
    let mut last_gene = genome[0].sense[0].clone();

    let mut windows = Windows::new((max_gene_length / window_step) as usize);

    let lines = io::BufReader::new(methylome_file).lines();
    for (i, line_result) in lines.enumerate() {
        if let Ok(line) = line_result {
            if i % 1000000 == 0 {
                println!("Done with methylation site {i}");
            }

            // If cg site could not be extracted, continue with the next line. Happens on header rows, for example.
            let Some(cg) = MethylationSite::from_methylome_file_line(&line) else {continue;};
            let gene = cg.find_gene(&genome, &last_gene, ignore_strand)?;
            last_gene = gene.clone();
            cg.place_in_windows(&gene, window_size, window_step, &mut windows)?;
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
