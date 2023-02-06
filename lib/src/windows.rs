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
    pub fn new(max_gene_length: u32, args: &Args) -> Self {
        let gene_window_count = if args.absolute {
            max_gene_length / args.window_step
        } else {
            100
        };
        let up_down_window_count = if args.absolute {
            args.cutoff / args.window_step
        } else {
            100
        };
        Windows {
            upstream: vec![Vec::new(); up_down_window_count as usize],
            gene: vec![Vec::new(); gene_window_count as usize],
            downstream: vec![Vec::new(); up_down_window_count as usize],
        }
    }
    pub fn get(&self, region: Region) -> &Vec<Window> {
        match region {
            Region::Upstream => &self.upstream,
            Region::Gene => &self.gene,
            Region::Downstream => &self.downstream,
        }
    }
    pub fn get_mut<'a>(&'a mut self, location: &Region) -> &'a mut Vec<Window> {
        match location {
            Region::Upstream => &mut self.upstream,
            Region::Gene => &mut self.gene,
            Region::Downstream => &mut self.downstream,
        }
    }

    pub fn iter_upstream(
        &self,
    ) -> std::slice::Iter<'_, std::vec::Vec<methylation_site::MethylationSite>> {
        self.upstream.iter()
    }

    pub fn iter_gene(
        &self,
    ) -> std::slice::Iter<'_, std::vec::Vec<methylation_site::MethylationSite>> {
        self.gene.iter()
    }

    pub fn iter_downstream(
        &self,
    ) -> std::slice::Iter<'_, std::vec::Vec<methylation_site::MethylationSite>> {
        self.downstream.iter()
    }

    pub fn combined(&self) -> Vec<Vec<MethylationSite>> {
        let combined = vec![
            self.upstream.clone(),
            self.gene.clone(),
            self.downstream.clone(),
        ];
        combined.concat()
    }

    pub fn inverse(mut self) -> Self {
        self.upstream = self.downstream.iter().rev().map(|a| a.to_owned()).collect();
        self.gene = self.gene.iter().rev().map(|a| a.to_owned()).collect();
        self.downstream = self.upstream.iter().rev().map(|a| a.to_owned()).collect();
        self
    }

    pub fn steady_state_methylation(&self) -> Vec<f64> {
        self.combined()
            .iter()
            .map(|w| w.iter().fold(0.0, |acc, cur| acc + cur.meth_lvl) / w.len() as f64)
            .collect()
    }

    pub fn print_steady_state_methylation(methylations: &[f64]) -> String {
        let mut output = String::new();
        for (i, average) in methylations.iter().enumerate() {
            output.push_str(&format!("{i};{average}\n"));
        }
        output
    }

    pub fn distribution(&self) -> Vec<i32> {
        self.combined().iter().map(|w| w.len() as i32).collect()
    }

    pub fn print_distribution(distribution: &[i32]) -> String {
        // In CSV format
        let mut output = String::new();

        for (i, count) in distribution.iter().enumerate() {
            output.push_str(&format!("{i};{count}\n"));
        }

        output
    }

    pub fn save(&self, args: &Args, filename: &OsString) -> Result<()> {
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
                    args.output_dir,
                    windows.1,
                    window * args.window_step as usize,
                    filename.to_str().unwrap()
                );
                let mut file = OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(&output_file)?;

                let metadata = file.metadata();
                if metadata.unwrap().len() == 0 {
                    // On first write to file, create header line
                    file.write_all("seqnames\tstart\tstrand\tcontext\tcounts.methylated\tcounts.total\tposteriorMax\tstatus\trc.meth.lvl\tcontext.trinucleotide\n".as_bytes())?;
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
    max_gene_length: u32,
    args: Args,
) -> Result<Windows> {
    let mut last_gene: Option<&Gene> = None;

    let mut windows = Windows::new(max_gene_length, &args);

    let lines = io::BufReader::new(methylome_file).lines();
    for (i, line_result) in lines.enumerate().skip(1) {
        // skip header row
        if let Ok(line) = line_result {
            if i % 1_000_000 == 0 {
                println!("Done with methylation site {i} ");
            }

            // If cg site could not be extracted from a file line, continue with the next line. Happens on header rows, for example.
            let Ok(cg) = MethylationSite::from_methylome_file_line(&line, args.invert) else {continue;};
            if last_gene.is_none() || !cg.is_in_gene(last_gene.unwrap(), args.cutoff) {
                last_gene = cg.find_gene(&genome, args.cutoff);
            }
            if let Some(gene) = last_gene {
                cg.place_in_windows(gene, &mut windows, &args);
                continue;
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

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::arguments::Args;

    #[test]
    fn new_absolute() {
        let args = Args {
            edges: PathBuf::new(),
            nodes: PathBuf::new(),
            db: false,
            invert: false,
            methylome: "/home/constantin/methylome/within_gbM_genes".to_string(),
            genome: "/home/constantin/methylome/gbM_gene_anotation_extract_Arabidopsis.bed"
                .to_string(),
            window_size: 512,
            window_step: 256,

            output_dir: "/home/constantin/windows".to_string(),
            absolute: true,
            cutoff: 2048,
        };
        let windows = super::Windows::new(4096, &args);
        assert_eq!(windows.upstream.len(), 8);
        assert_eq!(windows.gene.len(), 16);
        assert_eq!(windows.downstream.len(), 8);
    }
    #[test]
    fn new_relative() {
        let args = Args {
            db: false,
            edges: PathBuf::new(),
            nodes: PathBuf::new(),
            invert: false,
            methylome: "/home/constantin/methylome/within_gbM_genes".to_string(),
            genome: "/home/constantin/methylome/gbM_gene_anotation_extract_Arabidopsis.bed"
                .to_string(),
            window_size: 5,
            window_step: 1,

            output_dir: "/home/constantin/windows".to_string(),
            absolute: false,
            cutoff: 2048,
        };
        let windows = super::Windows::new(4096, &args);
        assert_eq!(windows.upstream.len(), 100);
        assert_eq!(windows.gene.len(), 100);
        assert_eq!(windows.downstream.len(), 100);
    }
}
