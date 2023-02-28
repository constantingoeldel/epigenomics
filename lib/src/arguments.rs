use std::path::PathBuf;

use clap::Parser;

/// simple tool to separate a methylome by position within a gene
#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Path of directory containing the methlyome files from which to extract the CG-sites
    #[arg(short, long)]
    pub methylome: PathBuf,

    /// Path of the annotation file containing information about beginning and end of gbM-genes
    #[arg(short, long)]
    pub genome: PathBuf,

    /// Size of the window in percent of the gbM-gene length or in basepair number if --absolute is supplied
    #[arg(short, long, default_value_t = 5)]
    pub window_size: u32,

    /// Size of the step between the start of each window. Default value is window-size, so no overlapp happens
    #[arg(long, short('s'), default_value_t = 0)]
    pub window_step: u32,

    /// Path of the directory where extracted segments shall be stored
    #[arg(short, long)]
    pub output_dir: PathBuf,

    /// Use absolute length in base-pairs for window size instead of percentage of gene length
    #[arg(short, long, default_value_t = false)]
    pub absolute: bool,

    /// Number of basepairs to include upstream and downstream of gene
    #[arg(short, long, default_value_t = 2048)]
    pub cutoff: u32,

    /// Invert strands, to switch from 5' to 3' and vice versa
    #[arg(short, long, default_value_t = false)]
    pub invert: bool,

    /// Use a Postgres database to do everything
    #[arg(long, default_value_t = true)]
    pub db: bool,

    /// Provide an edgefile
    #[arg(long, short)]
    pub edges: Option<std::path::PathBuf>,

    /// Provide a nodefile - paths will be updated to match the output directory
    #[arg(long, short)]
    pub nodes: Option<std::path::PathBuf>,

    /// Also run AlphaBeta on every window after extraction, results will be stored in the same directory as the segments
    #[arg(long, default_value_t = false)]
    pub alphabeta: bool,

    /// Name of the run to be used when storing the result in Postgres
    #[arg(long, default_value_t = format!("{:?}", std::time::Instant::now()))]
    pub name: String,

    /// Overwrite existing content in output directory? If false (default) it will reuse existing windows
    #[arg(long, short, default_value_t = false)]
    pub force: bool,
    /// Let the cutoff be the gene length instead of a fixed number.
    /// So if the gene is 1000 bp long, the cutoff will be 1000 bp instead of 2048 bp (the default).
    /// This option takes preference over the cutoff option.    
    #[arg(long, default_value_t = false)]
    pub cutoff_gene_length: bool,
}

impl Default for Args {
    fn default() -> Self {
        Args {
            db: false,
            invert: false,
            absolute: false,
            cutoff: 2048,
            genome: PathBuf::from("default"),
            methylome: PathBuf::from("also default"),
            output_dir: PathBuf::from("also default"),
            window_size: 5,
            window_step: 1,
            edges: None,
            nodes: None,
            alphabeta: false,
            name: String::new(),
            force: false,
            cutoff_gene_length: false,
        }
    }
}
