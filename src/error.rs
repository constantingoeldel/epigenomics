use std::io;

use thiserror::Error;

// use crate::methylation_site::MethylationSite;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Argument error")]
    Argument(#[from] clap::Error),

    #[error("Could not find the specified {0}! Does it exist? \nPath: {1}")]
    File(String, String),

    #[error("File error {0}")]
    FileSystem(#[from] io::Error),
    //#[error("No fitting gene found for cg site {0}!")]
    //NoCorrespondingGeneFound(MethylationSite),
}
