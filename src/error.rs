use std::io;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Argument error")]
    Argument(#[from] clap::Error),

    #[error("Could not find the specified {0}! Does it exist? \nPath: {1}")]
    File(String, String),

    #[error("File error {0}")]
    FileSystem(#[from] io::Error),

    #[error("Unable to extract CG site from line")]
    CGSite,

    #[error("Unable to convert: Are you passing a valid number? {0}")]
    NumberConversion(#[from] std::num::ParseIntError),
}
