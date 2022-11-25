use std::io;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Argument error")]
    ArgumentError(#[from] clap::Error),

    #[error("Could not find the specified {0}! Does it exist? \nPath: {1}")]
    FileError(String, String),

    #[error("File error")]
    FileSystemError(#[from] io::Error),
}
