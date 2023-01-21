use std::io;
#[derive(Debug, thiserror::Error)]
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

    #[error("Unable to convert: Are you passing a valid number? {0}")]
    FloatConversion(#[from] std::num::ParseFloatError),

    #[error("Error when querying the database: {0}")]
    DB(#[from] sqlx::Error),

    #[error("Error: {0}")]
    Simple(&'static str),
}
