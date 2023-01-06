use clap::Parser;
use extractor::{arguments::Args, extract};

fn main() {
    let args = Args::parse();
    match extract(args) {
        Ok(_) => println!("Done!"),
        Err(e) => println!("Error: {}", e),
    }
}
