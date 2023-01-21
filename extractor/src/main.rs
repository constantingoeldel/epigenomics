use clap::Parser;
use lib::{arguments::Args, *};

fn main() {
    let args = Args::parse();
    match extract(args) {
        Ok(_) => println!("Done!"),
        Err(e) => println!("Error: {}", e),
    }
}
