use std::path::Path;

use clap::Parser;
use lib::{arguments::Args, *};

fn main() {
    let args = Args::parse();
    match extract(args.clone()) {
        Err(e) => println!("Error: {e}"),
        Ok(max_gene_length) => {
            if args.alphabeta {
                let regions = vec![
                    ("gene", max_gene_length),
                    ("upstream", args.cutoff),
                    ("downstream", args.cutoff),
                ];

                for region in regions {
                    let max = if args.absolute { region.1 } else { 100 };
                    let side = region.0;

                    for window in (0..=max).step_by(args.window_step as usize) {
                        let alphabeta_result = alphabeta::run(
                            Path::new(&format!(
                                "{}/{}/{}/nodelist.txt",
                                &args.output_dir, side, window
                            )),
                            Path::new(&format!(
                                "{}/{}/{}/edgelist.txt",
                                &args.output_dir, side, window
                            )),
                            0.99,
                            1000,
                            format!("{}/{}/{}/", &args.output_dir, side, window),
                        );
                        match alphabeta_result {
                            Err(e) => println!("Error: {e}"),
                            Ok(_) => println!("Finished {side} {window}"),
                        }
                    }
                }
            }
        }
    }
}
