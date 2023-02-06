use std::{fmt::Write, path::Path, time::Duration};

use clap::Parser;
use db::import_results;
use indicatif::{HumanDuration, MultiProgress, ProgressBar, ProgressState, ProgressStyle};
use lib::{arguments::Args, structs::Region, *};
use tokio::time::Interval;

#[tokio::main]
async fn main() {
    let args = Args::parse();
    println!("Starting run {}", args.name);
    match extract(args.clone()) {
        Err(e) => println!("Error: {e}"),
        Ok(max_gene_length) => {
            if args.alphabeta {
                let regions = vec![
                    (Region::Gene, max_gene_length),
                    (Region::Upstream, args.cutoff),
                    (Region::Downstream, args.cutoff),
                ];
                let mut results = Vec::new();
                let total_steps = if args.absolute {
                    (max_gene_length + 2 * args.cutoff) / args.window_step
                } else {
                    (3 * 100) / args.window_step
                };
                let multi = MultiProgress::new();
                let pb = multi.add(ProgressBar::new(total_steps as u64));
                pb.set_message("Progress ");
                pb.enable_steady_tick(Duration::new(1, 0));
                pb.set_style(
                    ProgressStyle::with_template(
                        "{msg} {bar:40.magenta/blue} [{elapsed}] {pos:>7}/{len:7} ETA: {eta}",
                    )
                    .unwrap()
                    .with_key(
                        "eta",
                        |state: &ProgressState, w: &mut dyn Write| {
                            write!(
                                w,
                                "{}",
                                HumanDuration(Duration::from_secs(state.eta().as_secs()))
                            )
                            .unwrap();
                        },
                    ),
                );

                for region in regions {
                    let max = if args.absolute { region.1 } else { 100 };

                    for window in (0..max).step_by(args.window_step as usize) {
                        pb.inc(1);
                        let alphabeta_result = alphabeta::run(
                            Path::new(&format!(
                                "{}/{}/{}/nodelist.txt",
                                &args.output_dir, &region.0, window
                            )),
                            Path::new(&format!(
                                "{}/{}/{}/edgelist.txt",
                                &args.output_dir, &region.0, window
                            )),
                            0.99,
                            1000,
                            format!("{}/{}/{}/", &args.output_dir, region.0, window),
                            &multi,
                        );
                        match alphabeta_result {
                            Err(e) => println!("Error: {e}"),
                            Ok((model, errors)) => results.push((model, errors, region.0.clone())),
                        }
                    }
                }
                pb.finish();
                let db = db::connect()
                    .await
                    .expect("Could not connect to database: Did you provide a connection string?");
                import_results(&db, args.name, results).await.expect("Could not save results to a database. Your data is stored in files in each directory");
            }
        }
    }
}
