use std::{fmt::Write, fs::write, path::Path, time::Duration};

use clap::Parser;
use db::import_results;
use indicatif::{HumanDuration, MultiProgress, ProgressBar, ProgressState, ProgressStyle};
use lib::{arguments::Args, structs::Region, *};

#[tokio::main]
async fn main() {
    let args = Args::parse();
    println!("Starting run {}", args.name);
    match extract(args.clone()) {
        Err(e) => println!("Error: {e}"),
        Ok((max_gene_length, distribution)) => {
            if args.alphabeta {
                let regions = vec![
                    (Region::Upstream, args.cutoff),
                    (Region::Gene, max_gene_length),
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
                            10,
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

                let mut print = String::from("run;window;cg_count;region;alpha;beta;alpha_error;beta_error;1/2*(alpha+beta);pred_steady_state\n");

                let steady_state = |alpha: f64, beta: f64| {
                    (alpha * ((1.0 - alpha).powi(2) - (1.0 - beta).powi(2) - 1.0))
                        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0))
                };
                for (i, ((model, sd, region), d)) in
                    results.iter().zip(distribution.iter()).enumerate()
                {
                    print += &format!(
                        "{};{};{};{};{};{};{};{};{};{}\n",
                        args.name,
                        i,
                        d,
                        region,
                        model.alpha,
                        model.beta,
                        sd.alpha,
                        sd.beta,
                        0.5 * (model.alpha + model.beta),
                        steady_state(model.alpha, model.beta)
                    )
                }

                println!("{print}");

                write(args.output_dir + "/results.txt", print)
                    .expect("Could not save results to file.");
                let db = db::connect()
                    .await
                    .expect("Could not connect to database: Did you provide a connection string?");
                import_results(&db, args.name, results).await.expect("Could not save results to a database. Your data is stored in files in each directory");
            }
        }
    }
}
