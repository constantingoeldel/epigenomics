use alphabeta::*;
use clap::Parser;

#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    #[arg(short, long, default_value_t = 1000)]
    pub iterations: u32,

    /// Provide an edgefile
    #[arg(long, short)]
    pub edgelist: std::path::PathBuf,

    /// Provide a nodefile - paths will be updated to match the output directory
    #[arg(long, short)]
    pub nodelist: std::path::PathBuf,

    #[arg(long, short, default_value_t = 0.99)]
    pub posterior_max_filter: f64,
}

fn main() {
    let args = Args::parse();

    let (pedigree, p0uu) =
        Pedigree::build(&args.nodelist, &args.edgelist, args.posterior_max_filter)
            .expect("Error while building pedigree: ");
    pedigree.to_file("./data/pedigree_wt.txt");

    let model = ABneutral::run(&pedigree, p0uu, p0uu, 1.0, args.iterations).expect("Model failed");
    let result = BootModel::run(&pedigree, &model, p0uu, p0uu, 1.0, args.iterations)
        .expect("Bootstrap failed");

    model.to_file("./data/model_wt.txt", &result);

    println!("##########");
    println!("Results");
    println!("{model}");
    println!("##########");
    println!("Bootstrap");
    println!("{result}");
}
