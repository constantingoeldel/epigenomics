pub mod ABneutral;
pub mod BootModel;
pub mod divergence;
pub mod macros;
pub mod pedigree;
pub mod structs;
extern crate blas_src;

use std::path::Path;

use anyhow::{anyhow, Error};
pub use divergence::*;
pub use macros::*;
pub use pedigree::*;
pub use structs::*;
pub use ABneutral::*;
pub use BootModel::*;

pub fn run(
    nodelist: &Path,
    edgelist: &Path,
    posterior_max_filter: f64,
    iterations: u32,
    output: String,
) -> Result<(Model, StandardDeviations), Error> {
    let (pedigree, p0uu) = Pedigree::build(nodelist, edgelist, posterior_max_filter)
        .map_err(|e| anyhow!("Error while building pedigree: {}", e))?;
    pedigree.to_file(&format!("{output}/pedigree.txt"))?;

    let model = ABneutral::run(&pedigree, p0uu, p0uu, 1.0, iterations)
        .map_err(|e| anyhow!("Model failed: {}", e))?;
    let result = BootModel::run(&pedigree, &model, p0uu, p0uu, 1.0, iterations)
        .map_err(|e| anyhow!("Bootstrap failed: {}", e))?;

    println!("##########");
    println!("Results:\n");
    println!("{model}");
    println!("{result}");
    println!("##########");
    model.to_file(&format!("{output}/pedigree.txt"), &result)?;
    Ok((model, result))
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{assert_close, assert_close_10_percent, ABneutral, BootModel, Pedigree};

    #[test] // Recommended to run with --release
    fn end_to_end() {
        let (pedigree, p0uu) = Pedigree::build(
            Path::new("./data/desired_output/nodelist.txt"),
            Path::new("./data/desired_output/edgelist.txt"),
            0.99,
        )
        .expect("Error while building pedigree: ");

        let model = ABneutral::run(&pedigree, p0uu, p0uu, 1.0, 1000).expect("Model failed");
        let result =
            BootModel::run(&pedigree, &model, p0uu, p0uu, 1.0, 200).expect("Bootstrap failed");
        println!("{result}");
        assert_close_10_percent!(model.alpha, 5.7985750419976e-05);
        assert_close_10_percent!(model.beta, 0.00655710970515347);
        assert_close_10_percent!(p0uu, 0.991008120326199);
        assert_close!(model.alpha, 5.7985750419976e-05);
        assert_close!(model.beta, 0.00655710970515347);
        assert_close!(p0uu, 0.991008120326199);

        // assert_close_10_percent!(result.alpha, 1.19025049000535e-05);
        // assert_close_10_percent!(result.beta, 0.00138056339729951);
        // assert_close_10_percent!(result.alpha_beta, 0.594234575518036);
    }
}
