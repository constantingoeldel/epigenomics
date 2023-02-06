use std::sync::{
    atomic::{AtomicU32, Ordering},
    Mutex,
};

use crate::*;
use argmin::{core::Executor, solver::neldermead::NelderMead};
use rayon::prelude::*;

pub fn run(
    pedigree: &Pedigree,
    p0uu: f64,
    eqp: f64,
    eqp_weight: f64,
    n_starts: u32,
) -> Result<Model, Box<dyn std::error::Error>> {
    let p0mm = 1.0 - p0uu;
    let p0um = 0.0;
    let max_divergence = *pedigree
        .column(3)
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    assert_eq!(p0mm + p0uu + p0um, 1.0);

    let results: Mutex<Vec<Model>> = Mutex::new(Vec::new());
    let counter = AtomicU32::new(0);
    // Optimization loop
    (0..n_starts).into_par_iter().for_each(|_| {
        // Draw random starting values

        let problem = Problem {
            pedigree: pedigree.clone(),
            eqp_weight,
            eqp,
            p_mm: p0mm,
            p_um: p0um,
            p_uu: p0uu,
        };
        // Run Nelder-Mead optimization
        let nm = NelderMead::new(vec![
            Model::new(max_divergence).to_vec(),
            Model::new(max_divergence).to_vec(),
            Model::new(max_divergence).to_vec(),
            Model::new(max_divergence).to_vec(),
            Model::new(max_divergence).to_vec(),
        ]);

        let res = Executor::new(problem, nm)
            .configure(|state| {
                state
                    // .param(vec![alpha, beta, weight, intercept])
                    .max_iters(1000)
            })
            .run()
            .expect("Failed to run Nelder-Mead optimization");

        let m = Model::from_vec(&res.state.best_param.unwrap());

        // let predicted_mm = (m.alpha * ((1.0 - m.alpha).powi(2) - (1.0 - m.beta).powi(2) - 1.0))
        //     / ((m.alpha + m.beta) * ((m.alpha + m.beta - 1.0).powi(2) - 2.0));
        // let predicted_um = (4.0 * m.alpha * m.beta * (m.alpha + m.beta - 2.0))
        //     / ((m.alpha + m.beta) * ((m.alpha + m.beta - 1.0).powi(2) - 2.0));
        // let predicted_uu = (m.beta * ((1.0 - m.beta).powi(2) - (1.0 - m.alpha).powi(2) - 1.0))
        //     / ((m.alpha + m.beta) * ((m.alpha + m.beta - 1.0).powi(2) - 2.0));
        // dbg!(&m);
        results.lock().unwrap().push(m);
        let c = counter.fetch_add(1, Ordering::SeqCst);
        println!("Progress: {}%", ((c * 100) as f32 / (n_starts) as f32));
    });
    let mut results = results.into_inner().unwrap();
    // Calculating the least squares error for all results and selecting the best one
    results.sort_by(|a, b| {
        let pedigree = pedigree.clone();
        let divergence_a = divergence(&pedigree, p0mm, p0um, p0uu, a.alpha, a.beta, a.weight);
        let divergence_b = divergence(&pedigree, p0mm, p0um, p0uu, b.alpha, b.beta, b.weight);

        let lse_a = divergence_a
            .dt1t2
            .iter()
            .zip(pedigree.column(3))
            .map(|(div, ped)| (ped - a.intercept - div).powi(2))
            .sum::<f64>();
        let lse_b = divergence_b
            .dt1t2
            .iter()
            .zip(pedigree.column(3))
            .map(|(div, ped)| (ped - b.intercept - div).powi(2))
            .sum::<f64>();
        lse_a.partial_cmp(&lse_b).unwrap() // Sort ascending
    });

    // Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
    // Caution: Calculating predicted divergence based on lowest LSQ model: check the biology!", "\n")

    let best: &Model = &results[0];
    println!(
        "Done with optimizing: Best model {}, {}, {}, {}",
        best.alpha, best.beta, best.weight, best.intercept
    );
    let divergence = divergence(
        pedigree,
        p0mm,
        p0um,
        p0uu,
        best.alpha,
        best.beta,
        best.weight,
    );

    let mut delta_t = Vec::new();
    for row in pedigree.rows() {
        delta_t.push(row[2] - row[3] - 2.0 * row[1])
    }

    let mut predicted_divergence = Vec::new();

    for (i, _row) in pedigree.rows().into_iter().enumerate() {
        predicted_divergence.push(
            best.intercept + divergence.dt1t2[i], /*+ delta_t[i] * best.alpha*/
        );
    }

    let mut residual = Vec::new();

    for (i, row) in pedigree.rows().into_iter().enumerate() {
        residual.push(row[3] - predicted_divergence[i]);
    }

    // Generating theoretical fit

    // Not needed for now

    Ok(best.to_owned())
}

#[cfg(test)]
mod tests {
    // use super::*;
    // Takes quite a while to run
    //     #[test]
    //     fn it_calculates_model() {
    //         let pedigree = Pedigree::from_file("./data/pedigree.txt");

    //         let result = ABneutral::run(&pedigree, 0.75, 0.5, 0.7, 1).expect("Model failed");

    //         let r = Model::default();

    //         assert_close!(result.alpha, r.alpha);
    //         assert_close!(result.beta, r.beta);
    //         assert_close!(result.weight, r.weight);
    //         assert_close!(result.intercept, r.intercept);
    //     }
}
