// Search strategies: random + evolutionary.

use crate::catalog::{self, Entry};
use crate::classify::{self, Dynamics};
use crate::ode::{self, Coeffs, NCOEFFS};
use rand::Rng;
use std::path::Path;

const COEFF_RANGE: f64 = 2.0;
const TRAJ_SAMPLE: usize = 50000;
const TRAJ_DT: f64 = 0.01;

/// Random coefficient vector in [-range, range].
fn rand_coeffs(rng: &mut impl Rng) -> Coeffs {
    let mut c = [0.0; NCOEFFS];
    for v in c.iter_mut() {
        *v = rng.gen_range(
            -COEFF_RANGE..COEFF_RANGE,
        );
    }
    c
}

/// Report a chaotic find to stdout and catalog.
fn report(
    c: &Coeffs,
    data: &crate::lyapunov::LyapData,
    results_dir: &Path,
    found: usize,
    method: &str,
) {
    let hash = catalog::coeff_hash(c);
    println!(
        "  CHAOTIC #{}: hash={:016x} \
         λ1={:.4} dim={:.3}",
        found, hash, data.spectrum[0], data.ky_dim,
    );
    // Sample trajectory for catalog.
    let s0 = [0.1, 0.1, 0.1];
    let traj = ode::integrate_traj(
        c, &s0, TRAJ_DT, TRAJ_SAMPLE,
    );
    let entry = Entry::new(
        c,
        data.spectrum,
        data.ky_dim,
        &traj,
        method,
    );
    if let Err(e) = catalog::save(results_dir, &entry)
    {
        eprintln!("  save error: {e}");
    }
}

/// Random search: evaluate `count` random systems.
pub fn random_search(
    count: usize,
    results_dir: &Path,
) {
    let mut rng = rand::thread_rng();
    let mut found = 0usize;
    let mut evaluated = 0usize;
    println!("Random search: {count} systems");
    for i in 0..count {
        let c = rand_coeffs(&mut rng);
        evaluated += 1;
        match classify::classify(&c) {
            Dynamics::Chaotic(data) => {
                found += 1;
                report(
                    &c,
                    &data,
                    results_dir,
                    found,
                    "random",
                );
            }
            _ => {}
        }
        if (i + 1) % 1000 == 0 {
            println!(
                "  [{}/{}] chaotic: {} ({:.2}%)",
                i + 1,
                count,
                found,
                100.0 * found as f64
                    / evaluated as f64,
            );
        }
    }
    println!(
        "Done. {evaluated} evaluated, \
         {found} chaotic ({:.2}%)",
        100.0 * found as f64 / evaluated as f64,
    );
}

/// Evolutionary search.
pub fn evolve_search(
    generations: usize,
    pop_size: usize,
    results_dir: &Path,
) {
    let mut rng = rand::thread_rng();
    let mut pop: Vec<(Coeffs, f64)> = (0..pop_size)
        .map(|_| (rand_coeffs(&mut rng), f64::NEG_INFINITY))
        .collect();
    let mut found = 0usize;
    let mut total_eval = 0usize;
    let mutate_std = 0.3;

    println!(
        "Evolve: {generations} gens, pop {pop_size}"
    );

    for gen in 0..generations {
        // Evaluate fitness for new individuals.
        for item in pop.iter_mut() {
            if item.1 == f64::NEG_INFINITY {
                total_eval += 1;
                item.1 = fitness(&item.0);
                if let Dynamics::Chaotic(data) =
                    classify::classify(&item.0)
                {
                    found += 1;
                    report(
                        &item.0,
                        &data,
                        results_dir,
                        found,
                        "evolve",
                    );
                }
            }
        }

        // Sort by fitness descending.
        pop.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        if (gen + 1) % 50 == 0 {
            let best = pop[0].1;
            println!(
                "  gen {}: best_fit={:.4} \
                 chaotic={found} eval={total_eval}",
                gen + 1,
                best,
            );
        }

        // Keep top half, mutate to fill rest.
        let half = pop_size / 2;
        for i in half..pop_size {
            let parent = i % half;
            let mut child = pop[parent].0;
            mutate(
                &mut child,
                mutate_std,
                &mut rng,
            );
            pop[i] = (child, f64::NEG_INFINITY);
        }
    }
    println!(
        "Done. {total_eval} evaluated, \
         {found} chaotic",
    );
}

/// Fitness: higher = more interesting.
/// Positive Lyapunov + bounded = best.
fn fitness(c: &Coeffs) -> f64 {
    let s0 = [0.1f64, 0.1, 0.1];
    let (lyap, ok) = crate::lyapunov::max_lyapunov(
        c, &s0, 0.005, 1000, 5000,
    );
    if !ok || !lyap.is_finite() {
        return -100.0;
    }
    // Reward positive exponent, penalize huge.
    if lyap > 0.0 && lyap < 10.0 {
        lyap
    } else if lyap <= 0.0 {
        lyap - 10.0
    } else {
        -50.0
    }
}

/// Mutate coefficients with Gaussian noise.
fn mutate(
    c: &mut Coeffs,
    std: f64,
    rng: &mut impl Rng,
) {
    for v in c.iter_mut() {
        // Box-Muller for normal dist w/o dep.
        let u1: f64 = rng.gen_range(1e-10..1.0);
        let u2: f64 = rng.gen_range(0.0..
            std::f64::consts::TAU);
        let z = (-2.0 * u1.ln()).sqrt() * u2.cos();
        *v += std * z;
        *v = v.clamp(-COEFF_RANGE, COEFF_RANGE);
    }
}
