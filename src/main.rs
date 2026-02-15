// Attractor search: find novel strange attractors
// in 3D quadratic polynomial ODE systems.

mod catalog;
mod classify;
mod lyapunov;
mod ode;
mod search;

use std::path::PathBuf;

fn usage() {
    eprintln!(
        "Usage:\n  \
         attractor-search random --count N\n  \
         attractor-search evolve \
         --generations G --pop P"
    );
}

fn main() {
    let args: Vec<String> =
        std::env::args().collect();
    if args.len() < 2 {
        usage();
        std::process::exit(1);
    }

    let results = PathBuf::from("results");

    match args[1].as_str() {
        "random" => {
            let count = parse_flag(
                &args, "--count",
            )
            .unwrap_or(10000);
            search::random_search(count, &results);
        }
        "evolve" => {
            let gens = parse_flag(
                &args, "--generations",
            )
            .unwrap_or(500);
            let pop = parse_flag(
                &args, "--pop",
            )
            .unwrap_or(200);
            search::evolve_search(
                gens, pop, &results,
            );
        }
        _ => {
            usage();
            std::process::exit(1);
        }
    }
}

fn parse_flag(
    args: &[String],
    flag: &str,
) -> Option<usize> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
}
