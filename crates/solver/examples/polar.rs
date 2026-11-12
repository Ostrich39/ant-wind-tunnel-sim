//! Prints a table of lift, drag, and moment for a NACA 4-digit airfoil.
//! Usage: cargo run --release --example polar -- 2412

use solver::{naca4, solve};

fn main() {
    let code = std::env::args().nth(1).unwrap_or_else(|| "0012".into());
    let foil = match naca4(&code, 160) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };
    println!("{} (inviscid, {} panels)", foil.name, foil.points.len() - 1);
    println!("{:>7} {:>9} {:>9} {:>9}", "alpha", "cl", "cd", "cm");
    for a in (-4..=10).step_by(2) {
        let t = std::time::Instant::now();
        let s = solve(&foil, a as f64).expect("solve failed");
        println!(
            "{:>7} {:>9.4} {:>9.5} {:>9.4}   ({:.2} ms)",
            a, s.cl, s.cd, s.cm,
            t.elapsed().as_secs_f64() * 1000.0
        );
    }
}
