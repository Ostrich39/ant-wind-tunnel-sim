//! Validation tests. Each test compares the solver with a value from theory.

use solver::{joukowski, naca4, solve};
use std::f64::consts::PI;

/// Exact lift coefficient of a Joukowski airfoil:
/// cl = 8 * pi * (a / chord) * sin(alpha + beta).
fn joukowski_exact(a: f64, beta: f64, chord: f64, alpha_deg: f64) -> f64 {
    8.0 * PI * a / chord * (alpha_deg.to_radians() + beta).sin()
}

#[test]
fn joukowski_symmetric_matches_exact_lift() {
    let (foil, a, beta) = joukowski(0.1, 0.0, 200);
    for alpha in [2.0, 5.0, 8.0] {
        let exact = joukowski_exact(a, beta, foil.chord(), alpha);
        let cl = solve(&foil, alpha).unwrap().cl;
        let err = (cl - exact).abs() / exact;
        println!(
            "alpha {alpha}: cl {cl:.5}, exact {exact:.5}, error {:.3}%",
            err * 100.0
        );
        assert!(err < 0.01, "error {err}");
    }
}

#[test]
fn joukowski_cambered_matches_exact_lift() {
    // The sharp (cusp) trailing edge of this shape makes convergence slow,
    // thus this test uses more panels and a larger tolerance.
    let (foil, a, beta) = joukowski(0.1, 0.1, 400);
    for alpha in [0.0, 4.0] {
        let exact = joukowski_exact(a, beta, foil.chord(), alpha);
        let cl = solve(&foil, alpha).unwrap().cl;
        let err = (cl - exact).abs() / exact;
        println!(
            "alpha {alpha}: cl {cl:.5}, exact {exact:.5}, error {:.3}%",
            err * 100.0
        );
        assert!(err < 0.02, "error {err}");
    }
}

#[test]
fn error_decreases_when_panel_count_increases() {
    let error = |n: usize| {
        let (foil, a, beta) = joukowski(0.1, 0.1, n);
        let exact = joukowski_exact(a, beta, foil.chord(), 4.0);
        (solve(&foil, 4.0).unwrap().cl - exact).abs() / exact
    };
    let (e1, e2, e3) = (error(100), error(200), error(400));
    println!(
        "errors: {:.3}% {:.3}% {:.3}%",
        e1 * 100.0,
        e2 * 100.0,
        e3 * 100.0
    );
    assert!(e1 > e2 && e2 > e3);
}

#[test]
fn symmetric_airfoil_has_no_lift_at_zero_angle() {
    let s = solve(&naca4("0012", 160).unwrap(), 0.0).unwrap();
    assert!(s.cl.abs() < 1e-6, "cl {}", s.cl);
    assert!(s.cm.abs() < 1e-6, "cm {}", s.cm);
}

#[test]
fn inviscid_flow_has_no_drag() {
    // D'Alembert's paradox: drag is zero in inviscid 2D flow.
    for code in ["0012", "2412", "4415"] {
        let s = solve(&naca4(code, 160).unwrap(), 4.0).unwrap();
        println!(
            "NACA {code}: cl {:.4}, cd {:.5}, cm {:.4}",
            s.cl, s.cd, s.cm
        );
        assert!(s.cd.abs() < 3e-3, "cd {}", s.cd);
    }
}

/// Largest pressure coefficient error against the exact Joukowski solution,
/// for panels in the front 90% of the chord.
fn joukowski_cp_error(eps: f64, delta: f64, n: usize, alpha_deg: f64) -> f64 {
    let (foil, a, beta) = joukowski(eps, delta, n);
    let s = solve(&foil, alpha_deg).unwrap();
    let alpha = alpha_deg.to_radians();
    let (chord, x0) = (foil.chord(), foil.x_min());
    let mut worst: f64 = 0.0;
    for j in 0..n {
        // Point on the circle that maps to the middle of panel j.
        let phi = -beta - 2.0 * PI * (j as f64 + 0.5) / n as f64;
        let (zx, zy) = (-eps + a * phi.cos(), delta + a * phi.sin());
        let r2 = zx * zx + zy * zy;
        // |dz/dzeta| = |1 - 1/zeta^2|
        let (ix, iy) = ((zx * zx - zy * zy) / (r2 * r2), -2.0 * zx * zy / (r2 * r2));
        let stretch = (1.0 - ix).hypot(iy);
        // Exact surface speed divided by the free-stream speed.
        let v = 2.0 * ((phi - alpha).sin() + (alpha + beta).sin()).abs() / stretch;
        let exact_cp = 1.0 - v * v;
        if (s.mid[j][0] - x0) / chord < 0.9 {
            worst = worst.max((exact_cp - s.cp[j]).abs());
        }
    }
    worst
}

#[test]
fn surface_pressure_matches_exact_joukowski_solution() {
    for (eps, delta) in [(0.1, 0.0), (0.1, 0.1)] {
        for alpha in [0.0, 4.0] {
            let coarse = joukowski_cp_error(eps, delta, 200, alpha);
            let fine = joukowski_cp_error(eps, delta, 400, alpha);
            println!("delta {delta} alpha {alpha}: cp error {coarse:.4} -> {fine:.4}");
            assert!(fine < 0.05, "cp error {fine}");
            // The error must decrease when the panel count increases.
            assert!(fine < 0.6 * coarse, "no convergence: {coarse} -> {fine}");
        }
    }
}

#[test]
fn two_lift_calculations_agree() {
    // Lift from circulation and lift from pressure are independent results.
    for code in ["0012", "2412", "4415"] {
        for alpha in [-2.0, 4.0, 8.0] {
            let s = solve(&naca4(code, 200).unwrap(), alpha).unwrap();
            assert!(
                (s.cl - s.cl_pressure).abs() < 0.01,
                "NACA {code} alpha {alpha}: {} and {}",
                s.cl,
                s.cl_pressure
            );
        }
    }
}

#[test]
fn lift_is_antisymmetric_for_a_symmetric_airfoil() {
    let foil = naca4("0012", 160).unwrap();
    let (up, down) = (solve(&foil, 5.0).unwrap(), solve(&foil, -5.0).unwrap());
    assert!((up.cl + down.cl).abs() < 1e-9);
    assert!((up.cm + down.cm).abs() < 1e-9);
}

#[test]
fn naca_geometry_matches_the_defining_equation() {
    // NACA 0012 thickness at 30% chord is 12% (maximum thickness point).
    let foil = naca4("0012", 200).unwrap();
    let max_y = foil.points.iter().map(|p| p[1]).fold(0.0, f64::max);
    assert!((max_y - 0.06).abs() < 2e-4, "half thickness {max_y}");
    // NACA 2412: maximum camber is 2% of chord at 40% chord.
    let foil = naca4("2412", 200).unwrap();
    let n = foil.points.len();
    let camber_max = (0..n / 2)
        .map(|k| 0.5 * (foil.points[k][1] + foil.points[n - 1 - k][1]))
        .fold(0.0, f64::max);
    assert!((camber_max - 0.02).abs() < 2e-4, "camber {camber_max}");
}

#[test]
fn wrong_point_order_gives_an_error() {
    let mut foil = naca4("0012", 160).unwrap();
    foil.points.reverse();
    assert!(solve(&foil, 2.0).is_err());
    let mut open = naca4("0012", 160).unwrap();
    open.points.pop();
    assert!(solve(&open, 2.0).is_err());
    assert!(solve(&naca4("0012", 160).unwrap(), f64::NAN).is_err());
}

#[test]
fn stagnation_point_has_cp_near_one() {
    let s = solve(&naca4("0012", 160).unwrap(), 4.0).unwrap();
    let max = s.cp.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!((max - 1.0).abs() < 0.02, "max cp {max}");
}

#[test]
fn bad_input_gives_an_error() {
    assert!(naca4("12", 160).is_err());
    assert!(naca4("0012", 161).is_err());
    assert!(naca4("2012", 160).is_err());
}
