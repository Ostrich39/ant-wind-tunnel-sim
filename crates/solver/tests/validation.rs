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
        println!("alpha {alpha}: cl {cl:.5}, exact {exact:.5}, error {:.3}%", err * 100.0);
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
        println!("alpha {alpha}: cl {cl:.5}, exact {exact:.5}, error {:.3}%", err * 100.0);
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
    println!("errors: {:.3}% {:.3}% {:.3}%", e1 * 100.0, e2 * 100.0, e3 * 100.0);
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
        println!("NACA {code}: cl {:.4}, cd {:.5}, cm {:.4}", s.cl, s.cd, s.cm);
        assert!(s.cd.abs() < 3e-3, "cd {}", s.cd);
    }
}

#[test]
fn naca0012_lift_slope_is_near_thick_airfoil_estimate() {
    // Estimate for thick airfoils: slope = 2 * pi * (1 + 0.77 * t) per radian.
    let estimate = 2.0 * PI * (1.0 + 0.77 * 0.12) * 4.0_f64.to_radians();
    let cl = solve(&naca4("0012", 160).unwrap(), 4.0).unwrap().cl;
    println!("cl {cl:.4}, estimate {estimate:.4}");
    assert!((cl - estimate).abs() / estimate < 0.03);
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
