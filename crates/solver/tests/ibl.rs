//! Tests of the two-equation boundary layer method (Drela closure) against
//! exact and published solutions. Sources: REFERENCES.md, section 1.

use solver::closure::{self, Flow};
use solver::ibl::*;

fn linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| a + (b - a) * i as f64 / (n - 1) as f64)
        .collect()
}

#[test]
fn laminar_flat_plate_is_near_the_blasius_solution() {
    let re = 1.0e5;
    let xi = linspace(0.005, 1.0, 200);
    let ue = vec![1.0; xi.len()];
    let result = march(&xi, &ue, re, false, None, closure::N_CRIT).unwrap();
    assert!(result.transition.is_none() && result.separation.is_none());
    let last = result.states.last().unwrap();
    let d = last.derived(re);
    let (theta, cf, h) = (0.664 / re.sqrt(), 0.664 / re.sqrt(), 1.72 / 0.664);
    println!(
        "theta {:.4e} ({theta:.4e}), cf {:.4e} ({cf:.4e}), H {:.3} ({h:.3})",
        last.theta, d.cf, d.h
    );
    assert!((last.theta - theta).abs() / theta < 0.01);
    assert!((d.cf - cf).abs() / cf < 0.01);
    assert!((d.h - h).abs() / h < 0.01);
}

#[test]
fn stagnation_flow_is_near_the_hiemenz_solution() {
    // ue = k xi. Exact: theta = 0.292 / sqrt(re k), H = 2.216.
    let (re, k) = (1.0e5, 10.0);
    let xi = linspace(0.001, 0.05, 50);
    let ue: Vec<f64> = xi.iter().map(|x| k * x).collect();
    let result = march(&xi, &ue, re, true, None, closure::N_CRIT).unwrap();
    let exact = 0.292 / (re * k).sqrt();
    for state in [result.states[0], *result.states.last().unwrap()] {
        println!(
            "theta {:.4e} ({exact:.4e}), H {:.3}",
            state.theta,
            state.h()
        );
        assert!((state.theta - exact).abs() / exact < 0.02);
        assert!((state.h() - 2.216).abs() < 0.03);
    }
}

#[test]
fn howarth_flow_separates_near_the_exact_position() {
    // ue = 1 - x. The exact solution separates at x = 0.1199.
    let re = 1.0e4;
    let xi = linspace(0.0005, 0.2, 800);
    let ue: Vec<f64> = xi.iter().map(|x| 1.0 - x).collect();
    let result = march(&xi, &ue, re, false, None, closure::N_CRIT).unwrap();
    let (index, cause) = result.transition.expect("laminar separation");
    assert_eq!(cause, TransitionCause::LaminarSeparation);
    println!("H reaches {H_LIMIT_LAMINAR} at x = {:.4}", xi[index]);
    assert!((xi[index] - 0.1199).abs() < 0.006);
}

#[test]
fn turbulent_flat_plate_is_near_the_prandtl_schlichting_formula() {
    // Measured difference at Re 1e7 and 1e8: the method is 5% to 7% below
    // the formula.
    for re in [1.0e6_f64, 1.0e7, 1.0e8] {
        let xi = linspace(0.001, 1.0, 400);
        let ue = vec![1.0; xi.len()];
        let result = march(&xi, &ue, re, false, Some(0.0), closure::N_CRIT).unwrap();
        assert!(result.separation.is_none());
        // The trip is effective at the first station where the boundary
        // layer is thick enough for turbulent flow.
        let (index, cause) = result.transition.unwrap();
        assert_eq!(cause, TransitionCause::Forced);
        println!("Re {re:.0e}: trip effective at x = {:.4}", xi[index]);
        assert!(re * result.states[index].theta >= TRIP_MIN_RE_THETA);
        assert!(re * result.states[index - 1].theta < TRIP_MIN_RE_THETA);
        let last = result.states.last().unwrap();
        let cf_total = 2.0 * last.theta;
        let formula = 0.455 / re.log10().powf(2.58);
        println!(
            "Re {re:.0e}: CF {cf_total:.5}, formula {formula:.5}, H {:.3}",
            last.h()
        );
        // At Re 1e6 the laminar start covers 9% of the plate, thus the
        // formula for a fully turbulent plate does not apply.
        if re >= 1.0e7 {
            assert!((cf_total - formula).abs() / formula < 0.08);
        }
        assert_eq!(last.flow, Flow::Turbulent);
    }
}

#[test]
fn natural_transition_on_a_flat_plate() {
    // Envelope e^n method with n = 9 on a flat plate. This test records
    // the transition Reynolds number of the method.
    let re = 1.0e7;
    let xi = linspace(0.001, 1.0, 1000);
    let ue = vec![1.0; xi.len()];
    let result = march(&xi, &ue, re, false, None, closure::N_CRIT).unwrap();
    let (index, cause) = result.transition.unwrap();
    assert_eq!(cause, TransitionCause::Natural);
    let re_x = re * xi[index];
    // Measured: Re_x = 4.0e6. No reference value is recorded for this.
    println!("transition at Re_x = {re_x:.3e}");
    assert!(re_x > 1.0e6 && re_x < 1.0e7);
    // A lower critical value moves the transition forward.
    let early = march(&xi, &ue, re, false, None, 5.0).unwrap();
    assert!(early.transition.unwrap().0 < index);
}

#[test]
fn closure_values_at_reference_points() {
    // Laminar: the Blasius shape factor gives cf * Re_theta near 0.441.
    let cf = closure::cf_laminar(2.59, 1.0);
    assert!((cf - 0.441).abs() < 0.015, "{cf}");
    // The laminar relations are continuous at the branch points.
    assert!((closure::hs_laminar(4.35 - 1e-9) - closure::hs_laminar(4.35 + 1e-9)).abs() < 1e-6);
    assert!(
        (closure::cf_laminar(5.5 - 1e-9, 1.0) - closure::cf_laminar(5.5 + 1e-9, 1.0)).abs() < 1e-6
    );
    assert!(
        (closure::di_laminar(4.0 - 1e-9, 1.0) - closure::di_laminar(4.0 + 1e-9, 1.0)).abs() < 1e-6
    );
    // Turbulent H* is continuous at H0.
    let h0 = 3.0 + 400.0 / 2000.0;
    assert!(
        (closure::hs_turbulent(h0 - 1e-9, 2000.0) - closure::hs_turbulent(h0 + 1e-9, 2000.0)).abs()
            < 1e-6
    );
    // Turbulent skin friction decreases with the shape factor and is
    // negative in separated flow.
    assert!(closure::cf_turbulent(1.4, 5000.0) > closure::cf_turbulent(2.0, 5000.0));
    assert!(closure::cf_turbulent(4.0, 5000.0) < 0.0);
    // No amplification below the critical Reynolds number.
    assert_eq!(closure::amplification_rate(2.59, 1e-3, 50.0), 0.0);
    assert!(closure::amplification_rate(2.59, 1e-3, 2000.0) > 0.0);
}
