//! Tests of the boundary layer calculation against exact and published
//! solutions. Sources: REFERENCES.md, section 1.

use solver::boundary_layer::*;

fn linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| a + (b - a) * i as f64 / (n - 1) as f64)
        .collect()
}

#[test]
fn laminar_flat_plate_is_near_the_blasius_solution() {
    // Blasius (exact): theta = 0.664 x / sqrt(Re_x), cf = 0.664 / sqrt(Re_x),
    // H = 1.72 / 0.664 = 2.59.
    let re = 1.0e5;
    let s = linspace(0.0, 1.0, 201);
    let ue = vec![1.0; s.len()];
    let result = march(&s, &ue, re, None).unwrap();
    assert!(result.transition.is_none());
    let last = result.stations.last().unwrap();
    let (theta, cf, h) = (0.664 / re.sqrt(), 0.664 / re.sqrt(), 1.72 / 0.664);
    println!(
        "theta {:.3e} ({:.3e}), cf {:.3e} ({:.3e}), H {:.3} ({:.3})",
        last.theta, theta, last.cf, cf, last.h, h
    );
    assert!((last.theta - theta).abs() / theta < 0.015);
    assert!((last.cf - cf).abs() / cf < 0.015);
    assert!((last.h - h).abs() / h < 0.015);
}

#[test]
fn stagnation_flow_has_constant_thickness() {
    // For ue = k s, Thwaites gives theta^2 = 0.075 / (re k) at each station.
    // The exact (Hiemenz) value is theta = 0.292 / sqrt(re k), 6% larger.
    let (re, k) = (1.0e5, 10.0);
    let s = linspace(0.0, 0.05, 51);
    let ue: Vec<f64> = s.iter().map(|x| k * x).collect();
    let result = march(&s, &ue, re, None).unwrap();
    let thwaites = (0.075 / (re * k)).sqrt();
    for station in &result.stations {
        assert!((station.theta - thwaites).abs() / thwaites < 1e-9);
    }
    let hiemenz = 0.292 / (re * k).sqrt();
    assert!((thwaites - hiemenz).abs() / hiemenz < 0.07);
}

#[test]
fn howarth_flow_matches_the_closed_form_and_separates_at_the_published_position() {
    // Howarth flow: ue = 1 - x. Thwaites closed form:
    // theta^2 = (0.075 / re) * ((1 - x)^-6 - 1). Separation (lambda = -0.09)
    // at x = 0.123. The exact solution separates at x = 0.1199.
    let re = 1.0e4; // low value: no natural transition before separation
    let s = linspace(0.0, 0.2, 2001);
    let ue: Vec<f64> = s.iter().map(|x| 1.0 - x).collect();
    let result = march(&s, &ue, re, None).unwrap();
    let (index, cause) = result.transition.unwrap();
    assert_eq!(cause, TransitionCause::LaminarSeparation);
    println!("laminar separation at x = {:.4}", s[index]);
    assert!((s[index] - 0.123).abs() < 0.001);
    for station in &result.stations[1..index] {
        let exact = (0.075 / re * ((1.0 - station.s).powi(-6) - 1.0)).sqrt();
        assert!((station.theta - exact).abs() / exact < 1e-6);
    }
}

#[test]
fn turbulent_flat_plate_is_near_the_prandtl_schlichting_formula() {
    // Friction drag coefficient of one side of a smooth plate, turbulent
    // from the leading edge: CF = 0.455 / (log10 Re)^2.58.
    // Measured difference of Head's method: 7% at Re 1e6, 6% at Re 1e7,
    // 12% at Re 1e8. The program limits the Reynolds number to 1e7.
    for (re, tolerance) in [(1.0e6_f64, 0.10), (1.0e7, 0.10), (1.0e8, 0.15)] {
        let s = linspace(0.0, 1.0, 801);
        let ue = vec![1.0; s.len()];
        let result = march(&s, &ue, re, Some(0.0)).unwrap();
        assert!(result.separation.is_none());
        let last = result.stations.last().unwrap();
        let cf_total = 2.0 * last.theta;
        let formula = 0.455 / re.log10().powf(2.58);
        println!(
            "Re {re:.0e}: CF {cf_total:.5}, formula {formula:.5}, H {:.3}",
            last.h
        );
        assert!((cf_total - formula).abs() / formula < tolerance);
        assert!(last.h > 1.2 && last.h < 1.5);
    }
}

#[test]
fn natural_transition_obeys_the_michel_criterion() {
    let re = 1.0e7;
    let s = linspace(0.0, 1.0, 2001);
    let ue = vec![1.0; s.len()];
    let result = march(&s, &ue, re, None).unwrap();
    let (index, cause) = result.transition.unwrap();
    assert_eq!(cause, TransitionCause::Michel);
    let before = &result.stations[index - 1];
    assert!(re * before.theta <= michel_limit(re * before.s));
    let at = &result.stations[index];
    assert!(re * at.theta > michel_limit(re * at.s));
    assert_eq!(at.regime, Regime::Turbulent);
    assert_eq!(before.regime, Regime::Laminar);
}

#[test]
fn forced_transition_is_at_the_set_position() {
    let s = linspace(0.0, 1.0, 101);
    let ue = vec![1.0; s.len()];
    let result = march(&s, &ue, 1.0e5, Some(0.5)).unwrap();
    let (index, cause) = result.transition.unwrap();
    assert_eq!(cause, TransitionCause::Forced);
    assert!((s[index] - 0.5).abs() < 1e-9);
}

#[test]
fn correlations_are_continuous_and_invertible() {
    // The two branches of each Thwaites correlation agree at lambda = 0.
    assert!((thwaites_h(1e-9) - thwaites_h(-1e-9)).abs() < 1e-3);
    assert!((thwaites_l(1e-9) - thwaites_l(-1e-9)).abs() < 1e-6);
    // The slope of l is also continuous at lambda = 0.
    let (up, down) = (
        (thwaites_l(1e-4) - thwaites_l(0.0)) / 1e-4,
        (thwaites_l(0.0) - thwaites_l(-1e-4)) / 1e-4,
    );
    assert!((up - down).abs() < 0.01, "{up} {down}");
    // The inverse of Head's H1 relation. The published inverse has rounded
    // constants, thus the tolerance is 0.5%.
    for h in [1.2, 1.4, 1.59, 1.7, 2.0, 2.4] {
        assert!((head_h(head_h1(h)) - h).abs() / h < 0.005, "H {h}");
    }
}

#[test]
fn adverse_pressure_causes_turbulent_separation() {
    // Strong linear deceleration of a turbulent boundary layer.
    let s = linspace(0.0, 1.0, 401);
    let ue: Vec<f64> = s.iter().map(|x| 1.0 - 0.8 * x).collect();
    let result = march(&s, &ue, 1.0e6, Some(0.0)).unwrap();
    let index = result.separation.expect("separation");
    // The result list stops at the separation station.
    assert_eq!(result.stations.len(), index);
    // The shape factor increases toward separation.
    let last = result.stations.last().unwrap();
    assert!(
        last.h > 2.0 && last.h < TURBULENT_SEPARATION_H,
        "H {}",
        last.h
    );
    assert!(
        s[index] > 0.3 && s[index] < 0.8,
        "separation at {}",
        s[index]
    );
}

#[test]
fn bad_input_gives_an_error() {
    let s = linspace(0.0, 1.0, 11);
    let ue = vec![1.0; 11];
    assert!(march(&s, &ue, 0.0, None).is_err());
    assert!(march(&s, &ue, f64::NAN, None).is_err());
    assert!(march(&s[..2], &ue[..2], 1e5, None).is_err());
    assert!(march(&s, &ue[..10], 1e5, None).is_err());
    let mut reversed = s.clone();
    reversed.reverse();
    assert!(march(&reversed, &ue, 1e5, None).is_err());
    let mut negative = ue.clone();
    negative[5] = -1.0;
    assert!(march(&s, &negative, 1e5, None).is_err());
}
