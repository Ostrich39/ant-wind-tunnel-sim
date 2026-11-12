//! Tests of the airfoil friction drag calculation. No wind tunnel data is
//! in these tests. They check symmetry, trends, and numerical convergence.

use solver::{analyze, naca4, solve};

fn drag(code: &str, panels: usize, alpha: f64, re: f64, trip: Option<f64>) -> Option<f64> {
    let foil = naca4(code, panels).unwrap();
    let solution = solve(&foil, alpha).unwrap();
    analyze(&foil, &solution, re, trip).unwrap().cd
}

#[test]
fn symmetric_airfoil_at_zero_angle_has_equal_surfaces() {
    let foil = naca4("0012", 160).unwrap();
    let solution = solve(&foil, 0.0).unwrap();
    let result = analyze(&foil, &solution, 3.0e6, None).unwrap();
    let (up, low) = (result.upper.cd.unwrap(), result.lower.cd.unwrap());
    assert!((up - low).abs() / up < 1e-6, "{up} {low}");
    let (xu, xl) = (
        result.upper.transition.unwrap().0,
        result.lower.transition.unwrap().0,
    );
    assert!((xu - xl).abs() < 1e-9);
    assert!(
        result.stagnation_x.abs() < 1e-3,
        "stagnation at {}",
        result.stagnation_x
    );
}

#[test]
fn drag_decreases_when_reynolds_number_increases() {
    let values: Vec<f64> = [1.0e6, 3.0e6, 6.0e6, 1.0e7]
        .iter()
        .map(|re| drag("0012", 160, 0.0, *re, None).unwrap())
        .collect();
    println!("{values:?}");
    assert!(values.windows(2).all(|w| w[1] < w[0]));
}

#[test]
fn early_forced_transition_increases_drag() {
    let free = drag("0012", 160, 0.0, 6.0e6, None).unwrap();
    let tripped = drag("0012", 160, 0.0, 6.0e6, Some(0.05)).unwrap();
    println!("free {free:.5}, tripped {tripped:.5}");
    assert!(tripped > free);
}

#[test]
fn transition_moves_forward_on_the_upper_surface_with_angle() {
    let foil = naca4("0012", 160).unwrap();
    let mut previous = f64::INFINITY;
    for alpha in [0.0, 2.0, 4.0, 6.0] {
        let solution = solve(&foil, alpha).unwrap();
        let result = analyze(&foil, &solution, 3.0e6, None).unwrap();
        let x = result.upper.transition.unwrap().0;
        assert!(x < previous, "alpha {alpha}: {x}");
        previous = x;
    }
}

#[test]
fn drag_changes_little_when_the_panel_count_doubles() {
    for (code, alpha) in [("0012", 0.0), ("0012", 4.0), ("2412", 2.0)] {
        let coarse = drag(code, 160, alpha, 3.0e6, None).unwrap();
        let fine = drag(code, 320, alpha, 3.0e6, None).unwrap();
        println!("NACA {code} alpha {alpha}: {coarse:.5} -> {fine:.5}");
        assert!((coarse - fine).abs() / fine < 0.03);
    }
}

#[test]
fn separated_flow_gives_no_drag_value() {
    let foil = naca4("0012", 160).unwrap();
    let solution = solve(&foil, 14.0).unwrap();
    let result = analyze(&foil, &solution, 1.0e6, None).unwrap();
    assert!(result.cd.is_none());
    let x = result
        .upper
        .separation
        .expect("separation on the upper surface");
    assert!(x < solver::viscous::DRAG_REFERENCE_X);
}

#[test]
fn drag_is_a_valid_number_in_the_attached_flow_range() {
    for code in ["0006", "0012", "2412", "4415"] {
        for alpha in [-2.0, 0.0, 2.0, 4.0] {
            for re in [5.0e5, 3.0e6, 1.0e7] {
                let cd = drag(code, 160, alpha, re, None);
                if let Some(cd) = cd {
                    assert!(
                        cd.is_finite() && cd > 0.002 && cd < 0.03,
                        "{code} {alpha} {re}: {cd}"
                    );
                }
            }
        }
    }
    // The base case must have a result.
    assert!(drag("0012", 160, 0.0, 3.0e6, None).is_some());
}
