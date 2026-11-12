//! Tests of the two-way calculation (simultaneous solution of the outer
//! flow and the boundary layer), including a comparison with wind tunnel
//! data for the NACA 0012. Data: `data/naca0012/` (REFERENCES.md, section 4).

use solver::interaction::{FrictionSystem, Options};
use solver::naca4_for_friction;

const RE: f64 = 6.0e6;

fn system(code: &str, panels: usize) -> FrictionSystem {
    FrictionSystem::new(naca4_for_friction(code, panels).unwrap()).unwrap()
}

fn options(re: f64, trip: Option<f64>) -> Options {
    Options {
        re,
        forced_transition_x: trip,
        n_crit: 9.0,
    }
}

/// Reads a data file. Returns one list of rows for each zone.
fn read(name: &str) -> Vec<Vec<Vec<f64>>> {
    let path = format!("{}/../../data/naca0012/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut zones: Vec<Vec<Vec<f64>>> = vec![vec![]];
    for line in text.lines().map(str::trim) {
        if line.is_empty() || line.starts_with('#') || line.starts_with("variables") {
            continue;
        }
        if line.starts_with("zone") {
            zones.push(vec![]);
            continue;
        }
        zones.last_mut().unwrap().push(
            line.split_whitespace()
                .map(|v| v.parse().unwrap())
                .collect(),
        );
    }
    zones.retain(|zone| !zone.is_empty());
    zones
}

#[test]
fn symmetric_airfoil_at_zero_angle() {
    let (result, _) = system("0012", 120)
        .solve(0.0, &options(RE, None), None)
        .unwrap();
    assert!(result.converged);
    assert!(result.iterations <= 10, "{} iterations", result.iterations);
    assert!(result.cl.abs() < 1e-6, "cl {}", result.cl);
    assert_eq!(result.transition_upper, result.transition_lower);
    assert!(result.separation_upper.is_none() && result.separation_lower.is_none());
}

#[test]
fn opposite_angles_give_opposite_lift_and_equal_drag() {
    let friction = system("0012", 120);
    let (up, _) = friction.solve(4.0, &options(RE, None), None).unwrap();
    let (down, _) = friction.solve(-4.0, &options(RE, None), None).unwrap();
    assert!(up.converged && down.converged);
    assert!((up.cl + down.cl).abs() < 1e-6);
    assert!((up.cd - down.cd).abs() < 1e-7);
    assert_eq!(up.transition_upper, down.transition_lower);
}

#[test]
fn friction_decreases_the_lift() {
    let friction = system("2412", 120);
    for alpha in [0.0, 4.0, 8.0] {
        let (result, _) = friction.solve(alpha, &options(3.0e6, None), None).unwrap();
        assert!(result.converged, "alpha {alpha}");
        let ratio = result.cl / result.cl_inviscid;
        assert!(ratio > 0.80 && ratio < 0.97, "alpha {alpha}: ratio {ratio}");
    }
}

#[test]
fn result_changes_little_with_the_panel_count() {
    let (coarse, _) = system("0012", 120)
        .solve(4.0, &options(RE, Some(0.05)), None)
        .unwrap();
    let (fine, _) = system("0012", 200)
        .solve(4.0, &options(RE, Some(0.05)), None)
        .unwrap();
    println!(
        "cl {:.4} -> {:.4}, cd {:.5} -> {:.5}",
        coarse.cl, fine.cl, coarse.cd, fine.cd
    );
    assert!((coarse.cl - fine.cl).abs() / fine.cl < 0.005);
    assert!((coarse.cd - fine.cd).abs() / fine.cd < 0.03);
}

#[test]
fn start_from_an_earlier_result_gives_the_same_solution() {
    let friction = system("0012", 120);
    let settings = options(RE, Some(0.05));
    let (_, warm) = friction.solve(2.0, &settings, None).unwrap();
    let (from_warm, _) = friction.solve(3.0, &settings, Some(&warm)).unwrap();
    let (from_cold, _) = friction.solve(3.0, &settings, None).unwrap();
    assert!(from_warm.converged && from_cold.converged);
    assert!((from_warm.cl - from_cold.cl).abs() < 1e-4);
    assert!((from_warm.cd - from_cold.cd).abs() < 1e-5);
}

#[test]
fn drag_decreases_when_reynolds_number_increases() {
    let friction = system("0012", 120);
    let values: Vec<f64> = [1.0e6, 3.0e6, 6.0e6]
        .iter()
        .map(|re| {
            let (result, _) = friction.solve(0.0, &options(*re, None), None).unwrap();
            assert!(result.converged);
            result.cd
        })
        .collect();
    println!("{values:?}");
    assert!(values.windows(2).all(|w| w[1] < w[0]));
}

#[test]
fn lift_and_drag_with_trip_strip_are_near_the_ladson_data() {
    // Ladson: transition tripped, Re 6 million. Model: trip at 5% of chord
    // (assumption: the trip position is not in the data file).
    let friction = system("0012", 120);
    let settings = options(RE, Some(0.05));
    let (mut count, mut worst_lift, mut worst_drag): (usize, f64, f64) = (0, 0.0, 0.0);
    for row in read("CLCD_Ladson_expdata.dat").iter().flatten() {
        let (alpha, cl_data, cd_data) = (row[0], row[1], row[2]);
        if alpha.abs() > 14.5 {
            continue;
        }
        let (result, _) = friction.solve_robust(alpha, &settings, None).unwrap();
        assert!(result.converged, "alpha {alpha}");
        let drag_error = (result.cd - cd_data) / cd_data;
        println!("alpha {alpha:>6.2}: cl {:+.4} (data {cl_data:+.4}), cd {:.5} (data {cd_data:.5}, {:+.1}%)", result.cl, result.cd, drag_error * 100.0);
        // Lift: absolute difference, because the lift is near zero at 0 deg.
        let lift_error = (result.cl - cl_data).abs();
        worst_lift = worst_lift.max(lift_error / cl_data.abs().max(0.2));
        worst_drag = worst_drag.max(drag_error.abs());
        count += 1;
    }
    println!(
        "{count} points: largest lift difference {:.1}%, largest drag difference {:.1}%",
        worst_lift * 100.0,
        worst_drag * 100.0
    );
    assert!(count >= 36);
    // Measured: lift 8.2%, drag 8.7%.
    assert!(worst_lift < 0.10);
    assert!(worst_drag < 0.10);
}

#[test]
fn near_maximum_lift_the_drag_is_too_low() {
    // Ladson data from 14.5 deg to the maximum lift at approximately 17 deg.
    // Measured: lift 3.7% to 5.5% high, drag 3% to 17% low.
    let friction = system("0012", 120);
    let settings = options(RE, Some(0.05));
    let mut count = 0;
    for row in read("CLCD_Ladson_expdata.dat").iter().flatten() {
        let (alpha, cl_data, cd_data) = (row[0], row[1], row[2]);
        if !(14.5..=17.3).contains(&alpha) || cl_data < 1.4 {
            continue;
        }
        let (result, _) = friction.solve_robust(alpha, &settings, None).unwrap();
        assert!(result.converged, "alpha {alpha}");
        let (lift, drag) = (
            (result.cl - cl_data) / cl_data,
            (result.cd - cd_data) / cd_data,
        );
        println!(
            "alpha {alpha:.2}: lift {:+.1}%, drag {:+.1}%",
            lift * 100.0,
            drag * 100.0
        );
        assert!(lift > 0.0 && lift < 0.08, "alpha {alpha}: lift {lift}");
        assert!(drag < 0.0 && drag > -0.22, "alpha {alpha}: drag {drag}");
        count += 1;
    }
    assert!(count >= 6);
}

#[test]
fn abrupt_stall_is_not_predicted() {
    // Known limit. In the Ladson data, the airfoil stalls between 17 deg and
    // 18 deg: the lift decreases from 1.6 to 1.0. The model gives a slow
    // growth of the separation from the trailing edge, and its lift does not
    // decrease. This test records the limit. Change it when the model can
    // predict the stall.
    let friction = system("0012", 120);
    let (result, _) = friction
        .solve_robust(18.02, &options(RE, Some(0.05)), None)
        .unwrap();
    let cl_data = 0.9967;
    assert!(result.converged);
    assert!(result.cl > 1.5 * cl_data, "model lift {}", result.cl);
    // The model does show separation on the upper surface at this angle.
    assert!(result.separation_upper.is_some_and(|x| x < 0.9));
}

#[test]
fn solution_exists_for_the_normal_angle_range() {
    // Each condition must converge from the friction-free start, without
    // help from a solution at a different angle.
    for (code, re, last_angle) in [
        ("0012", 6.0e6, 18),
        ("2412", 3.0e6, 16),
        ("4412", 3.0e6, 14),
        ("4415", 1.0e6, 12),
        ("0012", 2.0e5, 12),
        ("0006", 5.0e5, 6),
    ] {
        let friction = system(code, 120);
        for step in -2..=last_angle / 2 {
            let alpha = step as f64 * 2.0;
            let (result, _) = friction.solve(alpha, &options(re, None), None).unwrap();
            assert!(result.converged, "NACA {code}, Re {re:.0e}, alpha {alpha}");
            assert!(
                result.cd > 0.003 && result.cd < 0.06,
                "NACA {code} alpha {alpha}: cd {}",
                result.cd
            );
        }
    }
}

#[test]
fn drag_with_free_transition_against_the_abbott_data() {
    // Abbott and von Doenhoff: no trip strip, digitized from a plot.
    let friction = system("0012", 120);
    let settings = options(RE, None);
    let mut polar = Vec::new();
    for step in -28..=28 {
        let (result, _) = friction
            .solve_robust(step as f64 * 0.5, &settings, None)
            .unwrap();
        if result.converged {
            polar.push((result.cl, result.cd));
        }
    }
    let model_cd = |cl: f64| {
        polar
            .windows(2)
            .find(|w| w[0].0 <= cl && cl <= w[1].0)
            .map(|w| {
                let t = (cl - w[0].0) / (w[1].0 - w[0].0);
                w[0].1 + t * (w[1].1 - w[0].1)
            })
    };
    let (mut count, mut worst): (usize, f64) = (0, 0.0);
    for row in read("0012.abbottdata.cd.dat").iter().flatten() {
        let (cl, cd_data) = (row[0], row[1]);
        let Some(cd) = model_cd(cl) else { continue };
        let error = (cd - cd_data) / cd_data;
        println!(
            "cl {cl:+.3}: model {cd:.5}, data {cd_data:.5}, {:+.1}%",
            error * 100.0
        );
        worst = worst.max(error.abs());
        count += 1;
    }
    println!(
        "{count} points: largest drag difference {:.1}%",
        worst * 100.0
    );
    assert!(count >= 14);
    assert!(worst < 0.25);
}

#[test]
fn deep_stall_has_no_converged_result() {
    // Thin airfoil at a high angle and a low Reynolds number: the flow
    // separates at the leading edge. The method has no solution here, and
    // it must report this.
    let friction = system("0006", 120);
    // An error result is also correct here: it also reports "no solution".
    if let Ok((result, _)) = friction.solve_robust(14.0, &options(5.0e5, None), None) {
        assert!(!result.converged);
    }
}

#[test]
fn bad_input_gives_an_error() {
    let friction = system("0012", 120);
    assert!(friction.solve(2.0, &options(0.0, None), None).is_err());
    assert!(friction.solve(2.0, &options(f64::NAN, None), None).is_err());
    assert!(friction.solve(f64::NAN, &options(RE, None), None).is_err());
}

#[test]
fn one_solution_is_fast_enough_for_live_update() {
    // Goal: new results in less than 0.1 second (optimized build).
    let friction = system("2412", 120);
    let start = std::time::Instant::now();
    let (result, _) = friction.solve(4.0, &options(3.0e6, None), None).unwrap();
    let elapsed = start.elapsed().as_secs_f64();
    assert!(result.converged);
    println!(
        "{:.0} ms, {} iterations",
        elapsed * 1000.0,
        result.iterations
    );
    if !cfg!(debug_assertions) {
        assert!(elapsed < 0.1, "{elapsed} s");
    }
}
