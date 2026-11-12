//! Tests of the experimental two-way calculation (displacement sources with
//! direct iteration). The program does not use this calculation, because the
//! lift result depends on the cut-off station `SOURCE_LIMIT_X` (see GOALS.md).
//! These tests record its behavior for the NACA 0012 with a trip strip.

use solver::{naca4, solve_coupled, PanelSystem};

const RE: f64 = 6.0e6;
const TRIP: Option<f64> = Some(0.05);

/// Rows (angle, lift, drag) of the Ladson data, all three grit sizes.
fn ladson() -> Vec<Vec<f64>> {
    let path = format!(
        "{}/../../data/naca0012/CLCD_Ladson_expdata.dat",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(str::trim)
        .filter(|l| l.starts_with(|c: char| c.is_ascii_digit() || c == '-' || c == '.'))
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect()
}

#[test]
fn iteration_converges_for_attached_flow() {
    let foil = naca4("0012", 200).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    for alpha in [0.0, 2.0, 4.0, 6.0, 8.0, 10.0] {
        let result = solve_coupled(&foil, &system, alpha, RE, TRIP).unwrap();
        assert!(result.converged, "alpha {alpha}");
        assert!(
            result.iterations <= 20,
            "alpha {alpha}: {} iterations",
            result.iterations
        );
    }
}

#[test]
fn symmetric_case_has_no_lift() {
    let foil = naca4("0012", 200).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let result = solve_coupled(&foil, &system, 0.0, RE, TRIP).unwrap();
    assert!(result.solution.cl.abs() < 1e-6);
}

#[test]
fn displacement_lowers_the_lift() {
    let foil = naca4("0012", 200).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    for alpha in [2.0, 4.0, 6.0, 8.0, 10.0] {
        let result = solve_coupled(&foil, &system, alpha, RE, TRIP).unwrap();
        let reduction = 1.0 - result.solution.cl / result.cl_inviscid;
        println!("alpha {alpha}: lift reduction {:.1}%", reduction * 100.0);
        assert!(
            reduction > 0.01 && reduction < 0.08,
            "alpha {alpha}: {reduction}"
        );
    }
}

#[test]
fn lift_is_nearer_to_the_ladson_data_but_still_too_high() {
    let foil = naca4("0012", 200).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let mut count = 0;
    for row in ladson().iter().filter(|row| (1.9..=11.2).contains(&row[0])) {
        let (alpha, cl_data) = (row[0], row[1]);
        let result = solve_coupled(&foil, &system, alpha, RE, TRIP).unwrap();
        let ratio = result.solution.cl / cl_data;
        assert!(
            ratio > 1.01 && ratio < 1.14,
            "alpha {alpha}: ratio {ratio:.3}"
        );
        assert!(
            (result.solution.cl - cl_data).abs() < (result.cl_inviscid - cl_data).abs(),
            "alpha {alpha}"
        );
        count += 1;
    }
    assert!(count >= 15);
}

#[test]
fn drag_is_near_the_ladson_data_to_eleven_degrees() {
    let foil = naca4("0012", 200).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let (mut count, mut worst): (usize, f64) = (0, 0.0);
    for row in ladson().iter().filter(|row| row[0].abs() <= 11.2) {
        let (alpha, cd_data) = (row[0], row[2]);
        let result = solve_coupled(&foil, &system, alpha, RE, TRIP).unwrap();
        let cd = result
            .viscous
            .cd
            .unwrap_or_else(|| panic!("no drag value at {alpha}"));
        let error = (cd - cd_data) / cd_data;
        worst = worst.max(error.abs());
        assert!(
            error.abs() < 0.08,
            "alpha {alpha}: model {cd:.5}, data {cd_data:.5}"
        );
        count += 1;
    }
    println!("{count} points, largest difference {:.1}%", worst * 100.0);
    assert!(count >= 26);
}
