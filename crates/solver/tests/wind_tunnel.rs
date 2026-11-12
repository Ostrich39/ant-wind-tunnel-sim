//! Comparison of the friction drag with wind tunnel data for the NACA 0012
//! at a Reynolds number of 6 million. Data: `data/naca0012/` (REFERENCES.md,
//! items 4.2, 4.4, 1.3).
//!
//! The tolerances are the measured differences plus a margin. They record
//! the accuracy of the model. They are not design targets.

use solver::{analyze, naca4, solve};

const RE: f64 = 6.0e6;

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
        let row: Vec<f64> = line
            .split_whitespace()
            .map(|v| v.parse().unwrap())
            .collect();
        zones.last_mut().unwrap().push(row);
    }
    zones.retain(|zone| !zone.is_empty());
    zones
}

#[test]
fn drag_with_trip_strip_is_near_the_ladson_data() {
    // Ladson: transition tripped, three grit sizes. Model: forced transition
    // at 5% of chord (assumption: the trip position is not in the data file).
    // Compared at the same angle of attack, for angles where the model has
    // attached flow.
    let foil = naca4("0012", 200).unwrap();
    let zones = read("CLCD_Ladson_expdata.dat");
    assert_eq!(zones.len(), 3);
    let (mut count, mut worst): (usize, f64) = (0, 0.0);
    for row in zones.iter().flatten().filter(|row| row[0].abs() <= 8.5) {
        let (alpha, cd_data) = (row[0], row[2]);
        let solution = solve(&foil, alpha).unwrap();
        let cd = analyze(&foil, &solution, RE, Some(0.05))
            .unwrap()
            .cd
            .expect("drag value");
        let error = (cd - cd_data) / cd_data;
        worst = worst.max(error.abs());
        count += 1;
        assert!(
            error.abs() < 0.07,
            "alpha {alpha}: model {cd:.5}, data {cd_data:.5}"
        );
    }
    println!("{count} points, largest difference {:.1}%", worst * 100.0);
    assert!(count >= 20);
}

#[test]
fn lift_without_friction_correction_is_higher_than_the_ladson_data() {
    // This test records a known limit: the model lift is too high.
    let foil = naca4("0012", 200).unwrap();
    for row in read("CLCD_Ladson_expdata.dat").iter().flatten() {
        let (alpha, cl_data) = (row[0], row[1]);
        if !(2.0..=8.5).contains(&alpha) {
            continue;
        }
        let cl = solve(&foil, alpha).unwrap().cl;
        let ratio = cl / cl_data;
        assert!(
            ratio > 1.05 && ratio < 1.20,
            "alpha {alpha}: ratio {ratio:.3}"
        );
    }
}

#[test]
fn drag_with_free_transition_is_near_the_abbott_data() {
    // Abbott and von Doenhoff: no trip strip. Compared at the same lift
    // coefficient. The data points were digitized from a plot.
    let foil = naca4("0012", 200).unwrap();
    let mut polar = Vec::new();
    for step in -48..=48 {
        let solution = solve(&foil, step as f64 * 0.25).unwrap();
        if let Some(cd) = analyze(&foil, &solution, RE, None).unwrap().cd {
            polar.push((solution.cl, cd));
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
    let zones = read("0012.abbottdata.cd.dat");
    let (mut low_lift, mut high_lift) = (0, 0);
    for row in zones.iter().flatten() {
        let (cl, cd_data) = (row[0], row[1]);
        let Some(cd) = model_cd(cl) else { continue };
        let error = (cd - cd_data) / cd_data;
        println!(
            "cl {cl:+.3}: model {cd:.5}, data {cd_data:.5}, {:+.1}%",
            error * 100.0
        );
        if cl.abs() <= 0.5 {
            // Low lift: the model is near the data.
            assert!(error.abs() < 0.07, "cl {cl}");
            low_lift += 1;
        } else {
            // Higher lift: the model drag is too low (known limit of a
            // one-way boundary layer calculation).
            assert!(error < 0.0 && error > -0.20, "cl {cl}");
            high_lift += 1;
        }
    }
    assert!(low_lift >= 8 && high_lift >= 6);
}
