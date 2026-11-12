//! Tests of the source panels that model the boundary layer displacement.

use solver::{naca4, PanelSystem};
use std::f64::consts::PI;

#[test]
fn zero_sources_give_the_same_result_as_no_sources() {
    let foil = naca4("2412", 120).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let a = system.solve(3.0, None).unwrap();
    let b = system.solve(3.0, Some(&vec![0.0; 120])).unwrap();
    assert!((a.cl - b.cl).abs() < 1e-12);
}

#[test]
fn flow_out_of_the_surface_is_equal_to_the_source_strength() {
    // Near the outer side of each panel, the normal velocity must be equal
    // to the source strength. Near the inner side it must be zero.
    let foil = naca4("0012", 120).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let m = system.panel_count();
    // A smooth source distribution with different values on each panel.
    let sigma: Vec<f64> = (0..m)
        .map(|j| 0.02 + 0.015 * (2.0 * PI * j as f64 / m as f64).sin())
        .collect();
    let alpha = 2.0;
    let solution = system.solve(alpha, Some(&sigma)).unwrap();
    let mut worst_out: f64 = 0.0;
    let mut worst_in: f64 = 0.0;
    // Panels away from the trailing edge, where the panels are very small,
    // and away from the leading edge. Near the stagnation point the normal
    // velocity changes quickly with the distance from the surface, thus a
    // check point at a small distance is not a good measure there.
    for j in 10..m - 10 {
        let (a, b) = (foil.points[j], foil.points[j + 1]);
        if 0.5 * (a[0] + b[0]) < 0.1 {
            continue;
        }
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let len = dx.hypot(dy);
        let (nx, ny) = (-dy / len, dx / len);
        let mid = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
        let offset = 0.05 * len;
        let outside = [mid[0] + nx * offset, mid[1] + ny * offset];
        let inside = [mid[0] - nx * offset, mid[1] - ny * offset];
        let vo = system.velocity_at(outside, alpha, &solution, Some(&sigma));
        let vi = system.velocity_at(inside, alpha, &solution, Some(&sigma));
        worst_out = worst_out.max((vo[0] * nx + vo[1] * ny - sigma[j]).abs());
        worst_in = worst_in.max((vi[0] * nx + vi[1] * ny).abs());
    }
    println!("largest error: outside {worst_out:.5}, inside {worst_in:.5}");
    assert!(worst_out < 0.002, "outside {worst_out}");
    assert!(worst_in < 0.001, "inside {worst_in}");
}

#[test]
fn total_flow_through_a_far_circle_is_equal_to_the_total_source() {
    let foil = naca4("2412", 120).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let m = system.panel_count();
    let sigma = vec![0.01; m];
    let solution = system.solve(1.0, Some(&sigma)).unwrap();
    let total: f64 = solution.panel_len.iter().map(|len| 0.01 * len).sum();
    let (radius, steps) = (5.0, 720);
    let mut flux = 0.0;
    for k in 0..steps {
        let angle = 2.0 * PI * (k as f64 + 0.5) / steps as f64;
        let (nx, ny) = (angle.cos(), angle.sin());
        let v = system.velocity_at(
            [0.5 + radius * nx, radius * ny],
            1.0,
            &solution,
            Some(&sigma),
        );
        flux += (v[0] * nx + v[1] * ny) * radius * 2.0 * PI / steps as f64;
    }
    println!("flux {flux:.6}, total source {total:.6}");
    assert!((flux - total).abs() / total < 1e-3);
}

#[test]
fn sources_do_not_change_the_surface_speed_of_a_friction_free_flow_without_sources() {
    // The velocity function agrees with the surface speed of the solution.
    let foil = naca4("0012", 120).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let solution = system.solve(4.0, None).unwrap();
    for j in [20, 40, 80, 100] {
        let (a, b) = (foil.points[j], foil.points[j + 1]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let len = dx.hypot(dy);
        let (tx, ty, nx, ny) = (dx / len, dy / len, -dy / len, dx / len);
        let mid = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
        let outside = [mid[0] + nx * 0.05 * len, mid[1] + ny * 0.05 * len];
        let v = system.velocity_at(outside, 4.0, &solution, None);
        let tangential = v[0] * tx + v[1] * ty;
        assert!(
            (tangential - solution.vt[j]).abs() < 0.02,
            "panel {j}: {tangential} {}",
            solution.vt[j]
        );
    }
}
