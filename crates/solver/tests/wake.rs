//! Tests of the wake line and of the source influence table.

use solver::panel::source_velocity;
use solver::wake::{SourceInfluence, Wake};
use solver::{naca4, PanelSystem};
use std::f64::consts::PI;

fn setup(alpha: f64) -> (PanelSystem, solver::Solution, Wake) {
    let foil = naca4("2412", 120).unwrap();
    let system = PanelSystem::new(&foil).unwrap();
    let inviscid = system.solve(alpha, None).unwrap();
    let wake = Wake::new(&system, alpha, &inviscid, 24, 1.0);
    (system, inviscid, wake)
}

#[test]
fn wake_follows_the_flow() {
    let alpha = 5.0;
    let (system, inviscid, wake) = setup(alpha);
    let total: f64 = wake.len.iter().sum();
    assert!((total - 1.0).abs() < 1e-6, "wake length {total}");
    assert_eq!(wake.points[0], system.points()[0]);
    // The flow direction at each station agrees with the wake direction.
    // The first panels are near the trailing edge stagnation point, where
    // the direction changes quickly.
    for w in 3..wake.panel_count() {
        let v = system.velocity_at(wake.mid[w], alpha, &inviscid, None);
        let (s, c) = wake.theta[w].sin_cos();
        let across = (-v[0] * s + v[1] * c) / v[0].hypot(v[1]);
        assert!(across.abs() < 0.01, "station {w}: {across}");
        assert!(wake.ue_inviscid[w] > 0.5 && wake.ue_inviscid[w] < 1.1);
    }
    // One chord behind the airfoil, the lift still turns the flow down.
    // Thus the wake angle is between zero and the angle of attack.
    let last = wake.theta.last().unwrap().to_degrees();
    assert!(last > 0.0 && last < alpha, "{last}");
}

#[test]
fn table_agrees_with_the_panel_solution_for_airfoil_sources() {
    let alpha = 3.0;
    let (system, inviscid, wake) = setup(alpha);
    let m = system.panel_count();
    let influence = SourceInfluence::new(&system, &system.source_response(), &wake);
    let sources: Vec<f64> = (0..m).map(|i| 0.01 * (0.3 * i as f64).sin()).collect();
    let with_sources = system.solve(alpha, Some(&sources)).unwrap();
    for i in 0..m {
        let predicted: f64 = inviscid.vt[i]
            + (0..m)
                .map(|p| influence.table[i][p] * sources[p])
                .sum::<f64>();
        assert!((predicted - with_sources.vt[i]).abs() < 1e-10, "panel {i}");
    }
    // Speed at the wake stations.
    for w in 0..wake.panel_count() {
        let v = system.velocity_at(wake.mid[w], alpha, &with_sources, Some(&sources));
        let (s, c) = wake.theta[w].sin_cos();
        let direct = v[0] * c + v[1] * s;
        let predicted: f64 = wake.ue_inviscid[w]
            + (0..m)
                .map(|p| influence.table[m + w][p] * sources[p])
                .sum::<f64>();
        assert!((predicted - direct).abs() < 1e-9, "wake station {w}");
    }
}

#[test]
fn wake_sources_keep_the_flow_tangent_and_conserve_the_flow() {
    let alpha = 3.0;
    let (system, inviscid, wake) = setup(alpha);
    let (m, nw) = (system.panel_count(), wake.panel_count());
    let influence = SourceInfluence::new(&system, &system.source_response(), &wake);
    let strength: Vec<f64> = (0..nw).map(|q| 0.02 / (1.0 + q as f64)).collect();
    // Vortex strengths with the wake sources.
    let mut gamma = inviscid.gamma.clone();
    for q in 0..nw {
        for (g, d) in gamma.iter_mut().zip(&influence.dgamma[m + q]) {
            *g += d * strength[q];
        }
    }
    let field = |point: [f64; 2]| {
        let mut solution = inviscid.clone();
        solution.gamma = gamma.clone();
        let mut v = system.velocity_at(point, alpha, &solution, None);
        for q in 0..nw {
            let [vx, vy] = source_velocity(point, wake.points[q], wake.points[q + 1]);
            v[0] += vx * strength[q];
            v[1] += vy * strength[q];
        }
        v
    };
    // 1. No flow goes through the airfoil surface.
    let (mid, angles) = (system.midpoints(), system.angles());
    for i in (0..m).step_by(7) {
        let (s, c) = angles[i].sin_cos();
        let outside = [mid[i][0] - 1e-6 * s, mid[i][1] + 1e-6 * c];
        let v = field(outside);
        let normal = -v[0] * s + v[1] * c;
        assert!(normal.abs() < 2e-3, "panel {i}: normal velocity {normal}");
    }
    // 2. The flow through a far circle is equal to the total source.
    let total: f64 = strength.iter().zip(&wake.len).map(|(s, l)| s * l).sum();
    let (radius, steps) = (20.0, 720);
    let mut flux = 0.0;
    for k in 0..steps {
        let phi = 2.0 * PI * (k as f64 + 0.5) / steps as f64;
        let v = field([0.5 + radius * phi.cos(), radius * phi.sin()]);
        flux += (v[0] * phi.cos() + v[1] * phi.sin()) * radius * 2.0 * PI / steps as f64;
    }
    println!("flux {flux:.6}, total source {total:.6}");
    assert!((flux - total).abs() / total < 1e-3);
    // 3. The table gives the same speed on the surface as the field.
    for i in (0..m).step_by(11) {
        let (s, c) = angles[i].sin_cos();
        let outside = [mid[i][0] - 1e-6 * s, mid[i][1] + 1e-6 * c];
        let v = field(outside);
        let direct = v[0] * c + v[1] * s;
        let predicted: f64 = inviscid.vt[i]
            + (0..nw)
                .map(|q| influence.table[i][m + q] * strength[q])
                .sum::<f64>();
        assert!(
            (predicted - direct).abs() < 2e-3,
            "panel {i}: {predicted} {direct}"
        );
    }
    // 4. The table gives the same speed at the wake stations as the field
    // (the field value is the average of the two sides of the wake panel).
    for w in 0..nw {
        let (s, c) = wake.theta[w].sin_cos();
        let offset = 1e-7;
        let above = field([wake.mid[w][0] - offset * s, wake.mid[w][1] + offset * c]);
        let below = field([wake.mid[w][0] + offset * s, wake.mid[w][1] - offset * c]);
        let direct = 0.5 * ((above[0] + below[0]) * c + (above[1] + below[1]) * s);
        let predicted: f64 = wake.ue_inviscid[w]
            + (0..nw)
                .map(|q| influence.table[m + w][m + q] * strength[q])
                .sum::<f64>();
        assert!(
            (predicted - direct).abs() < 1e-6,
            "wake station {w}: {predicted} {direct}"
        );
    }
}
