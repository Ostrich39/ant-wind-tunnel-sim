//! Friction drag of an airfoil: the boundary layer calculation applied to
//! the upper and lower surfaces of a panel solution.
//!
//! Limits of this method:
//! - One-way calculation. The boundary layer does not change the pressure,
//!   thus the lift has no friction correction.
//! - Valid only for attached flow. If the turbulent boundary layer
//!   separates before the reference station, there is no drag result.
//! - The laminar separation bubble is not modeled.

use crate::boundary_layer::{march, SurfaceResult, TransitionCause};
use crate::geometry::Airfoil;
use crate::panel::Solution;

/// Position (x/c) of the drag reference station. The inviscid flow has a
/// stagnation point at the trailing edge, which a real flow does not have.
/// Thus the drag formula uses the station at this position, not the last one.
pub const DRAG_REFERENCE_X: f64 = 0.98;

#[derive(Debug, Clone)]
pub struct SurfaceReport {
    pub boundary_layer: SurfaceResult,
    /// Chord position (x/c) of each input station.
    pub x_over_c: Vec<f64>,
    /// Surface distance of each input station, divided by the chord.
    pub s: Vec<f64>,
    /// Panel index of each input station. The first station is the
    /// stagnation point, which is between two panels.
    pub panel: Vec<usize>,
    /// Chord position (x/c) where turbulent flow starts, and the cause.
    pub transition: Option<(f64, TransitionCause)>,
    /// Chord position (x/c) of turbulent separation, if it occurs.
    pub separation: Option<f64>,
    /// Panel index at the transition position.
    pub transition_panel: Option<usize>,
    /// Panel index at the separation position.
    pub separation_panel: Option<usize>,
    /// Drag coefficient contribution of this surface.
    pub cd: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ViscousResult {
    pub upper: SurfaceReport,
    pub lower: SurfaceReport,
    /// Drag coefficient (Squire-Young formula). `None` if the flow separates
    /// before the reference station on one surface or on the two surfaces.
    pub cd: Option<f64>,
    /// Chord position (x/c) of the stagnation point.
    pub stagnation_x: f64,
}

/// Squire-Young formula: drag coefficient contribution of one surface from
/// the boundary layer at the trailing edge.
pub fn squire_young(theta: f64, ue: f64, h: f64) -> f64 {
    2.0 * theta * ue.powf(0.5 * (h + 5.0))
}

pub fn analyze(
    foil: &Airfoil,
    solution: &Solution,
    re: f64,
    forced_transition_x: Option<f64>,
) -> Result<ViscousResult, String> {
    let vt = &solution.vt;
    let len = &solution.panel_len;
    let m = vt.len();
    let (chord, x_min) = (foil.chord(), foil.x_min());
    let x_of = |j: usize| (solution.mid[j][0] - x_min) / chord;

    // Stagnation point: the surface speed changes from negative (flow to
    // the trailing edge along the lower surface) to positive.
    let j = (0..m - 1)
        .filter(|&j| vt[j] < 0.0 && vt[j + 1] >= 0.0)
        .min_by(|&a, &b| x_of(a).total_cmp(&x_of(b)))
        .ok_or("no stagnation point was found on the airfoil")?;
    let gap = 0.5 * (len[j] + len[j + 1]);
    let fraction = -vt[j] / (vt[j + 1] - vt[j]);
    let stagnation_x = x_of(j) + fraction * (x_of(j + 1) - x_of(j));

    let build = |indices: Vec<usize>, first_gap: f64| -> Result<SurfaceReport, String> {
        let mut s = vec![0.0];
        let mut ue = vec![0.0];
        let mut x = vec![stagnation_x];
        let mut panel = vec![indices[0]];
        let mut distance = first_gap;
        for (k, &i) in indices.iter().enumerate() {
            if k > 0 {
                distance += 0.5 * (len[indices[k - 1]] + len[i]);
            }
            // Stop if the surface speed changes sign again (not attached flow).
            if distance <= *s.last().unwrap() || vt[i].abs() <= 0.0 {
                break;
            }
            s.push(distance / chord);
            ue.push(vt[i].abs());
            x.push(x_of(i));
            panel.push(i);
        }
        // Convert a forced transition position (x/c) into a surface distance.
        let forced = forced_transition_x.map(|limit| {
            x.iter()
                .zip(&s)
                .find(|(xi, _)| **xi >= limit)
                .map_or(f64::INFINITY, |(_, si)| *si)
        });
        let boundary_layer = march(&s, &ue, re, forced)?;
        let transition = boundary_layer.transition.map(|(i, cause)| (x[i], cause));
        let separation = boundary_layer.separation.map(|i| x[i]);
        // Drag reference station: the last calculated station at or before
        // the reference position.
        let reference = boundary_layer
            .stations
            .iter()
            .zip(&x)
            .rfind(|(_, xi)| **xi <= DRAG_REFERENCE_X)
            .map(|(station, _)| *station);
        let reaches_reference = x[boundary_layer.stations.len() - 1] >= DRAG_REFERENCE_X
            || boundary_layer.stations.len() == x.len();
        let cd = match reference {
            Some(st) if reaches_reference => Some(squire_young(st.theta, st.ue, st.h)),
            _ => None,
        };
        let transition_panel = boundary_layer.transition.map(|(i, _)| panel[i]);
        let separation_panel = boundary_layer.separation.map(|i| panel[i]);
        Ok(SurfaceReport {
            boundary_layer,
            x_over_c: x,
            s,
            panel,
            transition,
            separation,
            transition_panel,
            separation_panel,
            cd,
        })
    };

    let upper = build((j + 1..m).collect(), (1.0 - fraction) * gap)?;
    let lower = build((0..=j).rev().collect(), fraction * gap)?;
    let cd = match (upper.cd, lower.cd) {
        (Some(a), Some(b)) => Some(a + b),
        _ => None,
    };
    Ok(ViscousResult {
        upper,
        lower,
        cd,
        stagnation_x,
    })
}
