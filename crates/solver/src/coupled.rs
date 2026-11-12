//! Two-way calculation: the boundary layer changes the outer flow.
//!
//! The boundary layer moves the outer flow away from the surface by the
//! displacement thickness. The model for this effect is a source on each
//! panel with the strength d(ue * delta_star)/ds. The calculation repeats
//! these steps until the lift does not change:
//! 1. Solve the outer flow with the sources.
//! 2. Calculate the boundary layer for the new surface speed.
//! 3. Calculate new sources from the boundary layer.
//!
//! Limits:
//! - This direct iteration is stable only for attached flow.
//! - The wake has no sources.
//! - After the limit station (`SOURCE_LIMIT_X`) and after a separation
//!   point, the displacement is a linear extension of the values before it.

use crate::geometry::Airfoil;
use crate::panel::{PanelSystem, Solution};
use crate::viscous::{analyze, SurfaceReport, ViscousResult};

/// The boundary layer result is used for the sources up to this chord
/// position. After it, the displacement is extended linearly.
pub const SOURCE_LIMIT_X: f64 = 0.95;
/// Part of the new sources that is used in each iteration.
pub const RELAXATION: f64 = 0.3;
pub const MAX_ITERATIONS: usize = 60;
/// The iteration stops when the lift coefficient changes by less than this.
pub const TOLERANCE: f64 = 1.0e-4;

#[derive(Debug, Clone)]
pub struct CoupledResult {
    /// Outer flow with the displacement effect.
    pub solution: Solution,
    /// Boundary layer for that outer flow.
    pub viscous: ViscousResult,
    /// Lift coefficient of the friction-free flow (no displacement effect).
    pub cl_inviscid: f64,
    pub iterations: usize,
    /// False if the lift did not settle in `MAX_ITERATIONS`.
    pub converged: bool,
}

/// Source strength on each panel from the boundary layer of one surface.
fn add_sources(surface: &SurfaceReport, sources: &mut [f64]) {
    let stations = &surface.boundary_layer.stations;
    let count = surface.s.len();
    // Mass defect m = ue * theta * H at the stations that have a result
    // before the limit station.
    let valid = stations
        .iter()
        .zip(&surface.x_over_c)
        .take_while(|(_, x)| **x <= SOURCE_LIMIT_X)
        .count();
    if valid < 3 {
        return;
    }
    let mut mass: Vec<f64> = stations[..valid]
        .iter()
        .map(|st| st.ue * st.theta * st.h)
        .collect();
    let slope = (mass[valid - 1] - mass[valid - 2]) / (surface.s[valid - 1] - surface.s[valid - 2]);
    for k in valid..count {
        mass.push(mass[valid - 1] + slope * (surface.s[k] - surface.s[valid - 1]));
    }
    for k in 1..count {
        let (lo, hi) = (k - 1, (k + 1).min(count - 1));
        let sigma = (mass[hi] - mass[lo]) / (surface.s[hi] - surface.s[lo]);
        sources[surface.panel[k]] = sigma;
    }
}

pub fn solve_coupled(
    foil: &Airfoil,
    system: &PanelSystem,
    alpha_deg: f64,
    re: f64,
    forced_transition_x: Option<f64>,
) -> Result<CoupledResult, String> {
    let m = system.panel_count();
    let mut sources = vec![0.0; m];
    let mut previous_cl = f64::NAN;
    let mut cl_inviscid = 0.0;
    for iteration in 0..MAX_ITERATIONS {
        let solution = system.solve(alpha_deg, Some(&sources))?;
        if iteration == 0 {
            cl_inviscid = solution.cl;
        }
        let viscous = analyze(foil, &solution, re, forced_transition_x)?;
        let converged = (solution.cl - previous_cl).abs() < TOLERANCE;
        if converged || iteration == MAX_ITERATIONS - 1 {
            return Ok(CoupledResult {
                solution,
                viscous,
                cl_inviscid,
                iterations: iteration + 1,
                converged,
            });
        }
        previous_cl = solution.cl;
        let mut target = vec![0.0; m];
        add_sources(&viscous.upper, &mut target);
        add_sources(&viscous.lower, &mut target);
        for (value, goal) in sources.iter_mut().zip(&target) {
            *value += RELAXATION * (goal - *value);
        }
    }
    unreachable!("the loop returns in its last iteration")
}
