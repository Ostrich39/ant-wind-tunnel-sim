//! Wake of the airfoil: the line that the flow follows after the trailing
//! edge. The boundary layers of the two surfaces join in the wake. The
//! wake has source panels that model its displacement effect.

use crate::panel::{source_velocity, PanelSystem, Solution};

#[derive(Debug, Clone)]
pub struct Wake {
    /// End points of the wake panels. The first point is the trailing edge.
    pub points: Vec<[f64; 2]>,
    /// Midpoint of each wake panel (the wake stations).
    pub mid: Vec<[f64; 2]>,
    /// Length of each wake panel.
    pub len: Vec<f64>,
    /// Direction angle of each wake panel.
    pub theta: Vec<f64>,
    /// Friction-free flow speed along the wake at each station.
    pub ue_inviscid: Vec<f64>,
}

/// Growth ratio r of a geometric series: first * (1 + r + ... + r^(n-1)) = total.
fn growth_ratio(first: f64, total: f64, n: usize) -> f64 {
    let sum = |r: f64| (0..n).map(|k| r.powi(k as i32)).sum::<f64>() * first;
    let (mut lo, mut hi) = (1.0, 10.0);
    if sum(lo) >= total {
        return 1.0;
    }
    for _ in 0..80 {
        let midpoint = 0.5 * (lo + hi);
        if sum(midpoint) < total {
            lo = midpoint;
        } else {
            hi = midpoint;
        }
    }
    0.5 * (lo + hi)
}

impl Wake {
    /// Builds the wake as a streamline of the friction-free flow.
    /// `length` is the wake length divided by the chord.
    pub fn new(
        system: &PanelSystem,
        alpha_deg: f64,
        inviscid: &Solution,
        panels: usize,
        length: f64,
    ) -> Self {
        let m = system.panel_count();
        let (angles, lengths) = (system.angles(), system.lengths());
        let trailing_edge = system.points()[0];
        let velocity = |p: [f64; 2]| system.velocity_at(p, alpha_deg, inviscid, None);
        let unit = |v: [f64; 2]| {
            let size = v[0].hypot(v[1]);
            [v[0] / size, v[1] / size]
        };
        // First direction: the bisector of the trailing edge angle.
        let upper = [angles[m - 1].cos(), angles[m - 1].sin()];
        let lower = [-angles[0].cos(), -angles[0].sin()];
        let mut direction = unit([upper[0] + lower[0], upper[1] + lower[1]]);

        let first = 0.5 * (lengths[0] + lengths[m - 1]);
        let ratio = growth_ratio(first, length * system.chord(), panels);
        let mut points = vec![trailing_edge];
        let mut step = first;
        for k in 0..panels {
            let start = points[k];
            if k > 0 {
                direction = unit(velocity(start));
            }
            // Predictor and corrector (Heun's method).
            let trial = [
                start[0] + step * direction[0],
                start[1] + step * direction[1],
            ];
            let second = unit(velocity(trial));
            let average = unit([direction[0] + second[0], direction[1] + second[1]]);
            points.push([start[0] + step * average[0], start[1] + step * average[1]]);
            step *= ratio;
        }

        let mut mid = Vec::with_capacity(panels);
        let mut len = Vec::with_capacity(panels);
        let mut theta = Vec::with_capacity(panels);
        let mut ue_inviscid = Vec::with_capacity(panels);
        for k in 0..panels {
            let (a, b) = (points[k], points[k + 1]);
            let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
            let size = dx.hypot(dy);
            let center = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
            let v = velocity(center);
            mid.push(center);
            len.push(size);
            theta.push(dy.atan2(dx));
            ue_inviscid.push((v[0] * dx + v[1] * dy) / size);
        }
        Self {
            points,
            mid,
            len,
            theta,
            ue_inviscid,
        }
    }

    pub fn panel_count(&self) -> usize {
        self.len.len()
    }
}

/// Effect of the sources on the flow speed at each station.
///
/// Stations: the airfoil panel midpoints (index 0 to m - 1, speed in the
/// direction of the point order), then the wake stations (speed along the
/// wake). `table[i][j]` is the change of the speed at station i for a unit
/// source on panel j (airfoil panels, then wake panels).
pub struct SourceInfluence {
    pub table: Vec<Vec<f64>>,
    /// `dgamma[j]` is the change of the vortex strengths for a unit source
    /// on panel j (airfoil panels, then wake panels).
    pub dgamma: Vec<Vec<f64>>,
}

impl SourceInfluence {
    /// `airfoil_response` is the result of `PanelSystem::source_response`.
    /// It depends only on the shape, thus the caller can keep it.
    pub fn new(
        system: &PanelSystem,
        airfoil_response: &(Vec<Vec<f64>>, Vec<Vec<f64>>),
        wake: &Wake,
    ) -> Self {
        let m = system.panel_count();
        let nw = wake.panel_count();
        let total = m + nw;
        let (dgamma_airfoil, dvt_airfoil) = airfoil_response;
        let mut table = vec![vec![0.0; total]; total];
        let mut dgamma = Vec::with_capacity(total);

        // Speed at each wake station for each unit of vortex strength and
        // each unit of airfoil source strength.
        let rows: Vec<(Vec<f64>, Vec<f64>)> = (0..nw)
            .map(|w| system.influence_at(wake.mid[w], wake.theta[w]))
            .collect();
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();

        // Airfoil sources.
        for p in 0..m {
            for i in 0..m {
                table[i][p] = dvt_airfoil[p][i];
            }
            for w in 0..nw {
                table[m + w][p] = dot(&rows[w].0, &dgamma_airfoil[p]) + rows[w].1[p];
            }
            dgamma.push(dgamma_airfoil[p].clone());
        }

        // Wake sources.
        let (mid, angles) = (system.midpoints(), system.angles());
        for q in 0..nw {
            let (a, b) = (wake.points[q], wake.points[q + 1]);
            let mut normal = vec![0.0; m];
            let mut tangential = vec![0.0; m];
            for i in 0..m {
                let [vx, vy] = source_velocity(mid[i], a, b);
                let (si, ci) = angles[i].sin_cos();
                normal[i] = -vx * si + vy * ci;
                tangential[i] = vx * ci + vy * si;
            }
            let gamma = system.gamma_response(&normal);
            let from_gamma = system.tangential_from_gamma(&gamma);
            for i in 0..m {
                table[i][m + q] = from_gamma[i] + tangential[i];
            }
            for w in 0..nw {
                let direct = if w == q {
                    0.0
                } else {
                    let [vx, vy] = source_velocity(wake.mid[w], a, b);
                    vx * wake.theta[w].cos() + vy * wake.theta[w].sin()
                };
                table[m + w][m + q] = dot(&rows[w].0, &gamma) + direct;
            }
            dgamma.push(gamma);
        }
        Self { table, dgamma }
    }
}
