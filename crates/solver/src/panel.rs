//! Inviscid linear-vortex panel method.
//!
//! The airfoil surface is divided into straight panels. Each panel has a
//! vortex sheet with a strength that changes linearly along the panel.
//! The solver finds the strengths that make the flow tangent to the
//! surface at each panel midpoint, with smooth flow at the trailing edge
//! (Kutta condition).
//!
//! Reference: Kuethe and Chow, "Foundations of Aerodynamics", 5th ed.,
//! section 5.10.
//!
//! The system can also have a source on each panel. A source pushes the
//! outer flow away from the surface. The coupled calculation uses sources
//! to model the displacement effect of the boundary layer.

use crate::geometry::Airfoil;
use crate::linalg;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Solution {
    /// Lift coefficient, from the total circulation (Kutta-Joukowski theorem).
    pub cl: f64,
    /// Lift coefficient, from the sum of the panel pressures. This is a
    /// second, less accurate calculation. Use it to check `cl`.
    pub cl_pressure: f64,
    /// Pressure drag coefficient. The theory value for inviscid flow is 0,
    /// thus this value shows the numerical error.
    pub cd: f64,
    /// Moment coefficient about the quarter-chord point (nose up = positive).
    pub cm: f64,
    /// Panel midpoints.
    pub mid: Vec<[f64; 2]>,
    /// Pressure coefficient at each panel midpoint.
    pub cp: Vec<f64>,
    /// Surface speed at each panel midpoint, divided by the free-stream speed.
    pub speed: Vec<f64>,
    /// Signed surface speed at each panel midpoint. Positive is in the
    /// direction of the point order (clockwise).
    pub vt: Vec<f64>,
    /// Length of each panel.
    pub panel_len: Vec<f64>,
    /// Vortex strength at each panel end point, divided by 2 pi and by the
    /// free-stream speed.
    pub gamma: Vec<f64>,
}

/// Influence of the linear vortex panel `j` (from `a` to `b`) on the
/// velocity at `point`, in the direction with angle `dir`. Returns
/// (normal 1, normal 2, tangential 1, tangential 2): the coefficients for
/// the vortex strengths at the two ends of the panel.
fn vortex_coefficients(point: [f64; 2], dir: f64, a: [f64; 2], b: [f64; 2]) -> [f64; 4] {
    let (px, py) = (b[0] - a[0], b[1] - a[1]);
    let s = px.hypot(py);
    let theta_j = py.atan2(px);
    let dx = point[0] - a[0];
    let dy = point[1] - a[1];
    let (sj, cj) = theta_j.sin_cos();
    let aa = -dx * cj - dy * sj;
    let bb = dx * dx + dy * dy;
    let (c, d) = (dir - theta_j).sin_cos();
    let e = dx * sj - dy * cj;
    let f = (1.0 + s * (s + 2.0 * aa) / bb).ln();
    let g = (e * s).atan2(bb + aa * s);
    let (s2, c2) = (dir - 2.0 * theta_j).sin_cos();
    let p = dx * s2 + dy * c2;
    let q = dx * c2 - dy * s2;
    let cn2 = d + 0.5 * q * f / s - (aa * c + d * e) * g / s;
    let cn1 = 0.5 * d * f + c * g - cn2;
    let ct2 = c + 0.5 * p * f / s + (aa * d - c * e) * g / s;
    let ct1 = 0.5 * c * f - d * g - ct2;
    [cn1, cn2, ct1, ct2]
}

/// Velocity at `point` caused by a source panel from `a` to `b` with a
/// strength of 1 for each unit of length. Returns (x part, y part).
fn source_velocity(point: [f64; 2], a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let (px, py) = (b[0] - a[0], b[1] - a[1]);
    let s = px.hypot(py);
    let (tx, ty) = (px / s, py / s);
    let (nx, ny) = (-ty, tx);
    let (dx, dy) = (point[0] - a[0], point[1] - a[1]);
    // Panel coordinates: x1 along the panel, y1 normal to the panel.
    let (x1, y1) = (dx * tx + dy * ty, dx * nx + dy * ny);
    let r1 = x1 * x1 + y1 * y1;
    let r2 = (x1 - s) * (x1 - s) + y1 * y1;
    let ut = (r1 / r2).ln() / (4.0 * PI);
    let un = (y1.atan2(x1 - s) - y1.atan2(x1)) / (2.0 * PI);
    [ut * tx + un * nx, ut * ty + un * ny]
}

/// The panel equations for one airfoil shape. The matrix is factored one
/// time, and then each solution (for an angle of attack and a set of
/// sources) is fast.
pub struct PanelSystem {
    points: Vec<[f64; 2]>,
    mid: Vec<[f64; 2]>,
    len: Vec<f64>,
    theta: Vec<f64>,
    /// Tangential influence of the vortex strengths.
    at: Vec<Vec<f64>>,
    lu: linalg::Lu,
    /// Normal and tangential velocity at midpoint i from a unit source on
    /// panel j. The panel's effect on itself is not included.
    source_n: Vec<Vec<f64>>,
    source_t: Vec<Vec<f64>>,
    chord: f64,
    x_ref: f64,
}

impl PanelSystem {
    pub fn new(foil: &Airfoil) -> Result<Self, String> {
        let pts = &foil.points;
        let m = pts.len().saturating_sub(1); // number of panels
        if m < 10 {
            return Err("the airfoil must have 10 panels or more".into());
        }
        if pts.iter().any(|p| !p[0].is_finite() || !p[1].is_finite()) {
            return Err("the airfoil has a coordinate that is not a finite number".into());
        }
        let (first, last) = (pts[0], pts[m]);
        if (first[0] - last[0]).hypot(first[1] - last[1]) > 1e-9 {
            return Err("the first point and the last point must be the same point".into());
        }
        // Shoelace formula: a negative area means clockwise point order.
        let area: f64 = (0..m)
            .map(|j| pts[j][0] * pts[j + 1][1] - pts[j + 1][0] * pts[j][1])
            .sum::<f64>()
            * 0.5;
        if area >= 0.0 {
            return Err("the points must be in clockwise order \
                        (trailing edge, lower surface, leading edge, upper surface)"
                .into());
        }

        let mut mid = Vec::with_capacity(m);
        let mut len = Vec::with_capacity(m);
        let mut theta = Vec::with_capacity(m);
        for j in 0..m {
            let (dx, dy) = (pts[j + 1][0] - pts[j][0], pts[j + 1][1] - pts[j][1]);
            let s = dx.hypot(dy);
            if s < 1e-12 {
                return Err(format!("panel {j} has zero length"));
            }
            mid.push([
                0.5 * (pts[j][0] + pts[j + 1][0]),
                0.5 * (pts[j][1] + pts[j + 1][1]),
            ]);
            len.push(s);
            theta.push(dy.atan2(dx));
        }

        // Influence coefficients: normal (an) and tangential (at) velocity at
        // midpoint i caused by the vortex strength at node j.
        let mut an = vec![vec![0.0; m + 1]; m + 1];
        let mut at = vec![vec![0.0; m + 1]; m];
        let mut source_n = vec![vec![0.0; m]; m];
        let mut source_t = vec![vec![0.0; m]; m];
        for i in 0..m {
            let (si, ci) = theta[i].sin_cos();
            for j in 0..m {
                let [cn1, cn2, ct1, ct2] = if i == j {
                    [-1.0, 1.0, 0.5 * PI, 0.5 * PI]
                } else {
                    vortex_coefficients(mid[i], theta[i], pts[j], pts[j + 1])
                };
                an[i][j] += cn1;
                an[i][j + 1] += cn2;
                at[i][j] += ct1;
                at[i][j + 1] += ct2;
                if i != j {
                    let [vx, vy] = source_velocity(mid[i], pts[j], pts[j + 1]);
                    source_n[i][j] = -vx * si + vy * ci;
                    source_t[i][j] = vx * ci + vy * si;
                }
            }
        }
        // Kutta condition: the vortex strengths at the trailing edge cancel.
        an[m][0] = 1.0;
        an[m][m] = 1.0;
        let lu = linalg::Lu::new(an).ok_or("the panel system has no unique solution")?;

        Ok(Self {
            points: pts.clone(),
            mid,
            len,
            theta,
            at,
            lu,
            source_n,
            source_t,
            chord: foil.chord(),
            x_ref: foil.x_min() + 0.25 * foil.chord(),
        })
    }

    /// Number of panels.
    pub fn panel_count(&self) -> usize {
        self.len.len()
    }

    /// Solves the flow for an angle of attack. `sources` is the source
    /// strength on each panel: the speed of the flow that comes out of the
    /// surface, divided by the free-stream speed. `None` means no sources.
    pub fn solve(&self, alpha_deg: f64, sources: Option<&[f64]>) -> Result<Solution, String> {
        let m = self.panel_count();
        if !alpha_deg.is_finite() {
            return Err("the angle of attack must be a finite number".into());
        }
        if let Some(sigma) = sources {
            if sigma.len() != m || sigma.iter().any(|v| !v.is_finite()) {
                return Err("the source list must have one finite value for each panel".into());
            }
        }
        let alpha = alpha_deg.to_radians();
        let dot = |row: &[f64], values: &[f64]| -> f64 {
            row.iter().zip(values).map(|(a, b)| a * b).sum()
        };

        // Flow tangency on the inner side of each panel. A source panel
        // makes a jump in normal velocity: +sigma/2 outside, -sigma/2 inside.
        // Thus the flow out of the surface is sigma and the inner flow is zero.
        // The test `flow_out_of_the_surface_is_equal_to_the_source_strength`
        // checks the signs and the scale of these terms.
        let mut rhs: Vec<f64> = (0..m)
            .map(|i| {
                let base = (self.theta[i] - alpha).sin();
                match sources {
                    Some(sigma) => base + 0.5 * sigma[i] - dot(&self.source_n[i], sigma),
                    None => base,
                }
            })
            .collect();
        rhs.push(0.0);
        let gamma = self.lu.solve(&rhs);

        let (mut fx, mut fy, mut mz) = (0.0, 0.0, 0.0);
        let mut cp = Vec::with_capacity(m);
        let mut speed = Vec::with_capacity(m);
        let mut vt = Vec::with_capacity(m);
        for i in 0..m {
            let mut v = (self.theta[i] - alpha).cos() + dot(&self.at[i], &gamma);
            if let Some(sigma) = sources {
                v += dot(&self.source_t[i], sigma);
            }
            let cpi = 1.0 - v * v;
            // Pressure pushes on the surface, opposite to the outward normal.
            let (nx, ny) = (-self.theta[i].sin(), self.theta[i].cos());
            let (dfx, dfy) = (-cpi * nx * self.len[i], -cpi * ny * self.len[i]);
            fx += dfx;
            fy += dfy;
            mz += (self.mid[i][0] - self.x_ref) * dfy - self.mid[i][1] * dfx;
            cp.push(cpi);
            speed.push(v.abs());
            vt.push(v);
        }
        let (sa, ca) = alpha.sin_cos();
        // Lift from the total circulation (Kutta-Joukowski theorem). This is
        // more accurate than the sum of the panel pressures.
        let circulation: f64 = (0..m)
            .map(|j| PI * (gamma[j] + gamma[j + 1]) * self.len[j])
            .sum();
        let chord = self.chord;
        Ok(Solution {
            cl: 2.0 * circulation / chord,
            cl_pressure: (fy * ca - fx * sa) / chord,
            cd: (fx * ca + fy * sa) / chord,
            cm: -mz / (chord * chord),
            mid: self.mid.clone(),
            cp,
            speed,
            vt,
            panel_len: self.len.clone(),
            gamma,
        })
    }

    /// Velocity (x part, y part) at a point that is not on the surface, for
    /// a solution of this system. Used to check the solution.
    pub fn velocity_at(
        &self,
        point: [f64; 2],
        alpha_deg: f64,
        solution: &Solution,
        sources: Option<&[f64]>,
    ) -> [f64; 2] {
        let alpha = alpha_deg.to_radians();
        let (mut vx, mut vy) = (alpha.cos(), alpha.sin());
        for j in 0..self.panel_count() {
            let (a, b) = (self.points[j], self.points[j + 1]);
            // Direction angle 0 gives: tangential = x part, normal = y part.
            let [cn1, cn2, ct1, ct2] = vortex_coefficients(point, 0.0, a, b);
            vx += ct1 * solution.gamma[j] + ct2 * solution.gamma[j + 1];
            vy += cn1 * solution.gamma[j] + cn2 * solution.gamma[j + 1];
            if let Some(sigma) = sources {
                let [sx, sy] = source_velocity(point, a, b);
                vx += sx * sigma[j];
                vy += sy * sigma[j];
            }
        }
        [vx, vy]
    }
}

/// Solves the friction-free flow for one airfoil and one angle of attack.
pub fn solve(foil: &Airfoil, alpha_deg: f64) -> Result<Solution, String> {
    PanelSystem::new(foil)?.solve(alpha_deg, None)
}
