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

use crate::geometry::Airfoil;
use crate::linalg;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Solution {
    /// Lift coefficient.
    pub cl: f64,
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
}

pub fn solve(foil: &Airfoil, alpha_deg: f64) -> Result<Solution, String> {
    let pts = &foil.points;
    let m = pts.len().saturating_sub(1); // number of panels
    if m < 10 {
        return Err("the airfoil must have 10 panels or more".into());
    }
    let alpha = alpha_deg.to_radians();

    let mut mid = Vec::with_capacity(m);
    let mut len = Vec::with_capacity(m);
    let mut theta = Vec::with_capacity(m);
    for j in 0..m {
        let (dx, dy) = (pts[j + 1][0] - pts[j][0], pts[j + 1][1] - pts[j][1]);
        let s = dx.hypot(dy);
        if s < 1e-12 {
            return Err(format!("panel {j} has zero length"));
        }
        mid.push([0.5 * (pts[j][0] + pts[j + 1][0]), 0.5 * (pts[j][1] + pts[j + 1][1])]);
        len.push(s);
        theta.push(dy.atan2(dx));
    }

    // Influence coefficients: normal (an) and tangential (at) velocity at
    // midpoint i caused by the vortex strength at node j.
    let mut an = vec![vec![0.0; m + 1]; m + 1];
    let mut at = vec![vec![0.0; m + 1]; m];
    for i in 0..m {
        for j in 0..m {
            let (cn1, cn2, ct1, ct2) = if i == j {
                (-1.0, 1.0, 0.5 * PI, 0.5 * PI)
            } else {
                let dx = mid[i][0] - pts[j][0];
                let dy = mid[i][1] - pts[j][1];
                let (sj, cj) = theta[j].sin_cos();
                let s = len[j];
                let a = -dx * cj - dy * sj;
                let b = dx * dx + dy * dy;
                let (c, d) = (theta[i] - theta[j]).sin_cos();
                let e = dx * sj - dy * cj;
                let f = (1.0 + s * (s + 2.0 * a) / b).ln();
                let g = (e * s).atan2(b + a * s);
                let (s2, c2) = (theta[i] - 2.0 * theta[j]).sin_cos();
                let p = dx * s2 + dy * c2;
                let q = dx * c2 - dy * s2;
                let cn2 = d + 0.5 * q * f / s - (a * c + d * e) * g / s;
                let cn1 = 0.5 * d * f + c * g - cn2;
                let ct2 = c + 0.5 * p * f / s + (a * d - c * e) * g / s;
                let ct1 = 0.5 * c * f - d * g - ct2;
                (cn1, cn2, ct1, ct2)
            };
            an[i][j] += cn1;
            an[i][j + 1] += cn2;
            at[i][j] += ct1;
            at[i][j + 1] += ct2;
        }
    }
    // Kutta condition: the vortex strengths at the trailing edge cancel.
    an[m][0] = 1.0;
    an[m][m] = 1.0;

    let mut rhs: Vec<f64> = theta.iter().map(|th| (th - alpha).sin()).collect();
    rhs.push(0.0);
    let gamma = linalg::solve(an, rhs).ok_or("the panel system has no unique solution")?;

    let chord = foil.chord();
    let x_ref = foil.x_min() + 0.25 * chord;
    let (mut fx, mut fy, mut mz) = (0.0, 0.0, 0.0);
    let mut cp = Vec::with_capacity(m);
    let mut speed = Vec::with_capacity(m);
    for i in 0..m {
        let v = (theta[i] - alpha).cos()
            + at[i].iter().zip(&gamma).map(|(a, g)| a * g).sum::<f64>();
        let cpi = 1.0 - v * v;
        // Pressure pushes on the surface, opposite to the outward normal.
        let (nx, ny) = (-theta[i].sin(), theta[i].cos());
        let (dfx, dfy) = (-cpi * nx * len[i], -cpi * ny * len[i]);
        fx += dfx;
        fy += dfy;
        mz += (mid[i][0] - x_ref) * dfy - mid[i][1] * dfx;
        cp.push(cpi);
        speed.push(v.abs());
    }
    let (sa, ca) = alpha.sin_cos();
    // Lift from the total circulation (Kutta-Joukowski theorem). This is
    // more accurate than the sum of the panel pressures.
    let circulation: f64 = (0..m).map(|j| PI * (gamma[j] + gamma[j + 1]) * len[j]).sum();
    Ok(Solution {
        cl: 2.0 * circulation / chord,
        cd: (fx * ca + fy * sa) / chord,
        cm: -mz / (chord * chord),
        mid,
        cp,
        speed,
    })
}
