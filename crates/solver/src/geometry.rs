//! Airfoil shapes.
//!
//! Point order: trailing edge -> lower surface -> leading edge ->
//! upper surface -> trailing edge (clockwise). The first point and the
//! last point are the same point.

use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Airfoil {
    pub name: String,
    pub points: Vec<[f64; 2]>,
}

impl Airfoil {
    /// Chord length: distance from the smallest x to the largest x.
    pub fn chord(&self) -> f64 {
        self.x_max() - self.x_min()
    }
    pub fn x_min(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p[0])
            .fold(f64::INFINITY, f64::min)
    }
    pub fn x_max(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p[0])
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

/// NACA 4-digit airfoil, for example "0012" or "2412".
/// `n_panels` must be an even number. The trailing edge is closed.
pub fn naca4(code: &str, n_panels: usize) -> Result<Airfoil, String> {
    naca4_with_spacing(code, n_panels, 1.0)
}

/// NACA 4-digit airfoil with panels for the two-way (friction) calculation.
/// The panels are small at the leading edge and medium at the trailing
/// edge (approximately 1% of the chord for 120 panels). Very small panels
/// at the trailing edge make the two-way calculation unstable.
pub fn naca4_for_friction(code: &str, n_panels: usize) -> Result<Airfoil, String> {
    naca4_with_spacing(code, n_panels, 0.85)
}

/// `spacing_end` sets the panel distribution: 1.0 gives cosine spacing
/// (small panels at the two ends); a smaller value gives larger panels at
/// the trailing edge.
fn naca4_with_spacing(code: &str, n_panels: usize, spacing_end: f64) -> Result<Airfoil, String> {
    let digits: Vec<u32> = code.chars().filter_map(|c| c.to_digit(10)).collect();
    if code.len() != 4 || digits.len() != 4 {
        return Err(format!("'{code}' is not a 4-digit NACA code"));
    }
    if n_panels < 10 || !n_panels.is_multiple_of(2) {
        return Err("n_panels must be an even number, 10 or more".into());
    }
    let m = digits[0] as f64 / 100.0; // maximum camber
    let p = digits[1] as f64 / 10.0; // position of maximum camber
    let t = (digits[2] * 10 + digits[3]) as f64 / 100.0; // thickness
    if t == 0.0 {
        return Err("thickness must not be zero".into());
    }
    if (m > 0.0) != (p > 0.0) {
        return Err("camber and camber position must both be zero or both be non-zero".into());
    }

    let n_side = n_panels / 2;
    let surface = |k: usize| -> ([f64; 2], [f64; 2]) {
        // Cosine spacing puts more points near the leading and trailing edges.
        let beta_end = PI * spacing_end;
        let x = (1.0 - (beta_end * k as f64 / n_side as f64).cos()) / (1.0 - beta_end.cos());
        let yt = 5.0
            * t
            * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3)
                - 0.1036 * x.powi(4));
        let (yc, dyc) = if m == 0.0 {
            (0.0, 0.0)
        } else if x < p {
            (
                m / (p * p) * (2.0 * p * x - x * x),
                2.0 * m / (p * p) * (p - x),
            )
        } else {
            let q = (1.0 - p) * (1.0 - p);
            (
                m / q * (1.0 - 2.0 * p + 2.0 * p * x - x * x),
                2.0 * m / q * (p - x),
            )
        };
        let th = dyc.atan();
        (
            [x - yt * th.sin(), yc + yt * th.cos()], // upper
            [x + yt * th.sin(), yc - yt * th.cos()], // lower
        )
    };

    let mut points = Vec::with_capacity(n_panels + 1);
    for k in (1..=n_side).rev() {
        points.push(surface(k).1);
    }
    points.push(surface(0).0);
    for k in 1..=n_side {
        points.push(surface(k).0);
    }
    // Make the first and last points exactly equal.
    let te = [1.0, 0.0];
    points[0] = te;
    *points.last_mut().unwrap() = te;
    Ok(Airfoil {
        name: format!("NACA {code}"),
        points,
    })
}

/// Joukowski airfoil. This shape has an exact mathematical solution,
/// thus it is used to check the solver.
///
/// `eps` sets the thickness and `delta` sets the camber.
/// Returns the airfoil and the exact lift coefficient function data:
/// (airfoil, circle radius `a`, zero-lift angle `beta` in radians).
pub fn joukowski(eps: f64, delta: f64, n_panels: usize) -> (Airfoil, f64, f64) {
    let a = ((1.0 + eps).powi(2) + delta * delta).sqrt();
    let beta = (delta / a).asin();
    let mut points = Vec::with_capacity(n_panels + 1);
    for k in 0..=n_panels {
        // Clockwise around the circle. Start and end at the trailing edge.
        let phi = -beta - 2.0 * PI * k as f64 / n_panels as f64;
        let (zx, zy) = (-eps + a * phi.cos(), delta + a * phi.sin());
        let r2 = zx * zx + zy * zy;
        points.push([zx + zx / r2, zy - zy / r2]); // z = zeta + 1/zeta
    }
    points[0] = [2.0, 0.0];
    points[n_panels] = [2.0, 0.0];
    (
        Airfoil {
            name: format!("Joukowski eps={eps} delta={delta}"),
            points,
        },
        a,
        beta,
    )
}
