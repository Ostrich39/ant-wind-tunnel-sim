//! Boundary layer (friction) calculation along one surface.
//!
//! Input: the flow speed at the edge of the boundary layer, at stations
//! along the surface. Output: the boundary layer thickness, shape, and
//! skin friction at each station.
//!
//! Methods (see REFERENCES.md, section 1):
//! - Laminar flow: Thwaites' method.
//! - Transition: Michel's criterion. Laminar separation also starts
//!   turbulent flow (the short separation bubble is not modeled).
//! - Turbulent flow: Head's entrainment method with the Ludwieg-Tillmann
//!   skin friction law.
//!
//! This is a one-way calculation. The boundary layer does not change the
//! pressure on the surface. Thus the method is valid only for attached flow.
//!
//! Units: lengths are divided by the chord, speeds are divided by the
//! free-stream speed. `re` is the Reynolds number for the chord.

/// Thwaites parameter at laminar separation.
pub const LAMINAR_SEPARATION_LAMBDA: f64 = -0.09;
/// Shape factor at turbulent separation. Published values are 1.8 to 2.4.
pub const TURBULENT_SEPARATION_H: f64 = 2.4;
/// Shape factor at the start of turbulent flow (assumption).
pub const TURBULENT_START_H: f64 = 1.4;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Regime {
    Laminar,
    Turbulent,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransitionCause {
    /// Michel's criterion (natural transition).
    Michel,
    /// The laminar boundary layer separated.
    LaminarSeparation,
    /// The user set the transition position.
    Forced,
}

#[derive(Debug, Clone, Copy)]
pub struct Station {
    /// Distance along the surface from the start.
    pub s: f64,
    /// Speed at the edge of the boundary layer.
    pub ue: f64,
    /// Momentum thickness.
    pub theta: f64,
    /// Shape factor (displacement thickness / momentum thickness).
    pub h: f64,
    /// Skin friction coefficient, referred to the local edge speed.
    pub cf: f64,
    pub regime: Regime,
}

#[derive(Debug, Clone)]
pub struct SurfaceResult {
    /// Stations with a result. This list stops at turbulent separation.
    pub stations: Vec<Station>,
    /// Index of the first turbulent station, and the cause.
    pub transition: Option<(usize, TransitionCause)>,
    /// Index of the station where the turbulent boundary layer separates.
    pub separation: Option<usize>,
}

/// Thwaites shape factor correlation, H(lambda).
pub fn thwaites_h(lambda: f64) -> f64 {
    let l = lambda.clamp(-0.1, 0.1);
    if l >= 0.0 {
        2.61 - 3.75 * l + 5.24 * l * l
    } else {
        2.088 + 0.0731 / (l + 0.14)
    }
}

/// Thwaites shear correlation, l(lambda).
pub fn thwaites_l(lambda: f64) -> f64 {
    let l = lambda.clamp(-0.1, 0.1);
    if l >= 0.0 {
        0.22 + 1.57 * l - 1.8 * l * l
    } else {
        0.22 + 1.402 * l + 0.018 * l / (l + 0.107)
    }
}

/// Michel's criterion: momentum-thickness Reynolds number at transition.
pub fn michel_limit(re_x: f64) -> f64 {
    1.174 * (1.0 + 22400.0 / re_x) * re_x.powf(0.46)
}

/// Head's entrainment shape factor H1 as a function of H.
pub fn head_h1(h: f64) -> f64 {
    let h = h.max(1.1001);
    if h <= 1.6 {
        3.3 + 0.8234 * (h - 1.1).powf(-1.287)
    } else {
        3.3 + 1.5501 * (h - 0.6778).powf(-3.064)
    }
}

/// Inverse of `head_h1`.
pub fn head_h(h1: f64) -> f64 {
    if h1 <= 3.3 {
        3.0
    } else if h1 < 5.3 {
        0.6778 + 1.1536 * (h1 - 3.3).powf(-0.326)
    } else {
        1.1 + 0.86 * (h1 - 3.3).powf(-0.777)
    }
}

/// Ludwieg-Tillmann skin friction law for turbulent flow.
pub fn ludwieg_tillmann_cf(h: f64, re_theta: f64) -> f64 {
    0.246 * 10f64.powf(-0.678 * h) * re_theta.powf(-0.268)
}

/// Integral of ue^5 over one interval, exact for a linear change of ue.
fn integral_ue5(ua: f64, ub: f64, ds: f64) -> f64 {
    if (ub - ua).abs() < 1e-12 * ua.abs().max(ub.abs()).max(1e-300) {
        ds * ua.powi(5)
    } else {
        ds * (ub.powi(6) - ua.powi(6)) / (6.0 * (ub - ua))
    }
}

/// Slope of ue at station i (difference formula for unequal spacing).
fn slope(s: &[f64], ue: &[f64], i: usize) -> f64 {
    let n = s.len();
    if i == 0 {
        (ue[1] - ue[0]) / (s[1] - s[0])
    } else if i == n - 1 {
        (ue[n - 1] - ue[n - 2]) / (s[n - 1] - s[n - 2])
    } else {
        let (h1, h2) = (s[i] - s[i - 1], s[i + 1] - s[i]);
        let (d1, d2) = ((ue[i] - ue[i - 1]) / h1, (ue[i + 1] - ue[i]) / h2);
        (d1 * h2 + d2 * h1) / (h1 + h2)
    }
}

/// Rates of change for Head's method: (d theta / ds, d H1 / ds).
fn head_rates(ue: f64, due: f64, theta: f64, h1: f64, re: f64) -> (f64, f64) {
    let h = head_h(h1);
    let cf = ludwieg_tillmann_cf(h, re * ue * theta);
    let dtheta = 0.5 * cf - (h + 2.0) * theta / ue * due;
    let entrainment = 0.0306 * (h1 - 3.0).max(1e-6).powf(-0.6169);
    let dh1 = (ue * entrainment - h1 * (ue * dtheta + theta * due)) / (ue * theta);
    (dtheta, dh1)
}

/// Calculates the boundary layer along one surface.
///
/// - `s`: distance along the surface, increasing, `s[0]` at the start.
/// - `ue`: edge speed at each station. `ue[0]` is 0 at a stagnation point,
///   or more than 0 at the sharp leading edge of a flat plate.
/// - `forced_transition`: if set, turbulent flow starts at this distance
///   (or before it, if natural transition occurs first).
pub fn march(
    s: &[f64],
    ue: &[f64],
    re: f64,
    forced_transition: Option<f64>,
) -> Result<SurfaceResult, String> {
    let n = s.len();
    if n < 3 || ue.len() != n {
        return Err("the surface must have 3 stations or more".into());
    }
    if !(re.is_finite() && re > 0.0) {
        return Err("the Reynolds number must be more than zero".into());
    }
    if s.windows(2).any(|w| w[1] <= w[0]) || s.iter().any(|v| !v.is_finite()) {
        return Err("the surface distance must increase at each station".into());
    }
    if ue.iter().any(|v| !v.is_finite() || *v < 0.0) || ue[1..].iter().any(|v| *v <= 0.0) {
        return Err("the edge speed must be more than zero after the first station".into());
    }

    let mut stations: Vec<Station> = Vec::with_capacity(n);
    let mut transition = None;

    // Laminar part: Thwaites' method.
    let mut integral = 0.0;
    for i in 0..n {
        let due = slope(s, ue, i);
        let theta = if i == 0 {
            if ue[0] == 0.0 {
                // Stagnation point: theta^2 = 0.075 / (re * due/ds).
                if due <= 0.0 {
                    return Err("the edge speed must increase from the stagnation point".into());
                }
                (0.075 / (re * due)).sqrt()
            } else {
                0.0
            }
        } else {
            integral += integral_ue5(ue[i - 1], ue[i], s[i] - s[i - 1]);
            (0.45 * integral / (re * ue[i].powi(6))).sqrt()
        };
        let lambda = theta * theta * re * due;
        let cf = if theta > 0.0 && ue[i] > 0.0 {
            2.0 * thwaites_l(lambda) / (re * ue[i] * theta)
        } else {
            f64::INFINITY
        };

        if i > 0 {
            let re_x = re * ue[i] * (s[i] - s[0]);
            let re_theta = re * ue[i] * theta;
            let cause = if forced_transition.is_some_and(|limit| s[i] >= limit) {
                Some(TransitionCause::Forced)
            } else if re_theta > michel_limit(re_x) {
                Some(TransitionCause::Michel)
            } else if lambda < LAMINAR_SEPARATION_LAMBDA {
                Some(TransitionCause::LaminarSeparation)
            } else {
                None
            };
            if let Some(cause) = cause {
                transition = Some((i, cause));
                let h = TURBULENT_START_H;
                stations.push(Station {
                    s: s[i],
                    ue: ue[i],
                    theta,
                    h,
                    cf: ludwieg_tillmann_cf(h, re_theta),
                    regime: Regime::Turbulent,
                });
                break;
            }
        }
        stations.push(Station {
            s: s[i],
            ue: ue[i],
            theta,
            h: thwaites_h(lambda),
            cf,
            regime: Regime::Laminar,
        });
    }

    // Turbulent part: Head's method, fourth-order Runge-Kutta steps.
    let mut separation = None;
    if let Some((start, _)) = transition {
        let mut theta = stations[start].theta;
        let mut h1 = head_h1(TURBULENT_START_H);
        const SUBSTEPS: usize = 8;
        'stations: for i in start + 1..n {
            let ds = (s[i] - s[i - 1]) / SUBSTEPS as f64;
            let due = (ue[i] - ue[i - 1]) / (s[i] - s[i - 1]);
            for k in 0..SUBSTEPS {
                let u = |frac: f64| ue[i - 1] + due * ds * (k as f64 + frac);
                let (a1, b1) = head_rates(u(0.0), due, theta, h1, re);
                let (a2, b2) =
                    head_rates(u(0.5), due, theta + 0.5 * ds * a1, h1 + 0.5 * ds * b1, re);
                let (a3, b3) =
                    head_rates(u(0.5), due, theta + 0.5 * ds * a2, h1 + 0.5 * ds * b2, re);
                let (a4, b4) = head_rates(u(1.0), due, theta + ds * a3, h1 + ds * b3, re);
                theta += ds * (a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0;
                h1 += ds * (b1 + 2.0 * b2 + 2.0 * b3 + b4) / 6.0;
                if !(theta.is_finite() && h1.is_finite()) || theta <= 0.0 {
                    separation = Some(i);
                    break 'stations;
                }
            }
            let h = head_h(h1);
            if h >= TURBULENT_SEPARATION_H {
                separation = Some(i);
                break;
            }
            stations.push(Station {
                s: s[i],
                ue: ue[i],
                theta,
                h,
                cf: ludwieg_tillmann_cf(h, re * ue[i] * theta),
                regime: Regime::Turbulent,
            });
        }
    }

    Ok(SurfaceResult {
        stations,
        transition,
        separation,
    })
}
