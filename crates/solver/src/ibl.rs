//! Two-equation integral boundary layer equations (Drela), written as
//! residuals between two stations, and a solver that marches along one
//! surface for a given edge speed ("direct" mode).
//!
//! Unknowns at each station: momentum thickness `theta`, displacement
//! thickness `dstar`, and a third variable `s`: the amplification factor
//! in laminar flow, or the square root of the shear stress coefficient in
//! turbulent flow and in the wake.
//!
//! Units: lengths are divided by the chord, speeds by the free-stream
//! speed. `re` is the Reynolds number for the chord.

use crate::closure::{self as cl, Flow};

/// Shape factor limits for the direct march. Above the limit, the direct
/// march is not valid (separation).
pub const H_LIMIT_LAMINAR: f64 = 3.8;
pub const H_LIMIT_TURBULENT: f64 = 2.5;
/// A trip strip starts turbulent flow only if the momentum-thickness
/// Reynolds number is at this value or more. A thinner boundary layer
/// cannot keep turbulent flow (assumption, see GOALS.md).
pub const TRIP_MIN_RE_THETA: f64 = 200.0;

#[derive(Debug, Clone, Copy)]
pub struct State {
    /// Distance along the surface.
    pub xi: f64,
    /// Edge speed.
    pub ue: f64,
    pub theta: f64,
    pub dstar: f64,
    /// Amplification factor (laminar) or square root of the shear stress
    /// coefficient (turbulent, wake).
    pub s: f64,
    pub flow: Flow,
}

/// Quantities that the closure relations give for one state.
#[derive(Debug, Clone, Copy)]
pub struct Derived {
    pub h: f64,
    pub re_theta: f64,
    pub hs: f64,
    pub cf: f64,
    /// Dissipation, 2 cD / H*.
    pub di: f64,
    pub us: f64,
    pub ctau_eq_sqrt: f64,
    pub delta: f64,
}

impl State {
    pub fn h(&self) -> f64 {
        let minimum = if self.flow == Flow::Wake {
            1.00005
        } else {
            1.05
        };
        (self.dstar / self.theta).max(minimum)
    }

    pub fn derived(&self, re: f64) -> Derived {
        let h = self.h();
        let re_theta = (re * self.ue * self.theta).max(1e-6);
        match self.flow {
            Flow::Laminar => Derived {
                h,
                re_theta,
                hs: cl::hs_laminar(h),
                cf: cl::cf_laminar(h, re_theta),
                di: cl::di_laminar(h, re_theta),
                us: 0.0,
                ctau_eq_sqrt: 0.0,
                delta: cl::delta(h, self.theta, self.dstar),
            },
            Flow::Turbulent | Flow::Wake => {
                let hs = cl::hs_turbulent(h, re_theta);
                let us = cl::slip_velocity(hs, h, self.flow);
                let wake = self.flow == Flow::Wake;
                let cf = if wake {
                    0.0
                } else {
                    cl::cf_turbulent(h, re_theta)
                };
                let ctau = self.s * self.s;
                let di = if wake {
                    cl::di_wake(h, hs, us, ctau, re_theta)
                } else {
                    cl::di_turbulent(h, hs, us, cf, ctau, re_theta)
                };
                Derived {
                    h,
                    re_theta,
                    hs,
                    cf,
                    di,
                    us,
                    ctau_eq_sqrt: cl::ctau_eq_sqrt(h, hs, us, re_theta, self.flow),
                    delta: cl::delta(h, self.theta, self.dstar),
                }
            }
        }
    }
}

/// Residuals of the three equations between station 1 and station 2.
/// The two states must have the same flow type.
pub fn residuals(a: &State, b: &State, re: f64) -> [f64; 3] {
    let (da, db) = (a.derived(re), b.derived(re));
    let x_log = (b.xi / a.xi).ln();
    let u_log = (b.ue / a.ue).ln();
    let dx = b.xi - a.xi;
    let h_avg = 0.5 * (da.h + db.h);

    // Upwind weight: more weight on station 2 when the shape factor
    // changes quickly.
    let ratio = (db.h - 1.0) / (da.h - 1.0);
    let c_up = if b.flow == Flow::Wake { 1.0 } else { 5.0 };
    let up = 1.0 - 0.5 * (-(ratio.ln()).powi(2) * c_up / (db.h * db.h)).exp();
    let blend = |v1: f64, v2: f64| (1.0 - up) * v1 + up * v2;

    // Momentum equation.
    let mid = State {
        xi: 0.5 * (a.xi + b.xi),
        ue: 0.5 * (a.ue + b.ue),
        theta: 0.5 * (a.theta + b.theta),
        dstar: 0.5 * (a.dstar + b.dstar),
        s: 0.5 * (a.s + b.s),
        flow: b.flow,
    };
    let cf_mid = mid.derived(re).cf;
    let cfx = 0.5 * cf_mid * mid.xi / mid.theta
        + 0.25 * (da.cf * a.xi / a.theta + db.cf * b.xi / b.theta);
    let momentum = (b.theta / a.theta).ln() + (2.0 + h_avg) * u_log - 0.5 * cfx * x_log;

    // Shape (kinetic energy) equation.
    let cfx_up = blend(da.cf * a.xi / a.theta, db.cf * b.xi / b.theta);
    let dix_up = blend(da.di * a.xi / a.theta, db.di * b.xi / b.theta);
    let shape = (db.hs / da.hs).ln() + (1.0 - h_avg) * u_log + x_log * (0.5 * cfx_up - dix_up);

    // Third equation.
    let third = match b.flow {
        Flow::Laminar => {
            let rate = 0.5
                * (cl::amplification_rate(da.h, a.theta, da.re_theta)
                    + cl::amplification_rate(db.h, b.theta, db.re_theta));
            b.s - a.s - rate * dx
        }
        Flow::Turbulent | Flow::Wake => {
            let eta_d = if b.flow == Flow::Wake {
                cl::ETA_D_WAKE
            } else {
                1.0
            };
            let s_up = blend(a.s, b.s);
            let eq_up = blend(da.ctau_eq_sqrt, db.ctau_eq_sqrt);
            let us_up = blend(da.us, db.us);
            let delta_avg = 0.5 * (da.delta + db.delta);
            let dstar_avg = 0.5 * (a.dstar + b.dstar);
            let h_up = blend(da.h, db.h);
            let cf_up = blend(da.cf, db.cf);
            let re_theta_avg = 0.5 * (da.re_theta + db.re_theta);
            let gc = if b.flow == Flow::Wake { 0.0 } else { cl::G_C };
            let hkc = (h_up - 1.0 - gc / re_theta_avg).max(0.01);
            let uq =
                (0.5 * cf_up - (hkc / (cl::G_A * eta_d * h_up)).powi(2)) / (cl::G_B * dstar_avg);
            2.0 * delta_avg * (b.s / a.s).ln()
                - cl::K_LAG / (cl::G_B * (1.0 + us_up)) * (eq_up - eta_d * s_up) * dx
                - 2.0 * delta_avg * (uq * dx - u_log)
        }
    };
    [momentum, shape, third]
}

/// Residuals at the first station after a stagnation point (similarity
/// station: ue is proportional to xi, and theta and H are constant).
pub fn similarity_residuals(b: &State, re: f64) -> [f64; 2] {
    let d = b.derived(re);
    let cfx = d.cf * b.xi / b.theta;
    let dix = d.di * b.xi / b.theta;
    [(2.0 + d.h) - 0.5 * cfx, (1.0 - d.h) + 0.5 * cfx - dix]
}

/// Solves a small system by Newton's method with a numerical Jacobian.
/// `unknowns` holds the start values and receives the result.
fn newton<const N: usize>(
    unknowns: &mut [f64; N],
    residual: impl Fn(&[f64; N]) -> [f64; N],
) -> bool {
    for _ in 0..40 {
        let r = residual(unknowns);
        if r.iter().any(|v| !v.is_finite()) {
            return false;
        }
        if r.iter().all(|v| v.abs() < 1e-10) {
            return true;
        }
        let mut jac = vec![vec![0.0; N]; N];
        for k in 0..N {
            let mut shifted = *unknowns;
            let step = 1e-6 * shifted[k].abs().max(1e-9);
            shifted[k] += step;
            let rs = residual(&shifted);
            for (row, jac_row) in jac.iter_mut().enumerate() {
                jac_row[k] = (rs[row] - r[row]) / step;
            }
        }
        let Some(delta) = crate::linalg::solve(jac, r.iter().map(|v| -v).collect()) else {
            return false;
        };
        // Limit the step so that theta and dstar stay positive.
        let mut scale: f64 = 1.0;
        for k in 0..N.min(2) {
            let relative = delta[k] / unknowns[k];
            if relative.abs() > 0.3 {
                scale = scale.min(0.3 / relative.abs());
            }
        }
        for k in 0..N {
            unknowns[k] += scale * delta[k];
        }
    }
    false
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransitionCause {
    /// The amplification factor reached the critical value.
    Natural,
    /// The laminar shape factor reached the limit (laminar separation).
    LaminarSeparation,
    /// The user set the transition position.
    Forced,
}

#[derive(Debug, Clone)]
pub struct MarchResult {
    /// Stations with a result. The list stops at turbulent separation.
    pub states: Vec<State>,
    /// Index of the first turbulent station, and the cause.
    pub transition: Option<(usize, TransitionCause)>,
    /// Index of the station where the direct march stopped (separation).
    pub separation: Option<usize>,
}

/// Marches the boundary layer along one surface for a given edge speed.
///
/// - `xi[0]` must be more than 0. If `stagnation` is true, the first
///   station is a similarity station after a stagnation point. If not, it
///   is the first station of a flat plate with a sharp leading edge.
/// - `forced_transition`: if set, turbulent flow starts at this distance.
pub fn march(
    xi: &[f64],
    ue: &[f64],
    re: f64,
    stagnation: bool,
    forced_transition: Option<f64>,
    n_crit: f64,
) -> Result<MarchResult, String> {
    let n = xi.len();
    if n < 3 || ue.len() != n {
        return Err("the surface must have 3 stations or more".into());
    }
    if !(re.is_finite() && re > 0.0) {
        return Err("the Reynolds number must be more than zero".into());
    }
    if xi[0] <= 0.0 || xi.windows(2).any(|w| w[1] <= w[0]) || xi.iter().any(|v| !v.is_finite()) {
        return Err("the surface distance must be more than zero and must increase".into());
    }
    if ue.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err("the edge speed must be more than zero".into());
    }

    // First station.
    let theta0 = if stagnation {
        (0.075 * xi[0] / (re * ue[0])).sqrt()
    } else {
        // Blasius thickness at the first station of a flat plate.
        0.664 * (xi[0] / (re * ue[0])).sqrt()
    };
    let mut first = State {
        xi: xi[0],
        ue: ue[0],
        theta: theta0,
        dstar: if stagnation {
            2.2 * theta0
        } else {
            2.59 * theta0
        },
        s: 0.0,
        flow: Flow::Laminar,
    };
    if stagnation {
        let mut unknowns = [first.theta, first.dstar];
        let ok = newton(&mut unknowns, |v| {
            similarity_residuals(
                &State {
                    theta: v[0],
                    dstar: v[1],
                    ..first
                },
                re,
            )
        });
        if !ok {
            return Err("the boundary layer start at the stagnation point has no solution".into());
        }
        first.theta = unknowns[0];
        first.dstar = unknowns[1];
    }

    let mut states = vec![first];
    let mut transition = None;
    let mut separation = None;
    for i in 1..n {
        let mut start = *states.last().unwrap();
        if start.flow == Flow::Laminar {
            let thick_enough = re * start.ue * start.theta >= TRIP_MIN_RE_THETA;
            let cause = if thick_enough && forced_transition.is_some_and(|limit| start.xi >= limit)
            {
                Some(TransitionCause::Forced)
            } else if start.s >= n_crit {
                Some(TransitionCause::Natural)
            } else if start.h() >= H_LIMIT_LAMINAR {
                Some(TransitionCause::LaminarSeparation)
            } else {
                None
            };
            if let Some(cause) = cause {
                start = turbulent_start(&start, re);
                transition = Some((i - 1, cause));
            }
        }
        let mut next = step(&start, xi[i], ue[i], re);
        if next.is_none() && start.flow == Flow::Laminar {
            // Laminar separation between two stations: start turbulent flow
            // at the previous station and solve this station again.
            start = turbulent_start(&start, re);
            transition = Some((i - 1, TransitionCause::LaminarSeparation));
            next = step(&start, xi[i], ue[i], re);
        }
        *states.last_mut().unwrap() = start;
        match next {
            Some(state) => states.push(state),
            None => {
                separation = Some(i);
                break;
            }
        }
    }
    Ok(MarchResult {
        states,
        transition,
        separation,
    })
}

/// Changes a laminar state into the first turbulent state at the same
/// station: same momentum thickness, shape factor not above the turbulent
/// limit, shear stress from the equilibrium value.
pub fn turbulent_start(state: &State, re: f64) -> State {
    let mut start = *state;
    start.flow = Flow::Turbulent;
    start.dstar = start.dstar.min(0.99 * H_LIMIT_TURBULENT * start.theta);
    start.s = 0.03;
    let d = start.derived(re);
    start.s = cl::ctau_start_sqrt(d.h, d.ctau_eq_sqrt).max(1e-4);
    start
}

/// Solves the next station. Returns `None` if there is no solution below
/// the shape factor limit (separation).
pub fn step(start: &State, xi: f64, ue: f64, re: f64) -> Option<State> {
    let guess = State { xi, ue, ..*start };
    let laminar = start.flow == Flow::Laminar;
    let mut unknowns = [guess.theta, guess.dstar, guess.s];
    let ok = newton(&mut unknowns, |v| {
        let s = if laminar { v[2] } else { v[2].max(1e-7) };
        residuals(
            start,
            &State {
                theta: v[0],
                dstar: v[1],
                s,
                ..guess
            },
            re,
        )
    });
    let next = State {
        theta: unknowns[0],
        dstar: unknowns[1],
        s: unknowns[2],
        ..guess
    };
    let limit = if laminar {
        H_LIMIT_LAMINAR + 0.5
    } else {
        H_LIMIT_TURBULENT
    };
    (ok && next.theta > 0.0 && next.dstar > 0.0 && next.h() < limit).then_some(next)
}
