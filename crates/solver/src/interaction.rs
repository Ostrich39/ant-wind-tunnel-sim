//! Two-way calculation by simultaneous solution (Drela's method).
//!
//! The outer flow and the boundary layer are solved together by Newton's
//! method. The unknowns at each station are the momentum thickness, the
//! mass defect (edge speed * displacement thickness), and the third
//! boundary layer variable. The edge speed is a linear function of the mass
//! defect at all stations, through the source influence table.
//!
//! Stations: the panel midpoints of the upper surface and of the lower
//! surface (from the stagnation point to the trailing edge), then the wake.
//!
//! Differences from XFOIL (see GOALS.md): transition is at a station, not
//! at a position between two stations. The Jacobian matrix is calculated
//! by finite differences and solved as a full matrix.

// Index loops are clearer than iterators for these matrix calculations.
#![allow(clippy::needless_range_loop)]

use crate::closure::{self as cl, Flow};
use crate::geometry::Airfoil;
use crate::ibl::{self, State};
use crate::linalg;
use crate::panel::PanelSystem;
use crate::wake::{SourceInfluence, Wake};

pub const WAKE_PANELS: usize = 20;
/// Wake length, divided by the chord.
pub const WAKE_LENGTH: f64 = 1.0;
pub const MAX_ITERATIONS: usize = 40;
/// The iteration stops if the residual does not decrease this many times.
const MAX_STALLED_STEPS: usize = 3;
/// The iteration stops when the largest relative change is below this.
pub const TOLERANCE: f64 = 1.0e-4;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Surface {
    Upper,
    Lower,
    Wake,
}

#[derive(Debug, Clone, Copy)]
struct Station {
    surface: Surface,
    /// Index in the list of airfoil panels followed by wake panels.
    panel: usize,
    /// Distance along the surface from the stagnation point.
    xi: f64,
    /// Chord position, x/c.
    x: f64,
    /// +1 if the flow is in the direction of the panel, -1 if not.
    sign: f64,
}

#[derive(Debug, Clone)]
pub struct StationResult {
    pub surface: Surface,
    /// Airfoil panel index. `None` for a wake station.
    pub panel: Option<usize>,
    pub x: f64,
    pub xi: f64,
    pub ue: f64,
    pub theta: f64,
    pub dstar: f64,
    pub h: f64,
    pub cf: f64,
    pub flow: Flow,
}

#[derive(Debug, Clone)]
pub struct InteractionResult {
    /// Lift coefficient with the effect of the boundary layer.
    pub cl: f64,
    /// Lift coefficient of the friction-free flow.
    pub cl_inviscid: f64,
    /// Drag coefficient, from the momentum thickness at the end of the wake.
    pub cd: f64,
    /// Chord position where turbulent flow starts, upper and lower surface.
    pub transition_upper: Option<f64>,
    pub transition_lower: Option<f64>,
    /// Chord position where the flow separates (skin friction below zero),
    /// upper and lower surface.
    pub separation_upper: Option<f64>,
    pub separation_lower: Option<f64>,
    /// Pressure coefficient at each airfoil panel midpoint.
    pub cp: Vec<f64>,
    pub stations: Vec<StationResult>,
    pub iterations: usize,
    /// False if the iteration did not reach the tolerance.
    pub converged: bool,
    /// Largest relative change in the last iteration.
    pub last_change: f64,
}

/// Start values from an earlier result, for a faster solution when the
/// inputs change by a small amount.
#[derive(Debug, Clone)]
pub struct WarmStart {
    /// Values for each airfoil panel, then for each wake panel.
    values: Vec<[f64; 3]>,
    flow: Vec<Flow>,
    surface: Vec<Surface>,
}

impl WarmStart {
    fn new(lay: &Layout, values: &[[f64; 3]], flow: &[Flow], panels: usize) -> Self {
        let mut start = Self {
            values: vec![[0.0; 3]; panels],
            flow: vec![Flow::Laminar; panels],
            surface: vec![Surface::Wake; panels],
        };
        for (k, station) in lay.stations.iter().enumerate() {
            start.values[station.panel] = values[k];
            start.flow[station.panel] = flow[k];
            start.surface[station.panel] = station.surface;
        }
        start
    }

    /// Start values for a station layout. A station that was on the other
    /// surface (the stagnation point moved) gets the values of the next
    /// station on its new surface.
    fn values_for(&self, lay: &Layout) -> Option<(Vec<[f64; 3]>, Vec<Flow>)> {
        if self.values.len() <= lay.stations.iter().map(|st| st.panel).max()? {
            return None;
        }
        let n = lay.stations.len();
        let mut values = vec![[0.0; 3]; n];
        let mut flow = vec![Flow::Laminar; n];
        for k in (0..n).rev() {
            let station = &lay.stations[k];
            if self.surface[station.panel] == station.surface {
                values[k] = self.values[station.panel];
                flow[k] = self.flow[station.panel];
            } else if k + 1 < n && lay.stations[k + 1].surface == station.surface {
                values[k] = values[k + 1];
                flow[k] = flow[k + 1];
            } else {
                return None;
            }
        }
        Some((values, flow))
    }
}

struct Layout {
    stations: Vec<Station>,
    upper: std::ops::Range<usize>,
    lower: std::ops::Range<usize>,
    wake: std::ops::Range<usize>,
    stagnation_panel: usize,
}

fn layout(system: &PanelSystem, foil: &Airfoil, wake: &Wake, vt: &[f64]) -> Result<Layout, String> {
    let m = system.panel_count();
    let (len, mid) = (system.lengths(), system.midpoints());
    let (chord, x_min) = (foil.chord(), foil.x_min());
    let x_of = |p: [f64; 2]| (p[0] - x_min) / chord;
    let j = (0..m - 1)
        .filter(|&j| vt[j] < 0.0 && vt[j + 1] >= 0.0)
        .min_by(|&a, &b| x_of(mid[a]).total_cmp(&x_of(mid[b])))
        .ok_or("no stagnation point was found on the airfoil")?;
    let gap = 0.5 * (len[j] + len[j + 1]);
    let fraction = (-vt[j] / (vt[j + 1] - vt[j])).clamp(0.05, 0.95);

    let mut stations = Vec::new();
    let mut distance = (1.0 - fraction) * gap;
    for p in j + 1..m {
        if p > j + 1 {
            distance += 0.5 * (len[p - 1] + len[p]);
        }
        stations.push(Station {
            surface: Surface::Upper,
            panel: p,
            xi: distance / chord,
            x: x_of(mid[p]),
            sign: 1.0,
        });
    }
    let upper = 0..stations.len();
    let upper_end = distance + 0.5 * len[m - 1];
    distance = fraction * gap;
    for p in (0..=j).rev() {
        if p < j {
            distance += 0.5 * (len[p + 1] + len[p]);
        }
        stations.push(Station {
            surface: Surface::Lower,
            panel: p,
            xi: distance / chord,
            x: x_of(mid[p]),
            sign: -1.0,
        });
    }
    let lower = upper.end..stations.len();
    distance = upper_end;
    for w in 0..wake.panel_count() {
        distance += if w == 0 {
            0.5 * wake.len[0]
        } else {
            0.5 * (wake.len[w - 1] + wake.len[w])
        };
        stations.push(Station {
            surface: Surface::Wake,
            panel: m + w,
            xi: distance / chord,
            x: x_of(wake.mid[w]),
            sign: 1.0,
        });
    }
    let wake_range = lower.end..stations.len();
    if upper.len() < 4 || lower.len() < 4 {
        return Err("a surface has less than 4 stations".into());
    }
    Ok(Layout {
        stations,
        upper,
        lower,
        wake: wake_range,
        stagnation_panel: j,
    })
}

/// Table that gives the source strength on each panel from the mass defect
/// at each station: `sources = table * mass`.
///
/// The flow out of a panel is the mass defect at its station minus the
/// mass defect at the station before it. This one-sided difference keeps
/// the total flow correct, and it does not permit a station-to-station
/// oscillation of the mass defect (a central difference permits it).
fn source_table(layout: &Layout, panel_len: &[f64]) -> Vec<Vec<f64>> {
    let n = layout.stations.len();
    let mut table = vec![vec![0.0; n]; panel_len.len()];
    let mut add = |range: &std::ops::Range<usize>, inflow: &[usize]| {
        for k in range.clone() {
            let row = &mut table[layout.stations[k].panel];
            let scale = 1.0 / panel_len[layout.stations[k].panel];
            row[k] += scale;
            if k == range.start {
                for &source in inflow {
                    row[source] -= scale;
                }
            } else {
                row[k - 1] -= scale;
            }
        }
    };
    add(&layout.upper, &[]);
    add(&layout.lower, &[]);
    add(&layout.wake, &[layout.upper.end - 1, layout.lower.end - 1]);
    table
}

/// The stations that one set of three equations uses.
enum Equation {
    /// First station after the stagnation point.
    Similarity(usize),
    /// Interval between two stations. `start_turbulent` is true if
    /// turbulent flow starts at the first station.
    Between {
        a: usize,
        b: usize,
        start_turbulent: bool,
    },
    /// First wake station: the two boundary layers join.
    Junction {
        upper: usize,
        lower: usize,
        wake: usize,
    },
}

impl Equation {
    fn stations(&self) -> Vec<usize> {
        match *self {
            Equation::Similarity(k) => vec![k],
            Equation::Between { a, b, .. } => vec![a, b],
            Equation::Junction { upper, lower, wake } => vec![upper, lower, wake],
        }
    }

    fn residual(&self, states: &[State], re: f64) -> [f64; 3] {
        match *self {
            Equation::Similarity(_) => {
                let r = ibl::similarity_residuals(&states[0], re);
                [r[0], r[1], states[0].s]
            }
            Equation::Between {
                start_turbulent, ..
            } => {
                let a = if start_turbulent {
                    turbulent_start(&states[0], re)
                } else {
                    states[0]
                };
                ibl::residuals(&a, &states[1], re)
            }
            Equation::Junction { .. } => {
                // The two boundary layers join into one wake state at the
                // trailing edge. The momentum defect and the mass defect of
                // the two layers are added. Then the wake equations apply
                // from this state to the first wake station.
                let (u, l, w) = (&states[0], &states[1], &states[2]);
                let shear = |state: &State| {
                    if state.flow == Flow::Laminar {
                        turbulent_start(state, re).s
                    } else {
                        state.s
                    }
                };
                let ue = 0.5 * (u.ue + l.ue);
                let theta = (u.ue * u.ue * u.theta + l.ue * l.ue * l.theta) / (ue * ue);
                let joined = State {
                    xi: u.xi,
                    ue,
                    theta,
                    dstar: (u.ue * u.dstar + l.ue * l.dstar) / ue,
                    s: (u.theta * shear(u) + l.theta * shear(l)) / (u.theta + l.theta),
                    flow: Flow::Wake,
                };
                ibl::residuals(&joined, w, re)
            }
        }
    }
}

/// First turbulent state at a laminar station (no change of thickness).
fn turbulent_start(state: &State, re: f64) -> State {
    let mut start = State {
        flow: Flow::Turbulent,
        s: 0.03,
        ..*state
    };
    let d = start.derived(re);
    start.s = cl::ctau_start_sqrt(d.h, d.ctau_eq_sqrt).max(1e-4);
    start
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    pub re: f64,
    /// Forced transition position (x/c), if a trip strip is used.
    pub forced_transition_x: Option<f64>,
    pub n_crit: f64,
}

pub fn solve_interaction(
    foil: &Airfoil,
    system: &PanelSystem,
    airfoil_response: &(Vec<Vec<f64>>, Vec<Vec<f64>>),
    alpha_deg: f64,
    options: &Options,
    warm_start: Option<&WarmStart>,
) -> Result<(InteractionResult, WarmStart), String> {
    let re = options.re;
    if !(re.is_finite() && re > 0.0) {
        return Err("the Reynolds number must be more than zero".into());
    }
    let m = system.panel_count();
    let inviscid = system.solve(alpha_deg, None)?;
    let wake = Wake::new(system, alpha_deg, &inviscid, WAKE_PANELS, WAKE_LENGTH);
    let influence = SourceInfluence::new(system, airfoil_response, &wake);
    let panel_len: Vec<f64> = system.lengths().iter().chain(&wake.len).copied().collect();
    let speed_inviscid: Vec<f64> = inviscid
        .vt
        .iter()
        .chain(&wake.ue_inviscid)
        .copied()
        .collect();

    let mut vt_current = inviscid.vt.clone();
    let mut carried: Option<WarmStart> = warm_start.cloned();
    let mut outcome = None;
    // Outer loop: the stagnation point can move to a different panel.
    for _ in 0..4 {
        let lay = layout(system, foil, &wake, &vt_current)?;
        let n = lay.stations.len();
        let sources_from_mass = source_table(&lay, &panel_len);
        // d[k][j]: change of the edge speed at station k for a unit of mass
        // defect at station j.
        let mut d = vec![vec![0.0; n]; n];
        for k in 0..n {
            let row = &influence.table[lay.stations[k].panel];
            for p in 0..panel_len.len() {
                let factor = lay.stations[k].sign * row[p];
                if factor != 0.0 {
                    for j in 0..n {
                        d[k][j] += factor * sources_from_mass[p][j];
                    }
                }
            }
        }
        let mut ue_inviscid: Vec<f64> = lay
            .stations
            .iter()
            .map(|st| st.sign * speed_inviscid[st.panel])
            .collect();
        // First station of each surface: the edge speed is proportional to
        // the distance from the stagnation point, with the slope from the
        // second station. The speed at the first station itself is very
        // small, and a small move of the stagnation point changes its sign.
        for range in [&lay.upper, &lay.lower] {
            let (first, second) = (range.start, range.start + 1);
            let ratio = lay.stations[first].xi / lay.stations[second].xi;
            ue_inviscid[first] = ratio * ue_inviscid[second];
            for j in 0..n {
                d[first][j] = ratio * d[second][j];
            }
        }

        // Start values.
        let (mut values, mut flow) = match carried.as_ref().and_then(|start| start.values_for(&lay))
        {
            Some(start) => start,
            None => start_values(&lay, &ue_inviscid, options)?,
        };

        let edge_speed = |values: &[[f64; 3]]| -> Vec<f64> {
            (0..n)
                .map(|k| {
                    let extra: f64 = (0..n).map(|j| d[k][j] * values[j][1]).sum();
                    (ue_inviscid[k] + extra).max(1e-3)
                })
                .collect()
        };
        let state_of = |k: usize, values: &[[f64; 3]], ue: &[f64], flow: &[Flow]| State {
            xi: lay.stations[k].xi,
            ue: ue[k],
            theta: values[k][0],
            dstar: values[k][1] / ue[k],
            s: values[k][2],
            flow: flow[k],
        };

        // Residuals of all equations, and the Jacobian matrix if wanted.
        let assemble = |values: &[[f64; 3]], flow: &[Flow], with_jacobian: bool| {
            let ue = edge_speed(values);
            let equations = equations(&lay, flow);
            let mut residual = vec![0.0; 3 * n];
            let mut jacobian = if with_jacobian {
                vec![vec![0.0; 3 * n]; 3 * n]
            } else {
                Vec::new()
            };
            for (k, equation) in equations.iter().enumerate() {
                let members = equation.stations();
                let states: Vec<State> = members
                    .iter()
                    .map(|&i| state_of(i, values, &ue, flow))
                    .collect();
                let base = equation.residual(&states, re);
                for row in 0..3 {
                    residual[3 * k + row] = base[row];
                }
                if !with_jacobian {
                    continue;
                }
                for (local, &i) in members.iter().enumerate() {
                    // Partial derivatives for theta, dstar, s, and ue.
                    let mut partial = [[0.0; 3]; 4];
                    for variable in 0..4 {
                        let mut shifted = states.clone();
                        let target = &mut shifted[local];
                        let value = match variable {
                            0 => &mut target.theta,
                            1 => &mut target.dstar,
                            2 => &mut target.s,
                            _ => &mut target.ue,
                        };
                        let step = 1e-6 * value.abs().max(1e-8);
                        *value += step;
                        let moved = equation.residual(&shifted, re);
                        for row in 0..3 {
                            partial[variable][row] = (moved[row] - base[row]) / step;
                        }
                    }
                    let (mass, speed) = (values[i][1], ue[i]);
                    for row in 0..3 {
                        let line = &mut jacobian[3 * k + row];
                        line[3 * i] += partial[0][row];
                        line[3 * i + 1] += partial[1][row] / speed;
                        line[3 * i + 2] += partial[2][row];
                        // Effect of the mass defect at all stations through
                        // the edge speed at this station.
                        let through_speed =
                            partial[3][row] - partial[1][row] * mass / (speed * speed);
                        if through_speed != 0.0 {
                            for j in 0..n {
                                line[3 * j + 1] += through_speed * d[i][j];
                            }
                        }
                    }
                }
            }
            (residual, jacobian)
        };
        let norm = |r: &[f64]| r.iter().map(|v| v * v).sum::<f64>().sqrt();

        let mut converged = false;
        let mut iterations = 0;
        let mut last_change = f64::INFINITY;
        let mut stalled = 0;
        for iteration in 0..MAX_ITERATIONS {
            iterations = iteration + 1;
            let (residual, jacobian) = assemble(&values, &flow, true);
            if residual.iter().any(|v| !v.is_finite()) {
                break;
            }
            let rhs: Vec<f64> = residual.iter().map(|v| -v).collect();
            // A matrix with no unique solution ends this attempt without a
            // converged result. A different start can still succeed.
            let Some(delta) = linalg::solve(jacobian, rhs) else {
                break;
            };

            // Limit the step: no thickness changes by more than a set part.
            let mut relax: f64 = 1.0;
            let mut change: f64 = 0.0;
            for k in 0..n {
                let scaled = [
                    delta[3 * k] / values[k][0],
                    delta[3 * k + 1] / values[k][1],
                    if flow[k] == Flow::Laminar {
                        delta[3 * k + 2] / 10.0
                    } else {
                        delta[3 * k + 2] / (10.0 * values[k][2])
                    },
                ];
                for value in scaled {
                    change = change.max(value.abs());
                    if value > 1.5 {
                        relax = relax.min(1.5 / value);
                    } else if value < -0.5 {
                        relax = relax.min(-0.5 / value);
                    }
                }
            }
            // Line search: make the step smaller until the residual decreases.
            let before = norm(&residual);
            let apply = |relax: f64| {
                let mut next = values.clone();
                for k in 0..n {
                    next[k][0] += relax * delta[3 * k];
                    next[k][1] += relax * delta[3 * k + 1];
                    next[k][2] += relax * delta[3 * k + 2];
                    if flow[k] == Flow::Laminar {
                        next[k][2] = next[k][2].max(0.0);
                    } else {
                        next[k][2] = next[k][2].clamp(1e-7, 0.3);
                    }
                }
                next
            };
            let mut next = apply(relax);
            let mut after = norm(&assemble(&next, &flow, false).0);
            let mut halvings = 0;
            while (after >= before || after.is_nan()) && halvings < 8 {
                relax *= 0.5;
                halvings += 1;
                next = apply(relax);
                after = norm(&assemble(&next, &flow, false).0);
            }
            values = next;
            last_change = change;
            // No progress: the iteration is in a region without a solution
            // path (for example deep stall). Stop early.
            stalled = if after < 0.999 * before || before < 1e-9 {
                0
            } else {
                stalled + 1
            };
            if stalled >= MAX_STALLED_STEPS {
                break;
            }
            let ue_new = edge_speed(&values);
            let changed = update_flow(&lay, &mut values, &mut flow, &ue_new, options);
            if change < TOLERANCE && !changed {
                converged = true;
                break;
            }
        }

        // Surface speed with the sources, to check the stagnation point.
        let ue = edge_speed(&values);
        for k in lay.upper.clone().chain(lay.lower.clone()) {
            vt_current[lay.stations[k].panel] = lay.stations[k].sign * ue[k];
        }
        let moved = layout(system, foil, &wake, &vt_current)
            .map(|next| next.stagnation_panel != lay.stagnation_panel)
            .unwrap_or(false);
        carried = Some(WarmStart::new(&lay, &values, &flow, panel_len.len()));
        outcome = Some((
            lay,
            values,
            flow,
            ue,
            sources_from_mass,
            iterations,
            converged,
            last_change,
        ));
        if !moved || !converged {
            break;
        }
        // The stagnation point moved to a different panel: solve again
        // with the new station layout, from the values of this solution.
    }

    let (lay, values, flow, ue, sources_from_mass, iterations, converged, last_change) =
        outcome.ok_or("the two-way calculation did not run")?;
    let n = lay.stations.len();
    // Vortex strengths with the effect of all sources.
    let mut gamma = inviscid.gamma.clone();
    for p in 0..panel_len.len() {
        let sigma: f64 = (0..n).map(|j| sources_from_mass[p][j] * values[j][1]).sum();
        for (g, change) in gamma.iter_mut().zip(&influence.dgamma[p]) {
            *g += change * sigma;
        }
    }
    let mut cp = vec![0.0; m];
    let mut stations = Vec::with_capacity(n);
    for k in 0..n {
        let state = State {
            xi: lay.stations[k].xi,
            ue: ue[k],
            theta: values[k][0],
            dstar: values[k][1] / ue[k],
            s: values[k][2],
            flow: flow[k],
        };
        let derived = state.derived(re);
        if lay.stations[k].surface != Surface::Wake {
            cp[lay.stations[k].panel] = 1.0 - ue[k] * ue[k];
        }
        stations.push(StationResult {
            surface: lay.stations[k].surface,
            panel: (lay.stations[k].surface != Surface::Wake).then_some(lay.stations[k].panel),
            x: lay.stations[k].x,
            xi: state.xi,
            ue: state.ue,
            theta: state.theta,
            dstar: state.dstar,
            h: derived.h,
            cf: derived.cf,
            flow: state.flow,
        });
    }
    let first = |range: &std::ops::Range<usize>, test: &dyn Fn(&StationResult) -> bool| {
        stations[range.clone()]
            .iter()
            .find(|st| test(st))
            .map(|st| st.x)
    };
    let turbulent = |st: &StationResult| st.flow == Flow::Turbulent;
    let separated = |st: &StationResult| st.cf < 0.0;
    let end = stations.last().unwrap();
    let cd = 2.0 * end.theta * end.ue.powf(0.5 * (end.h + 5.0));
    let result = InteractionResult {
        cl: system.lift_from_gamma(&gamma),
        cl_inviscid: inviscid.cl,
        cd,
        transition_upper: first(&lay.upper, &turbulent),
        transition_lower: first(&lay.lower, &turbulent),
        separation_upper: first(&lay.upper, &separated),
        separation_lower: first(&lay.lower, &separated),
        cp,
        stations,
        iterations,
        converged,
        last_change,
    };
    let warm = carried.ok_or("the two-way calculation did not run")?;
    Ok((result, warm))
}

/// The equations for each station, in station order.
fn equations(lay: &Layout, flow: &[Flow]) -> Vec<Equation> {
    let mut list = Vec::with_capacity(lay.stations.len());
    for range in [&lay.upper, &lay.lower] {
        for k in range.clone() {
            if k == range.start {
                list.push(Equation::Similarity(k));
            } else {
                let start_turbulent = flow[k] == Flow::Turbulent && flow[k - 1] == Flow::Laminar;
                list.push(Equation::Between {
                    a: k - 1,
                    b: k,
                    start_turbulent,
                });
            }
        }
    }
    for k in lay.wake.clone() {
        if k == lay.wake.start {
            list.push(Equation::Junction {
                upper: lay.upper.end - 1,
                lower: lay.lower.end - 1,
                wake: k,
            });
        } else {
            list.push(Equation::Between {
                a: k - 1,
                b: k,
                start_turbulent: false,
            });
        }
    }
    list
}

/// Sets the flow type of each station from the amplification factor and
/// the trip strip. Returns true if a flow type changed.
///
/// Turbulent flow starts at the first laminar station that reaches the
/// critical amplification factor or the trip strip. In one solution, the
/// start of turbulent flow can only move toward the leading edge. This
/// rule prevents a jump of the transition point between two stations
/// from one iteration to the next.
fn update_flow(
    lay: &Layout,
    values: &mut [[f64; 3]],
    flow: &mut [Flow],
    ue: &[f64],
    options: &Options,
) -> bool {
    let mut changed = false;
    for range in [&lay.upper, &lay.lower] {
        let mut start = None;
        for k in range.clone() {
            if flow[k] != Flow::Laminar {
                start = Some(k.saturating_sub(1).max(range.start));
                break;
            }
            let thick = options.re * ue[k] * values[k][0] >= ibl::TRIP_MIN_RE_THETA;
            let tripped = thick
                && options
                    .forced_transition_x
                    .is_some_and(|x| lay.stations[k].x >= x);
            if tripped || values[k][2] >= options.n_crit {
                start = Some(k);
                break;
            }
        }
        for k in range.clone() {
            let turbulent = start.is_some_and(|s| k > s);
            let wanted = if turbulent {
                Flow::Turbulent
            } else {
                Flow::Laminar
            };
            if flow[k] != wanted {
                changed = true;
                flow[k] = wanted;
                values[k][2] = if turbulent { 0.03 } else { values[k - 1][2] };
            }
        }
    }
    changed
}

/// Start values: a direct boundary layer march with the friction-free
/// edge speed. After a separation point, the last values are repeated.
fn start_values(
    lay: &Layout,
    ue: &[f64],
    options: &Options,
) -> Result<(Vec<[f64; 3]>, Vec<Flow>), String> {
    let n = lay.stations.len();
    let mut values = vec![[0.0; 3]; n];
    let mut flow = vec![Flow::Laminar; n];
    for range in [&lay.upper, &lay.lower] {
        let xi: Vec<f64> = range.clone().map(|k| lay.stations[k].xi).collect();
        let speed: Vec<f64> = range.clone().map(|k| ue[k].max(1e-3)).collect();
        let forced = options.forced_transition_x.map(|limit| {
            range
                .clone()
                .find(|&k| lay.stations[k].x >= limit)
                .map_or(f64::INFINITY, |k| lay.stations[k].xi)
        });
        let result = ibl::march(&xi, &speed, options.re, true, forced, options.n_crit)?;
        let mut last = result.states[0];
        let mut amplification = 0.0;
        for (offset, k) in range.clone().enumerate() {
            if let Some(state) = result.states.get(offset) {
                last = *state;
            }
            let turbulent = result.transition.is_some_and(|(index, _)| offset > index);
            flow[k] = if turbulent {
                Flow::Turbulent
            } else {
                Flow::Laminar
            };
            if last.flow == Flow::Laminar {
                amplification = last.s;
            }
            let s = match (turbulent, last.flow) {
                (true, Flow::Laminar) => 0.03,
                (true, _) => last.s,
                (false, _) => amplification,
            };
            // After the march stops, the mass defect stays constant.
            let mass = if offset < result.states.len() {
                speed[offset] * last.dstar
            } else {
                values[k - 1][1]
            };
            values[k] = [last.theta, mass, s];
        }
    }
    // Wake: the sum of the two boundary layers at the trailing edge, then
    // a direct march. If the march has no solution at a station, the
    // momentum thickness stays constant and the shape factor decreases.
    let (u, l) = (lay.upper.end - 1, lay.lower.end - 1);
    let shear = |k: usize| {
        if flow[k] == Flow::Turbulent {
            values[k][2]
        } else {
            0.03
        }
    };
    let theta = values[u][0] + values[l][0];
    let dstar = values[u][1] / ue[u].max(1e-3) + values[l][1] / ue[l].max(1e-3);
    let first = lay.wake.start;
    let mut state = State {
        xi: lay.stations[first].xi,
        ue: ue[first].max(1e-3),
        theta,
        dstar: theta * (dstar / theta).clamp(1.2, 2.4),
        s: (values[u][0] * shear(u) + values[l][0] * shear(l)) / theta,
        flow: Flow::Wake,
    };
    for k in lay.wake.clone() {
        if k > first {
            let (xi, speed) = (lay.stations[k].xi, ue[k].max(1e-3));
            state = match ibl::step(&state, xi, speed, options.re) {
                Some(next) if next.h() > 1.001 => next,
                _ => {
                    let h = 1.0 + (state.h() - 1.0) / (1.0 + 10.0 * (xi - state.xi));
                    State {
                        xi,
                        ue: speed,
                        dstar: h * state.theta,
                        ..state
                    }
                }
            };
        }
        flow[k] = Flow::Wake;
        values[k] = [state.theta, state.ue * state.dstar, state.s];
    }
    Ok((values, flow))
}

/// Solves the two-way calculation and tries more than one start:
/// 1. the start values from an earlier result, if given;
/// 2. a start from the friction-free flow;
/// 3. a sequence of angles from zero to the target angle, where each
///    solution is the start for the next angle.
///
/// The result can still have `converged == false` (for example in stall).
pub fn solve_robust(
    foil: &Airfoil,
    system: &PanelSystem,
    airfoil_response: &(Vec<Vec<f64>>, Vec<Vec<f64>>),
    alpha_deg: f64,
    options: &Options,
    warm_start: Option<&WarmStart>,
) -> Result<(InteractionResult, WarmStart), String> {
    let solve = |alpha: f64, warm: Option<&WarmStart>| {
        solve_interaction(foil, system, airfoil_response, alpha, options, warm)
    };
    let mut fallback = None;
    for warm in [warm_start, None] {
        if warm.is_none() && fallback.is_some() && warm_start.is_none() {
            break;
        }
        match solve(alpha_deg, warm) {
            Ok((result, start)) if result.converged => return Ok((result, start)),
            Ok(pair) => fallback = Some(Ok(pair)),
            Err(error) => fallback = Some(Err(error)),
        }
    }
    // Sequence of angles.
    let mut alpha = 0.0;
    let mut step = 2.0_f64.copysign(alpha_deg);
    let mut warm = match solve(0.0, None) {
        Ok((result, start)) if result.converged => start,
        _ => return fallback.unwrap_or(Err("the two-way calculation has no result".into())),
    };
    let mut tries = 0;
    while (alpha_deg - alpha).abs() > 1e-9 && tries < 14 {
        tries += 1;
        let next = if (alpha_deg - alpha).abs() <= step.abs() {
            alpha_deg
        } else {
            alpha + step
        };
        match solve(next, Some(&warm)) {
            Ok((result, start)) if result.converged => {
                if next == alpha_deg {
                    return Ok((result, start));
                }
                alpha = next;
                warm = start;
            }
            other => {
                if next == alpha_deg {
                    if let Ok(pair) = other {
                        fallback = Some(Ok(pair));
                    }
                }
                step *= 0.5;
                if step.abs() < 0.2 {
                    break;
                }
            }
        }
    }
    fallback.unwrap_or(Err("the two-way calculation has no result".into()))
}

/// An airfoil with the data that the two-way calculation needs. Make it
/// one time for each shape, then use it for each angle and Reynolds number.
pub struct FrictionSystem {
    pub foil: Airfoil,
    pub system: PanelSystem,
    response: (Vec<Vec<f64>>, Vec<Vec<f64>>),
}

impl FrictionSystem {
    pub fn new(foil: Airfoil) -> Result<Self, String> {
        let system = PanelSystem::new(&foil)?;
        let response = system.source_response();
        Ok(Self {
            foil,
            system,
            response,
        })
    }

    /// One attempt, from the given start values or from the friction-free flow.
    pub fn solve(
        &self,
        alpha_deg: f64,
        options: &Options,
        warm_start: Option<&WarmStart>,
    ) -> Result<(InteractionResult, WarmStart), String> {
        solve_interaction(
            &self.foil,
            &self.system,
            &self.response,
            alpha_deg,
            options,
            warm_start,
        )
    }

    /// Two attempts: from the given start values, then from the
    /// friction-free flow. This is faster than `solve_robust` when there is
    /// no solution.
    pub fn solve_quick(
        &self,
        alpha_deg: f64,
        options: &Options,
        warm_start: Option<&WarmStart>,
    ) -> Result<(InteractionResult, WarmStart), String> {
        if warm_start.is_some() {
            if let Ok((result, start)) = self.solve(alpha_deg, options, warm_start) {
                if result.converged {
                    return Ok((result, start));
                }
            }
        }
        self.solve(alpha_deg, options, None)
    }

    /// More than one attempt. See `solve_robust`.
    pub fn solve_robust(
        &self,
        alpha_deg: f64,
        options: &Options,
        warm_start: Option<&WarmStart>,
    ) -> Result<(InteractionResult, WarmStart), String> {
        solve_robust(
            &self.foil,
            &self.system,
            &self.response,
            alpha_deg,
            options,
            warm_start,
        )
    }
}
