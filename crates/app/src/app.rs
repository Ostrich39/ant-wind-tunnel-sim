//! Program state and window layout.

use crate::view::{self, Mapper};
use eframe::egui::{
    self, vec2, Align2, Color32, FontId, Pos2, Rect, RichText, Sense, Shape, Stroke, Ui,
};
use solver::closure::Flow;
use solver::interaction::{FrictionSystem, InteractionResult, Options, Surface, WarmStart};
use solver::{naca4_for_friction, Airfoil, Solution};

const BACKGROUND: Color32 = Color32::from_rgb(15, 19, 24);
const PANEL: Color32 = Color32::from_rgb(22, 27, 34);
const GRID: Color32 = Color32::from_rgb(40, 48, 58);
const AXIS: Color32 = Color32::from_rgb(92, 104, 118);
const TEXT_DIM: Color32 = Color32::from_rgb(139, 150, 163);
const TEXT: Color32 = Color32::from_rgb(226, 232, 238);
const BODY_FILL: Color32 = Color32::from_rgb(52, 62, 74);
const BODY_LINE: Color32 = Color32::from_rgb(214, 222, 230);
const LOWER_LINE: Color32 = Color32::from_rgb(150, 160, 172);
const TRANSITION: Color32 = Color32::from_rgb(123, 216, 143);
const SEPARATION: Color32 = Color32::from_rgb(255, 209, 102);

/// The inputs that set the result. A new solution is necessary only when
/// these values change.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Inputs {
    /// Maximum camber, percent of chord (first NACA digit).
    pub camber: u32,
    /// Position of maximum camber, tenths of chord (second NACA digit).
    pub camber_pos: u32,
    /// Thickness, percent of chord (last two NACA digits).
    pub thickness: u32,
    /// Angle of attack, degrees.
    pub alpha: f64,
    /// Number of panels (even).
    pub panels: usize,
    /// Reynolds number for the chord.
    pub reynolds: f64,
    /// If true, turbulent flow starts at `trip_x` or before it.
    pub trip: bool,
    /// Forced transition position, x/c.
    pub trip_x: f64,
    /// Critical amplification factor for natural transition. A lower
    /// value is for a flow with more turbulence.
    pub n_crit: f64,
}

impl Default for Inputs {
    fn default() -> Self {
        Self {
            camber: 2,
            camber_pos: 4,
            thickness: 12,
            alpha: 4.0,
            panels: 120,
            reynolds: 3.0e6,
            trip: false,
            trip_x: 0.05,
            n_crit: 9.0,
        }
    }
}

impl Inputs {
    /// Makes the inputs valid. A NACA 4-digit shape with camber must have a
    /// camber position, and a shape without camber must not have one.
    pub fn normalized(mut self) -> Self {
        self.camber = self.camber.min(9);
        self.thickness = self.thickness.clamp(1, 40);
        self.camber_pos = if self.camber == 0 {
            0
        } else {
            self.camber_pos.clamp(1, 9)
        };
        self.panels = (self.panels.clamp(60, 200) / 2) * 2;
        self.alpha = if self.alpha.is_finite() {
            self.alpha.clamp(-20.0, 20.0)
        } else {
            0.0
        };
        self.reynolds = if self.reynolds.is_finite() {
            self.reynolds.clamp(1.0e5, 1.0e7)
        } else {
            3.0e6
        };
        self.trip_x = if self.trip_x.is_finite() {
            self.trip_x.clamp(0.01, 0.9)
        } else {
            0.05
        };
        self.n_crit = if self.n_crit.is_finite() {
            self.n_crit.clamp(3.0, 12.0)
        } else {
            9.0
        };
        self
    }

    pub fn code(&self) -> String {
        format!("{}{}{:02}", self.camber, self.camber_pos, self.thickness)
    }
}

pub struct Computed {
    pub foil: Airfoil,
    /// Friction-free flow.
    pub solution: Solution,
    /// Two-way result (flow and boundary layer). `None` if the calculation
    /// has no converged result; `friction_note` then gives the cause.
    pub friction: Option<InteractionResult>,
    pub friction_note: Option<String>,
    pub solve_ms: f64,
}

impl Computed {
    /// Pressure coefficient on each panel: with friction if that result
    /// is available, friction-free if not.
    pub fn cp(&self) -> &[f64] {
        match &self.friction {
            Some(result) => &result.cp,
            None => &self.solution.cp,
        }
    }
}

/// Keeps the data of the current airfoil shape and the last solution, so
/// that a new result after a small input change is fast.
#[derive(Default)]
pub struct Analyzer {
    shape: Option<(String, usize)>,
    friction: Option<FrictionSystem>,
    warm: Option<WarmStart>,
}

impl Analyzer {
    pub fn compute(&mut self, inputs: &Inputs) -> Result<Computed, String> {
        let start = std::time::Instant::now();
        let shape = (inputs.code(), inputs.panels);
        if self.shape.as_ref() != Some(&shape) || self.friction.is_none() {
            let foil = naca4_for_friction(&shape.0, shape.1)?;
            self.friction = Some(FrictionSystem::new(foil)?);
            self.shape = Some(shape);
        }
        let system = self
            .friction
            .as_ref()
            .ok_or("the airfoil data is not available")?;
        let solution = system.system.solve(inputs.alpha, None)?;
        let options = Options {
            re: inputs.reynolds,
            forced_transition_x: inputs.trip.then_some(inputs.trip_x),
            n_crit: inputs.n_crit,
        };
        let (friction, friction_note) =
            match system.solve_quick(inputs.alpha, &options, self.warm.as_ref()) {
                Ok((result, warm)) if result.converged => {
                    self.warm = Some(warm);
                    (Some(result), None)
                }
                Ok(_) => {
                    self.warm = None;
                    (None, Some(NO_SOLUTION_TEXT.to_owned()))
                }
                Err(message) => {
                    self.warm = None;
                    (None, Some(message))
                }
            };
        Ok(Computed {
            foil: system.foil.clone(),
            solution,
            friction,
            friction_note,
            solve_ms: start.elapsed().as_secs_f64() * 1000.0,
        })
    }
}

const NO_SOLUTION_TEXT: &str = "The friction calculation has no solution at this condition. \
     The cause is usually stall (the flow separates from a large part of the surface). \
     Decrease the angle of attack, or move the angle slowly from a lower value.";

pub struct TunnelApp {
    inputs: Inputs,
    solved_for: Option<Inputs>,
    analyzer: Analyzer,
    result: Result<Computed, String>,
    show_pressure_lines: bool,
}

impl TunnelApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut visuals = egui::Visuals::dark();
        visuals.panel_fill = PANEL;
        visuals.window_fill = PANEL;
        visuals.override_text_color = Some(TEXT);
        cc.egui_ctx.set_visuals(visuals);

        let inputs = Inputs::default();
        Self {
            inputs,
            solved_for: None,
            analyzer: Analyzer::default(),
            result: Err("no result yet".to_owned()),
            show_pressure_lines: true,
        }
    }

    fn controls(&mut self, ui: &mut Ui) {
        ui.add_space(10.0);
        ui.label(
            RichText::new("ANT WIND TUNNEL")
                .size(11.0)
                .color(TEXT_DIM)
                .strong(),
        );
        ui.label(
            RichText::new(format!("NACA {}", self.inputs.code()))
                .size(30.0)
                .monospace()
                .color(TEXT),
        );
        ui.add_space(12.0);

        section(ui, "FLOW");
        ui.add(
            egui::Slider::new(&mut self.inputs.alpha, -10.0..=15.0)
                .step_by(0.1)
                .fixed_decimals(1)
                .suffix("°")
                .text("Angle of attack"),
        );
        ui.add(
            egui::Slider::new(&mut self.inputs.reynolds, 1.0e5..=1.0e7)
                .logarithmic(true)
                .custom_formatter(|v, _| format!("{:.2} M", v / 1.0e6))
                .custom_parser(|text| {
                    text.trim()
                        .trim_end_matches('M')
                        .trim()
                        .parse::<f64>()
                        .ok()
                        .map(|v| v * 1.0e6)
                })
                .text("Reynolds number"),
        );
        ui.add_space(12.0);

        section(ui, "SHAPE");
        ui.add(
            egui::Slider::new(&mut self.inputs.thickness, 4..=30)
                .suffix(" %")
                .text("Thickness"),
        );
        ui.add(
            egui::Slider::new(&mut self.inputs.camber, 0..=9)
                .suffix(" %")
                .text("Camber"),
        );
        let has_camber = self.inputs.camber > 0;
        if has_camber && self.inputs.camber_pos == 0 {
            self.inputs.camber_pos = 4;
        }
        ui.add_enabled(
            has_camber,
            egui::Slider::new(&mut self.inputs.camber_pos, 1..=9)
                .custom_formatter(|v, _| format!("{:.0} %", v * 10.0))
                .text("Camber position"),
        );
        ui.add_space(12.0);

        section(ui, "BOUNDARY LAYER");
        ui.checkbox(&mut self.inputs.trip, "Force transition (trip strip)");
        ui.add_enabled(
            self.inputs.trip,
            egui::Slider::new(&mut self.inputs.trip_x, 0.01..=0.5)
                .custom_formatter(|v, _| format!("{:.0} %", v * 100.0))
                .text("Trip position"),
        );
        ui.add_space(12.0);

        ui.add(
            egui::Slider::new(&mut self.inputs.n_crit, 3.0..=12.0)
                .step_by(0.5)
                .fixed_decimals(1)
                .text("Transition factor N"),
        );
        ui.label(
            RichText::new("N = 9: clean wind tunnel. Lower N: more turbulence in the flow.")
                .size(11.0)
                .color(TEXT_DIM),
        );
        ui.add_space(12.0);

        section(ui, "SOLVER");
        ui.add(
            egui::Slider::new(&mut self.inputs.panels, 60..=200)
                .step_by(20.0)
                .text("Panels"),
        );
        ui.checkbox(&mut self.show_pressure_lines, "Show pressure lines");
        ui.add_space(12.0);

        section(ui, "RESULTS");
        match &self.result {
            Ok(done) => {
                match &done.friction {
                    Some(friction) => {
                        readout(ui, "Lift coefficient", &format!("{:+.4}", friction.cl));
                        ui.label(
                            RichText::new(format!(
                                "Without friction: {:+.4}",
                                friction.cl_inviscid
                            ))
                            .size(11.0)
                            .color(TEXT_DIM),
                        );
                        readout(ui, "Drag coefficient", &format!("{:.5}", friction.cd));
                        readout(
                            ui,
                            "Lift / drag",
                            &format!("{:.1}", friction.cl / friction.cd),
                        );
                        ui.add_space(4.0);
                        ui.label(
                            RichText::new("Start of turbulent flow")
                                .size(12.0)
                                .color(TEXT_DIM),
                        );
                        ui.label(
                            RichText::new(format!(
                                "Upper {}   Lower {}",
                                position_text(friction.transition_upper),
                                position_text(friction.transition_lower)
                            ))
                            .size(13.0)
                            .color(TEXT),
                        );
                        if let Some(text) = separation_text(friction) {
                            ui.label(RichText::new(text).size(12.0).color(SEPARATION));
                        }
                    }
                    None => {
                        readout(
                            ui,
                            "Lift coefficient (no friction)",
                            &format!("{:+.4}", done.solution.cl),
                        );
                        ui.label(RichText::new("Drag coefficient").size(12.0).color(TEXT_DIM));
                        ui.label(RichText::new("No value").size(20.0).color(SEPARATION));
                        if let Some(note) = &done.friction_note {
                            ui.label(RichText::new(note).size(12.0).color(SEPARATION));
                        }
                    }
                }
                ui.label(
                    RichText::new(format!("Solved in {:.0} ms", done.solve_ms))
                        .size(11.0)
                        .color(TEXT_DIM),
                );
            }
            Err(message) => {
                ui.colored_label(view::PRESSURE, format!("No result: {message}"));
            }
        }

        ui.add_space(12.0);
        ui.label(
            RichText::new(
                "Compared with wind tunnel data for the NACA 0012 only. \
                 The method is not reliable near stall and after stall. \
                 See GOALS.md for the measured differences.",
            )
            .size(11.0)
            .color(TEXT_DIM),
        );
        ui.add_space(10.0);
    }

    fn airfoil_view(&self, ui: &mut Ui, size: egui::Vec2) {
        let (response, painter) = ui.allocate_painter(size, Sense::hover());
        let rect = response.rect;
        painter.rect_filled(rect, 0.0, BACKGROUND);
        let Ok(done) = &self.result else { return };

        let pts = &done.foil.points;
        let cp = done.cp();
        let n = pts.len() - 1;
        let map = view::fit_uniform(rect.shrink(16.0), -0.12, 1.12, -0.30, 0.30);
        // Show the airfoil at its angle of attack. The flow is horizontal.
        let angle = -self.inputs.alpha.to_radians();
        let pivot = [0.25, 0.0];
        let world: Vec<[f64; 2]> = pts.iter().map(|p| view::rotate(*p, pivot, angle)).collect();
        let screen: Vec<Pos2> = world.iter().map(|p| map.to_screen(p[0], p[1])).collect();

        // Flow direction marker.
        let y = rect.top() + 30.0;
        let (x0, x1) = (rect.left() + 24.0, rect.left() + 84.0);
        let stroke = Stroke::new(1.5, AXIS);
        painter.line_segment([Pos2::new(x0, y), Pos2::new(x1, y)], stroke);
        painter.line_segment([Pos2::new(x1, y), Pos2::new(x1 - 7.0, y - 4.0)], stroke);
        painter.line_segment([Pos2::new(x1, y), Pos2::new(x1 - 7.0, y + 4.0)], stroke);
        label(
            &painter,
            Pos2::new(x1 + 8.0, y),
            Align2::LEFT_CENTER,
            "FLOW",
            TEXT_DIM,
        );

        // Pressure lines: length shows the size of the pressure coefficient.
        if self.show_pressure_lines {
            let length_per_cp = 0.05 * map.x_scale();
            for j in 0..n {
                let (a, b) = (world[j], world[j + 1]);
                let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                let len = dx.hypot(dy);
                // Outward normal for clockwise point order.
                let (nx, ny) = (-dy / len, dx / len);
                let cp = cp[j];
                let reach = (cp.abs().min(8.0) as f32) * length_per_cp;
                let base = map.to_screen(0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]));
                // Screen y points down, thus the sign of ny changes.
                let tip = base + vec2(nx as f32, -ny as f32) * reach;
                painter.line_segment(
                    [base, tip],
                    Stroke::new(1.0, view::cp_color(cp).gamma_multiply(0.75)),
                );
            }
            legend(&painter, rect);
        }

        // Body fill: one triangle mesh between the lower and upper surfaces.
        // A mesh has no seams between its triangles.
        let mut mesh = egui::Mesh::default();
        for point in &screen[..n] {
            mesh.colored_vertex(*point, BODY_FILL);
        }
        let index = |i: usize| (i % n) as u32;
        for k in 0..n / 2 {
            let (a, b, c, d) = (index(k), index(k + 1), index(n - k - 1), index(n - k));
            mesh.add_triangle(a, b, c);
            mesh.add_triangle(a, c, d);
        }
        painter.add(Shape::mesh(mesh));
        painter.add(Shape::closed_line(
            screen[..n].to_vec(),
            Stroke::new(1.5, BODY_LINE),
        ));

        // Boundary layer markers: transition (T) and separation (S).
        if let Some(friction) = &done.friction {
            let marker = |panel: usize, text: &str, color: Color32| {
                let (a, b) = (world[panel], world[panel + 1]);
                let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                let len = dx.hypot(dy);
                // Outward normal in screen coordinates (screen y points down).
                let normal = vec2((-dy / len) as f32, -(dx / len) as f32);
                let base = map.to_screen(0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]));
                painter.line_segment([base, base + normal * 14.0], Stroke::new(2.0, color));
                painter.circle_filled(base, 3.0, color);
                painter.text(
                    base + normal * 24.0,
                    Align2::CENTER_CENTER,
                    text,
                    FontId::proportional(13.0),
                    color,
                );
            };
            for surface in [Surface::Upper, Surface::Lower] {
                let stations = friction.stations.iter().filter(|st| st.surface == surface);
                if let Some(panel) = stations
                    .clone()
                    .find(|st| st.flow == Flow::Turbulent)
                    .and_then(|st| st.panel)
                {
                    marker(panel, "T", TRANSITION);
                }
                if let Some(panel) = stations
                    .clone()
                    .find(|st| st.cf < 0.0)
                    .and_then(|st| st.panel)
                {
                    marker(panel, "S", SEPARATION);
                }
            }
        }
    }

    fn pressure_plot(&self, ui: &mut Ui, size: egui::Vec2) {
        let (response, painter) = ui.allocate_painter(size, Sense::hover());
        let rect = response.rect;
        painter.rect_filled(rect, 0.0, BACKGROUND);
        let Ok(done) = &self.result else { return };

        let plot = Rect::from_min_max(rect.min + vec2(64.0, 30.0), rect.max - vec2(24.0, 44.0));
        let cp_min = done.cp().iter().cloned().fold(-1.0, f64::min);
        let top_value = (cp_min - 0.25).floor();
        // Negative pressure coefficient is at the top (aerodynamics convention).
        let map = Mapper {
            rect: plot,
            x0: 0.0,
            x1: 1.0,
            y0: 1.0,
            y1: top_value,
        };

        for value in view::nice_ticks(top_value, 1.0, 6) {
            let y = map.to_screen(0.0, value).y;
            let color = if value == 0.0 { AXIS } else { GRID };
            painter.line_segment(
                [Pos2::new(plot.left(), y), Pos2::new(plot.right(), y)],
                Stroke::new(1.0, color),
            );
            label(
                &painter,
                Pos2::new(plot.left() - 8.0, y),
                Align2::RIGHT_CENTER,
                &format!("{value:.1}"),
                TEXT_DIM,
            );
        }
        for value in view::nice_ticks(0.0, 1.0, 5) {
            let x = map.to_screen(value, 0.0).x;
            painter.line_segment(
                [Pos2::new(x, plot.top()), Pos2::new(x, plot.bottom())],
                Stroke::new(1.0, GRID),
            );
            label(
                &painter,
                Pos2::new(x, plot.bottom() + 6.0),
                Align2::CENTER_TOP,
                &format!("{value:.1}"),
                TEXT_DIM,
            );
        }
        label(
            &painter,
            Pos2::new(plot.center().x, rect.bottom() - 8.0),
            Align2::CENTER_BOTTOM,
            "Position along chord, x/c",
            TEXT_DIM,
        );
        label(
            &painter,
            Pos2::new(rect.left() + 12.0, rect.top() + 8.0),
            Align2::LEFT_TOP,
            "PRESSURE COEFFICIENT  (negative values at the top)",
            TEXT_DIM,
        );

        let foil = &done.foil;
        let (x_min, chord) = (foil.x_min(), foil.chord());
        let cp = done.cp();
        let n = cp.len();
        let point = |j: usize| map.to_screen((done.solution.mid[j][0] - x_min) / chord, cp[j]);
        let clip = painter.with_clip_rect(plot.expand(1.0));
        // Panels 0..n/2 are on the lower surface. The others are on the upper surface.
        clip.line(
            (0..n / 2).map(point).collect(),
            Stroke::new(1.5, LOWER_LINE),
        );
        clip.line((n / 2..n).map(point).collect(), Stroke::new(2.0, BODY_LINE));

        let key = Pos2::new(plot.right() - 150.0, plot.top() + 12.0);
        for (row, (name, color, width)) in [
            ("Upper surface", BODY_LINE, 2.0),
            ("Lower surface", LOWER_LINE, 1.5),
        ]
        .iter()
        .enumerate()
        {
            let y = key.y + row as f32 * 18.0;
            painter.line_segment(
                [Pos2::new(key.x, y), Pos2::new(key.x + 24.0, y)],
                Stroke::new(*width, *color),
            );
            label(
                &painter,
                Pos2::new(key.x + 32.0, y),
                Align2::LEFT_CENTER,
                name,
                TEXT,
            );
        }
    }
}

impl eframe::App for TunnelApp {
    fn ui(&mut self, ui: &mut Ui, _frame: &mut eframe::Frame) {
        egui::Panel::left("controls")
            .exact_size(300.0)
            .resizable(false)
            .show(ui, |ui| {
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        self.controls(ui);
                    });
            });

        // Live update: solve again when an input changes.
        self.inputs = self.inputs.normalized();
        if self.solved_for != Some(self.inputs) {
            self.result = self.analyzer.compute(&self.inputs);
            self.solved_for = Some(self.inputs);
        }

        egui::CentralPanel::default()
            .frame(egui::Frame::NONE.fill(BACKGROUND))
            .show(ui, |ui| {
                ui.spacing_mut().item_spacing = vec2(0.0, 0.0);
                let size = ui.available_size();
                let top = (size.y * 0.52).round();
                self.airfoil_view(ui, vec2(size.x, top));
                let (line, _) = ui.allocate_exact_size(vec2(size.x, 1.0), Sense::hover());
                ui.painter().rect_filled(line, 0.0, GRID);
                self.pressure_plot(ui, vec2(size.x, size.y - top - 1.0));
            });
    }
}

fn section(ui: &mut Ui, title: &str) {
    ui.label(RichText::new(title).size(11.0).color(TEXT_DIM).strong());
    ui.add_space(4.0);
}

fn readout(ui: &mut Ui, name: &str, value: &str) {
    ui.label(RichText::new(name).size(12.0).color(TEXT_DIM));
    ui.label(RichText::new(value).size(26.0).monospace().color(TEXT));
}

fn label(painter: &egui::Painter, pos: Pos2, anchor: Align2, text: &str, color: Color32) {
    painter.text(pos, anchor, text, FontId::proportional(12.0), color);
}

/// Text for a chord position, in percent.
fn position_text(x: Option<f64>) -> String {
    match x {
        Some(x) => format!("{:.0} %", x * 100.0),
        None => "none".to_owned(),
    }
}

/// Text that tells where the flow separates, if it does.
fn separation_text(friction: &InteractionResult) -> Option<String> {
    let mut parts = Vec::new();
    for (name, x) in [
        ("upper", friction.separation_upper),
        ("lower", friction.separation_lower),
    ] {
        if let Some(x) = x {
            parts.push(format!("{:.0} % of chord on the {name} surface", x * 100.0));
        }
    }
    (!parts.is_empty()).then(|| format!("Flow separation from {}.", parts.join(" and from ")))
}

fn legend(painter: &egui::Painter, rect: Rect) {
    let x = rect.right() - 250.0;
    let y = rect.top() + 30.0;
    for (row, (name, color)) in [
        ("Suction (below free-stream pressure)", view::SUCTION),
        ("Pressure (above free-stream pressure)", view::PRESSURE),
        ("T: start of turbulent flow", TRANSITION),
        ("S: flow separation", SEPARATION),
    ]
    .iter()
    .enumerate()
    {
        let y = y + row as f32 * 18.0;
        painter.line_segment(
            [Pos2::new(x, y), Pos2::new(x + 16.0, y)],
            Stroke::new(2.0, *color),
        );
        label(
            painter,
            Pos2::new(x + 24.0, y),
            Align2::LEFT_CENTER,
            name,
            TEXT_DIM,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn inputs() -> Inputs {
        Inputs::default()
    }

    #[test]
    fn code_has_four_digits() {
        assert_eq!(inputs().code(), "2412");
        let thin = Inputs {
            camber: 0,
            camber_pos: 0,
            thickness: 6,
            ..inputs()
        };
        assert_eq!(thin.code(), "0006");
    }

    #[test]
    fn normalized_inputs_give_a_flow_result() {
        // Each shape that the controls can make must give a friction-free
        // result. The friction result can be absent (no solution).
        let mut analyzer = Analyzer::default();
        for camber in 0..=9 {
            for camber_pos in [0, 2, 5, 9] {
                for thickness in [4, 12, 30] {
                    let raw = Inputs {
                        camber,
                        camber_pos,
                        thickness,
                        alpha: 3.0,
                        panels: 60,
                        ..inputs()
                    };
                    let done = analyzer.compute(&raw.normalized());
                    assert!(done.is_ok(), "{raw:?}: {:?}", done.err());
                    let done = done.unwrap();
                    assert!(done.friction.is_some() || done.friction_note.is_some());
                    assert_eq!(done.cp().len(), 60);
                }
            }
        }
    }

    #[test]
    fn normalization_repairs_bad_values() {
        let bad = Inputs {
            camber: 50,
            camber_pos: 0,
            thickness: 0,
            alpha: f64::NAN,
            panels: 41,
            reynolds: f64::NAN,
            trip: true,
            trip_x: 7.0,
            n_crit: f64::NAN,
        };
        let fixed = bad.normalized();
        assert_eq!((fixed.camber, fixed.camber_pos, fixed.thickness), (9, 1, 1));
        assert_eq!(fixed.panels, 60);
        assert_eq!(fixed.alpha, 0.0);
        assert_eq!(fixed.reynolds, 3.0e6);
        assert_eq!(fixed.trip_x, 0.9);
        assert_eq!(fixed.n_crit, 9.0);
        let flat = Inputs {
            camber: 0,
            camber_pos: 7,
            ..inputs()
        }
        .normalized();
        assert_eq!(flat.camber_pos, 0);
    }

    #[test]
    fn start_values_give_a_friction_result() {
        let mut analyzer = Analyzer::default();
        let done = analyzer.compute(&inputs()).unwrap();
        let friction = done.friction.expect("friction result");
        assert!(friction.cd > 0.0 && friction.cl < friction.cl_inviscid);
        assert!(done.friction_note.is_none());
        assert!(position_text(friction.transition_upper).ends_with('%'));
        // A trip strip moves the start of turbulent flow forward.
        let tripped = analyzer
            .compute(&Inputs {
                trip: true,
                trip_x: 0.05,
                ..inputs()
            })
            .unwrap();
        let tripped = tripped.friction.expect("friction result with trip strip");
        assert!(tripped.transition_lower.unwrap() < friction.transition_lower.unwrap());
        assert!(tripped.cd > friction.cd);
    }

    #[test]
    fn no_solution_gives_a_message_and_the_friction_free_lift() {
        let mut analyzer = Analyzer::default();
        let stalled = Inputs {
            camber: 0,
            camber_pos: 0,
            thickness: 6,
            alpha: 14.0,
            reynolds: 5.0e5,
            ..inputs()
        };
        let done = analyzer.compute(&stalled).unwrap();
        assert!(done.friction.is_none());
        assert!(done.friction_note.unwrap().contains("no solution"));
        assert!(done.solution.cl > 1.0);
        // The next result at a normal condition is available again.
        let done = analyzer
            .compute(&Inputs {
                thickness: 12,
                alpha: 2.0,
                reynolds: 3.0e6,
                ..stalled
            })
            .unwrap();
        assert!(done.friction.is_some());
    }

    #[test]
    fn separation_text_names_the_surface() {
        let mut analyzer = Analyzer::default();
        let done = analyzer
            .compute(&Inputs {
                camber: 0,
                camber_pos: 0,
                alpha: 16.0,
                reynolds: 6.0e6,
                ..inputs()
            })
            .unwrap();
        if let Some(friction) = &done.friction {
            if friction.separation_upper.is_some() {
                assert!(separation_text(friction).unwrap().contains("upper surface"));
            }
        }
        assert_eq!(position_text(None), "none");
        assert_eq!(position_text(Some(0.25)), "25 %");
    }

    #[test]
    fn a_small_input_change_is_fast_enough_for_live_update() {
        // Goal: new results in less than 0.1 second (optimized build).
        let mut analyzer = Analyzer::default();
        analyzer.compute(&inputs()).unwrap();
        let done = analyzer
            .compute(&Inputs {
                alpha: 4.1,
                ..inputs()
            })
            .unwrap();
        assert!(done.friction.is_some());
        println!("{:.0} ms", done.solve_ms);
        if !cfg!(debug_assertions) {
            assert!(done.solve_ms < 100.0, "{} ms", done.solve_ms);
        }
    }
}
