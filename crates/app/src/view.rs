//! Drawing math that does not depend on the window: coordinate mapping,
//! axis ticks, colors, and rotation. These functions have unit tests.

use eframe::egui::{Color32, Pos2, Rect};

/// Linear map from world coordinates to screen coordinates.
/// `y0` is at the bottom of the rectangle and `y1` is at the top.
#[derive(Debug, Clone, Copy)]
pub struct Mapper {
    pub rect: Rect,
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
}

impl Mapper {
    pub fn to_screen(self, x: f64, y: f64) -> Pos2 {
        let u = (x - self.x0) / (self.x1 - self.x0);
        let v = (y - self.y0) / (self.y1 - self.y0);
        Pos2::new(
            self.rect.left() + u as f32 * self.rect.width(),
            self.rect.bottom() - v as f32 * self.rect.height(),
        )
    }

    /// Screen pixels for one world unit in x.
    pub fn x_scale(self) -> f32 {
        self.rect.width() / (self.x1 - self.x0) as f32
    }
}

/// Makes a map with the same scale in x and y. The world box
/// `[x0, x1] x [y0, y1]` fits in `rect` and is in the center.
pub fn fit_uniform(rect: Rect, x0: f64, x1: f64, y0: f64, y1: f64) -> Mapper {
    let (w, h) = (x1 - x0, y1 - y0);
    let scale = (rect.width() as f64 / w).min(rect.height() as f64 / h);
    let (half_w, half_h) = (
        0.5 * rect.width() as f64 / scale,
        0.5 * rect.height() as f64 / scale,
    );
    let (cx, cy) = (0.5 * (x0 + x1), 0.5 * (y0 + y1));
    Mapper {
        rect,
        x0: cx - half_w,
        x1: cx + half_w,
        y0: cy - half_h,
        y1: cy + half_h,
    }
}

/// Axis tick values with a step of 1, 2, or 5 times a power of ten.
pub fn nice_ticks(lo: f64, hi: f64, target_count: usize) -> Vec<f64> {
    let (lo, hi) = if lo <= hi { (lo, hi) } else { (hi, lo) };
    let span = hi - lo;
    if !span.is_finite() || span <= 0.0 || target_count == 0 {
        return vec![lo];
    }
    let raw = span / target_count as f64;
    let power = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|m| m * power)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * power);
    let first = (lo / step).ceil() as i64;
    let last = (hi / step + 1e-9).floor() as i64;
    (first..=last).map(|i| i as f64 * step).collect()
}

/// Rotates a point about a pivot. A positive angle is counter-clockwise.
pub fn rotate(p: [f64; 2], pivot: [f64; 2], angle_rad: f64) -> [f64; 2] {
    let (s, c) = angle_rad.sin_cos();
    let (dx, dy) = (p[0] - pivot[0], p[1] - pivot[1]);
    [pivot[0] + c * dx - s * dy, pivot[1] + s * dx + c * dy]
}

pub const SUCTION: Color32 = Color32::from_rgb(86, 180, 255);
pub const PRESSURE: Color32 = Color32::from_rgb(255, 122, 89);

/// Color for a pressure coefficient: blue for suction (below the
/// free-stream pressure), orange-red for pressure (above it).
pub fn cp_color(cp: f64) -> Color32 {
    if cp < 0.0 {
        SUCTION
    } else {
        PRESSURE
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eframe::egui::pos2;

    fn rect() -> Rect {
        Rect::from_min_max(pos2(100.0, 50.0), pos2(500.0, 250.0))
    }

    #[test]
    fn mapper_maps_corners() {
        let m = Mapper {
            rect: rect(),
            x0: 0.0,
            x1: 1.0,
            y0: -2.0,
            y1: 2.0,
        };
        assert_eq!(m.to_screen(0.0, -2.0), pos2(100.0, 250.0));
        assert_eq!(m.to_screen(1.0, 2.0), pos2(500.0, 50.0));
        assert_eq!(m.to_screen(0.5, 0.0), pos2(300.0, 150.0));
    }

    #[test]
    fn mapper_supports_an_inverted_axis() {
        // Pressure plots put negative values at the top.
        let m = Mapper {
            rect: rect(),
            x0: 0.0,
            x1: 1.0,
            y0: 1.0,
            y1: -3.0,
        };
        assert_eq!(m.to_screen(0.0, 1.0).y, 250.0);
        assert_eq!(m.to_screen(0.0, -3.0).y, 50.0);
    }

    #[test]
    fn uniform_fit_has_equal_scales_and_is_centered() {
        let m = fit_uniform(rect(), 0.0, 1.0, -0.25, 0.25);
        let a = m.to_screen(0.0, 0.0);
        let b = m.to_screen(1.0, 0.0);
        let c = m.to_screen(0.0, 1.0);
        assert!(((b.x - a.x) - (a.y - c.y)).abs() < 1e-3);
        assert_eq!(m.to_screen(0.5, 0.0), rect().center());
        // The full world box is inside the rectangle.
        assert!(rect().expand(0.01).contains(m.to_screen(0.0, -0.25)));
        assert!(rect().expand(0.01).contains(m.to_screen(1.0, 0.25)));
    }

    #[test]
    fn ticks_use_round_steps() {
        assert_eq!(
            nice_ticks(0.0, 1.0, 5),
            vec![0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0]
        );
        assert_eq!(nice_ticks(-3.0, 1.0, 4), vec![-3.0, -2.0, -1.0, 0.0, 1.0]);
        assert_eq!(nice_ticks(1.0, -3.0, 4), vec![-3.0, -2.0, -1.0, 0.0, 1.0]);
        assert_eq!(nice_ticks(2.0, 2.0, 4), vec![2.0]);
    }

    #[test]
    fn rotation_by_a_quarter_turn() {
        let p = rotate([2.0, 1.0], [1.0, 1.0], std::f64::consts::FRAC_PI_2);
        assert!((p[0] - 1.0).abs() < 1e-12 && (p[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn color_shows_the_sign_of_the_pressure() {
        assert_eq!(cp_color(-0.5), SUCTION);
        assert_eq!(cp_color(0.5), PRESSURE);
    }
}
