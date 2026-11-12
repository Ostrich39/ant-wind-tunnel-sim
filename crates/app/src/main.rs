//! Desktop program: airfoil view, pressure plot, and live controls.

// Do not open a console window on Windows in release builds.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod app;
mod view;

use eframe::egui;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Ant Wind Tunnel")
            .with_inner_size([1200.0, 760.0])
            .with_min_inner_size([900.0, 600.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Ant Wind Tunnel",
        options,
        Box::new(|cc| Ok(Box::new(app::TunnelApp::new(cc)))),
    )
}
