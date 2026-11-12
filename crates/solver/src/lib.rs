//! 2D airfoil flow solver.
//!
//! Stage 1: inviscid (no friction) linear-vortex panel method.
//! Stage 3a: one-way boundary layer calculation (friction drag, transition).
//! Stage 3b (experimental): two-way calculation with displacement sources.

pub mod boundary_layer;
pub mod closure;
pub mod coupled;
pub mod geometry;
pub mod ibl;
pub mod linalg;
pub mod panel;
pub mod viscous;

pub use coupled::{solve_coupled, CoupledResult};
pub use geometry::{joukowski, naca4, Airfoil};
pub use panel::{solve, PanelSystem, Solution};
pub use viscous::{analyze, ViscousResult};
