//! 2D airfoil flow solver.
//!
//! Stage 1: inviscid (no friction) linear-vortex panel method.
//! Stage 3a: one-way boundary layer calculation (friction drag, transition).

pub mod boundary_layer;
pub mod geometry;
pub mod linalg;
pub mod panel;
pub mod viscous;

pub use geometry::{joukowski, naca4, Airfoil};
pub use panel::{solve, Solution};
pub use viscous::{analyze, ViscousResult};
