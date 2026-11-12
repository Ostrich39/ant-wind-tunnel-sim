//! 2D airfoil flow solver.
//!
//! Stage 1: inviscid (no friction) linear-vortex panel method.
//! Later stages add the boundary layer (friction, drag, transition).

pub mod geometry;
pub mod linalg;
pub mod panel;

pub use geometry::{joukowski, naca4, Airfoil};
pub use panel::{solve, Solution};
