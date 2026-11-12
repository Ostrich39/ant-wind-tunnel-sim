//! 2D airfoil flow solver.
//!
//! Stage 1: inviscid (no friction) linear-vortex panel method.
//! Stage 3a: one-way boundary layer calculation (friction drag, transition).
//! Stage 3b: two-way calculation (Drela's method): the boundary layer and
//! the outer flow are solved together.

pub mod boundary_layer;
pub mod closure;
pub mod geometry;
pub mod ibl;
pub mod interaction;
pub mod linalg;
pub mod panel;
pub mod viscous;
pub mod wake;

pub use geometry::{joukowski, naca4, naca4_for_friction, Airfoil};
pub use interaction::{FrictionSystem, InteractionResult, Options};
pub use panel::{solve, PanelSystem, Solution};
pub use viscous::{analyze, ViscousResult};
