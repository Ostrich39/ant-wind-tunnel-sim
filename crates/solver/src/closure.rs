//! Closure relations for the two-equation integral boundary layer method
//! of Drela (incompressible form).
//!
//! The integral equations have more unknowns than equations. These
//! relations supply the missing quantities as functions of the shape
//! factor H and the momentum-thickness Reynolds number.
//!
//! Sources: REFERENCES.md, items 2.1, 2.2, 2.4. The equations were written
//! from the papers and compared with the pages of item 2.4 (equations 45
//! to 79 and Table 2). No program code was copied.

/// Constants of the turbulent closure.
pub const G_A: f64 = 6.7;
pub const G_B: f64 = 0.75;
pub const G_C: f64 = 18.0;
/// Shear lag constant.
pub const K_LAG: f64 = 5.6;
/// Constants for the shear stress at the start of turbulent flow.
pub const C_TAU_START: f64 = 1.8;
pub const E_TAU_START: f64 = 3.3;
/// Dissipation length ratio in the wake.
pub const ETA_D_WAKE: f64 = 0.9;
/// Critical amplification factor for transition (clean wind tunnel).
pub const N_CRIT: f64 = 9.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Flow {
    Laminar,
    Turbulent,
    Wake,
}

/// Laminar energy shape factor H*(H).
pub fn hs_laminar(h: f64) -> f64 {
    let t = h - 4.35;
    if h < 4.35 {
        0.0111 * t * t / (h + 1.0) - 0.0278 * t.powi(3) / (h + 1.0) + 1.528
            - 0.0002 * (t * h).powi(2)
    } else {
        0.015 * t * t / h + 1.528
    }
}

/// Laminar skin friction coefficient.
pub fn cf_laminar(h: f64, re_theta: f64) -> f64 {
    if h < 5.5 {
        (0.0727 * (5.5 - h).powi(3) / (h + 1.0) - 0.07) / re_theta
    } else {
        (0.015 * (1.0 - 1.0 / (h - 4.5)).powi(2) - 0.07) / re_theta
    }
}

/// Laminar dissipation, in the form 2 cD / H*.
pub fn di_laminar(h: f64, re_theta: f64) -> f64 {
    if h < 4.0 {
        (0.00205 * (4.0 - h).powf(5.5) + 0.207) / re_theta
    } else {
        let t = h - 4.0;
        (-0.0016 * t * t / (1.0 + 0.02 * t * t) + 0.207) / re_theta
    }
}

/// Turbulent energy shape factor H*(H, Re_theta).
pub fn hs_turbulent(h: f64, re_theta: f64) -> f64 {
    let h0 = (3.0 + 400.0 / re_theta).min(4.0);
    let re = re_theta.max(200.0);
    if h < h0 {
        let hr = (h0 - h) / (h0 - 1.0);
        1.5 + 4.0 / re + 1.5 * (0.5 - 4.0 / re) * hr * hr / (h + 0.5)
    } else {
        let log_re = re.ln();
        let t = h - h0 + 4.0 / log_re;
        1.5 + 4.0 / re + (h - h0).powi(2) * (0.007 * log_re / (t * t) + 0.015 / h)
    }
}

/// Turbulent skin friction coefficient.
pub fn cf_turbulent(h: f64, re_theta: f64) -> f64 {
    let b = re_theta.log10().max(1.303);
    // Smooth lower limit of the exponent for very large shape factors
    // (Fidkowski, equation 58).
    let mut a = -1.33 * h;
    if a < -17.0 {
        a = -20.0 + 3.0 * ((a + 17.0) / 3.0).exp();
    }
    0.3 * a.exp() * b.powf(-1.74 - 0.31 * h) + 0.00011 * ((4.0 - h / 0.875).tanh() - 1.0)
}

/// Slip velocity ratio at the edge of the wall layer.
pub fn slip_velocity(hs: f64, h: f64, flow: Flow) -> f64 {
    let us = 0.5 * hs * (1.0 - (h - 1.0) / (G_B * h));
    us.min(if flow == Flow::Wake { 0.99995 } else { 0.98 })
}

/// Shape factor term of the equilibrium relations.
fn hkc(h: f64, re_theta: f64, flow: Flow) -> f64 {
    let gc = if flow == Flow::Wake { 0.0 } else { G_C };
    (h - 1.0 - gc / re_theta).max(0.01)
}

/// Square root of the equilibrium shear stress coefficient.
pub fn ctau_eq_sqrt(h: f64, hs: f64, us: f64, re_theta: f64, flow: Flow) -> f64 {
    let c = hkc(h, re_theta, flow);
    (hs * (h - 1.0) * c * c / (2.0 * G_A * G_A * G_B * (1.0 - us) * h * h * h)).sqrt()
}

/// Boundary layer thickness.
pub fn delta(h: f64, theta: f64, delta_star: f64) -> f64 {
    ((3.15 + 1.72 / (h - 1.0)) * theta + delta_star).min(12.0 * theta)
}

/// Turbulent dissipation, in the form 2 cD / H*. `ctau` is the shear
/// stress coefficient (not its square root).
pub fn di_turbulent(h: f64, hs: f64, us: f64, cf: f64, ctau: f64, re_theta: f64) -> f64 {
    let wall = 0.5 * cf * us * (2.0 / hs) * 0.5 * (1.0 + ((h - 1.0) * re_theta.ln() / 2.1).tanh());
    let outer = ctau * (0.995 - us) * 2.0 / hs;
    let stress = 0.3 * (0.995 - us).powi(2) / (hs * re_theta);
    (wall + outer + stress).max(di_laminar(h, re_theta))
}

/// Wake dissipation, in the form 2 cD / H* (Fidkowski, equations 66 and 67).
///
/// The paper gives no extra factor for the two halves of the wake, and
/// this code follows the paper. A factor of 2 changes the drag by less than
/// 2% and the lift by less than 0.5% for the NACA 0012 (see GOALS.md).
pub fn di_wake(h: f64, hs: f64, us: f64, ctau: f64, re_theta: f64) -> f64 {
    let outer = ctau * (0.995 - us) * 2.0 / hs;
    let stress = 0.3 * (0.995 - us).powi(2) / (hs * re_theta);
    let laminar = 2.2 * (1.0 - 1.0 / h).powi(2) / (h * hs * re_theta);
    (outer + stress).max(laminar)
}

/// Growth rate of the amplification factor, d n / d xi (envelope e^n method).
pub fn amplification_rate(h: f64, theta: f64, re_theta: f64) -> f64 {
    let hm = 1.0 / (h - 1.0);
    let log_critical = 2.492 * hm.powf(0.43) + 0.7 * (1.0 + (14.0 * hm - 9.24).tanh());
    let s = (re_theta.log10() - (log_critical - 0.1)) / 0.2;
    let ramp = if s <= 0.0 {
        return 0.0;
    } else if s >= 1.0 {
        1.0
    } else {
        3.0 * s * s - 2.0 * s.powi(3)
    };
    let f = -0.05 + 2.7 * hm - 5.5 * hm * hm + 3.0 * hm.powi(3) + 0.1 * (-20.0 * hm).exp();
    let g = 0.028 * (h - 1.0) - 0.0345 * (-(3.87 * hm - 2.52).powi(2)).exp();
    ramp * f * g / theta
}

/// Square root of the shear stress coefficient at the start of turbulent flow.
pub fn ctau_start_sqrt(h: f64, ctau_eq_sqrt: f64) -> f64 {
    C_TAU_START * ctau_eq_sqrt * (-E_TAU_START / (h - 1.0)).exp()
}

/// Small extra growth rate of the amplification factor near the critical
/// value (Fidkowski, equation 76). It keeps the factor increasing after it
/// reaches the critical value.
pub fn amplification_extra(n: f64, n_crit: f64, theta: f64) -> f64 {
    0.001 * (1.0 + (5.0 * (n - n_crit)).tanh()) / theta
}
