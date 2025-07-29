use crate::{CriticalState, R};

/// The default and recommended equation of state of this library.
pub type DefaultEos = PengRobinson;

pub trait EquationOfState {
    /// Compute the molecular attraction parameter of the EoS, aka. the A parameter.
    ///
    /// # Arguments
    ///  * `cs` - The critical state of the molecule
    ///  * `w`  - The acentric factor of the molecule (no dimension)
    ///  * `t`  - The temperature of the gas, in K
    fn a(cs: &CriticalState, w: f64, t: f64) -> f64;

    /// Compute the molecular volume parameter of the EoS, aka. the B parameter.
    ///
    /// # Arguments
    ///  * `cs` - The critical state of the molecule
    fn b(cs: &CriticalState) -> f64;

    /// Modification of the molecular attraction parameter of the EoS, aka. the C parameter.
    ///
    /// # Arguments
    ///  * `cs` - The critical state of the molecule
    fn c(cs: &CriticalState) -> f64;

    /// Compute the gas pressure for given parameters and state.
    ///
    /// # Arguments
    ///  * `a`  - The molecular attraction parameter
    ///  * `b`  - The molecular volume parameter
    ///  * `c`  - The modified molecular attraction parameter
    ///  * `vm` - The molar volume of the gas, in m^3/mol
    ///  * `t`  - The temperature of the gas, in K
    fn pressure(a: f64, b: f64, c: f64, vm: f64, t: f64) -> f64;

    /// The Z polyn [a3, a2, a1, a0] such as `a3*Z^3 + a2*Z^2 + a1*Z + a0 = 0`
    ///
    /// # Arguments
    ///  * `a` - The molecular attraction parameter
    ///  * `b` - The molecular volume parameter
    ///  * `c` - The modified molecular attraction parameter
    ///  * `p` - The pressure of the gas, in Pa
    ///  * `t` - The temperature of the gas, in K
    fn z_polyn(a: f64, b: f64, c: f64, p: f64, t: f64) -> [f64; 4];
}

/// The ideal gas law
pub enum IdealGas {}

impl EquationOfState for IdealGas {
    fn a(_cs: &CriticalState, _w: f64, _t: f64) -> f64 {
        0.0
    }

    fn b(_cs: &CriticalState) -> f64 {
        0.0
    }

    fn c(_cs: &CriticalState) -> f64 {
        0.0
    }

    fn pressure(_a: f64, _b: f64, _c: f64, vm: f64, t: f64) -> f64 {
        R * t / vm
    }

    fn z_polyn(_a: f64, _b: f64, _c: f64, _p: f64, _t: f64) -> [f64; 4] {
        [0.0, 0.0, 1.0, -1.0]
    }
}

/// The Van der Waals equation of state
pub enum VanDerWaals {}

impl EquationOfState for VanDerWaals {
    fn a(cs: &CriticalState, _w: f64, _t: f64) -> f64 {
        return 27.0 * R * R * cs.t * cs.t / (64.0 * cs.p);
    }

    fn b(cs: &CriticalState) -> f64 {
        return R * cs.t / (8.0 * cs.p);
    }

    fn c(_cs: &CriticalState) -> f64 {
        0.0
    }

    fn pressure(a: f64, b: f64, _c: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (vm * vm)
    }

    fn z_polyn(a: f64, b: f64, _c: f64, p: f64, t: f64) -> [f64; 4] {
        let a = a * p / (R * R * t * t);
        let b = b * p / (R * t);

        let a3 = 1f64;
        let a2 = -b - 1f64;
        let a1 = a;
        let a0 = -a * b;

        [a3, a2, a1, a0]
    }
}

/// The Redlich-Kwong equation of state
pub enum RedlichKwong {}

impl EquationOfState for RedlichKwong {
    fn a(cs: &CriticalState, _w: f64, _t: f64) -> f64 {
        0.42748023 * R * R * cs.t.powf(2.5) / cs.p
    }

    fn b(cs: &CriticalState) -> f64 {
        0.08664035 * R * cs.t / cs.p
    }

    fn c(_cs: &CriticalState) -> f64 {
        0.0
    }

    fn pressure(a: f64, b: f64, _c: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (t.sqrt() * vm * (vm + b))
    }

    fn z_polyn(a: f64, b: f64, _c: f64, p: f64, t: f64) -> [f64; 4] {
        let a = a * p / (R * R * t.powf(2.5));
        let b = b * p / (R * t);

        let a3 = 1f64;
        let a2 = -1f64;
        let a1 = a - b * b - b;
        let a0 = -a * b;

        [a3, a2, a1, a0]
    }
}

/// The Soave-Redlich-Kwong equation of state
pub enum SoaveRedlichKwong {}

impl EquationOfState for SoaveRedlichKwong {
    fn a(cs: &CriticalState, w: f64, t: f64) -> f64 {
        let m = 0.48 + 1.574 * w - 0.176 * w * w;
        let sq_a = 1f64 + m * (1f64 - (t / cs.t).sqrt());
        let alpha = sq_a * sq_a;

        alpha * 0.42748023 * R * R * cs.t * cs.t / cs.p
    }

    fn b(cs: &CriticalState) -> f64 {
        return 0.08664035 * R * cs.t / cs.p;
    }

    fn c(_cs: &CriticalState) -> f64 {
        0.0
    }

    fn pressure(a: f64, b: f64, _c: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (vm * (vm + b))
    }

    fn z_polyn(a: f64, b: f64, _c: f64, p: f64, t: f64) -> [f64; 4] {
        let a = a * p / (R * R * t * t);
        let b = b * p / (R * t);

        let a3 = 1f64;
        let a2 = -1f64;
        let a1 = a - b * b - b;
        let a0 = -a * b;

        [a3, a2, a1, a0]
    }
}

/// The Peng-Robinson equation of state
pub enum PengRobinson {}

impl EquationOfState for PengRobinson {
    fn a(cs: &CriticalState, w: f64, t: f64) -> f64 {
        let m = if w <= 0.491 {
            0.37464 + 1.56226 * w - 0.26992 * w * w
        } else {
            0.379642 + 1.487503 * w - 0.164423 * w * w - 0.016666 * w * w * w
        };
        let sq_a = 1f64 + m * (1f64 - (t / cs.t).sqrt());
        let alpha = sq_a * sq_a;
        
        alpha * 0.4572355289213821 * R * R * cs.t * cs.t / cs.p
    }

    fn b(cs: &CriticalState) -> f64 {
        0.07779607390388844 * R * cs.t / cs.p
    }

    fn c(_cs: &CriticalState) -> f64 {
        0.0
    }

    fn pressure(a: f64, b: f64, _c: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (vm * vm + 2.0 * b * vm - b * b)
    }

    fn z_polyn(a: f64, b: f64, _c: f64, p: f64, t: f64) -> [f64; 4] {
        let a = a * p / (R * R * t * t);
        let b = b * p / (R * t);

        let a3 = 1f64;
        let a2 = b - 1f64;
        let a1 = -3f64 * b * b - 2f64 * b + a;
        let a0 = b * b * b + b * b - a * b;

        [a3, a2, a1, a0]
    }
}

pub enum PatelTejaValderrama {} 

impl EquationOfState for PatelTejaValderrama {
    fn a(cs: &CriticalState, w: f64, t: f64) -> f64 {
        let m = 0.46283 + 3.58230 * w * cs.z() + 8.19417 * w * w * cs.z() * cs.z();
        let sq_a = 1f64 + m * (1f64 - (t / cs.t).sqrt());
        let alpha = sq_a * sq_a;
        let omega_a = 0.66121 - 0.76105 * cs.z();
        omega_a * alpha * R * R * cs.t * cs.t / cs.p
    }

    fn b(cs: &CriticalState) -> f64 {
        let omega_b = 0.02207 + 0.20868 * cs.z();
        omega_b * R * cs.t / cs.p
    }

    fn c(cs: &CriticalState) -> f64 {
        let omega_c = 0.57765 - 1.78080 * cs.z();
        omega_c * R * cs.t / cs.p
    }

    fn pressure(a: f64, b: f64, c: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (vm * (vm + b) + c * (vm - b))
    }

    fn z_polyn(a: f64, b: f64, c: f64, p: f64, t: f64) -> [f64; 4] {
        let a = a * p / (R * R * t * t);
        let b = b * p / (R * t);
        let c = c * p / (R * t);

        let a3 = 1f64;
        let a2 = c - 1f64;
        let a1 = -2f64 * b * c - b*b -b -c + a;
        let a0 = b * b * c + b * c - a * b;

        [a3, a2, a1, a0]
    }
}

/// An equation of state determined at runtime
#[derive(Debug, Clone, Copy)]
pub enum Eos {
    /// The ideal gas law
    IdealGas,
    /// The Van der Waals equation of state
    VanDerWaals,
    /// The Redlich-Kwong equation of state
    RedlichKwong,
    /// The Soave-Redlich-Kwong equation of state
    SoaveRedlichKwong,
    /// The Peng-Robinson equation of state
    PengRobinson,
    /// The Patel-Teja-Valderrama equation of state
    PatelTejaValderrama,
}

impl Default for Eos {
    fn default() -> Self {
        Eos::PengRobinson
    }
}
