use crate::R;

pub trait EquationOfState {
    /// Compute the molecular attraction parameter of the EoS, aka. the A parameter.
    ///
    /// # Arguments
    ///  * `pc` - The critical pressure of the molecule, in Pa
    ///  * `tc` - The critical temperature of the molecule, in K
    ///  * `w`  - The acentric factor of the molecule (no dimension)
    ///  * `t`  - The temperature of the gas, in K
    fn a(pc: f64, tc: f64, w: f64, t: f64) -> f64;

    /// Compute the molecular volume parameter of the EoS, aka. the B parameter.
    ///
    /// # Arguments
    ///  * `pc` - The critical pressure of the molecule, in Pa
    ///  * `tc` - The critical temperature of the molecule, in K
    fn b(pc: f64, tc: f64) -> f64;

    /// Compute the gas pressure for given parameters.
    fn pressure(a: f64, b: f64, vm: f64, t: f64) -> f64;

    /// The Z polyn [a3, a2, a1, a0] such as `a3*Z^3 + a2*Z^2 + a1*Z + a0 = 0`
    fn z_polyn(a: f64, b: f64, p: f64, t: f64) -> [f64; 4];
}

/// The ideal gas law
pub enum IdealGas {}

impl EquationOfState for IdealGas {
    fn a(_pc: f64, _tc: f64, _w: f64, _t: f64) -> f64 {
        0.0
    }

    fn b(_pc: f64, _tc: f64) -> f64 {
        0.0
    }

    fn pressure(_a: f64, _b: f64, vm: f64, t: f64) -> f64 {
        R * t / vm
    }

    fn z_polyn(_a: f64, _b: f64, _p: f64, _t: f64) -> [f64; 4] {
        [0.0, 0.0, 1.0, -1.0]
    }
}

/// The Van der Waals equation of state
pub enum VanDerWaals {}

impl EquationOfState for VanDerWaals {
    fn a(pc: f64, tc: f64, _w: f64, _t: f64) -> f64 {
        return 27.0 * R * R * tc * tc / (64.0 * pc);
    }

    fn b(pc: f64, tc: f64) -> f64 {
        return R * tc / (8.0 * pc);
    }

    fn pressure(a: f64, b: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (vm * vm)
    }

    fn z_polyn(a: f64, b: f64, p: f64, t: f64) -> [f64; 4] {
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
    fn a(pc: f64, tc: f64, _w: f64, _t: f64) -> f64 {
        return 0.42748023 * R * R * tc.powf(2.5) / pc;
    }

    fn b(pc: f64, tc: f64) -> f64 {
        return 0.08664035 * R * tc / pc;
    }

    fn pressure(a: f64, b: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (t.sqrt() * vm * (vm + b))
    }

    fn z_polyn(a: f64, b: f64, p: f64, t: f64) -> [f64; 4] {
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
    fn a(pc: f64, tc: f64, w: f64, t: f64) -> f64 {
        let m = 0.48 + 1.574 * w - 0.176 * w * w;
        let alpha = 1f64 + m * (1f64 - (t / tc).sqrt());
        let alpha = alpha * alpha;
        return alpha * 0.42748023 * R * R * tc * tc / pc;
    }

    fn b(pc: f64, tc: f64) -> f64 {
        return 0.08664035 * R * tc / pc;
    }

    fn pressure(a: f64, b: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (vm * (vm + b))
    }

    fn z_polyn(a: f64, b: f64, p: f64, t: f64) -> [f64; 4] {
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
    fn a(pc: f64, tc: f64, w: f64, t: f64) -> f64 {
        let m = if w <= 0.491 {
            0.37464 + 1.56226 * w - 0.26992 * w * w
        } else {
            0.379642 + 1.487503 * w - 0.164423 * w * w - 0.016666 * w * w * w
        };
        let alpha = 1f64 + m * (1f64 - (t / tc).sqrt());
        let alpha = alpha * alpha;
        return alpha * 0.4572355289213821 * R * R * tc * tc / pc;
    }

    fn b(pc: f64, tc: f64) -> f64 {
        return 0.07779607390388844 * R * tc / pc;
    }

    fn pressure(a: f64, b: f64, vm: f64, t: f64) -> f64 {
        R * t / (vm - b) - a / (vm * vm + 2.0 * b * vm - b * b)
    }

    fn z_polyn(a: f64, b: f64, p: f64, t: f64) -> [f64; 4] {
        let a = a * p / (R * R * t * t);
        let b = b * p / (R * t);

        let a3 = 1f64;
        let a2 = b - 1f64;
        let a1 = -3f64 * b * b - 2f64 * b + a;
        let a0 = b * b * b + b * b - a * b;

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
}
