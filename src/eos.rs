use std::borrow::Borrow;

use crate::{CriticalState, R};

/// The default and recommended equation of state of this library.
pub type DefaultEos = PengRobinson;

/// The A and B parameters of an equation of state.
/// Most equations of state use these parameters to compute the pressure of a gas.
#[derive(Debug, Clone, Copy)]
pub struct AbParams {
    /// The molecular attraction parameter
    pub a: f64,
    /// The molecular volume parameter
    pub b: f64,
}

/// The A, B and C parameters of an equation of state.
/// Some equations of state use this additional C parameters to provide better accuracy in some cases.
#[derive(Debug, Clone, Copy)]
pub struct AbcParams {
    /// The molecular attraction parameter
    pub a: f64,
    /// The molecular volume parameter
    pub b: f64,
    /// The additional parameter
    pub c: f64,
}

/// Mixing rules for equations of state parameters.
pub trait MixingRules {
    fn mix<P>(mixture_params: P) -> Self
    where
        P: IntoIterator + Clone,
        P::Item: Borrow<(f64, Self)>;
}

impl MixingRules for () {
    fn mix<P>(_mixture_params: P) -> Self
    where
        P: IntoIterator + Clone,
        P::Item: Borrow<(f64, Self)>,
    {
        ()
    }
}

/// Mixing rules for equations of state parameters that use the A and B parameters.
impl MixingRules for AbParams {
    fn mix<P>(mixture_params: P) -> Self
    where
        P: IntoIterator + Clone,
        P::Item: Borrow<(f64, Self)>,
    {
        let mut a = 0.0;
        let mut b = 0.0;
        for params in mixture_params.clone() {
            let (fi, pi) = params.borrow();
            for params in mixture_params.clone() {
                let (fj, pj) = params.borrow();
                a += fi * fj * (pi.a * pj.a).sqrt();
            }
            b += fi * pi.b;
        }
        AbParams { a, b }
    }
}

/// Mixing rules for equations of state parameters that use the A, B and C parameters.
impl MixingRules for AbcParams {
    fn mix<P>(mixture_params: P) -> Self
    where
        P: IntoIterator + Clone,
        P::Item: Borrow<(f64, Self)>,
    {
        let mut a = 0.0;
        let mut b = 0.0;
        let mut c = 0.0;
        for params in mixture_params.clone() {
            let (fi, pi) = params.borrow();
            for params in mixture_params.clone() {
                let (fj, pj) = params.borrow();
                a += fi * fj * (pi.a * pj.a).sqrt();
            }
            b += fi * pi.b;
            c += fi * pi.c;
        }
        AbcParams { a, b, c }
    }
}

pub trait EquationOfState {
    /// The parameters of the equation of state
    type Params: MixingRules;

    /// Compute the parameters of the equation of state.
    ///
    /// # Arguments
    ///  * `cs` - The critical state of the molecule
    ///  * `w`  - The acentric factor of the molecule (no dimension)
    ///  * `t`  - The temperature of the gas, in K
    fn params(cs: &CriticalState, w: f64, t: f64) -> Self::Params;

    /// Compute the gas pressure for given parameters and state.
    ///
    /// # Arguments
    ///  * `params` - The equation parameters
    ///  * `vm`     - The molar volume of the gas, in m^3/mol
    ///  * `t`      - The temperature of the gas, in K
    fn pressure(params: &Self::Params, vm: f64, t: f64) -> f64;

    /// The Z polyn [a3, a2, a1, a0] such as `a3*Z^3 + a2*Z^2 + a1*Z + a0 = 0`
    ///
    /// # Arguments
    ///  * `params` - The equation parameters
    ///  * `p`      - The pressure of the gas, in Pa
    ///  * `t`      - The temperature of the gas, in K
    fn z_polyn(params: &Self::Params, p: f64, t: f64) -> [f64; 4];
}

/// The ideal gas law
pub enum IdealGas {}

impl EquationOfState for IdealGas {
    type Params = ();
    fn params(_cs: &CriticalState, _w: f64, _t: f64) -> Self::Params {
        // No parameters needed for the ideal gas law
        ()
    }

    fn pressure(_params: &Self::Params, vm: f64, t: f64) -> f64 {
        R * t / vm
    }

    fn z_polyn(_params: &Self::Params, _p: f64, _t: f64) -> [f64; 4] {
        // Z = 1
        [0.0, 0.0, 1.0, -1.0]
    }
}

/// The Van der Waals equation of state
pub enum VanDerWaals {}

impl EquationOfState for VanDerWaals {
    type Params = AbParams;

    fn params(cs: &CriticalState, _w: f64, _t: f64) -> Self::Params {
        let a = 27.0 * R * R * cs.t * cs.t / (64.0 * cs.p);
        let b = R * cs.t / (8.0 * cs.p);
        AbParams { a, b }
    }

    fn pressure(params: &Self::Params, vm: f64, t: f64) -> f64 {
        let AbParams { a, b } = *params;
        R * t / (vm - b) - a / (vm * vm)
    }

    fn z_polyn(params: &Self::Params, p: f64, t: f64) -> [f64; 4] {
        let a = params.a * p / (R * R * t * t);
        let b = params.b * p / (R * t);

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
    type Params = AbParams;

    fn params(cs: &CriticalState, _w: f64, _t: f64) -> Self::Params {
        let a = 0.42748023 * R * R * cs.t.powf(2.5) / cs.p;
        let b = 0.08664035 * R * cs.t / cs.p;

        AbParams { a, b }
    }

    fn pressure(params: &Self::Params, vm: f64, t: f64) -> f64 {
        let AbParams { a, b } = *params;
        R * t / (vm - b) - a / (t.sqrt() * vm * (vm + b))
    }

    fn z_polyn(params: &Self::Params, p: f64, t: f64) -> [f64; 4] {
        let a = params.a * p / (R * R * t.powf(2.5));
        let b = params.b * p / (R * t);

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
    type Params = AbParams;

    fn params(cs: &CriticalState, w: f64, t: f64) -> Self::Params {
        let m = 0.48 + 1.574 * w - 0.176 * w * w;
        let sq_a = 1f64 + m * (1f64 - (t / cs.t).sqrt());
        let alpha = sq_a * sq_a;

        let a = alpha * 0.42748023 * R * R * cs.t * cs.t / cs.p;
        let b = 0.08664035 * R * cs.t / cs.p;

        AbParams { a, b }
    }

    fn pressure(params: &Self::Params, vm: f64, t: f64) -> f64 {
        let AbParams { a, b } = *params;
        R * t / (vm - b) - a / (vm * (vm + b))
    }

    fn z_polyn(params: &Self::Params, p: f64, t: f64) -> [f64; 4] {
        let a = params.a * p / (R * R * t * t);
        let b = params.b * p / (R * t);

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
    type Params = AbParams;

    fn params(cs: &CriticalState, w: f64, t: f64) -> Self::Params {
        let m = if w <= 0.491 {
            0.37464 + 1.56226 * w - 0.26992 * w * w
        } else {
            0.379642 + 1.487503 * w - 0.164423 * w * w - 0.016666 * w * w * w
        };
        let sq_a = 1f64 + m * (1f64 - (t / cs.t).sqrt());
        let alpha = sq_a * sq_a;

        let a = alpha * 0.4572355289213821 * R * R * cs.t * cs.t / cs.p;
        let b = 0.07779607390388844 * R * cs.t / cs.p;

        AbParams { a, b }
    }

    fn pressure(params: &Self::Params, vm: f64, t: f64) -> f64 {
        let AbParams { a, b } = *params;
        R * t / (vm - b) - a / (vm * vm + 2.0 * b * vm - b * b)
    }

    fn z_polyn(params: &Self::Params, p: f64, t: f64) -> [f64; 4] {
        let a = params.a * p / (R * R * t * t);
        let b = params.b * p / (R * t);

        let a3 = 1f64;
        let a2 = b - 1f64;
        let a1 = -3f64 * b * b - 2f64 * b + a;
        let a0 = b * b * b + b * b - a * b;

        [a3, a2, a1, a0]
    }
}

pub enum PatelTejaValderrama {}

impl EquationOfState for PatelTejaValderrama {
    type Params = AbcParams;

    fn params(cs: &CriticalState, w: f64, t: f64) -> Self::Params {
        let zc = cs.z();

        let m = 0.46283 + 3.58230 * w * zc + 8.19417 * w * w * zc * zc;
        let sq_a = 1f64 + m * (1f64 - (t / cs.t).sqrt());
        let alpha = sq_a * sq_a;
        let omega_a = 0.66121 - 0.76105 * zc;
        let a = omega_a * alpha * R * R * cs.t * cs.t / cs.p;

        let omega_b = 0.02207 + 0.20868 * zc;
        let b = omega_b * R * cs.t / cs.p;

        let omega_c = 0.57765 - 1.78080 * zc;
        let c = omega_c * R * cs.t / cs.p;

        AbcParams { a, b, c }
    }

    fn pressure(params: &Self::Params, vm: f64, t: f64) -> f64 {
        let AbcParams { a, b, c } = *params;
        R * t / (vm - b) - a / (vm * (vm + b) + c * (vm - b))
    }

    fn z_polyn(params: &Self::Params, p: f64, t: f64) -> [f64; 4] {
        let a = params.a * p / (R * R * t * t);
        let b = params.b * p / (R * t);
        let c = params.c * p / (R * t);

        let a3 = 1f64;
        let a2 = c - 1f64;
        let a1 = -2f64 * b * c - b * b - b - c + a;
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
