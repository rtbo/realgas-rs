pub mod eos;
mod gas;
pub mod molecules;

use eos::{Eos, EquationOfState};
pub use gas::{Gas, Mixture, Molecule};

/// Universal gas constant in J/mol.K
pub const R: f64 = 8.31446262;

/// State trait of a gas.
/// All values here are intensive.
pub trait State {
    /// The molar mass of the gas, in kg/mol
    fn molar_mass(&self) -> f64;

    /// The molecular attraction parameter
    fn a<E: EquationOfState>(&self, t: f64) -> f64;

    /// The molecular volume parameter
    fn b<E: EquationOfState>(&self) -> f64;

    /// The modified molecular volume parameter
    fn c<E: EquationOfState>(&self) -> f64;

    /// Compute the pressure of the gas for the molar volume and temperature
    fn pressure<E: EquationOfState>(&self, vm: f64, t: f64) -> f64 {
        let a = self.a::<E>(t);
        let b = self.b::<E>();
        let c = self.c::<E>();
        E::pressure(a, b, c, vm, t)
    }

    /// Compute the compression factor Z such as Z = PV/RT
    ///
    /// Effectively resolves the cubic equation of state as a function of `p` and `t`.
    ///
    /// # Arguments
    ///  * `p` - The pressure of the gas, in Pa
    ///  * `t` - The temperature of the gas, in K
    ///
    /// # Panics
    /// This function will panic of no positive real root can be found, which is generally
    /// an indication that the parameters have physical non-sense.
    fn z<E: EquationOfState>(&self, p: f64, t: f64) -> f64 {
        use roots::Roots;

        let a = self.a::<E>(t);
        let b = self.b::<E>();
        let c = self.b::<E>();
        let [a3, a2, a1, a0] = E::z_polyn(a, b, c, p, t);
        let roots = roots::find_roots_cubic(a3, a2, a1, a0);
        let z = match roots {
            Roots::No([]) => None,
            Roots::One([r]) => Some(r),
            Roots::Two([r1, r2]) => Some(r1.max(r2)),
            Roots::Three([r1, r2, r3]) => Some(r1.max(r2).max(r3)),
            _ => unreachable!(),
        };
        z.filter(|&z| z > 0.0)
            .expect("Should have a found a positive real root")
    }

    /// Compute the molar volume the gas in m^3/mol
    fn molar_volume<E: EquationOfState>(&self, p: f64, t: f64) -> f64 {
        let z = self.z::<E>(p, t);
        z * R * t / p
    }

    /// Compute the specific mass of the gas in kg/m^3
    fn specific_mass<E: EquationOfState>(&self, p: f64, t: f64) -> f64 {
        let z = self.z::<E>(p, t);
        self.molar_mass() * p / (z * R * t)
    }
}

/// An helper trait to compute extensive state
pub trait ExtensiveState: State {
    /// Compute the amount of mols for given pressure, volume and temperature.
    ///
    /// # Panics
    /// This function can panic if the parameters have physical non-sense
    fn mols<E: EquationOfState>(&self, p: f64, v: f64, t: f64) -> f64 {
        let z = self.z::<E>(p, t);
        p * v / z / R / t
    }

    /// Compute the volume of the gas for given pressure, mols and temperature.
    ///
    /// # Panics
    /// This function can panic if the parameters have physical non-sense
    fn volume<E: EquationOfState>(&self, p: f64, n: f64, t: f64) -> f64 {
        let z = self.z::<E>(p, t);
        n * z * R * t / p
    }

    /// Compute the mass of the gas for given pressure, volume and temperature.
    fn mass<E: EquationOfState>(&self, p: f64, v: f64, t: f64) -> f64 {
        let n = self.mols::<E>(p, v, t);
        self.molar_mass() * n
    }
}

/// State trait of a gas for equation of state known at runtime.
/// All values here are intensive.
pub trait StateEos: State {
    /// The molecular attraction parameter
    fn a_eos(&self, eos: Eos, t: f64) -> f64 {
        match eos {
            Eos::IdealGas => self.a::<eos::IdealGas>(t),
            Eos::VanDerWaals => self.a::<eos::VanDerWaals>(t),
            Eos::RedlichKwong => self.a::<eos::RedlichKwong>(t),
            Eos::SoaveRedlichKwong => self.a::<eos::SoaveRedlichKwong>(t),
            Eos::PengRobinson => self.a::<eos::PengRobinson>(t),
            Eos::PatelTejaValderrama => self.a::<eos::PatelTejaValderrama>(t),
        }
    }

    /// The molecular volume parameter
    fn b_eos(&self, eos: Eos) -> f64 {
        match eos {
            Eos::IdealGas => self.b::<eos::IdealGas>(),
            Eos::VanDerWaals => self.b::<eos::VanDerWaals>(),
            Eos::RedlichKwong => self.b::<eos::RedlichKwong>(),
            Eos::SoaveRedlichKwong => self.b::<eos::SoaveRedlichKwong>(),
            Eos::PengRobinson => self.b::<eos::PengRobinson>(),
            Eos::PatelTejaValderrama => self.b::<eos::PatelTejaValderrama>(),
        }
    }

    /// Compute the pressure of the gas for the molar volume and temperature.
    fn pressure_eos(&self, eos: Eos, vm: f64, t: f64) -> f64 {
        match eos {
            Eos::IdealGas => self.pressure::<eos::IdealGas>(vm, t),
            Eos::VanDerWaals => self.pressure::<eos::VanDerWaals>(vm, t),
            Eos::RedlichKwong => self.pressure::<eos::RedlichKwong>(vm, t),
            Eos::SoaveRedlichKwong => self.pressure::<eos::SoaveRedlichKwong>(vm, t),
            Eos::PengRobinson => self.pressure::<eos::PengRobinson>(vm, t),
            Eos::PatelTejaValderrama => self.pressure::<eos::PatelTejaValderrama>(vm, t),
        }
    }

    /// Compute the compression factor Z such as Z = PV/RT
    ///
    /// Effectively resolves the cubic equation of state as a function of `p` and `t`.
    ///
    /// # Arguments
    ///  * `p` - The pressure of the gas, in Pa
    ///  * `t` - The temperature of the gas, in K
    ///
    /// # Panics
    /// This function will panic of no positive real root can be found, which is generally
    /// an indication that the parameters have physical non-sense.
    fn z_eos(&self, eos: Eos, p: f64, t: f64) -> f64 {
        match eos {
            Eos::IdealGas => self.z::<eos::IdealGas>(p, t),
            Eos::VanDerWaals => self.z::<eos::VanDerWaals>(p, t),
            Eos::RedlichKwong => self.z::<eos::RedlichKwong>(p, t),
            Eos::SoaveRedlichKwong => self.z::<eos::SoaveRedlichKwong>(p, t),
            Eos::PengRobinson => self.z::<eos::PengRobinson>(p, t),
            Eos::PatelTejaValderrama => self.z::<eos::PatelTejaValderrama>(p, t),
        }
    }
}

/// An helper trait to compute extensive state for equation of state known at runtime.
pub trait ExtensiveStateEos: StateEos {
    /// Compute the amount of mols for given pressure, volume and temperature.
    ///
    /// # Panics
    /// This function can panic if the parameters have physical non-sense
    fn mols_eos(&self, eos: Eos, p: f64, v: f64, t: f64) -> f64 {
        let z = self.z_eos(eos, p, t);
        p * v / z / R / t
    }

    /// Compute the volume of the gas for the pressure, mols and temperature.
    ///
    /// # Panics
    /// This function can panic if the parameters have physical non-sense
    fn volume_eos(&self, eos: Eos, p: f64, n: f64, t: f64) -> f64 {
        let z = self.z_eos(eos, p, t);
        n * z * R * t / p
    }
}

impl State for Molecule {
    fn a<E: EquationOfState>(&self, t: f64) -> f64 {
        E::a(self.pc, self.tc, self.zc(), self.w, t)
    }

    fn b<E: EquationOfState>(&self) -> f64 {
        E::b(self.pc, self.tc, self.zc())
    }

    fn c<E: EquationOfState>(&self) -> f64 {
        E::c(self.pc, self.tc, self.zc())
    }

    fn molar_mass(&self) -> f64 {
        self.m
    }
}

impl ExtensiveState for Molecule {}
impl StateEos for Molecule {}
impl ExtensiveStateEos for Molecule {}

impl State for Mixture {
    fn a<E: EquationOfState>(&self, t: f64) -> f64 {
        let mut res = 0f64;
        for (fi, pi) in self.comps.iter() {
            let ai = E::a(pi.pc, pi.tc, pi.zc(), pi.w, t);
            for (fj, pj) in self.comps.iter() {
                let aj = E::a(pj.pc, pj.tc, pj.zc(), pj.w, t);
                res += fi * fj * (ai * aj).sqrt();
            }
        }
        res
    }

    fn b<E: EquationOfState>(&self) -> f64 {
        self.comps
            .iter()
            .fold(0.0, |s, (f, p)| s + f * E::b(p.pc, p.tc, p.zc()))
    }

    fn c<E: EquationOfState>(&self) -> f64 {
        self.comps
            .iter()
            .fold(0.0, |s, (f, p)| s + f * E::c(p.pc, p.tc, p.zc()))
    }

    fn molar_mass(&self) -> f64 {
        self.comps
            .iter()
            .fold(0.0, |s, (f, p)| s + f * p.m)
    }
}

impl ExtensiveState for Mixture {}
impl StateEos for Mixture {}
impl ExtensiveStateEos for Mixture {}

impl State for Gas {
    fn a<E: EquationOfState>(&self, t: f64) -> f64 {
        match self {
            Gas::Molecule(props) => props.a::<E>(t),
            Gas::Mixture(mix) => mix.a::<E>(t),
        }
    }

    fn b<E: EquationOfState>(&self) -> f64 {
        match self {
            Gas::Molecule(props) => props.b::<E>(),
            Gas::Mixture(mix) => mix.b::<E>(),
        }
    }

    fn c<E: EquationOfState>(&self) -> f64 {
        match self {
            Gas::Molecule(props) => props.c::<E>(),
            Gas::Mixture(mix) => mix.c::<E>(),
        }
    }

    fn molar_mass(&self) -> f64 {
        match self {
            Gas::Molecule(props) => props.molar_mass(),
            Gas::Mixture(mix) => mix.molar_mass(),
        }
    }
}

impl ExtensiveState for Gas {}
impl StateEos for Gas {}
impl ExtensiveStateEos for Gas {}

#[cfg(test)]
mod tests {
    use super::{State};
    use crate::{eos, molecules};
    use float_eq::assert_float_eq;

    #[test]
    fn h2_mobility() {
        // H2 in mobility storage is reputed at 39.75 kg/m3
        let h2 = molecules::H2;
        let h2_storage_mass = 39.75; // kg/m3
        type E = eos::PengRobinson;

        // H2 mobility storage conditions (70 MPa, 15°C)
        let p = 70.0 * 1e6 + 101325.0;
        let t = 20.0 + 273.15;
        let mass = h2.specific_mass::<E>(p, t);
        assert_float_eq!(mass, h2_storage_mass, r2nd <= 0.05);

        // H2 mobility storage fueling conditions (87.5 MPa, 80°C)
        let p = 87.5 * 1e6 + 101325.0;
        let t = 85.0 + 273.15;
        let mass = h2.specific_mass::<E>(p, t);
        assert_float_eq!(mass, h2_storage_mass, r2nd <= 0.07);
    }
}
