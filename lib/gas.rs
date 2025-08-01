use crate::{Pvt, compounds};
use std::{borrow::Borrow, cmp::Reverse, fmt, num::ParseFloatError, str::FromStr};

/// A gas molecule, represented by its physical properties.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Molecule {
    /// The molar mass in kg/mol
    pub m: f64,
    /// The critical state of this molecule
    pub critical_state: Pvt,
    /// The acentric factor
    pub w: f64,
}

impl PartialOrd for Molecule {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.m
            .partial_cmp(&other.m)
            .or_else(|| self.critical_state.partial_cmp(&other.critical_state))
            .or_else(|| self.w.partial_cmp(&other.w))
    }
}

/// A mixture of several gases
#[derive(Debug, Clone, PartialEq)]
pub struct Mixture {
    pub(crate) comps: Vec<(f64, Molecule)>,
}

/// A mixture error
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MixtureError {
    MixtureNotWhole,
    InvalidFraction(f64),
}

impl fmt::Display for MixtureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MixtureError::MixtureNotWhole => write!(f, "The sum of fractions does not equal to 100%"),
            MixtureError::InvalidFraction(fraction) => write!(f, "{:.1}% isn't a valid molar fraction", fraction),
        }
        
    }
}

impl std::error::Error for MixtureError {}

/// A component to build a mixture
#[derive(Debug, Clone)]
pub enum Comp {
    Factor(f64, Gas),
    Remainder(Gas),
}

impl Mixture {
    pub fn new<I>(comps: I) -> Result<Mixture, MixtureError>
    where
        I: IntoIterator,
        I::Item: Borrow<Comp>,
    {
        let mut tmp: Vec<(bool, f64, Molecule)> = Vec::new(); // first tuple field means remainder
        let mut fill = 0f64;
        let mut num_voids = 0;

        for c in comps {
            let c = c.borrow();

            let (f, g) = match c {
                Comp::Factor(f, g) => (*f, g),
                Comp::Remainder(g) => (f64::NAN, g),
            };
            if f.is_nan() {
                num_voids += 1;
            } else {
                if f <= 0f64 || f >= 1f64 {
                    return Err(MixtureError::InvalidFraction(f));
                }
                fill += f;
            }
            match g {
                Gas::Molecule(m) => {
                    tmp.push((f.is_nan(), f, *m));
                }
                Gas::Mixture(Mixture { comps }) => {
                    for c in comps {
                        if f.is_nan() {
                            tmp.push((true, c.0, c.1));
                        } else {
                            tmp.push((false, f * c.0, c.1));
                        }
                    }
                }
            }
        }

        if fill > 1.0 {
            return Err(MixtureError::MixtureNotWhole);
        }
        if fill != 1.0 && num_voids == 0 {
            return Err(MixtureError::MixtureNotWhole);
        }

        if num_voids > 0 {
            let void_attrib = (1.0 - fill) / num_voids as f64;
            for c in &mut tmp {
                if c.0 {
                    if c.1.is_nan() {
                        c.1 = void_attrib;
                    } else {
                        c.1 *= void_attrib;
                    }
                }
            }
        }

        let mut comps: Vec<(f64, Molecule)> = tmp.into_iter().map(|(_, f, m)| (f, m)).collect();

        // Following sort and merge make the components always the same for a given mixture.
        // e.g. mixing air with O2 will result with a single O2 component instead of 2,
        // and components will always be in the same order.
        // This makes mixtures trivially comparable

        // sort with decreasing order of ratio, followed by decreasing order of molar mass
        // followed by decreasing order of critical parameters
        comps.sort_by(|(fa, ma), (fb, mb)| {
            Reverse((*fa, ma))
                .partial_cmp(&Reverse((*fb, mb)))
                .unwrap()
        });

        // merge gases that have identical properties
        let mut i1 = 0;
        let mut i2 = 1;
        while i2 < comps.len() {
            if comps[i1].1 == comps[i2].1 {
                comps[i1].0 += comps[i2].0;
                comps.remove(i2);
            } else {
                i1 += 1;
                i2 += 1;
            }
        }

        debug_assert!(comps.iter().map(|(f, _)| *f).sum::<f64>() > 0.9999999);
        debug_assert!(comps.iter().map(|(f, _)| *f).sum::<f64>() < 1.0000001);

        Ok(Mixture { comps })
    }
}

/// A generic gas, that can be either a molecule or a mixture.
#[derive(Debug, Clone, PartialEq)]
pub enum Gas {
    Molecule(Molecule),
    Mixture(Mixture),
}

impl From<Molecule> for Gas {
    fn from(value: Molecule) -> Self {
        Gas::Molecule(value)
    }
}

impl From<Mixture> for Gas {
    fn from(value: Mixture) -> Self {
        Gas::Mixture(value)
    }
}

#[cfg(test)]
impl From<&Gas> for Gas {
    fn from(value: &Gas) -> Self {
        value.clone()
    }
}

#[cfg(test)]
impl From<&Molecule> for Gas {
    fn from(value: &Molecule) -> Self {
        Gas::Molecule(*value)
    }
}

#[cfg(test)]
impl From<&Mixture> for Gas {
    fn from(value: &Mixture) -> Self {
        Gas::Mixture(value.clone())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum GasParseError {
    UnknownMolecule(String),
    Mixture(MixtureError),
    Float(ParseFloatError),
    Other(String),
}

impl From<MixtureError> for GasParseError {
    fn from(value: MixtureError) -> Self {
        GasParseError::Mixture(value)
    }
}
impl From<ParseFloatError> for GasParseError {
    fn from(value: ParseFloatError) -> Self {
        GasParseError::Float(value)
    }
}

impl fmt::Display for GasParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GasParseError::UnknownMolecule(m) => write!(f, "Can't lookup {m} as a known molecule"),
            GasParseError::Mixture(m) => m.fmt(f),
            GasParseError::Float(err) => err.fmt(f),
            GasParseError::Other(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for GasParseError {}

impl FromStr for Gas {
    type Err = GasParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let scomps: Vec<&str> = s.split("+").collect();

        if scomps.is_empty() {
            Err(GasParseError::Mixture(MixtureError::MixtureNotWhole))
        } else if scomps.len() == 1 {
            compounds::lookup(&scomps[0])
                .ok_or_else(|| GasParseError::UnknownMolecule(scomps[0].to_string()))
        } else {
            let mut mcomps = Vec::<Comp>::new();
            for sc in scomps {
                let sfrac: Vec<&str> = sc.split("%").collect();
                if sfrac.len() > 2 {
                    return Err(GasParseError::Other(format!("Can't parse {sc} as a compound fraction")));
                }
                let symbol = *sfrac.iter().last().unwrap();
                let g = compounds::lookup(symbol)
                    .ok_or_else(|| GasParseError::UnknownMolecule(symbol.to_string()))?;
                if sfrac.len() == 1 {
                    mcomps.push(Comp::Remainder(g));
                } else {
                    let frac = sfrac[0]
                        .parse::<f64>()?;
                    mcomps.push(Comp::Factor(frac / 100.0, g));
                }
            }

            Ok(Gas::Mixture(Mixture::new(mcomps)?))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Comp, Gas, Mixture};
    use crate::{Molecule, State, compounds, eos::PengRobinson, gas::MixtureError};
    use float_eq::assert_float_eq;

    fn assert_molecule_eq(lhs: &Molecule, rhs: &Molecule, rtol: f64) {
        assert_float_eq!(lhs.m, rhs.m, r1st <= rtol);
        assert_float_eq!(lhs.critical_state.p, rhs.critical_state.p, r1st <= rtol);
        assert_float_eq!(lhs.critical_state.v, rhs.critical_state.v, r1st <= rtol);
        assert_float_eq!(lhs.critical_state.t, rhs.critical_state.t, r1st <= rtol);
        assert_float_eq!(lhs.w, rhs.w, r1st <= rtol);
    }

    fn assert_mixture_eq(lhs: &Mixture, rhs: &Mixture, rtol: f64) {
        if lhs.comps.len() != rhs.comps.len() {
            panic!("assertion failed: lhs and rhs are mixtures with different components count");
        }
        for idx in 0..lhs.comps.len() {
            let cl = &lhs.comps[idx];
            let cr = &rhs.comps[idx];
            assert_float_eq!(cl.0, cr.0, r1st <= rtol);
            assert_molecule_eq(&cl.1, &cr.1, rtol);
        }
    }

    fn assert_gas_eq(lhs: &Gas, rhs: &Gas, rtol: f64) {
        match (lhs, rhs) {
            (Gas::Molecule(_), Gas::Mixture(_)) => {
                panic!("assertion failed, lhs is a molecule, rhs is a mixture")
            }
            (Gas::Mixture(_), Gas::Molecule(_)) => {
                panic!("assertion failed, lhs is a mixture, rhs is a molecule")
            }
            (Gas::Molecule(lhs), Gas::Molecule(rhs)) => {
                assert_molecule_eq(lhs, rhs, rtol);
            }
            (Gas::Mixture(lhs), Gas::Mixture(rhs)) => {
                assert_mixture_eq(lhs, rhs, rtol);
            }
        }
    }

    #[test]
    fn parse_molecule_works() {
        let gas: Gas = "N2".parse().expect("should parse N2");
        assert_eq!(gas, Gas::from(compounds::N2));
    }

    #[test]
    fn parse_dry_air_works() {
        let parsed_air: Gas = "78.08%N2+20.95%O2+0.93%Ar+CO2"
            .parse()
            .expect("should parse dry air composition");
        let built_air = Gas::Mixture(
            Mixture::new(vec![
                Comp::Factor(0.7808, compounds::N2.into()),
                Comp::Factor(0.2095, compounds::O2.into()),
                Comp::Factor(0.0093, compounds::AR.into()),
                Comp::Remainder(compounds::CO2.into()),
            ])
            .unwrap(),
        );

        let p = 200.0 * 1e5;
        let t = 273.15 - 80.0;
        let z1 = parsed_air.z::<PengRobinson>(p, t);
        let z2 = built_air.z::<PengRobinson>(p, t);
        assert_float_eq!(z1, z2, ulps <= 4);
        assert_gas_eq(&parsed_air, &built_air, 0.00001);
    }

    #[test]
    fn mixture_new_reports_mixture_not_whole() {
        fn assert(res: Result<Mixture, MixtureError>) {
            assert!(res.is_err());
            assert_eq!(res.unwrap_err(), MixtureError::MixtureNotWhole);
        }

        assert(Mixture::new(&[]));

        assert(Mixture::new(&[
            Comp::Factor(0.5, compounds::N2.into()),
            Comp::Factor(0.3, compounds::O2.into()),
            Comp::Factor(0.1, compounds::AR.into()),
        ]));

        assert(Mixture::new(&[
            Comp::Factor(0.5, compounds::N2.into()),
            Comp::Factor(0.5, compounds::O2.into()),
            Comp::Factor(0.1, compounds::AR.into()),
        ]));
    }

    #[test]
    fn can_compare_identical_mixtures_built_in_any_order() {
        let air_n2 = 0.7808;
        let air_o2 = 0.2095;
        let air_ar = 0.0093;
        let air_co2 = 0.0004;
        let air = Mixture::new(&[
            Comp::Factor(air_o2, compounds::O2.into()),
            Comp::Factor(air_n2, compounds::N2.into()),
            Comp::Factor(air_ar, compounds::AR.into()),
            Comp::Factor(air_co2, compounds::CO2.into()),
        ])
        .unwrap();
        let mix1 = Mixture::new(&[
            Comp::Factor(0.1, compounds::O2.into()),
            Comp::Remainder(air.clone().into()),
        ])
        .unwrap();
        let mix2 = Mixture::new(&[
            Comp::Factor(0.9, air.into()),
            Comp::Remainder(compounds::O2.into()),
        ])
        .unwrap();
        let mix3 = Mixture::new(&[
            Comp::Factor(air_n2 * 0.9, compounds::N2.into()),
            Comp::Factor(air_o2 * 0.9 + 0.1, compounds::O2.into()),
            Comp::Factor(air_ar * 0.9, compounds::AR.into()),
            Comp::Factor(air_co2 * 0.9, compounds::CO2.into()),
        ])
        .unwrap();
        let mix4 = Mixture::new(&[
            Comp::Factor(air_n2 * 0.9, compounds::N2.into()),
            Comp::Factor(air_o2 * 0.9, compounds::O2.into()),
            Comp::Factor(air_ar * 0.9, compounds::AR.into()),
            Comp::Factor(air_co2 * 0.9, compounds::CO2.into()),
            Comp::Factor(0.1, compounds::O2.into()),
        ])
        .unwrap();

        assert_eq!(mix1.comps.len(), 4);
        assert_eq!(mix1.comps[0].1, compounds::N2);
        assert_eq!(mix1.comps[1].1, compounds::O2);
        assert_eq!(mix1.comps[2].1, compounds::AR);
        assert_eq!(mix1.comps[3].1, compounds::CO2);
        assert_mixture_eq(&mix1, &mix2, 0.00001);
        assert_mixture_eq(&mix2, &mix3, 0.00001);
        assert_mixture_eq(&mix3, &mix4, 0.00001);
    }
}
