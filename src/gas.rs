use std::str::FromStr;

/// A gas molecule, represented by its physical properties.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Molecule {
    /// The critical pressure in Pa
    pub pc: f64,
    /// The critical temperature in K
    pub tc: f64,
    /// The acentric factor
    pub w: f64,
}

/// A mixture of several gases
#[derive(Debug, Clone)]
pub struct Mixture {
    pub(crate) comps: Vec<(f64, Molecule)>,
}

/// A mixture error
pub enum MixtureError {
    MixtureNotWhole,
    InvalidFactor,
}

/// A component to build a mixture
#[derive(Debug, Clone)]
pub enum Comp {
    Factor(f64, Gas),
    Remainder(Gas),
}

impl Mixture {
    pub fn new<I>(comps: I) -> Result<Mixture, MixtureError>
    where
        I: IntoIterator<Item = Comp>,
    {
        let mut tmp: Vec<(bool, f64, Molecule)> = Vec::new(); // first tuple field means remainder
        let mut fill = 0f64;
        let mut num_voids = 0;

        for c in comps {
            let (f, g) = match c {
                Comp::Factor(f, g) => (f, g),
                Comp::Remainder(g) => (f64::NAN, g),
            };
            if f.is_nan() {
                num_voids += 1;
            } else {
                if f <= 0f64 || f >= 1f64 {
                    return Err(MixtureError::InvalidFactor);
                }
                fill += f;
            }
            match g {
                Gas::Molecule(m) => {
                    tmp.push((f.is_nan(), f, m));
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

        let comps: Vec<(f64, Molecule)> = tmp.into_iter().map(|(_, f, m)| (f, m)).collect();
        
        debug_assert!(comps.iter().map(|(f, _)| *f).sum::<f64>() > 0.9999999);
        debug_assert!(comps.iter().map(|(f, _)| *f).sum::<f64>() < 1.0000001);
        
        Ok(Mixture { comps })
    }
}

/// A generic gas, that can be either a molecule or a mixture.
#[derive(Debug, Clone)]
pub enum Gas {
    Molecule(Molecule),
    Mixture(Mixture),
}
