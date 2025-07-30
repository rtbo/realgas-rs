/// Physical constants of gas molecules
/// source: http://www.kaylaiacovino.com/Petrology_Tools/Critical_Constants_and_Acentric_Factors.htm
use crate::{Gas, Mixture, Molecule, Pvt};

pub fn lookup<S>(name: S) -> Option<Gas>
where
    S: AsRef<str>,
{
    match name.as_ref() {
        "dry_air" => Some(dry_air().into()),
        "Ar" => Some(AR.into()),
        "Br2" => Some(BR2.into()),
        "Cl2" => Some(CL2.into()),
        "F2" => Some(F2.into()),
        "He" => Some(HE.into()),
        "H2" => Some(H2.into()),
        "I2" => Some(I2.into()),
        "Kr" => Some(KR.into()),
        "Ne" => Some(NE.into()),
        "N2" => Some(N2.into()),
        "O2" => Some(O2.into()),
        "Xe" => Some(XE.into()),
        "C2H2" => Some(C2H2.into()),
        "C6H6" => Some(C6H6.into()),
        "C4H10" => Some(C4H10.into()),
        "C4H8" => Some(C4H8.into()),
        "C6H12" => Some(C6H12.into()),
        "C3H6" => Some(C3H6.into()),
        "C2H6" => Some(C2H6.into()),
        "C2H4" => Some(C2H4.into()),
        "NH3" => Some(NH3.into()),
        "CO2" => Some(CO2.into()),
        "CO" => Some(CO.into()),
        "NO" => Some(NO.into()),
        "SO2" => Some(SO2.into()),
        "SO3" => Some(SO3.into()),
        "H2O" => Some(H2O.into()),
        "CH3COOH" => Some(CH3COOH.into()),
        "C3H6O" => Some(C3H6O.into()),
        "C2H5OH" => Some(C2H5OH.into()),
        "CH3OH" => Some(CH3OH.into()),
        "CH3CL" => Some(CH3CL.into()),
        _ => None,
    }
}

/// Air mixture
pub fn dry_air() -> Mixture {
    use crate::gas::Comp;
    Mixture::new(&[
        Comp::Factor(0.7808, N2.into()),
        Comp::Factor(0.2095, O2.into()),
        Comp::Factor(0.0093, AR.into()),
        Comp::Factor(0.0004, CO2.into()),
    ])
    .unwrap()
}

/// Argon
pub const AR: Molecule = Molecule {
    critical_state: Pvt {
        p: 48.7 * 1e5,
        v: 74.9 * 1e-6,
        t: 150.8,
    },
    w: 0.001,
    m: 0.039948,
};

/// Bromine
pub const BR2: Molecule = Molecule {
    critical_state: Pvt {
        p: 103.4 * 1e5,
        v: 127.2 * 1e-6,
        t: 588.0,
    },
    w: 0.108,
    m: 0.159808,
};

/// Chlore
pub const CL2: Molecule = Molecule {
    critical_state: Pvt {
        p: 79.8 * 1e5,
        v: 123.8 * 1e-6,
        t: 416.9,
    },
    w: 0.09,
    m: 0.070906,
};

/// Fluor
pub const F2: Molecule = Molecule {
    critical_state: Pvt {
        p: 52.2 * 1e5,
        v: 66.3 * 1e-6,
        t: 144.3,
    },
    w: 0.054,
    m: 0.0379968,
};

/// Helium
pub const HE: Molecule = Molecule {
    critical_state: Pvt {
        p: 2.27 * 1e5,
        v: 57.4 * 1e-6,
        t: 5.19,
    },
    w: -0.365,
    m: 0.004002602,
};

/// Hydrogen
pub const H2: Molecule = Molecule {
    critical_state: Pvt {
        p: 12.9 * 1e5,
        v: 64.3 * 1e-6,
        t: 33.0,
    },
    w: -0.216,
    m: 0.00201588,
};

/// Iode
pub const I2: Molecule = Molecule {
    critical_state: Pvt {
        p: 116.5 * 1e5,
        v: 155.0 * 1e-6,
        t: 819.0,
    },
    w: 0.229,
    m: 0.25380894,
};

/// Krypton
pub const KR: Molecule = Molecule {
    critical_state: Pvt {
        p: 55.0 * 1e5,
        v: 91.2 * 1e-6,
        t: 209.4,
    },
    w: 0.005,
    m: 0.083798,
};

/// Neon
pub const NE: Molecule = Molecule {
    critical_state: Pvt {
        p: 27.6 * 1e5,
        v: 41.6 * 1e-6,
        t: 44.4,
    },
    w: -0.029,
    m: 0.0201797,
};

/// Nitrogen
pub const N2: Molecule = Molecule {
    critical_state: Pvt {
        p: 33.9 * 1e5,
        v: 89.8 * 1e-6,
        t: 126.2,
    },
    w: 0.039,
    m: 0.0280134,
};

/// Oxygen
pub const O2: Molecule = Molecule {
    critical_state: Pvt {
        p: 50.4 * 1e5,
        v: 73.4 * 1e-6,
        t: 154.6,
    },
    w: 0.025,
    m: 0.0319988,
};

/// Xenon
pub const XE: Molecule = Molecule {
    critical_state: Pvt {
        p: 58.4 * 1e5,
        v: 66.3 * 1e-6,
        t: 289.7,
    },
    w: 0.008,
    m: 0.131293,
};

/// Acetylene
pub const C2H2: Molecule = Molecule {
    critical_state: Pvt {
        p: 61.4 * 1e5,
        v: 112.7 * 1e-6,
        t: 308.3,
    },
    w: 0.19,
    m: 0.0260373,
};

/// Benzene
pub const C6H6: Molecule = Molecule {
    critical_state: Pvt {
        p: 48.9 * 1e5,
        v: 259.0 * 1e-6,
        t: 562.1,
    },
    w: 0.212,
    m: 0.0781118,
};

/// Butane
pub const C4H10: Molecule = Molecule {
    critical_state: Pvt {
        p: 38.0 * 1e5,
        v: 255.0 * 1e-6,
        t: 425.2,
    },
    w: 0.199,
    m: 0.0581222,
};

/// Cyclobutane
pub const C4H8: Molecule = Molecule {
    critical_state: Pvt {
        p: 49.9 * 1e5,
        v: 210.0 * 1e-6,
        t: 460.0,
    },
    w: 0.181,
    m: 0.0561063,
};

/// Cyclohexane
pub const C6H12: Molecule = Molecule {
    critical_state: Pvt {
        p: 40.7 * 1e5,
        v: 308. * 1e-6,
        t: 553.8,
    },
    w: 0.212,
    m: 0.0841595,
};

/// Cyclopropane
pub const C3H6: Molecule = Molecule {
    critical_state: Pvt {
        p: 54.9 * 1e5,
        v: 163.0 * 1e-6,
        t: 397.8,
    },
    w: 0.130,
    m: 0.0420797,
};

/// Ethane
pub const C2H6: Molecule = Molecule {
    critical_state: Pvt {
        p: 48.8 * 1e5,
        v: 148.3 * 1e-6,
        t: 305.4,
    },
    w: 0.099,
    m: 0.030069,
};

/// Ethylene
pub const C2H4: Molecule = Molecule {
    critical_state: Pvt {
        p: 50.4 * 1e5,
        v: 130.4 * 1e-6,
        t: 282.4,
    },
    w: 0.089,
    m: 0.0280532,
};

/// Ammonia
pub const NH3: Molecule = Molecule {
    critical_state: Pvt {
        p: 113.5 * 1e5,
        v: 72.5 * 1e-6,
        t: 405.5,
    },
    w: 0.250,
    m: 0.01703052,
};

/// Carbon dioxide
pub const CO2: Molecule = Molecule {
    critical_state: Pvt {
        p: 73.8 * 1e5,
        v: 93.9 * 1e-6,
        t: 304.1,
    },
    w: 0.239,
    m: 0.0440095,
};

/// Carbon monoxide
pub const CO: Molecule = Molecule {
    critical_state: Pvt {
        p: 35.0 * 1e5,
        v: 93.2 * 1e-6,
        t: 132.9,
    },
    w: 0.066,
    m: 0.0280101,
};

/// Nitric oxide
pub const NO: Molecule = Molecule {
    critical_state: Pvt {
        p: 64.8 * 1e5,
        v: 57.7 * 1e-6,
        t: 180.0,
    },
    w: 0.588,
    m: 0.0300061,
};

/// Sulfur dioxide
pub const SO2: Molecule = Molecule {
    critical_state: Pvt {
        p: 78.8 * 1e5,
        v: 122.2 * 1e-6,
        t: 430.8,
    },
    w: 0.256,
    m: 0.064066,
};

/// Sulfur trioxide
pub const SO3: Molecule = Molecule {
    critical_state: Pvt {
        p: 82.1 * 1e5,
        v: 127.3 * 1e-6,
        t: 491.0,
    },
    w: 0.481,
    m: 0.080066,
};

/// Water
pub const H2O: Molecule = Molecule {
    critical_state: Pvt {
        p: 221.2 * 1e5,
        v: 57.1 * 1e-6,
        t: 647.3,
    },
    w: 0.344,
    m: 0.01801528,
};

/// Acetic acid
pub const CH3COOH: Molecule = Molecule {
    critical_state: Pvt {
        p: 57.9 * 1e5,
        v: 66.3 * 1e-6,
        t: 592.7,
    },
    w: 0.09,
    m: 0.060052,
};

/// Acetone
pub const C3H6O: Molecule = Molecule {
    critical_state: Pvt {
        p: 47.0 * 1e5,
        v: 209.0 * 1e-6,
        t: 508.1,
    },
    w: 0.304,
    m: 0.0580791,
};

/// Ethanol
pub const C2H5OH: Molecule = Molecule {
    critical_state: Pvt {
        p: 61.4 * 1e5,
        v: 167.1 * 1e-6,
        t: 513.9,
    },
    w: 0.644,
    m: 0.04606844,
};

/// Methanol
pub const CH3OH: Molecule = Molecule {
    critical_state: Pvt {
        p: 80.9 * 1e5,
        v: 118.0 * 1e-6,
        t: 512.6,
    },
    w: 0.556,
    m: 0.03204294,
};

/// Methyl Chloride
pub const CH3CL: Molecule = Molecule {
    critical_state: Pvt {
        p: 67.0 * 1e5,
        v: 138.9 * 1e-6,
        t: 416.3,
    },
    w: 0.153,
    m: 0.0504905,
};
