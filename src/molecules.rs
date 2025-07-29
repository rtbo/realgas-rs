/// Physical constants of gas molecules
/// source: http://www.kaylaiacovino.com/Petrology_Tools/Critical_Constants_and_Acentric_Factors.htm
use crate::Molecule;

pub fn lookup_molecule<S>(name: S) -> Option<Molecule> 
where S: AsRef<str>
{
    match name.as_ref() {
        "Ar" => Some(AR),
        "Br2" => Some(BR2),
        "Cl2" => Some(CL2),
        "F2" => Some(F2),
        "He" => Some(HE),
        "H2" => Some(H2),
        "I2" => Some(I2),
        "Kr" => Some(KR),
        "Ne" => Some(NE),
        "N2" => Some(N2),
        "O2" => Some(O2),
        "Xe" => Some(XE),
        "C2H2" => Some(C2H2),
        "C6H6" => Some(C6H6),
        "C4H10" => Some(C4H10),
        "C4H8" => Some(C4H8),
        "C6H12" => Some(C6H12),
        "C3H6" => Some(C3H6),
        "C2H6" => Some(C2H6),
        "C2H4" => Some(C2H4),
        "NH3" => Some(NH3),
        "CO2" => Some(CO2),
        "CO" => Some(CO),
        "NO" => Some(NO),
        "SO2" => Some(SO2),
        "SO3" => Some(SO3),
        "H2O" => Some(H2O),
        "CH3COOH" => Some(CH3COOH),
        "C3H6O" => Some(C3H6O),
        "C2H5OH" => Some(C2H5OH),
        "CH3OH" => Some(CH3OH),
        "CH3CL" => Some(CH3CL),
        _ => None,
    }
}

/// Argon
pub const AR: Molecule = Molecule {
    tc: 150.8,
    pc: 48.7 * 1e5,
    vc: 74.9 * 1e-6,
    w: 0.001,
    m: 0.039948,
};

/// Bromine
pub const BR2: Molecule = Molecule {
    tc: 588.0,
    pc: 103.4 * 1e5,
    vc: 127.2 * 1e-6,
    w: 0.108,
    m: 0.159808,
};

/// Chlore
pub const CL2: Molecule = Molecule {
    tc: 416.9,
    pc: 79.8 * 1e5,
    vc: 123.8 * 1e-6,
    w: 0.09,
    m: 0.070906,
};

/// Fluor
pub const F2: Molecule = Molecule {
    tc: 144.3,
    pc: 52.2 * 1e5,
    vc: 66.3 * 1e-6,
    w: 0.054,
    m: 0.0379968,
};

/// Helium
pub const HE: Molecule = Molecule {
    tc: 5.19,
    pc: 2.27 * 1e5,
    vc: 57.4 * 1e-6,
    w: -0.365,
    m: 0.004002602,
};

/// Hydrogen
pub const H2: Molecule = Molecule {
    tc: 33.0,
    pc: 12.9 * 1e5,
    vc: 64.3 * 1e-6,
    w: -0.216,
    m: 0.00201588,
};

/// Iode
pub const I2: Molecule = Molecule {
    tc: 819.0,
    pc: 116.5 * 1e5,
    vc: 155.0 * 1e-6,
    w: 0.229,
    m: 0.25380894,
};

/// Krypton
pub const KR: Molecule = Molecule {
    tc: 209.4,
    pc: 55.0 * 1e5,
    vc: 91.2 * 1e-6,
    w: 0.005,
    m: 0.083798,
};

/// Neon
pub const NE: Molecule = Molecule {
    tc: 44.4,
    pc: 27.6 * 1e5,
    vc: 41.6 * 1e-6,
    w: -0.029,
    m: 0.0201797,
};

/// Nitrogen
pub const N2: Molecule = Molecule {
    tc: 126.2,
    pc: 33.9 * 1e5,
    vc: 89.8 * 1e-6,
    w: 0.039,
    m: 0.0280134,
};

/// Oxygen
pub const O2: Molecule = Molecule {
    tc: 154.6,
    pc: 50.4 * 1e5,
    vc: 73.4 * 1e-6,
    w: 0.025,
    m: 0.0319988,
};

/// Xenon
pub const XE: Molecule = Molecule {
    tc: 289.7,
    pc: 58.4 * 1e5,
    vc: 66.3 * 1e-6,
    w: 0.008,
    m: 0.131293,
};

/// Acetylene
pub const C2H2: Molecule = Molecule {
    tc: 308.3,
    pc: 61.4 * 1e5,
    vc: 112.7 * 1e-6,
    w: 0.19,
    m: 0.0260373,
};

/// Benzene
pub const C6H6: Molecule = Molecule {
    tc: 562.1,
    pc: 48.9 * 1e5,
    vc: 259.0 * 1e-6,
    w: 0.212,
    m: 0.0781118,
};

/// Butane
pub const C4H10: Molecule = Molecule {
    tc: 425.2,
    pc: 38.0 * 1e5,
    vc: 255.0 * 1e-6,
    w: 0.199,
    m: 0.0581222,
};

/// Cyclobutane
pub const C4H8: Molecule = Molecule {
    tc: 460.0,
    pc: 49.9 * 1e5,
    vc: 210.0 * 1e-6,
    w: 0.181,
    m: 0.0561063,
};

/// Cyclohexane
pub const C6H12: Molecule = Molecule {
    tc: 553.8,
    pc: 40.7 * 1e5,
    vc: 308. * 1e-6,
    w: 0.212,
    m: 0.0841595,
};

/// Cyclopropane
pub const C3H6: Molecule = Molecule {
    tc: 397.8,
    pc: 54.9 * 1e5,
    vc: 163.0 * 1e-6,
    w: 0.130,
    m: 0.0420797,
};

/// Ethane
pub const C2H6: Molecule = Molecule {
    tc: 305.4,
    pc: 48.8 * 1e5,
    vc: 148.3 * 1e-6,
    w: 0.099,
    m: 0.030069,
};

/// Ethylene
pub const C2H4: Molecule = Molecule {
    tc: 282.4,
    pc: 50.4 * 1e5,
    vc: 130.4 * 1e-6,
    w: 0.089,
    m: 0.0280532,
};

/// Ammonia
pub const NH3: Molecule = Molecule {
    tc: 405.5,
    pc: 113.5 * 1e5,
    vc: 72.5 * 1e-6,
    w: 0.250,
    m: 0.01703052,
};

/// Carbon dioxide
pub const CO2: Molecule = Molecule {
    tc: 304.1,
    pc: 73.8 * 1e5,
    vc: 93.9 * 1e-6,
    w: 0.239,
    m: 0.0440095,
};

/// Carbon monoxide
pub const CO: Molecule = Molecule {
    tc: 132.9,
    pc: 35.0 * 1e5,
    vc: 93.2 * 1e-6,
    w: 0.066,
    m: 0.0280101,
};

/// Nitric oxide
pub const NO: Molecule = Molecule {
    tc: 180.0,
    pc: 64.8 * 1e5,
    vc: 57.7 * 1e-6,
    w: 0.588,
    m: 0.0300061,
};

/// Sulfur dioxide
pub const SO2: Molecule = Molecule {
    tc: 430.8,
    pc: 78.8 * 1e5,
    vc: 122.2 * 1e-6,
    w: 0.256,
    m: 0.064066,
};

/// Sulfur trioxide
pub const SO3: Molecule = Molecule {
    tc: 491.0,
    pc: 82.1 * 1e5,
    vc: 127.3 * 1e-6,
    w: 0.481,
    m: 0.080066,
};

/// Water
pub const H2O: Molecule = Molecule {
    tc: 647.3,
    pc: 221.2 * 1e5,
    vc: 57.1 * 1e-6,
    w: 0.344,
    m: 0.01801528,
};

/// Acetic acid
pub const CH3COOH: Molecule = Molecule {
    tc: 592.7,
    pc: 57.9 * 1e5,
    vc: 66.3 * 1e-6,
    w: 0.09,
    m: 0.060052,
};

/// Acetone
pub const C3H6O: Molecule = Molecule {
    tc: 508.1,
    pc: 47.0 * 1e5,
    vc: 209.0 * 1e-6,
    w: 0.304,
    m: 0.0580791,
};

/// Ethanol
pub const C2H5OH: Molecule = Molecule {
    tc: 513.9,
    pc: 61.4 * 1e5,
    vc: 167.1 * 1e-6,
    w: 0.644,
    m: 0.04606844,
};

/// Methanol
pub const CH3OH: Molecule = Molecule {
    tc: 512.6,
    pc: 80.9 * 1e5,
    vc: 118.0 * 1e-6,
    w: 0.556,
    m: 0.03204294,
};

/// Methyl Chloride
pub const CH3CL: Molecule = Molecule {
    tc: 416.3,
    pc: 67.0 * 1e5,
    vc: 138.9 * 1e-6,
    w: 0.153,
    m: 0.0504905,
};
