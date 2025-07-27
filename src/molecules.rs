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
        _ => None,
    }
}

/// Argon
pub const AR: Molecule = Molecule {
    tc: 150.8f64,
    pc: 4_870_000f64,
    w: 0.001f64,
};

/// Brome
pub const BR2: Molecule = Molecule {
    tc: 588f64,
    pc: 10_340_000f64,
    w: 0.108f64,
};

/// Chlore
pub const CL2: Molecule = Molecule {
    tc: 416.9f64,
    pc: 7_980_000f64,
    w: 0.09f64,
};

/// Fluor
pub const F2: Molecule = Molecule {
    tc: 144.3f64,
    pc: 5_220_000f64,
    w: 0.054f64,
};

/// Helium
pub const HE: Molecule = Molecule {
    tc: 5.19f64,
    pc: 227_000f64,
    w: -0.365f64,
};

/// Hydrogen
pub const H2: Molecule = Molecule {
    tc: 33f64,
    pc: 1_290_000f64,
    w: -0.216f64,
};

/// Iode
pub const I2: Molecule = Molecule {
    tc: 819f64,
    pc: 11_650_000f64,
    w: 0.229f64,
};

/// Krypton
pub const KR: Molecule = Molecule {
    tc: 209.4f64,
    pc: 5_500_000f64,
    w: 0.005f64,
};

/// Neon
pub const NE: Molecule = Molecule {
    tc: 44.4f64,
    pc: 2_760_000f64,
    w: -0.029f64,
};

/// Nitrogen
pub const N2: Molecule = Molecule {
    tc: 126.2f64,
    pc: 3_390_000f64,
    w: 0.039f64,
};

/// Oxygen
pub const O2: Molecule = Molecule {
    tc: 154.6f64,
    pc: 5_040_000f64,
    w: 0.025f64,
};

/// Xenon
pub const XE: Molecule = Molecule {
    tc: 289.7f64,
    pc: 5_840_000f64,
    w: 0.008f64,
};

/// Acetylene
pub const C2H2: Molecule = Molecule {
    tc: 308.3f64,
    pc: 6_140_000f64,
    w: 0.19f64,
};

/// Benzene
pub const C6H6: Molecule = Molecule {
    tc: 562.1f64,
    pc: 4_890_000f64,
    w: 0.212f64,
};

/// Butane
pub const C4H10: Molecule = Molecule {
    tc: 425.2f64,
    pc: 3_800_000f64,
    w: 0.199f64,
};

/// Cyclobutane
pub const C4H8: Molecule = Molecule {
    tc: 460f64,
    pc: 4_990_000f64,
    w: 0.181f64,
};

/// Cyclohexane
pub const C6H12: Molecule = Molecule {
    tc: 553.8f64,
    pc: 4_070_000f64,
    w: 0.212f64,
};

/// Cyclopropane
pub const C3H6: Molecule = Molecule {
    tc: 397.8f64,
    pc: 5_490_000f64,
    w: 0.130f64,
};

/// Ethane
pub const C2H6: Molecule = Molecule {
    tc: 305.4f64,
    pc: 4_880_000f64,
    w: 0.099f64,
};

/// Ethylene
pub const C2H4: Molecule = Molecule {
    tc: 282.4f64,
    pc: 5_040_000f64,
    w: 0.089f64,
};

/// Ammonia
pub const NH3: Molecule = Molecule {
    tc: 405.5f64,
    pc: 11_350_000f64,
    w: 0.250f64,
};

/// Carbon dioxide
pub const CO2: Molecule = Molecule {
    tc: 304.1f64,
    pc: 7_380_000f64,
    w: 0.239f64,
};

/// Carbon monoxide
pub const CO: Molecule = Molecule {
    tc: 132.9f64,
    pc: 3_500_000f64,
    w: 0.066f64,
};

/// Nitric oxide
pub const NO: Molecule = Molecule {
    tc: 180f64,
    pc: 6_480_000f64,
    w: 0.588f64,
};

/// Sulfur dioxide
pub const SO2: Molecule = Molecule {
    tc: 430.8f64,
    pc: 7_880_000f64,
    w: 0.256f64,
};

/// Sulfur trioxide
pub const SO3: Molecule = Molecule {
    tc: 491f64,
    pc: 8_210_000f64,
    w: 0.481f64,
};

/// Water
pub const H2O: Molecule = Molecule {
    tc: 647.3f64,
    pc: 22_120_000f64,
    w: 0.344f64,
};

/// Acetic acid
pub const CH3COOH: Molecule = Molecule {
    tc: 592.7f64,
    pc: 5_790_000f64,
    w: 0.09f64,
};

/// Acetone
pub const C3H6O: Molecule = Molecule {
    tc: 508.1f64,
    pc: 4_700_000f64,
    w: 0.304f64,
};

/// Ethanol
pub const C2H5OH: Molecule = Molecule {
    tc: 513.9f64,
    pc: 6_140_000f64,
    w: 0.644f64,
};

/// Methanol
pub const CH3OH: Molecule = Molecule {
    tc: 512.6f64,
    pc: 8_090_000f64,
    w: 0.556f64,
};
