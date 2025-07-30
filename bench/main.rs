use realgas::compounds;

mod z;

const EXP_Z_H2_CSV: &str = include_str!("exp/z_h2.csv");
const EXP_Z_N2_CSV: &str = include_str!("exp/z_n2.csv");
const EXP_Z_WATER_CSV: &str = include_str!("exp/z_water.csv");
const EXP_Z_AIR_CSV: &str = include_str!("exp/z_air.csv");

fn main() {
    bench_z();
}

fn bench_z() {
    let h2 = compounds::H2.into();
    z::do_gas(EXP_Z_H2_CSV, "H2", &h2, &[40.0, 300.0, 2000.0]);

    let n2 = compounds::N2.into();
    z::do_gas(EXP_Z_N2_CSV, "N2", &n2, &[80.0, 300.0, 1000.0]);

    let water = compounds::H2O.into();
    z::do_gas(EXP_Z_WATER_CSV, "water", &water, &[400.0, 800.0, 2000.0]);

    let air = compounds::dry_air().into();
    z::do_gas(EXP_Z_AIR_CSV, "air", &air, &[100.0, 300.0, 1000.0]);
}
