use realgas::compounds;


mod zbench;

const EXP_Z_AIR_CSV: &str = include_str!("exp/z_air.csv");

fn main() {
    let air = compounds::dry_air().into();

    zbench::do_gas(EXP_Z_AIR_CSV, "air", &air, &[100.0, 300.0, 1000.0]);
}
