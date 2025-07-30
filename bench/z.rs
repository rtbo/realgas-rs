use plotters::{element::DashedPathElement, style::{Color, ShapeStyle, BLUE, CYAN, GREEN, MAGENTA, RED, YELLOW}};
use realgas::{eos::{self, EquationOfState}, Gas, State};

#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub t: f64,
    pub z: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Data {
    pub pcols: Vec<f64>,
    pub zrows: Vec<Row>,
}

impl Data {
    pub fn new() -> Self {
        Data {
            pcols: Vec::new(),
            zrows: Vec::new(),
        }
    }

    pub fn row(&self, t: f64) -> Option<&Row> {
        self.zrows.iter().find(|row| (row.t - t).abs() < f64::EPSILON)
    }

    pub fn pressures(&self) -> &[f64] {
        &self.pcols
    }

    pub fn temperatures(&self) -> Vec<f64> {
        self.zrows.iter().map(|row| row.t).collect()
    }

    pub fn gen_eos<E: EquationOfState>(gas: &Gas, pressures: &[f64], temperatures: &[f64]) -> Data {
        let mut data = Data {
            pcols: pressures.to_vec(),
            zrows: Vec::new(),
        };

        for &t in temperatures {
            let mut z_row = Row { t, z: Vec::new() };
            for &p in pressures {
                let z = gas.z::<E>(p, t);
                z_row.z.push(z);
            }
            data.zrows.push(z_row);
        }

        data
    }

    pub fn from_csv(csv_data: &str) -> Self {
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(csv_data.as_bytes());

        let mut data = Data::new();

        let head = rdr.headers().expect("Failed to read headers");
        for header in head.iter().skip(1) {
            let p = header.parse::<f64>().unwrap() * 1e5;
            data.pcols.push(p);
        }

        for result in rdr.records() {
            match result {
                Ok(record) => {
                    let t = record.get(0).unwrap().parse().unwrap();
                    let mut z = Vec::new();
                    for field in record.iter().skip(1) {
                        let value = if field.is_empty() {
                            f64::NAN
                        } else {
                            field.trim().parse().expect("Failed to parse field")
                        };
                        z.push(value);
                    }
                    data.zrows.push(Row { t, z });
                }
                Err(e) => eprintln!("Error reading record: {}", e),
            }
        }

        data
    }

    pub fn _to_csv(&self) -> String {
        let mut wtr = csv::Writer::from_writer(vec![]);
        wtr.write_record(std::iter::once("T".to_string()).chain(self.pcols.iter().map(|p| (p * 1e-5).to_string())))
            .expect("Failed to write header");

        for row in &self.zrows {
            let record: Vec<String> = std::iter::once(row.t.to_string())
                .chain(row.z.iter().map(|z| z.to_string()))
                .collect();
            wtr.write_record(record).expect("Failed to write record");
        }

        String::from_utf8(wtr.into_inner().expect("Failed to get inner writer")).expect("Failed to convert to string")
    }

    pub fn _write_csv(&self, filename: &str) {
        let csv_data = self._to_csv();
        std::fs::write(filename, csv_data)
            .expect("Failed to write CSV file");
    }
}

struct Series<'a> {
    name: &'a str,
    data: &'a Data,
    style: ShapeStyle,
    dashed: bool,
}

pub fn do_gas(exp_csv: &str, gas_name: &str, gas: &Gas, plot_temps: &[f64]) {
    
    let exp = Data::from_csv(exp_csv);
    let pressures = exp.pressures();
    let temperatures = exp.temperatures();

    let vdw = Data::gen_eos::<eos::VanDerWaals>(gas, pressures, &temperatures);
    let rk = Data::gen_eos::<eos::RedlichKwong>(gas, pressures, &temperatures);
    let srk = Data::gen_eos::<eos::SoaveRedlichKwong>(gas, pressures, &temperatures);
    let pr = Data::gen_eos::<eos::PengRobinson>(gas, pressures, &temperatures);
    let ptv = Data::gen_eos::<eos::PatelTejaValderrama>(gas, pressures, &temperatures);

    let series = &[
        Series {
            name: "Experimental",
            data: &exp,
            style: BLUE.stroke_width(2),
            dashed: true,
        },
        Series {
            name: "Van der Waals",
            data: &vdw,
            style: RED.into(),
            dashed: false,
        },
        Series {
            name: "Redlich-Kwong",
            data: &rk,
            style: YELLOW.into(),
            dashed: false,
        },
        Series {
            name: "Soave-Redlich-Kwong",
            data: &srk,
            style: GREEN.into(),
            dashed: false,
        },
        Series {
            name: "Peng-Robinson",
            data: &pr,
            style: CYAN.into(),
            dashed: false,
        },
        Series {
            name: "Patel-Teja-Valderrama",
            data: &ptv,
            style: MAGENTA.into(),
            dashed: false,
        },
    ];

    for t in plot_temps {
        let path = format!("bench/gen/z_{}_{}.png", gas_name, t);
        plot_bench(gas_name, series, *t, &path);
    }
}

fn plot_bench(gas_name: &str, series: &[Series], temperature: f64, path: &str) {
    use plotters::prelude::*;

    let caption = format!("Z factor of {} ({}K, experimental vs EoS)", gas_name, temperature);

    let pressures = series[0].data.pressures();
    for Series { data, .. } in series {
        debug_assert_eq!(data.pressures(), pressures);
    }

    let (z_min, z_max) = series
        .iter()
        .map(|Series { data, .. } | {
            let row = data.row(temperature).expect("No data for this temperature");
            (
                row.z.iter().cloned().fold(f64::INFINITY, f64::min),
                row.z.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            )
        })
        .fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(min1, max1), (min2, max2)| (min1.min(min2), max1.max(max2)),
        );

    let z_min = if z_min > 0.5 { 0.5 } else { 0.0 };
    let z_max = {
        let mut z = 1.0;
        while z < z_max {
            z += 0.5;
        }
        z
    };

    let pressures_range = 0.0..*pressures.last().unwrap() * 1e-5;
    let z_range = z_min..z_max; // Adjusted range for better visibility

    let root = BitMapBackend::new(path, (1200, 900)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption(caption.as_str(), ("sans-serif", 36))
        .margin(30)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(pressures_range, z_range)
        .unwrap();

    chart
        .configure_mesh()
        .x_desc("Pressure [bar]")
        .y_desc("Z factor")
        .label_style(("sans-serif", 20))
        .x_label_formatter(&|x| format!("{:.0}", x))
        .draw()
        .unwrap();

    for s in series {
        let row = s.data.row(temperature).expect("No data for this temperature");
        if s.dashed {
            chart
                .draw_series(DashedLineSeries::new(
                    pressures.iter().zip(row.z.iter()).map(|(p, z)| (p * 1e-5, *z)),
                    5, 5,
                    s.style,
                ))
                .unwrap()
                .label(s.name)
                .legend(|(x, y)| DashedPathElement::new(vec![(x, y), (x + 15, y)], 5, 5, s.style));
        } else {
            chart
                .draw_series(LineSeries::new(
                    pressures.iter().zip(row.z.iter()).map(|(p, z)| (p * 1e-5, *z)),
                    s.style.clone(),
                ))
                .unwrap()
                .label(s.name)
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], s.style));
        }
    }

    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::LowerRight)
        .label_font(("sans-serif", 20))
        .background_style(&WHITE)
        .border_style(&BLACK)
        .draw()
        .unwrap();

    root.present().expect("Failed to present the drawing area");

}
