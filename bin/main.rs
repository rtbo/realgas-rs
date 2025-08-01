use std::{fmt, num::ParseFloatError, process::ExitCode, str::FromStr};

use clap::{Parser, Subcommand};
use realgas::{Gas, StateEos, eos::Eos};

/// Utility that performs real gas physics calculations.
#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Compute and print compressibility factor to stdout
    Z {
        /// Specify the gas to be used.
        #[arg(short = 'g', long)]
        gas: String,

        /// Equation of state used for computation
        #[arg(short='e', long, default_value_t=String::from("PR"))]
        eos: String,

        /// Specify the pressure or range of abs. pressure in bar
        #[arg(short = 'p', long)]
        pressure: String,

        /// Specify the pressure or range of temperature in °C
        #[arg(short = 't', long)]
        #[clap(allow_hyphen_values = true)]
        temperature: String,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{}", err);
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> anyhow::Result<()> {
    match &cli.command {
        Command::Z {
            gas,
            eos,
            pressure,
            temperature,
        } => {
            let gas: Gas = gas.parse()?;
            let eos: Eos = eos.parse()?;
            let pressure: Var = pressure.parse()?;
            let temperature: Var = temperature.parse()?;
            match (pressure, temperature) {
                (Var::Scalar(p), Var::Scalar(t)) => {
                    let p = p * 1e5;
                    let z = gas.z_eos(eos, p, t);
                    println!("{z}");
                }
                (p, t) => {
                    let p = p.to_vec();
                    let t = t.to_vec();
                    if t.first()
                        .copied()
                        .expect("Should have at least one temperature")
                        < -273.15
                    {
                        anyhow::bail!("Temperature below zero K !");
                    }
                    // print CSV header
                    print!("Temp");
                    for p in p.iter() {
                        print!(",{p}");
                    }
                    print!("\n");
                    for t in t.iter().copied() {
                        print!("{t}");
                        let t = t + 273.15;
                        for p in p.iter() {
                            let p = p * 1e5;
                            let z = gas.z_eos(eos, p, t);
                            print!(",{z}");
                        }
                        print!("\n");
                    }
                }
            }
        }
    }
    Ok(())
}

enum Var {
    Scalar(f64),
    Range {
        start: f64,
        end: f64,
        step: Option<f64>,
    },
}

impl Var {
    fn to_vec(&self) -> Vec<f64> {
        match self {
            &Var::Scalar(v) => vec![v],
            &Var::Range { start, end, step } => {
                let step = step.unwrap_or(1.0);
                let cap = ((end - start) / step) as usize;
                let mut res = Vec::with_capacity(cap);
                let mut v = start;
                while v <= (end + 2.0 * f64::EPSILON) {
                    res.push(v);
                    v += step;
                }
                res
            }
        }
    }
}

#[derive(Debug)]
enum ParseVarError {
    Empty,
    Float(ParseFloatError),
    Range(String),
}

impl fmt::Display for ParseVarError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseVarError::Empty => write!(f, "Can't parse variable from an empty string"),
            ParseVarError::Float(err) => err.fmt(f),
            ParseVarError::Range(msg) => msg.fmt(f),
        }
    }
}

impl From<ParseFloatError> for ParseVarError {
    fn from(value: ParseFloatError) -> Self {
        ParseVarError::Float(value)
    }
}

impl std::error::Error for ParseVarError {}

impl FromStr for Var {
    type Err = ParseVarError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseVarError::Empty);
        }

        let v = {
            let mut v: Vec<f64> = Vec::new();
            for s in s.split(':') {
                let n: f64 = s.parse()?;
                v.push(n);
            }
            v
        };

        match v.len() {
            1 => {
                let val = v[0];
                Ok(Var::Scalar(val))
            }
            2 => {
                let start = v[0];
                let end = v[1];
                if end <= start {
                    Err(ParseVarError::Range(
                        "Range end must be higher than start".into(),
                    ))
                } else {
                    Ok(Var::Range {
                        start,
                        end,
                        step: None,
                    })
                }
            }
            3 => {
                let start = v[0];
                let end = v[1];
                let step = v[2];
                if end <= start {
                    Err(ParseVarError::Range(
                        "Range stop must be higher than start".into(),
                    ))
                } else if step <= 0f64 {
                    Err(ParseVarError::Range("Range step must be positive".into()))
                } else {
                    Ok(Var::Range {
                        start,
                        end,
                        step: Some(step),
                    })
                }
            }
            _ => Err(ParseVarError::Range(format!(
                "Can't parse \"{}\" as a range",
                s
            ))),
        }
    }
}
