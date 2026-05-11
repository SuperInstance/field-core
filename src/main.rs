//! field-core CLI
//!
//! Command-line interface for interacting with a continuous constraint field.

use clap::{Parser, Subcommand};
use field_core::{
    ConstraintField, FieldPosition, FieldQuery, FieldStatus,
};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "field", version, about = "Continuous constraint field CLI")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Embed a constraint position in the field
    Embed {
        /// Position coordinates as comma-separated values (e.g., "1.0,2.0")
        #[arg(long)]
        position: String,
        /// Weight (confidence) of the position [0, 1]
        #[arg(long, default_value = "0.95")]
        weight: f64,
        /// Stiffness (spread) of the Gaussian
        #[arg(long, default_value = "200.0")]
        stiffness: f64,
        /// Path to field state file
        #[arg(long, default_value = ".field-state.json")]
        state: PathBuf,
    },
    /// Read the field at a query point
    Read {
        /// Query coordinates as comma-separated values
        #[arg(long)]
        query: String,
        /// Time offset for confidence decay
        #[arg(long, default_value = "0")]
        time: f64,
        /// Tolerance for nearby position matching
        #[arg(long, default_value = "3.0")]
        tolerance: f64,
        /// Path to field state file
        #[arg(long, default_value = ".field-state.json")]
        state: PathBuf,
    },
    /// Run one propagation step
    Propagate {
        /// Path to field state file
        #[arg(long, default_value = ".field-state.json")]
        state: PathBuf,
    },
    /// Show field topology state
    Topology {
        /// Path to field state file
        #[arg(long, default_value = ".field-state.json")]
        state: PathBuf,
    },
    /// List all field positions
    Nails {
        /// Path to field state file
        #[arg(long, default_value = ".field-state.json")]
        state: PathBuf,
    },
    /// Show field health status
    Status {
        /// Path to field state file
        #[arg(long, default_value = ".field-state.json")]
        state: PathBuf,
    },
}

fn load_field(path: &PathBuf) -> ConstraintField {
    if path.exists() {
        let data = std::fs::read_to_string(path)
            .unwrap_or_else(|e| {
                eprintln!("Warning: could not read state file: {e}");
                String::new()
            });
        if !data.is_empty() {
            serde_json::from_str(&data).unwrap_or_else(|e| {
                eprintln!("Warning: could not parse state file: {e}");
                ConstraintField::new()
            })
        } else {
            ConstraintField::new()
        }
    } else {
        ConstraintField::new()
    }
}

fn save_field(field: &ConstraintField, path: &PathBuf) {
    let data = serde_json::to_string_pretty(field).expect("serialization failed");
    std::fs::write(path, data).unwrap_or_else(|e| {
        eprintln!("Warning: could not write state file: {e}");
    });
}

fn parse_coords(s: &str) -> Vec<f64> {
    s.split(',')
        .map(|x| x.trim().parse().unwrap_or_else(|e| {
            eprintln!("Error: invalid coordinate '{x}': {e}");
            std::process::exit(1);
        }))
        .collect()
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Embed { position, weight, stiffness, state } => {
            let mut field = load_field(&state);
            let coords = parse_coords(&position);
            let tile = field.embed(coords, weight, stiffness);
            save_field(&field, &state);
            println!("{}", serde_json::to_string_pretty(&tile).unwrap());
        }
        Commands::Read { query, time, tolerance, state } => {
            let field = load_field(&state);
            let point = parse_coords(&query);
            let q = FieldQuery { point, time, tolerance };
            let reading = field.read(&q);
            println!("{}", serde_json::to_string_pretty(&reading).unwrap());
        }
        Commands::Propagate { state } => {
            let mut field = load_field(&state);
            let result = field.propagate();
            save_field(&field, &state);
            println!("{}", serde_json::to_string_pretty(&result).unwrap());
        }
        Commands::Topology { state } => {
            let field = load_field(&state);
            let topo = field.topology();
            println!("{}", serde_json::to_string_pretty(&topo).unwrap());
        }
        Commands::Nails { state } => {
            let field = load_field(&state);
            let nails: Vec<&FieldPosition> = field.positions().iter().collect();
            println!("{}", serde_json::to_string_pretty(&nails).unwrap());
        }
        Commands::Status { state } => {
            let field = load_field(&state);
            let status: FieldStatus = (&field).into();
            println!("{}", serde_json::to_string_pretty(&status).unwrap());
        }
    }
}
