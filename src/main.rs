mod estimate;
mod geode;
mod math;
mod noise;
mod search;
mod version;

use crate::{
    estimate::estimate,
    search::search_budding,
    version::{MC17, MC18, MC19},
};
use clap::{Parser, ValueEnum, ValueHint, value_parser};
use std::{io, path};

const WORLD_LIMIT: i64 = 30_000_000 / 16;
const WORLD_RANGE: std::ops::RangeInclusive<i64> = -WORLD_LIMIT..=WORLD_LIMIT;

#[derive(Debug, Copy, Clone, ValueEnum)]
pub enum VersionArgument {
    /// 1.17
    #[clap(name = "1.17")]
    MC17,

    /// 1.18
    #[clap(name = "1.18")]
    MC18,

    /// 1.19+
    #[clap(name = "1.19")]
    MC19,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about=None)]
struct Args {
    /// Minecraft version
    #[arg(short, long, value_enum, default_value_t = VersionArgument::MC19)]
    minecraft_version: VersionArgument,

    /// World seed
    #[arg(short, long, allow_hyphen_values = true, default_value_t = 0)]
    seed: i64,

    /// Search radius
    #[arg(short = 'r', long, default_value_t = 1000, value_parser = value_parser!(u32).range(1..=WORLD_LIMIT))]
    search_radius: u32,

    /// Minimum number of geodes per loaded area
    #[arg(short, long, default_value_t = 20, value_parser = value_parser!(u32).range(1..))]
    geode_threshold: u32,

    /// Minimum number of budding amethyst per loaded area
    #[arg(short, long, default_value_t = 800, value_parser = value_parser!(u32).range(1..))]
    budding_threshold: u32,

    /// Random tickable radius
    #[arg(long, default_value_t = 6, value_parser = value_parser!(u8).range(3..=64))]
    loaded_radius: u8,

    /// Search center chunk x
    #[arg(long, allow_negative_numbers = true, default_value_t = 0, value_parser = value_parser!(i64).range(WORLD_RANGE))]
    center_x: i64,

    /// Search center chunk z
    #[arg(long, allow_negative_numbers = true, default_value_t = 0, value_parser = value_parser!(i64).range(WORLD_RANGE))]
    center_z: i64,

    /// Where to save results
    #[arg(long, default_value = "output.json", default_missing_value = None, num_args=0..=1, value_hint = ValueHint::FilePath)]
    output_path: Option<path::PathBuf>,

    /// Estimate number of clusters without running search
    #[arg(long, default_value = "false", default_missing_value = "true")]
    estimate: bool,
}

fn main() -> io::Result<()> {
    let args = &Args::parse();
    let version = args.minecraft_version;
    if args.estimate {
        match version {
            VersionArgument::MC17 => estimate::<MC17>(args),
            VersionArgument::MC18 => estimate::<MC18>(args),
            VersionArgument::MC19 => estimate::<MC19>(args),
        };
    } else {
        match version {
            VersionArgument::MC17 => search_budding::<MC17>(args),
            VersionArgument::MC18 => search_budding::<MC18>(args),
            VersionArgument::MC19 => search_budding::<MC19>(args),
        }?;
    }

    Ok(())
}
