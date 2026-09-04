mod estimate;
mod geode;
mod math;
mod noise;
mod search;
mod version;

use std::{ops::RangeInclusive, path::PathBuf, thread};

use anyhow::Result;
use clap::{Parser, ValueEnum, ValueHint, value_parser};

use crate::{
    estimate::estimate_clusters,
    search::search,
    version::{MC17, MC18, MC182, MC194, MC263},
};

const WORLD_BORDER: i64 = 30_000_000 / 16;
const BORDER_RANGE: RangeInclusive<i64> = -WORLD_BORDER..=WORLD_BORDER;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum VersionArgument {
    /// 1.17-1.17.1
    #[clap(name = "1.17")]
    MC17,

    /// 1.18-1.18.1
    #[clap(name = "1.18")]
    MC18,

    /// 1.18.2-1.19.3
    #[clap(name = "1.18.2")]
    MC182,

    /// 1.19.4-26.2
    #[clap(name = "1.19.4")]
    MC194,

    /// 26.3+
    #[clap(name = "26.3")]
    MC263,
}

fn thread_parser(threads_str: &str) -> Result<usize> {
    let threads: usize = threads_str.parse()?;
    if threads == 0 {
        return Ok(thread::available_parallelism()?.get());
    }
    Ok(threads)
}

#[derive(Debug, Parser)]
#[command(author, version, about, long_about=None)]
struct Args {
    /// Minecraft version
    #[arg(short, long, value_enum, default_value_t = Self::default().minecraft_version)]
    minecraft_version: VersionArgument,

    /// World seed
    #[arg(short, long, allow_hyphen_values = true, default_value_t = Self::default().seed)]
    seed: i64,

    /// Search radius
    #[arg(short = 'r', long, default_value_t = Self::default().search_radius, value_parser = value_parser!(u32).range(1..=WORLD_BORDER))]
    search_radius: u32,

    /// Minimum number of geodes per loaded area
    #[arg(short, long, default_value_t = Self::default().geode_threshold, value_parser = value_parser!(i16).range(1..))]
    geode_threshold: i16,

    /// Minimum number of budding amethyst per loaded area
    #[arg(short, long, default_value_t = Self::default().budding_threshold, value_parser = value_parser!(u32).range(1..))]
    budding_threshold: u32,

    /// Number of threads to use (0 will use all cores)
    #[arg(long, default_value_t = Self::default().threads, value_parser = thread_parser)]
    threads: usize,

    /// Random tickable radius
    #[arg(long, default_value_t = Self::default().loaded_radius, value_parser = value_parser!(u8).range(0..=64))]
    loaded_radius: u8,

    /// Search center chunk x
    #[arg(long, allow_negative_numbers = true, default_value_t = Self::default().center_x, value_parser = value_parser!(i64).range(BORDER_RANGE))]
    center_x: i64,

    /// Search center chunk z
    #[arg(long, allow_negative_numbers = true, default_value_t = Self::default().center_z, value_parser = value_parser!(i64).range(BORDER_RANGE))]
    center_z: i64,

    /// Where to save results
    #[arg(long, default_value_os = Self::DEFAULT_OUTPUT_PATH, default_missing_value = None, num_args=0..=1, value_hint = ValueHint::FilePath)]
    output_path: Option<PathBuf>,

    /// Estimate number of clusters without running search
    #[arg(long, default_value_t = Self::default().estimate, default_missing_value = "true")]
    estimate: bool,
}

impl Args {
    const DEFAULT_OUTPUT_PATH: &str = "output.json";
}

impl Default for Args {
    fn default() -> Self {
        Self {
            minecraft_version: VersionArgument::MC263,
            seed: 0,
            search_radius: 1000,
            geode_threshold: 20,
            budding_threshold: 800,
            threads: 1,
            loaded_radius: 6,
            center_x: 0,
            center_z: 0,
            output_path: Some(PathBuf::from(Args::DEFAULT_OUTPUT_PATH)),
            estimate: false,
        }
    }
}

fn main() -> Result<()> {
    let args = &Args::parse();
    let version = args.minecraft_version;

    if args.estimate {
        match version {
            VersionArgument::MC17 => estimate_clusters::<MC17>(args),
            VersionArgument::MC18
            | VersionArgument::MC182
            | VersionArgument::MC194
            | VersionArgument::MC263 => estimate_clusters::<MC18>(args),
        }
    } else {
        match version {
            VersionArgument::MC17 => search::<MC17>(args),
            VersionArgument::MC18 => search::<MC18>(args),
            VersionArgument::MC182 => search::<MC182>(args),
            VersionArgument::MC194 => search::<MC194>(args),
            VersionArgument::MC263 => search::<MC263>(args),
        }?;
    }

    Ok(())
}
