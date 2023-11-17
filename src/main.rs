#![allow(warnings)]
use clap::{Parser, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashMap;

mod geode;
mod noise;
mod random;
mod search;

use geode::Geode;

#[derive(Debug, Copy, Clone, ValueEnum)]
pub enum GameVersion {
    /// 1.17
    #[clap(name = "17")]
    MC17,

    /// 1.18+
    #[clap(name = "18")]
    MC18,

    /// 1.20+
    #[clap(name = "20")]
    MC20,

    /// 1.17 top and 1.18+ bottom
    #[clap(name = "merged")]
    MCMerged,
}

// CLI arguments with clap
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Game version
    #[arg(value_enum, short = 'v', long, default_value_t = GameVersion::MC18)]
    game_version: GameVersion,

    /// World seed
    #[arg(short, long, allow_hyphen_values = true, default_value_t = 0)]
    seed: i64,

    /// Search radius
    #[arg(short = 'r', long, default_value_t = 1000)]
    search_radius: usize,

    /// Minimum number of geodes per area
    #[arg(short, long, default_value_t = 25)]
    geode_threshold: u8,

    /// Minimum number of budding amethyst per area
    #[arg(short, long, default_value_t = 900)]
    amethyst_threshold: u32,

    /// x coordinate of center chunk
    #[arg(long, default_value_t = 0)]
    center_x: i32,

    /// z coordinate of center chunk
    #[arg(long, default_value_t = 0)]
    center_z: i32,

    /// Number of threads to use
    #[arg(long, default_value_t = 1)]
    threads: u8,

    /// Search Mode
    #[arg(long, default_value_t = 1)]
    mode: u8,
}

fn initialize_progress_bar(length: u64) -> ProgressBar {
    let progress_style = ProgressStyle::default_spinner()
        .progress_chars("*-")
        .template("{spinner:.green} [{elapsed}] [{bar:.green/white}] ({eta_precise})")
        .unwrap();
    ProgressBar::new(length).with_style(progress_style)
}

fn main() {
    let args = Args::parse();

    let seed = args.seed;
    let search_radius = args.search_radius as i64;
    let geode_threshold = args.geode_threshold;
    let budding_threshold = args.amethyst_threshold;

    let mut cached_geodes = HashMap::new();

    let mut finder = Geode::new(seed, args.game_version);
    let mut locations = search(&mut finder, args);

    for loc in locations {
        let min_x = loc.0 - 12;
        let max_x = loc.0 + 1;
        let min_z = loc.1 - 12;
        let max_z = loc.1 + 1;

        let mut area_budding_count = 0;
        for i in min_x..max_x {
            for j in min_z..max_z {
                let count = match cached_geodes.get(&(i, j)) {
                    Some(cached_count) => {
                        *cached_count
                    }

                    None => {
                        let geode_budding_count = finder.generate(i, j);
                        cached_geodes.insert((i, j), geode_budding_count);
                        geode_budding_count
                    }
                };

                area_budding_count += count;
            }
        }

        if area_budding_count >= budding_threshold as i32 {
            println!(
                "Geode cluster with {} budding amethyst centered at {} {}",
                area_budding_count,
                (loc.0 - 6) * 16,
                (loc.1 - 6) * 16
            );
        }
    }
}

fn search(finder: &mut Geode, args: Args) -> Vec<(i64, i64)> {
    let seed = args.seed;
    let search_radius = args.search_radius as i64;
    let geode_threshold = args.geode_threshold;


    let search_length = search_radius as usize * 2 + 1;
    let progress_bar = initialize_progress_bar(search_length as u64);

    let mut locations: Vec<(i64, i64)> = vec![];

    let mut sum_index: usize = 0;
    let mut previous_sums = vec![vec![0; search_length]; 13];
    let mut current_sums = vec![0u8; search_length];

    for i in -search_radius..=search_radius {
        let mut slice = [0u8; 13];
        let mut slice_index = 0;
        let mut slice_sum = 0u8;
        let sum_slice = &mut previous_sums[sum_index];

        for j in -search_radius..=search_radius {
            let is_geode = finder.check_chunk(i, j) as u8;

            slice_sum += is_geode;
            slice_sum -= slice[slice_index];

            slice[slice_index] = is_geode;
            slice_index = (slice_index + 1) % 13;

            let k = (j + search_radius) as usize;
            let slice_u16 = slice_sum as u8;
            current_sums[k] += slice_u16 - sum_slice[k];
            sum_slice[k] = slice_u16;

            if current_sums[k] >= geode_threshold {
                locations.push((i, j));
                println!("Found region with {} geodes", current_sums[k]);
            }
        }

        sum_index = (sum_index + 1) % 13;
        progress_bar.inc(1);
    }

    progress_bar.finish();
    println!("{} potential locations found.", locations.len());

    locations
}
