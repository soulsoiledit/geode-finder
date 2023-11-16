use clap::{Parser, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashMap;
//
mod geode;
mod noise;
mod random;
use geode::GeodeGenerator;

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
    geode_threshold: u16,

    /// Minimum number of budding amethyst per area
    #[arg(short, long, default_value_t = 1000)]
    amethyst_threshold: u16,

    /// x coordinate of center chunk
    #[arg(long, default_value_t = 0)]
    center_x: i32,

    /// z coordinate of center chunk
    #[arg(long, default_value_t = 0)]
    center_z: i32,

    /// Number of threads to use
    #[arg(long, default_value_t = 1)]
    threads: u8,

    // /// Search Mode
    // #[arg(long, default_value_t = 1)]
    // threads: u8,
}

fn main() {
    let args = Args::parse();
    let seed = args.seed;
    let search_radius = args.search_radius;
    let geode_threshold = args.geode_threshold;
    let budding_threshold = args.budding_threshold;
    let is_17 = args.game_version.as_str() == "1.17";
    let salt: i64 = if is_17 { 20000 } else { 20002 };

    if is_17 {
        println!("Running geode search for version 1.17...");
    } else {
        println!("Running geode search for versions 1.18+...");
    }

    let mut finder = GeodeGenerator::new(seed, is_17);
    let mut locations: Vec<(i64, i64)> = vec![];

    {
        let progress_style = ProgressStyle::with_template(
            "{spinner:.green} [{elapsed}] [{bar:.green/white}] ({eta_precise})",
        )
        .unwrap()
        .progress_chars("ùwú");
        let progress_bar = ProgressBar::new(search_radius as u64 * 2).with_style(progress_style);

        let ignored_columns = -search_radius + 13;
        let mut sum_index: usize = 0;
        let mut previous_sums = vec![vec![0; (search_radius * 2) as usize]; 15];
        let mut current_sums = vec![0; (search_radius * 2) as usize];

        for i in -search_radius..search_radius {
            let mut slice = [0; 15];
            let mut slice_index = 0;
            let mut slice_sum = 0;

            for j in -search_radius..search_radius {
                let k = (j + search_radius) as usize;

                slice_sum -= slice[slice_index];
                slice[slice_index] = 0;

                finder.set_decorator_seed(i as i64, j as i64, salt);

                slice[slice_index] += finder.check_chunk() as i32;
                slice_sum += slice[slice_index];
                slice_index = (slice_index + 1) % 15;

                current_sums[k] -= previous_sums[sum_index][k];
                current_sums[k] += slice_sum;
                previous_sums[sum_index][k] = slice_sum;

                if i > ignored_columns && j > ignored_columns && current_sums[k] >= geode_threshold
                {
                    locations.push((i as i64, j as i64));
                }
            }

            sum_index = (sum_index + 1) % 15;
            progress_bar.inc(1);
        }
        progress_bar.finish();
    }
    println!("{}", locations.len());

    let mut geodes = HashMap::new();
    for loc in locations {
        let min_x = loc.0 - 14;
        let min_z = loc.1 - 14;

        let mut area_budding_count = 0;
        for i in min_x..min_x + 15 {
            for j in min_z..min_z + 15 {
                if geodes.contains_key(&(i, j)) {
                    area_budding_count += geodes[&(i, j)];
                } else {
                    finder.set_decorator_seed(i, j, salt);
                    if finder.check_chunk() {
                        let geode_budding_count = finder.generate(i, j);
                        area_budding_count += geode_budding_count;
                        geodes.insert((i, j), geode_budding_count);
                    }
                }
            }
        }

        if area_budding_count >= budding_threshold {
            println!(
                "Geode cluster with {} budding amethyst centered at {} {}",
                area_budding_count,
                (loc.0 - 7) * 16,
                (loc.1 - 7) * 16
            );
        }
    }
}
