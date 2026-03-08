mod geode;
mod math;
mod noise;
mod version;

use crate::{
    geode::Geode,
    version::{MC17, MC18, MC19, Version},
};
use clap::{Parser, ValueEnum, ValueHint::FilePath, builder::RangedI64ValueParser, value_parser};
use indicatif::{ProgressBar, ProgressStyle};
use serde::Serialize;
use serde_json::to_writer_pretty;
use statrs::distribution::{Binomial, DiscreteCDF};
use std::{collections::HashMap, fmt, fs, io, path};

const WORLD_LIMIT: i64 = 30_000_000 / 16;

#[derive(Debug, Copy, Clone, ValueEnum)]
pub enum MinecraftVersion {
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

macro_rules! versioned {
    ($ver:expr, $fn:ident, $($x:expr),*) => {
        match $ver {
            MinecraftVersion::MC17 => $fn::<MC17>($($x),*),
            MinecraftVersion::MC18 => $fn::<MC18>($($x),*),
            MinecraftVersion::MC19 => $fn::<MC19>($($x),*),
        }
    };
}

fn parse_chunk_coord() -> RangedI64ValueParser {
    value_parser!(i64).range(-WORLD_LIMIT..=WORLD_LIMIT)
}

fn parse_search_radius() -> RangedI64ValueParser<u32> {
    value_parser!(u32).range(1..=WORLD_LIMIT)
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about=None)]
struct Args {
    /// Minecraft version
    #[arg(short, long, value_enum, default_value_t = MinecraftVersion::MC19)]
    minecraft_version: MinecraftVersion,

    /// World seed
    #[arg(short, long, allow_hyphen_values = true, default_value_t = 0)]
    seed: i64,

    /// Search radius
    #[arg(short = 'r', long, default_value_t = 1000, value_parser = parse_search_radius())]
    search_radius: u32,

    /// Minimum number of geodes per loaded area
    #[arg(short, long, default_value_t = 20)]
    geode_threshold: u32,

    /// Minimum number of budding amethyst per loaded area
    #[arg(short, long, default_value_t = 800)]
    budding_threshold: u32,

    /// Random tickable radius
    #[arg(long, default_value_t = 6)]
    loaded_radius: u16,

    /// Search center chunk x
    #[arg(long, allow_negative_numbers = true, default_value_t = 0, value_parser = parse_chunk_coord())]
    center_x: i64,

    /// Search center chunk z
    #[arg(long, allow_negative_numbers = true, default_value_t = 0, value_parser = parse_chunk_coord())]
    center_z: i64,

    /// Where to save results
    #[arg(long, default_value = "output.json", default_missing_value = None, num_args=0..=1, value_hint = FilePath)]
    output_path: Option<path::PathBuf>,

    /// Estimate number of clusters without running search
    #[arg(long, default_value = "false", default_missing_value = "true")]
    estimate: bool,
}

#[derive(Debug, Serialize, PartialEq, PartialOrd, Eq, Ord)]
struct Cluster {
    budding_count: u32,
    geode_count: u32,
    center_x: i64,
    center_z: i64,
}

impl fmt::Display for Cluster {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} budding amethyst centered at ({}, {})",
            self.budding_count, self.center_x, self.center_z,
        )
    }
}

fn main() {
    let args = Args::parse();
    if args.estimate {
        versioned!(args.minecraft_version, estimate, &args);
    } else {
        let _ = versioned!(args.minecraft_version, search_budding, &args);
    }
}

fn estimate<V: Version>(args: &Args) {
    // Determined empirically
    const AVERAGE_BUDDING: f64 = 35.875;

    let loaded_area = (u64::from(args.loaded_radius) * 2 + 1).pow(2);
    let total_search_area = (f64::from(args.search_radius) * 2.0 + 1.0).powi(2);

    let binomial = Binomial::new(f64::from(V::CHANCE), loaded_area).unwrap();

    let geode_cluster_chance = 1.0 - binomial.cdf(u64::from(args.geode_threshold) - 1);
    let expected_geode_clusters = (total_search_area * geode_cluster_chance).floor();

    let required_geodes = (f64::from(args.budding_threshold) / AVERAGE_BUDDING).floor() as u64;
    let budding_cluster_chance = 1.0 - binomial.cdf(required_geodes - 1);
    let expected_budding_clusters = (total_search_area * budding_cluster_chance)
        .min(expected_geode_clusters)
        .floor();

    println!("Estimated geode clusters: {expected_geode_clusters}");
    println!("Estimated budding amethyst clusters: {expected_budding_clusters}");
}

fn search_geodes<V: Version>(args: &Args) -> Vec<Cluster> {
    let search_diameter = args.search_radius as usize * 2 + 1;
    let loaded_radius = i64::from(args.loaded_radius);
    let loaded_diameter = usize::from(args.loaded_radius) * 2 + 1;

    let start_x = args.center_x - i64::from(args.search_radius);
    let start_z = args.center_z - i64::from(args.search_radius);

    let mut geode = Geode::<V>::new(args.seed);
    let mut history: Vec<u32> = vec![0; search_diameter * loaded_diameter];
    let mut column_sum: Vec<u32> = vec![0; search_diameter];
    let mut clusters: Vec<Cluster> = vec![];

    let progress_bar = if cfg!(test) {
        ProgressBar::hidden()
    } else {
        ProgressBar::new(search_diameter as u64).with_style(
            ProgressStyle::with_template(
                "[{bar}] {msg} potential clusters found ({eta_precise} left)",
            )
            .unwrap(),
        )
    };
    progress_bar.set_message(clusters.len().to_string());

    for idz in progress_bar.wrap_iter(0..search_diameter) {
        let z = start_z + idz as i64;

        let mut geode_count = 0;
        let slice: usize = (idz % loaded_diameter) * search_diameter;
        let current: &mut [u32] = &mut history[slice..slice + search_diameter];

        for idx in 0..search_diameter {
            let x = start_x + idx as i64;

            let old_geode = current[idx];
            let column = &mut column_sum[idx];
            let is_geode = geode.check_fast(x, z);

            *column += is_geode;
            *column -= old_geode;
            geode_count += *column;
            current[idx] = is_geode;

            if let Some(old_column) = idx.checked_sub(loaded_diameter) {
                geode_count -= column_sum[old_column];
                if geode_count >= args.geode_threshold {
                    clusters.push(Cluster {
                        center_x: x - loaded_radius,
                        center_z: z - loaded_radius,
                        geode_count,
                        budding_count: 0,
                    });
                }
            }
        }

        progress_bar.set_message(clusters.len().to_string());
    }

    progress_bar.finish();

    clusters
}

fn save_results(path: &Option<&path::Path>, data: &[Cluster]) -> io::Result<()> {
    if let Some(path) = path {
        let file = fs::File::create(path)?;
        let writer = io::BufWriter::new(file);
        to_writer_pretty(writer, data)?;
    }
    Ok(())
}

fn search_budding<V: Version>(args: &Args) -> std::io::Result<Vec<Cluster>> {
    let loaded_radius = i64::from(args.loaded_radius);
    let path = &args.output_path.as_deref();

    let mut geode = Geode::<V>::new(args.seed);
    let geode_clusters = search_geodes::<V>(args);
    let save_interval = (geode_clusters.len() / 10).max(1);

    let mut cached_budding: HashMap<(i64, i64), u32> = HashMap::with_capacity(1000);
    let mut budding_clusters: Vec<Cluster> = Vec::with_capacity(100);

    for (i, cluster) in geode_clusters.iter().enumerate() {
        let center_x = cluster.center_x;
        let center_y = cluster.center_z;

        let min_x = center_x - loaded_radius;
        let max_x = center_x + loaded_radius;
        let min_z = center_y - loaded_radius;
        let max_z = center_y + loaded_radius;

        let mut budding_count = 0;
        for x in min_x..=max_x {
            for z in min_z..=max_z {
                budding_count += *cached_budding
                    .entry((x, z))
                    .or_insert_with(|| geode.generate(x, z));
            }
        }

        if budding_count >= args.budding_threshold {
            let center_x_blocks = center_x * 16;
            let center_z_blocks = center_y * 16;
            let cluster = Cluster {
                center_x: center_x_blocks,
                center_z: center_z_blocks,
                geode_count: cluster.geode_count,
                budding_count,
            };
            eprintln!("{}", cluster);
            budding_clusters.push(cluster)
        }

        if i > 0 && i % save_interval == 0 {
            save_results(path, &budding_clusters)?;
        }
    }

    save_results(path, &budding_clusters)?;
    Ok(budding_clusters)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_args() -> Args {
        Args {
            minecraft_version: MinecraftVersion::MC19,
            seed: 0,
            search_radius: 1000,
            geode_threshold: 20,
            budding_threshold: 800,
            output_path: None,
            loaded_radius: 6,
            center_x: 0,
            center_z: 0,
            estimate: false,
        }
    }

    #[test]
    fn test_search_geodes() {
        let clusters = search_geodes::<MC19>(&test_args());
        println!("{}", clusters.iter().max().unwrap());
        assert_eq!(clusters.len(), 189);
        assert_eq!(clusters.iter().max().unwrap().geode_count, 23);
    }

    #[test]
    fn test_search_budding() {
        if let Ok(clusters) = search_budding::<MC19>(&test_args()) {
            assert_eq!(clusters.len(), 40);
            assert_eq!(clusters.iter().max().unwrap().budding_count, 911);
        }
    }
}
