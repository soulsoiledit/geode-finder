use std::{
    cmp,
    collections::HashMap,
    fmt,
    fs::File,
    io::BufWriter,
    ops,
    sync::{
        Arc,
        atomic::{self, AtomicU64},
    },
    thread, time,
};

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::{
    ThreadPool, ThreadPoolBuilder,
    iter::{IntoParallelIterator, ParallelIterator},
    slice::{ParallelSlice, ParallelSliceMut},
};
use serde::Serialize;

use crate::{Args, geode::Geode, version::Version};

#[derive(Debug, Serialize, PartialEq, PartialOrd, Eq, Ord)]
pub struct GeodeCluster {
    geode_count: u32,
    center_x: i64,
    center_z: i64,
}

#[derive(Debug, Serialize, PartialEq, PartialOrd, Eq, Ord)]
pub struct BuddingCluster {
    budding_count: u32,
    geode_count: u32,
    center_x: i64,
    center_z: i64,
}

impl fmt::Display for BuddingCluster {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} budding amethyst centered at ({}, {})",
            self.budding_count, self.center_x, self.center_z,
        )
    }
}

fn build_configured_pool(threads: usize) -> Result<ThreadPool> {
    Ok(ThreadPoolBuilder::new().num_threads(threads).build()?)
}

#[derive(Clone)]
struct SharedSearchConfig {
    seed: i64,
    loaded_radius: u8,
    geode_threshold: u32,
    progress_position: Arc<AtomicU64>,
    total_clusters: Arc<AtomicU64>,
    z_range: ops::Range<i64>,
}

const PROGRESS_INTERVAL: u16 = 256;
fn search_geodes_tile<V: Version>(
    shared: SharedSearchConfig,
    start_x: i64,
    end_x: i64,
) -> Vec<GeodeCluster> {
    let geode = Geode::<V>::new(shared.seed);
    let mut clusters: Vec<GeodeCluster> = Vec::with_capacity(256);

    let loaded_radius = i64::from(shared.loaded_radius);
    let loaded_diameter = usize::from(shared.loaded_radius) * 2 + 1;

    let tile_width = usize::try_from((end_x - start_x).abs()).unwrap();
    let mut chunk_history: Vec<u8> = vec![0; tile_width * loaded_diameter];
    // Add loaded_diameter to accomodate buffer zone and index without branch later
    let mut column_history: Vec<u32> = vec![0; tile_width + loaded_diameter];

    let mut chunk_history_index = 0;
    let chunk_history_reset = chunk_history.len();
    for (idz, z) in shared.z_range.enumerate() {
        let center_z = z - loaded_radius;
        let scaled_z = geode.scale_z(z);

        // Slice here to keep current row in cache
        let chunk_history_slice =
            &mut chunk_history[chunk_history_index..chunk_history_index + tile_width];

        let mut geode_count = 0;
        for (idx, x) in (start_x..end_x).enumerate() {
            let is_geode = u8::from(geode.check_fast(x, scaled_z));

            let old_geode = chunk_history_slice[idx];
            chunk_history_slice[idx] = is_geode;

            let new_column =
                column_history[idx + loaded_diameter] + u32::from(is_geode) - u32::from(old_geode);
            column_history[idx + loaded_diameter] = new_column;

            geode_count += new_column;
            geode_count -= column_history[idx];
            if geode_count >= shared.geode_threshold {
                clusters.push(GeodeCluster {
                    center_x: x - loaded_radius,
                    center_z,
                    geode_count,
                });
                shared
                    .total_clusters
                    .fetch_add(1, atomic::Ordering::Relaxed);
            }
        }

        chunk_history_index += tile_width;
        if chunk_history_index == chunk_history_reset {
            chunk_history_index = 0;
        }

        if idz % usize::from(PROGRESS_INTERVAL) == 0 {
            shared
                .progress_position
                .fetch_add(u64::from(PROGRESS_INTERVAL), atomic::Ordering::Relaxed);
        }
    }

    clusters
}

fn search_geodes<V: Version>(args: &crate::Args) -> Result<Vec<GeodeCluster>> {
    let search_radius = usize::try_from(args.search_radius)?;
    let loaded_radius = args.loaded_radius;
    let geode_threshold = args.geode_threshold;

    let search_diameter = search_radius * 2 + 1;
    let loaded_diameter = i64::from(loaded_radius) * 2 + 1;

    let start_x = args.center_x - i64::from(args.search_radius);
    let start_z = args.center_z - i64::from(args.search_radius);

    let end_x = start_x + search_diameter as i64;
    let end_z = start_z + search_diameter as i64;

    let progress_position = Arc::new(AtomicU64::new(0));
    let total_clusters = Arc::new(AtomicU64::new(0));
    let progress_bar = if cfg!(test) {
        ProgressBar::hidden()
    } else {
        let style = ProgressStyle::with_template(
            "[{bar}] {percent_precise}% ({eta_precise} left), {msg} potential clusters found ",
        )?;
        ProgressBar::new(u64::try_from(search_diameter * args.threads).unwrap())
            .with_style(style)
            .with_message("0")
    };

    // Run UI on a different thread
    let handle = {
        let progress_position = progress_position.clone();
        let total_clusters = total_clusters.clone();
        let progress_bar = progress_bar.clone();

        thread::spawn(move || {
            loop {
                progress_bar.set_position(progress_position.load(atomic::Ordering::Relaxed));
                progress_bar
                    .set_message(total_clusters.load(atomic::Ordering::Relaxed).to_string());
                thread::sleep(time::Duration::from_millis(100));

                if progress_bar.is_finished() {
                    break;
                }
            }
        })
    };

    let shared = SharedSearchConfig {
        seed: args.seed,
        loaded_radius,
        geode_threshold,
        progress_position,
        total_clusters,
        z_range: start_z..end_z,
    };

    let chunk_size = search_diameter.div_ceil(args.threads) as i64;
    let geodes: Vec<GeodeCluster> = (0..args.threads)
        .into_par_iter()
        .flat_map(|i| {
            let tile_start_x = start_x + (i as i64 * chunk_size);
            if tile_start_x >= end_x {
                return vec![];
            }
            let tile_end_x = (tile_start_x + chunk_size + loaded_diameter).min(end_x);

            search_geodes_tile::<V>(shared.clone(), tile_start_x, tile_end_x)
        })
        .collect();

    progress_bar.finish();
    handle
        .join()
        .map_err(|e| anyhow::anyhow!("UI thread panicked: {:?}", e))?;

    Ok(geodes)
}

fn process_clusters<V: Version>(
    geode_clusters: &[GeodeCluster],
    seed: i64,
    loaded_radius: i64,
    budding_threshold: u32,
) -> Vec<BuddingCluster> {
    let mut geode = Geode::<V>::new(seed);
    let mut budding_clusters: Vec<BuddingCluster> = Vec::with_capacity(256);
    let loaded_area = usize::try_from(loaded_radius).unwrap().pow(2);
    let mut cached_budding: HashMap<(i64, i64), u32> = HashMap::with_capacity(loaded_area * 2);

    for cluster in geode_clusters {
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

        if budding_count >= budding_threshold {
            let center_x_blocks = center_x * 16;
            let center_z_blocks = center_y * 16;
            let cluster = BuddingCluster {
                center_x: center_x_blocks,
                center_z: center_z_blocks,
                geode_count: cluster.geode_count,
                budding_count,
            };
            println!("{}", cluster);
            budding_clusters.push(cluster)
        }
    }

    budding_clusters
}

fn search_budding<V: Version>(args: &Args) -> Result<Vec<BuddingCluster>> {
    let seed = args.seed;
    let loaded_radius = i64::from(args.loaded_radius);
    let budding_threshold = args.budding_threshold;

    let geode_clusters = search_geodes::<V>(args)?;

    let chunk_size = geode_clusters.len().div_ceil(args.threads).max(1);
    let mut budding_clusters: Vec<BuddingCluster> = geode_clusters
        .par_chunks(chunk_size)
        .flat_map(|clusters| {
            process_clusters::<V>(clusters, seed, loaded_radius, budding_threshold)
        })
        .collect();

    budding_clusters
        .par_sort_unstable_by_key(|a| (cmp::Reverse(a.budding_count), a.center_z, a.center_x));

    Ok(budding_clusters)
}

pub fn search<V: Version>(args: &Args) -> Result<Vec<BuddingCluster>> {
    let writer = args
        .output_path
        .as_ref()
        .map(File::create)
        .transpose()?
        .map(BufWriter::new);

    let budding_clusters =
        build_configured_pool(args.threads)?.install(|| search_budding::<V>(args))?;

    if let Some(w) = writer {
        serde_json::to_writer_pretty(w, &budding_clusters)?;
    }

    Ok(budding_clusters)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{VersionArgument, version::MC19};

    fn test_args() -> Args {
        Args {
            minecraft_version: VersionArgument::MC19,
            output_path: None,
            ..Default::default()
        }
    }

    #[test]
    fn search_geodes() -> Result<()> {
        let args = &test_args();
        let pool = build_configured_pool(args.threads)?;
        let clusters = pool.install(|| super::search_geodes::<MC19>(args))?;

        assert_eq!(clusters.len(), 189);
        assert_eq!(clusters.iter().max().map_or(0, |c| c.geode_count), 23);
        Ok(())
    }

    #[test]
    fn search() -> Result<()> {
        let clusters = super::search::<MC19>(&test_args())?;
        assert_eq!(clusters.len(), 40);
        assert_eq!(clusters.iter().max().map_or(0, |c| c.budding_count), 911);
        Ok(())
    }
}
