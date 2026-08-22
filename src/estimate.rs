use crate::{Args, WORLD_BORDER, version::Version};

/// P(X >= threshold)
// Poisson provides a close enough approximation (n >= 49 > 20, p = 1/24 or 1/53 <= 0.05)
fn poisson_tail(trials: f64, probability: f64, threshold: u32) -> f64 {
    if threshold == 0 {
        return 1.0;
    }

    let lambda = trials * probability;
    let mut current = (-lambda).exp();
    let mut cdf = current;

    for i in (1..threshold).map(f64::from) {
        current *= lambda / i;
        cdf += current;
    }

    (1.0 - cdf).max(0.0)
}

pub fn estimate_clusters<V: Version>(args: &Args) {
    // Determined empirically
    const AVERAGE_BUDDING: f64 = 35.875;

    let loaded_area = f64::from(args.loaded_radius).mul_add(2.0, 1.0).powi(2);
    let search_area = {
        let start_x = (args.center_x - i64::from(args.search_radius)).max(-WORLD_BORDER);
        let start_z = (args.center_z - i64::from(args.search_radius)).max(-WORLD_BORDER);
        let end_x = (args.center_x + i64::from(args.search_radius)).min(WORLD_BORDER);
        let end_z = (args.center_z + i64::from(args.search_radius)).min(WORLD_BORDER);
        ((start_x).abs_diff(end_x) * (start_z).abs_diff(end_z)) as f64
    };
    let probability = f64::from(V::CHANCE);

    let geode_cluster_chance = poisson_tail(loaded_area, probability, args.geode_threshold);
    let expected_geode_clusters = (search_area * geode_cluster_chance).floor();

    let required_geodes = (f64::from(args.budding_threshold) / AVERAGE_BUDDING).ceil() as u32;
    let budding_cluster_chance = poisson_tail(loaded_area, probability, required_geodes);
    let expected_budding_clusters = (search_area * budding_cluster_chance)
        .floor()
        .min(expected_geode_clusters);

    println!("Estimated geode clusters: {expected_geode_clusters}");
    println!("Estimated budding amethyst clusters: {expected_budding_clusters}");
}
