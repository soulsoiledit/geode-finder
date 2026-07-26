use crate::{Args, version::Version};

pub fn poisson_tail(trials: f64, probability: f64, threshold: u32) -> f64 {
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

    let loaded_diameter = f64::from(args.loaded_radius) * 2.0 + 1.0;
    let loaded_area = loaded_diameter.powi(2);

    let search_diameter = f64::from(args.search_radius) * 2.0 + 1.0;
    let search_area = search_diameter.powi(2);

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
