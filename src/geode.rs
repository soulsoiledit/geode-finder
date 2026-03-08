use crate::{
    math::{Block, Random, inv_sqrt},
    noise::{Noise, NormalNoise},
    version::Version,
};

pub struct Geode<V: Version> {
    seed: i64,
    x_scale: i64,
    z_scale: i64,
    random: V::RANDOM,
    noise: NormalNoise,
}

impl<V: Version> Geode<V> {
    // reversed from next_float
    const CHANCE_INT: i32 = (V::CHANCE * V::RANDOM::FLOAT_MULTIPLIER.recip()) as i32;

    pub fn new(seed: i64) -> Self {
        let mut random = V::RANDOM::new(seed);
        let noise = V::new_normal_noise(seed);
        Self {
            seed,
            x_scale: (random.next_long() | 1).wrapping_mul(16),
            z_scale: (random.next_long() | 1).wrapping_mul(16),
            random,
            noise,
        }
    }

    fn set_feature_seed(&mut self, chunk_x: i64, chunk_z: i64) {
        let scaled_x = chunk_x.wrapping_mul(self.x_scale);
        let scaled_z = chunk_z.wrapping_mul(self.z_scale);
        let decoration_seed = scaled_x.wrapping_add(scaled_z) ^ self.seed;
        let feature_seed = decoration_seed.wrapping_add(V::SALT);
        self.random.set_seed(feature_seed);
    }

    pub fn check(&mut self, chunk_x: i64, chunk_z: i64) -> bool {
        self.set_feature_seed(chunk_x, chunk_z);
        self.random.next_float() < V::CHANCE
    }

    pub fn check_fast(&mut self, chunk_x: i64, chunk_z: i64) -> u32 {
        self.set_feature_seed(chunk_x, chunk_z);
        u32::from(self.random.next_bits(24) < Self::CHANCE_INT)
    }

    pub fn generate(&mut self, chunk_x: i64, chunk_z: i64) -> u32 {
        if !self.check(chunk_x, chunk_z) {
            return 0;
        }

        let origin = {
            let x = self.random.next_int(16) + chunk_x as i32 * 16;
            let z = self.random.next_int(16) + chunk_z as i32 * 16;
            let y = V::Y_RANGE.sample(&mut self.random);
            Block::new(x, y, z)
        };

        let num_points = V::POINTS.sample(&mut self.random);
        let num_points_f = f64::from(num_points);
        let crack_size_adjustment = num_points_f / f64::from(V::RADIUS.max);

        let air_dist = inv_sqrt(V::AIR_LAYER);
        let amethyst_dist = inv_sqrt(V::AMETHYST_LAYER + crack_size_adjustment);
        let basalt_dist = inv_sqrt(V::BASALT_LAYER + crack_size_adjustment);

        let crack_size = inv_sqrt(
            V::CRACK_SIZE
                + self.random.next_double() / 2.0
                + crack_size_adjustment * f64::from(num_points > 3),
        );
        let should_generate_crack = f64::from(self.random.next_float()) < V::CRACK_CHANCE;

        let mut points = Vec::with_capacity(num_points as usize);
        for _ in 0..num_points {
            let mut next_coord = || V::RADIUS.sample(&mut self.random);
            let point = origin.add(next_coord(), next_coord(), next_coord());
            let offset = f64::from(V::POINT_OFFSET.sample(&mut self.random));
            points.push((point, offset));
        }

        let crack_points = should_generate_crack.then(|| {
            let crack = num_points * 2 + 1;
            let cracks = [(crack, 0), (0, crack), (crack, crack), (0, 0)];
            let (dx, dz) = cracks[self.random.next_int(4) as usize];
            [7, 5, 1].map(|dy| origin.add(dx, dy, dz))
        });

        let mut budding_count = 0;

        let range = |o: i32| o - V::OFFSET..=o + V::OFFSET;
        for z in range(origin.z) {
            let zf = f64::from(z);
            for y in range(origin.y) {
                let yf = f64::from(y);
                for x in range(origin.x) {
                    let xf = f64::from(x);
                    let block = Block::new(x, y, z);

                    let noise_offset = self.noise.get_value(xf, yf, zf) * V::NOISE_MULTIPLIER;

                    let mut shell_sum = 0.0;
                    for (point, offset) in &points {
                        shell_sum += V::inv_sqrt(V::distance_sq(&block, point) + (*offset));
                    }
                    shell_sum += noise_offset * num_points_f;

                    let in_solid_geode = shell_sum >= basalt_dist && shell_sum < air_dist;
                    if in_solid_geode {
                        let mut crack_sum = 0.0;
                        if let Some(crack_points) = crack_points {
                            for point in crack_points {
                                crack_sum +=
                                    V::inv_sqrt(V::distance_sq(&block, &point) + V::CRACK_OFFSET);
                            }
                            crack_sum += noise_offset * 3.0;
                        }

                        let in_crack = should_generate_crack && crack_sum >= crack_size;
                        if !in_crack && shell_sum >= amethyst_dist {
                            let place_budding =
                                f64::from(self.random.next_float()) < V::BUDDING_CHANCE;
                            if place_budding {
                                budding_count += 1;
                                self.random.skip(1);
                            }
                        }
                    }
                }
            }
        }

        budding_count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::version;

    const SEED: i64 = 0;
    const RANGE: i64 = 16;
    const FAR: i64 = 1_875_000;

    #[derive(Debug, PartialEq)]
    struct TestGeodeResult {
        geode_count: u32,
        budding_count: u32,
    }

    fn generate<V: Version>(center: i64, geode: &mut Geode<V>) -> TestGeodeResult {
        let start = center - RANGE;
        let end = center + RANGE;

        let mut geode_count = 0;
        let mut fast_geode_count = 0;
        let mut budding_count = 0;
        for cx in start..=end {
            for cz in start..=end {
                geode_count += geode.check(cx, cz) as u32;
                fast_geode_count += geode.check_fast(cx, cz);
                budding_count += geode.generate(cx, cz);
            }
        }

        assert!(geode_count == fast_geode_count);
        TestGeodeResult {
            geode_count,
            budding_count,
        }
    }

    #[test]
    fn generate_17() {
        let mut finder = Geode::<version::MC17>::new(SEED);
        assert_eq!(
            generate(0, &mut finder),
            TestGeodeResult {
                geode_count: 22,
                budding_count: 804
            }
        );
        assert_eq!(
            generate(FAR, &mut finder),
            TestGeodeResult {
                geode_count: 15,
                budding_count: 560
            }
        );
    }

    #[test]
    fn generate_18() {
        let mut finder = Geode::<version::MC18>::new(SEED);
        assert_eq!(
            generate(0, &mut finder),
            TestGeodeResult {
                geode_count: 37,
                budding_count: 1167
            }
        );
        assert_eq!(
            generate(FAR, &mut finder),
            TestGeodeResult {
                geode_count: 39,
                budding_count: 1518,
            }
        );
    }

    #[test]
    fn generate_19() {
        let mut finder = Geode::<version::MC19>::new(SEED);
        assert_eq!(
            generate(0, &mut finder),
            TestGeodeResult {
                geode_count: 37,
                budding_count: 1166
            }
        );
        assert_eq!(
            generate(FAR, &mut finder),
            TestGeodeResult {
                geode_count: 39,
                budding_count: 1521
            }
        );
    }
}
