use crate::noise::DoublePerlinNoiseSampler;
use crate::random::JavaRandom;
use crate::random::Random;
use crate::random::Xoroshiro128PlusPlusRandom;
use crate::GameVersion;

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // TODO: decouple this from random
    // fn fast_inv_sqrt() {
    // let mut random = Xoroshiro128PlusPlusRandom::with_seed(1);
    // assert_eq!(fast_inverse_sqrt(random.next_double()), 1.0274868578944794);
    // assert_eq!(fast_inverse_sqrt(random.next_double()), 1.0526759381538786);
    // assert_eq!(fast_inverse_sqrt(random.next_double()), 1.5075261791112862);
    // }

    // TODO: get real values for these too
    // #[test]
    // fn generate_17() {
    //     let mut finder = Geode::new(1, GameVersion::MC17);
    //     // finder.generate(10, 10);
    // }
    //
    // #[test]
    // fn generate_18() {
    //     let mut finder = Geode::new(1, GameVersion::MC18);
    //     // finder.generate(10, 10);
    // }

    #[test]
    fn generate_20() {
        let mut finder = Geode::new(1, GameVersion::MC20);

        let range = 128;
        let mut geode_count = 0;
        for i in -range..=range {
            for j in -range..=range {
                if finder.check_chunk(i, j) {
                    geode_count += 1;
                }
            }
        }

        assert_eq!(geode_count, 2759);

        let range = 32;
        let mut amethyst_count = 0;
        for i in -range..=range {
            for j in -range..=range {
                if finder.check_chunk(i, j) {
                    amethyst_count += finder.generate(i, j);
                }
            }
        }

        assert_eq!(amethyst_count, 6704);
    }

    #[test]
    fn generate_20_below_y0() {
        let mut finder = Geode::new(1, GameVersion::MC20);

        let range = 128;
        let mut geode_count = 0;
        for i in -range..=range {
            for j in -range..=range {
                if finder.check_chunk_below_y0(i, j) {
                    geode_count += 1;
                }
            }
        }

        assert_eq!(geode_count, 1803);

        let range = 32;
        let mut amethyst_count = 0;
        for i in -range..=range {
            for j in -range..=range {
                if finder.check_chunk(i, j) {
                    amethyst_count += finder.generate_below_y0(i, j);
                }
            }
        }

        assert_eq!(amethyst_count, 4047);
    }
}

const FILLER_THICKNESS: f64 = 1.7;
const INNER_THICKNESS: f64 = 2.2;
const OUTER_THICKNESS: f64 = 4.2;

const CRACK_CHANCE: f64 = 0.95;
const CRACK_SIZE: f64 = 2.0;
const CRACK_OFFSET: i32 = 2;

const BUDDING_CHANCE: f64 = 0.083;
const OUTER_WALL_DIST: (i32, i32) = (4, 6);

const DISTRIBUTION_POINTS: (i32, i32) = (3, 4);
const POINT_OFFSET: (i32, i32) = (1, 2);

const OFFSET: i32 = 16;
const NOISE_MULTIPLIER: f64 = 0.05;

struct BlockPos {
    x: i32,
    y: i32,
    z: i32,
}

fn find_squared_distance_17(pos1: &BlockPos, pos2: &BlockPos) -> f64 {
    let d = (pos1.x - pos2.x) as f64 + 0.5;
    let e = (pos1.y - pos2.y) as f64 + 0.5;
    let f = (pos1.z - pos2.z) as f64 + 0.5;

    d * d + e * e + f * f
}

fn find_squared_distance_18(pos1: &BlockPos, pos2: &BlockPos) -> f64 {
    let d = (pos1.x - pos2.x) as f64;
    let e = (pos1.y - pos2.y) as f64;
    let f = (pos1.z - pos2.z) as f64;

    d * d + e * e + f * f
}

impl BlockPos {
    fn add(&self, dx: i32, dy: i32, dz: i32) -> BlockPos {
        BlockPos {
            x: self.x + dx,
            y: self.y + dy,
            z: self.z + dz,
        }
    }
}

fn fast_inverse_sqrt(x: f64) -> f64 {
    let y = f64::from_bits(0x5FE6EB50C7B537AA - (&x.to_bits() >> 1));
    y * (1.5 - (0.5 * &x) * y * y)
}

fn real_inverse_sqrt(x: f64) -> f64 {
    1.0 / x.sqrt()
}

pub struct Geode {
    random: Box<dyn Random>,
    noise: DoublePerlinNoiseSampler,
    chance: f32,
    seed: i64,
    salt: i64,
    seed_high: i64,
    seed_low: i64,
    y_min: i32,
    y_max: i32,
    inverse_sqrt: fn(val: f64) -> f64,
    find_squared_distance: fn(pos1: &BlockPos, pos2: &BlockPos) -> f64,
}

impl Geode {
    pub fn new(seed: i64, game_version: GameVersion) -> Self {
        let mut random: Box<dyn Random> = match game_version {
            GameVersion::MC17 => Box::new(JavaRandom::with_seed(seed)),
            _ => Box::new(Xoroshiro128PlusPlusRandom::with_seed(seed)),
        };

        let chance = match game_version {
            GameVersion::MC17 => 1.0 / 53.0,
            _ => 1.0 / 24.0,
        };

        let salt = match game_version {
            GameVersion::MC17 => 20000,
            _ => 20002,
        };

        let (y_min, y_max) = match game_version {
            GameVersion::MC17 => (6, 46),
            _ => (-58, 30),
        };

        let find_squared_distance = match game_version {
            GameVersion::MC17 => find_squared_distance_17,
            _ => find_squared_distance_18,
        };

        let inverse_sqrt = match game_version {
            GameVersion::MC20 => real_inverse_sqrt,
            _ => fast_inverse_sqrt,
        };

        let mut noise_random = JavaRandom::with_seed(seed);

        Geode {
            seed_high: random.next_long() | 1,
            seed_low: random.next_long() | 1,
            chance,
            seed,
            salt,
            y_min,
            y_max,
            inverse_sqrt,
            find_squared_distance,
            random,
            noise: DoublePerlinNoiseSampler::new(&mut noise_random, game_version),
        }
    }

    pub fn set_decorator_seed(&mut self, chunk_x: i64, chunk_z: i64) {
        let bx = chunk_x.wrapping_shl(4).wrapping_mul(self.seed_high);
        let bz = chunk_z.wrapping_shl(4).wrapping_mul(self.seed_low);
        let pop_seed = bx.wrapping_add(bz) ^ self.seed;
        let dec_seed = pop_seed.wrapping_add(self.salt);
        self.random.set_seed(dec_seed);
    }

    pub fn check_chunk(&mut self, chunk_x: i64, chunk_z: i64) -> bool {
        self.set_decorator_seed(chunk_x, chunk_z);
        self.random.next_float() < self.chance
    }

    pub fn generate(&mut self, chunk_x: i64, chunk_z: i64) -> i32 {
        if !self.check_chunk(chunk_x, chunk_z) {
            return 0;
        }

        let origin = BlockPos {
            x: self.random.next_int(16) + 16 * chunk_x as i32,
            z: self.random.next_int(16) + 16 * chunk_z as i32,
            y: self.random.next_between(self.y_min, self.y_max),
        };

        let distribution_points = self
            .random
            .next_between(DISTRIBUTION_POINTS.0, DISTRIBUTION_POINTS.1);
        let d: f64 = (distribution_points as f64) / OUTER_WALL_DIST.1 as f64;

        let inv_filling_thickness = real_inverse_sqrt(FILLER_THICKNESS);
        let inv_inner_thickness = real_inverse_sqrt(INNER_THICKNESS + d);
        let inv_outer_thickness = real_inverse_sqrt(OUTER_THICKNESS + d);

        let l = real_inverse_sqrt(
            CRACK_SIZE
                + self.random.next_double() / 2.0
                + if distribution_points > 3 { d } else { 0.0 },
        );

        let generate_crack = (self.random.next_float() as f64) < CRACK_CHANCE;

        let mut block_list1: Vec<(BlockPos, i32)> = vec![];
        for _ in 0..distribution_points {
            let dx = self
                .random
                .next_between(OUTER_WALL_DIST.0, OUTER_WALL_DIST.1);
            let dy = self
                .random
                .next_between(OUTER_WALL_DIST.0, OUTER_WALL_DIST.1);
            let dz = self
                .random
                .next_between(OUTER_WALL_DIST.0, OUTER_WALL_DIST.1);

            let point_offset = self.random.next_between(POINT_OFFSET.0, POINT_OFFSET.1);
            block_list1.push((origin.add(dx, dy, dz), point_offset));
        }

        let mut block_list2: Vec<BlockPos> = vec![];
        if generate_crack {
            let n = self.random.next_int(4);
            let o = distribution_points * 2 + 1;

            let dx = if n == 0 || n == 2 { o } else { 0 };
            let dz = if n == 1 || n == 2 { o } else { 0 };

            block_list2 = vec![
                origin.add(dx, 7, dz),
                origin.add(dx, 5, dz),
                origin.add(dx, 1, dz),
            ];
        }

        let mut budding_count = 0;

        for z in (origin.z - OFFSET)..=(origin.z + OFFSET) {
            for y in (origin.y - OFFSET)..=(origin.y + OFFSET) {
                for x in (origin.x - OFFSET)..=(origin.x + OFFSET) {
                    let block3 = BlockPos { x, y, z };

                    let r = self
                        .noise
                        .sample(block3.x as f64, block3.y as f64, block3.z as f64)
                        * NOISE_MULTIPLIER;
                    let mut s: f64 = 0.0;
                    let mut t: f64 = 0.0;

                    for (coord, val) in &block_list1 {
                        s += (self.inverse_sqrt)(
                            (self.find_squared_distance)(&block3, coord) + (*val as f64),
                        ) + r;
                    }

                    for block4 in &block_list2 {
                        t += (self.inverse_sqrt)(
                            (self.find_squared_distance)(&block3, block4) + CRACK_OFFSET as f64,
                        ) + r;
                    }

                    if s < inv_outer_thickness {
                        continue;
                    }

                    if generate_crack && t >= l && s < inv_filling_thickness {
                        continue;
                    }

                    if s >= inv_filling_thickness {
                        continue;
                    }

                    if s >= inv_inner_thickness {
                        let place_budding = (self.random.next_float() as f64) < BUDDING_CHANCE;

                        if place_budding {
                            budding_count += 1;
                            self.random.next_float();
                        }
                    }
                }
            }
        }

        budding_count
    }
}

// How should geode struct work
// - default method for generate() that uses fields
// - needs separate for block position squared distance
// - needs separate for inverse square root
// - salt is different
// - random is different
// - noise is different
// need different outputs depending on situation
// - face count
// - face count + budding list for processing
