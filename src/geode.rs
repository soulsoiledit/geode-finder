use crate::noise::DoublePerlinNoiseSampler;
use crate::random::JavaRandom;
use crate::random::Random;
use crate::random::Xoroshiro128PlusPlus;

const GEODE_FILLER_THICKNESS: f64 = 1.7;
const GEODE_INNER_THICKNESS: f64 = 2.2;
const GEODE_OUTER_THICKNESS: f64 = 4.2;

const GEODE_CRACK_CHANCE: f64 = 0.95;
const GEODE_CRACK_SIZE: f64 = 2.0;
const GEODE_CRACK_OFFSET: i32 = 2;

const GEODE_BUDDING_CHANCE: f64 = 0.083;
const GEODE_OUTER_WALL_DIST: (i32, i32) = (4, 6);

const GEODE_DISTRIBUTION_POINTS: (i32, i32) = (3, 4);
const GEODE_POINT_OFFSET: (i32, i32) = (1, 2);

const GEODE_MIN_OFFSET: i32 = -16;
const GEODE_MAX_OFFSET: i32 = 16;

const GEODE_NOISE_MULTIPLIER: f64 = 0.05;

struct BlockPos {
    x: i32,
    y: i32,
    z: i32,
}

pub struct GeodeGenerator {
    random: Box<dyn Random>,
    world_seed: i64,
    a: i64,
    b: i64,
    is_17: bool,
    chance: f32,
    y_min: i32,
    y_max: i32,
    noise_source: DoublePerlinNoiseSampler,
}

impl BlockPos {
    fn get_squared_dist(&self, pos: &BlockPos, is_17: bool) -> f64 {
        let d: f64;
        let e: f64;
        let f: f64;

        if is_17 {
            d = (self.x - pos.x) as f64 + 0.5;
            e = (self.y - pos.y) as f64 + 0.5;
            f = (self.z - pos.z) as f64 + 0.5;
        } else {
            d = (self.x - pos.x) as f64;
            e = (self.y - pos.y) as f64;
            f = (self.z - pos.z) as f64;
        }

        d * d + e * e + f * f
    }

    fn add(&self, dx: i32, dy: i32, dz: i32) -> BlockPos {
        BlockPos {
            x: self.x + dx,
            y: self.y + dy,
            z: self.z + dz,
        }
    }
}

impl GeodeGenerator {
    pub fn new(world_seed: i64, is_17: bool) -> Self {
        let mut random: Box<dyn Random>;
        let chance: f32;
        let y_min: i32;
        let y_max: i32;

        if is_17 {
            random = Box::new(JavaRandom::with_seed(world_seed));
            chance = 1.0 / 53.0;
            y_min = 6;
            y_max = 46;
        } else {
            random = Box::new(Xoroshiro128PlusPlus::with_seed(world_seed));
            chance = 1.0 / 24.0;
            y_min = -58;
            y_max = 30;
        }

        GeodeGenerator {
            a: random.next_long() | 1,
            b: random.next_long() | 1,
            random,
            world_seed,
            is_17,
            chance,
            y_min,
            y_max,
            noise_source: DoublePerlinNoiseSampler::new(JavaRandom::with_seed(world_seed), is_17),
        }
    }

    pub fn set_decorator_seed(&mut self, chunk_x: i64, chunk_z: i64, salt: i64) -> () {
        self.random.set_seed(
            (chunk_x
                .wrapping_shl(4)
                .wrapping_mul(self.a)
                .wrapping_add(chunk_z.wrapping_shl(4).wrapping_mul(self.b))
                ^ self.world_seed)
                .wrapping_add(salt),
        );
    }

    pub fn check_chunk(&mut self) -> bool {
        self.random.next_float() < self.chance
    }

    pub fn fast_inv_sqrt(&mut self, x: f64) -> f64 {
        let y = f64::from_bits(6910469410427058090 - (x.to_bits() >> 1));
        y * (1.5 - 0.5 * x * y * y)
    }

    pub fn generate(&mut self, chunk_x: i64, chunk_z: i64) -> i32 {
        let origin = BlockPos {
            x: self.random.next_int(16) + 16 * chunk_x as i32,
            z: self.random.next_int(16) + 16 * chunk_z as i32,
            y: self.random.next_between(self.y_min, self.y_max),
        };

        let distribution_points = self
            .random
            .next_between(GEODE_DISTRIBUTION_POINTS.0, GEODE_DISTRIBUTION_POINTS.1);
        let d: f64 = (distribution_points as f64) / GEODE_OUTER_WALL_DIST.1 as f64;

        let inv_filling_thickness = 1.0 / GEODE_FILLER_THICKNESS.sqrt();
        let inv_inner_thickness = 1.0 / (GEODE_INNER_THICKNESS + d).sqrt();
        let inv_outer_thickness = 1.0 / (GEODE_OUTER_THICKNESS + d).sqrt();

        let l = 1.0
            / (GEODE_CRACK_SIZE
                + self.random.next_double() / 2.0
                + if distribution_points > 3 { d } else { 0.0 })
            .sqrt();

        let generate_crack = (self.random.next_float() as f64) < GEODE_CRACK_CHANCE;

        let mut block_list1: Vec<(BlockPos, i32)> = vec![];
        for _ in 0..distribution_points {
            let dx = self
                .random
                .next_between(GEODE_OUTER_WALL_DIST.0, GEODE_OUTER_WALL_DIST.1);
            let dy = self
                .random
                .next_between(GEODE_OUTER_WALL_DIST.0, GEODE_OUTER_WALL_DIST.1);
            let dz = self
                .random
                .next_between(GEODE_OUTER_WALL_DIST.0, GEODE_OUTER_WALL_DIST.1);

            let point_offset = self
                .random
                .next_between(GEODE_POINT_OFFSET.0, GEODE_POINT_OFFSET.1);
            block_list1.push((origin.add(dx, dy, dz), point_offset));
        }

        let mut block_list2: Vec<BlockPos> = vec![];
        if generate_crack {
            let n = self.random.next_int(4);
            let o = distribution_points * 2 + 1;

            match n {
                0 => {
                    block_list2 = vec![
                        origin.add(o, 7, 0),
                        origin.add(o, 5, 0),
                        origin.add(o, 1, 0)
                    ]
                },

                1 => {
                    block_list2 = vec![
                        origin.add(0, 7, o),
                        origin.add(0, 5, o),
                        origin.add(0, 1, o)
                    ]
                }

                2 => {
                    block_list2 = vec![
                        origin.add(o, 7, o),
                        origin.add(o, 5, o),
                        origin.add(o, 1, o)
                    ]
                }

                _ => {
                    block_list2 = vec![
                        origin.add(0, 7, 0),
                        origin.add(0, 5, 0),
                        origin.add(0, 1, 0)
                    ]
                }
            }
        }

        let mut budding_count = 0;
        for z in GEODE_MIN_OFFSET..=GEODE_MAX_OFFSET {
            for y in GEODE_MIN_OFFSET..=GEODE_MAX_OFFSET {
                for x in GEODE_MIN_OFFSET..=GEODE_MAX_OFFSET {
                    let block3 = origin.add(x, y, z);
                    let r =
                        self.noise_source
                            .sample(block3.x as f64, block3.y as f64, block3.z as f64)
                            * GEODE_NOISE_MULTIPLIER;
                    let mut s = 0f64;
                    let mut t = 0f64;

                    for (coord, val) in &block_list1 {
                        s += self.fast_inv_sqrt(
                            block3.get_squared_dist(coord, self.is_17) + (val.clone() as f64),
                        ) + r;
                    }

                    for block4 in &block_list2 {
                        t += self.fast_inv_sqrt(
                            block3.get_squared_dist(block4, self.is_17) + GEODE_CRACK_OFFSET as f64,
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
                        let place_budding =
                            (self.random.next_float() as f64) < GEODE_BUDDING_CHANCE;

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
