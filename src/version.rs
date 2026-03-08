use crate::{
    math::{Block, JavaRandom, Random, UniformInt, Xoroshiro128PlusPlusRandom, inv_sqrt},
    noise::{ImprovedNoise, NormalNoise, PerlinNoise},
};

pub trait Version {
    type RANDOM: Random;

    const AIR_LAYER: f64 = 1.7;
    const AMETHYST_LAYER: f64 = 2.2;
    const BASALT_LAYER: f64 = 4.2;

    const CRACK_CHANCE: f64 = 0.95;
    const CRACK_SIZE: f64 = 2.0;
    const CRACK_OFFSET: f64 = 2.0;

    const BUDDING_CHANCE: f64 = 0.083;
    const RADIUS: UniformInt = UniformInt::new(4, 6);
    const POINTS: UniformInt = UniformInt::new(3, 4);
    const POINT_OFFSET: UniformInt = UniformInt::new(1, 2);

    const OFFSET: i32 = 16;
    const NOISE_MULTIPLIER: f64 = 0.05;

    const SALT: i64 = 20002;
    const CHANCE: f32 = 1.0 / 24.0;
    const Y_RANGE: UniformInt = UniformInt::new(-58, 30);

    fn new_normal_noise(seed: i64) -> NormalNoise {
        let mut noise_random = JavaRandom::new(seed);
        let first = Self::new_noise(&mut noise_random);
        let second = Self::new_noise(&mut noise_random);
        NormalNoise::new(PerlinNoise::new(first), PerlinNoise::new(second))
    }

    fn new_noise(random: &mut JavaRandom) -> ImprovedNoise {
        ImprovedNoise::new(&mut random.fork_from_hash())
    }

    fn distance_sq(pos1: &Block, pos2: &Block) -> f64 {
        let dx = f64::from(pos1.x - pos2.x);
        let dy = f64::from(pos1.y - pos2.y);
        let dz = f64::from(pos1.z - pos2.z);
        dx * dx + dy * dy + dz * dz
    }

    // fast inv sqrt
    fn inv_sqrt(x: f64) -> f64 {
        let i = f64::from_bits(0x5FE6_EB50_C7B5_37AA - (x.to_bits() >> 1));
        i * (1.5 - (0.5 * x) * i * i)
    }
}

pub struct MC17;
impl Version for MC17 {
    type RANDOM = JavaRandom;

    const SALT: i64 = 20000;
    const CHANCE: f32 = 1.0 / 53.0;
    const Y_RANGE: UniformInt = UniformInt::new(6, 46);

    fn new_noise(random: &mut JavaRandom) -> ImprovedNoise {
        random.skip(262 * 4);
        ImprovedNoise::new(random)
    }

    // 1.17 measures from centers instead of corners
    fn distance_sq(pos1: &Block, pos2: &Block) -> f64 {
        let dx = f64::from(pos1.x - pos2.x) + 0.5;
        let dy = f64::from(pos1.y - pos2.y) + 0.5;
        let dz = f64::from(pos1.z - pos2.z) + 0.5;
        dx * dx + dy * dy + dz * dz
    }
}

pub struct MC18;
impl Version for MC18 {
    type RANDOM = Xoroshiro128PlusPlusRandom;
}

pub struct MC19;
impl Version for MC19 {
    type RANDOM = Xoroshiro128PlusPlusRandom;

    fn inv_sqrt(val: f64) -> f64 {
        inv_sqrt(val)
    }
}
