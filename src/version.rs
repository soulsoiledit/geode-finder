use std::ops::RangeInclusive;

use crate::{
    math::{Block, JavaRandom, Random, Xoroshiro128PlusPlusRandom},
    noise::{Float, NormalNoise, PerlinNoise},
};

pub trait Version: Sized {
    type Base: Version;
    type Random: Random;
    type Float: Float;

    const AIR_LAYER: f64 = 1.7;
    const AMETHYST_LAYER: f64 = 2.2;
    const BASALT_LAYER: f64 = 4.2;

    const CRACK_CHANCE: f64 = 0.95;
    const CRACK_SIZE: f64 = 2.0;
    const CRACK_OFFSET: f64 = 2.0;

    const BUDDING_CHANCE: f64 = 0.083;
    const RADIUS: RangeInclusive<i32> = 4..=6;
    const POINTS: RangeInclusive<i32> = 3..=4;
    const POINT_OFFSET: RangeInclusive<i32> = 1..=2;

    const OFFSET: i32 = 16;
    const NOISE_MULTIPLIER: f64 = 0.05;

    const OCTAVE: i32 = -4;
    const AMPLITUDE: f64 = 1.0;

    const SALT: i64 = Self::Base::SALT;
    const CHANCE: f32 = Self::Base::CHANCE;
    const Y_RANGE: RangeInclusive<i32> = Self::Base::Y_RANGE;

    fn inv_sqrt(x: f64) -> f64 {
        Self::Base::inv_sqrt(x)
    }

    fn new_perlin_pair<F: Float>(
        random: &mut JavaRandom,
        octave: i32,
    ) -> (PerlinNoise<F>, PerlinNoise<F>) {
        Self::Base::new_perlin_pair(random, octave)
    }

    fn distance_sq(pos1: &Block, pos2: &Block) -> f64 {
        Self::Base::distance_sq(pos1, pos2)
    }

    fn new_normal_noise<V: Version>(
        random: &mut JavaRandom,
        octave: i32,
        amplitude: f64,
    ) -> NormalNoise<V> {
        Self::Base::new_normal_noise(random, octave, amplitude)
    }

    fn apply_amplitude<F: Float>(amplitude: F, layer1: F, layer2: F) -> F {
        Self::Base::apply_amplitude(amplitude, layer1, layer2)
    }
}

pub struct MC17;
impl Version for MC17 {
    type Base = Self;
    type Random = JavaRandom;
    type Float = f64;

    const SALT: i64 = 20000;
    const CHANCE: f32 = 1.0 / 53.0;
    const Y_RANGE: RangeInclusive<i32> = 6..=46;

    fn new_perlin_pair<F: Float>(
        random: &mut JavaRandom,
        octave: i32,
    ) -> (PerlinNoise<F>, PerlinNoise<F>) {
        let steps = 262 * octave.unsigned_abs() as usize;
        random.skip(steps);
        let first = PerlinNoise::new(random);
        random.skip(steps);
        let second = PerlinNoise::new(random);
        (first, second)
    }

    // measures from block centers
    fn distance_sq(pos1: &Block, pos2: &Block) -> f64 {
        let dx = f64::from(pos1.x - pos2.x) + 0.5;
        let dy = f64::from(pos1.y - pos2.y) + 0.5;
        let dz = f64::from(pos1.z - pos2.z) + 0.5;
        dx * dx + dy * dy + dz * dz
    }

    fn inv_sqrt(x: f64) -> f64 {
        const QUAKE_MAGIC_64: u64 = 0x5FE6_EB50_C7B5_37AA;
        let i = f64::from_bits(QUAKE_MAGIC_64 - (x.to_bits() >> 1));
        i * (1.5 - (0.5 * x) * i * i)
    }

    fn new_normal_noise<V: Version>(
        random: &mut JavaRandom,
        octave: i32,
        amplitude: f64,
    ) -> NormalNoise<V> {
        NormalNoise::new(random, octave, amplitude).unnormalized(amplitude)
    }

    fn apply_amplitude<F: Float>(amplitude: F, layer1: F, layer2: F) -> F {
        amplitude * (layer1 + layer2)
    }
}

pub struct MC18;
impl Version for MC18 {
    type Base = MC17;
    type Random = Xoroshiro128PlusPlusRandom;
    type Float = <Self::Base as Version>::Float;

    const SALT: i64 = 20002;
    const CHANCE: f32 = 1.0 / 24.0;
    const Y_RANGE: RangeInclusive<i32> = -58..=30;

    fn new_perlin_pair<F: Float>(
        random: &mut JavaRandom,
        _octave: i32,
    ) -> (PerlinNoise<F>, PerlinNoise<F>) {
        let first = PerlinNoise::new(&mut random.fork_from_hash());
        let second = PerlinNoise::new(&mut random.fork_from_hash());
        (first, second)
    }
}

pub struct MC182;
impl Version for MC182 {
    type Base = MC18;
    type Random = <Self::Base as Version>::Random;
    type Float = <Self::Base as Version>::Float;

    fn distance_sq(pos1: &Block, pos2: &Block) -> f64 {
        let dx = f64::from(pos1.x - pos2.x);
        let dy = f64::from(pos1.y - pos2.y);
        let dz = f64::from(pos1.z - pos2.z);
        dx * dx + dy * dy + dz * dz
    }
}

pub struct MC194;
impl Version for MC194 {
    type Base = MC182;
    type Random = <Self::Base as Version>::Random;
    type Float = <Self::Base as Version>::Float;

    fn inv_sqrt(x: f64) -> f64 {
        x.sqrt().recip()
    }
}

pub struct MC263;
impl Version for MC263 {
    type Base = MC194;
    type Random = <Self::Base as Version>::Random;
    type Float = f32;

    fn new_normal_noise<V: Version>(
        random: &mut JavaRandom,
        octave: i32,
        amplitude: f64,
    ) -> NormalNoise<V> {
        NormalNoise::new(random, octave, amplitude)
    }

    fn apply_amplitude<F: Float>(amplitude: F, layer1: F, layer2: F) -> F {
        // distribution doesn't maintain bit parity
        amplitude * layer1 + amplitude * layer2
    }
}
