use std::{array, f64::consts::SQRT_2, fmt, marker::PhantomData};

use num_traits::Float as NumTraitsFloat;

use crate::{
    math::{JavaRandom, Random},
    version::Version,
};

pub trait Float: NumTraitsFloat + From<f32> + Into<f64> + fmt::Debug {
    const ZERO: Self;
    const ONE: Self;
    const NEG_ONE: Self;
    fn from_f64(val: f64) -> Self;
}

impl Float for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const NEG_ONE: Self = -1.0;
    fn from_f64(val: f64) -> Self {
        val
    }
}

impl Float for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const NEG_ONE: Self = -1.0;
    fn from_f64(val: f64) -> Self {
        val as f32
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PerlinNoise<F: Float> {
    x_offset: f64,
    y_offset: f64,
    z_offset: f64,
    permutation: [u8; 256],
    _output: PhantomData<F>,
}

impl<F: Float> PerlinNoise<F> {
    const NOISE_OFFSET_SCALE: f64 = 256.0;
    const ROUND_OFF: f64 = (1 << 25) as f64;
    const HALF_ROUND_OFF: f64 = ((1 << 24) as f64).next_down();
    const GRADIENT: [[F; 3]; 16] = [
        [F::ONE, F::ONE, F::ZERO],
        [F::NEG_ONE, F::ONE, F::ZERO],
        [F::ONE, F::NEG_ONE, F::ZERO],
        [F::NEG_ONE, F::NEG_ONE, F::ZERO],
        [F::ONE, F::ZERO, F::ONE],
        [F::NEG_ONE, F::ZERO, F::ONE],
        [F::ONE, F::ZERO, F::NEG_ONE],
        [F::NEG_ONE, F::ZERO, F::NEG_ONE],
        [F::ZERO, F::ONE, F::ONE],
        [F::ZERO, F::NEG_ONE, F::ONE],
        [F::ZERO, F::ONE, F::NEG_ONE],
        [F::ZERO, F::NEG_ONE, F::NEG_ONE],
        [F::ONE, F::ONE, F::ZERO],
        [F::ZERO, F::NEG_ONE, F::ONE],
        [F::NEG_ONE, F::ONE, F::ZERO],
        [F::ZERO, F::NEG_ONE, F::NEG_ONE],
    ];

    pub fn new(random: &mut JavaRandom) -> Self {
        let x_offset = random.next_double() * Self::NOISE_OFFSET_SCALE;
        let y_offset = random.next_double() * Self::NOISE_OFFSET_SCALE;
        let z_offset = random.next_double() * Self::NOISE_OFFSET_SCALE;

        let mut permutation: [u8; 256] = array::from_fn(|i| i as u8);
        for i in 0..256 {
            let offset = random.next_int(256 - i as i32) as usize;
            permutation.swap(i, i + offset);
        }

        Self {
            x_offset,
            y_offset,
            z_offset,
            permutation,
            _output: PhantomData,
        }
    }

    fn permute(&self, i: i32) -> i32 {
        i32::from(self.permutation[(i & 0xFF) as usize])
    }

    fn grad_dot(&self, hash: i32, x: F, y: F, z: F) -> F {
        let [gx, gy, gz] = Self::GRADIENT[(hash & 0xF) as usize];
        gx * x + gy * y + gz * z
    }

    fn wrap(&self, value: f64) -> f64 {
        if (-Self::HALF_ROUND_OFF..Self::HALF_ROUND_OFF).contains(&value) {
            value
        } else {
            value - (value / Self::ROUND_OFF + 0.5).floor() * Self::ROUND_OFF
        }
    }

    fn smoothstep(&self, x: F) -> F {
        x * x * x * (x * (x * 6.0.into() - 15.0.into()) + 10.0.into())
    }

    fn lerp(&self, tx: F, v0: F, v1: F) -> F {
        v0 + tx * (v1 - v0)
    }

    fn bilerp(&self, tx: F, ty: F, corners: [F; 4]) -> F {
        let [v00, v10, v01, v11] = corners;
        self.lerp(ty, self.lerp(tx, v00, v10), self.lerp(tx, v01, v11))
    }

    fn trilerp(&self, tx: F, ty: F, tz: F, corners: [F; 8]) -> F {
        let [v000, v100, v010, v110, v001, v101, v011, v111] = corners;
        self.lerp(
            tz,
            self.bilerp(tx, ty, [v000, v100, v010, v110]),
            self.bilerp(tx, ty, [v001, v101, v011, v111]),
        )
    }

    fn sample_lerp(&self, whole_pos: [i32; 3], frac_pos: [F; 3], raw_dy: F) -> F {
        let [x, y, z] = whole_pos;
        let [dx, dy, dz] = frac_pos;

        let v0 = self.permute(x);
        let v1 = self.permute(x + 1);

        let v00 = self.permute(v0 + y);
        let v10 = self.permute(v1 + y);
        let v01 = self.permute(v0 + y + 1);
        let v11 = self.permute(v1 + y + 1);

        let corners = [
            self.grad_dot(self.permute(v00 + z), dx, dy, dz),
            self.grad_dot(self.permute(v10 + z), dx - 1.0.into(), dy, dz),
            self.grad_dot(self.permute(v01 + z), dx, dy - 1.0.into(), dz),
            self.grad_dot(self.permute(v11 + z), dx - 1.0.into(), dy - 1.0.into(), dz),
            self.grad_dot(self.permute(v00 + z + 1), dx, dy, dz - 1.0.into()),
            self.grad_dot(
                self.permute(v10 + z + 1),
                dx - 1.0.into(),
                dy,
                dz - 1.0.into(),
            ),
            self.grad_dot(
                self.permute(v01 + z + 1),
                dx,
                dy - 1.0.into(),
                dz - 1.0.into(),
            ),
            self.grad_dot(
                self.permute(v11 + z + 1),
                dx - 1.0.into(),
                dy - 1.0.into(),
                dz - 1.0.into(),
            ),
        ];

        self.trilerp(
            self.smoothstep(dx),
            self.smoothstep(raw_dy),
            self.smoothstep(dz),
            corners,
        )
    }

    fn get(&self, x: f64, y: f64, z: f64) -> F {
        let offset_x = self.wrap(x) + self.x_offset;
        let offset_y = self.wrap(y) + self.y_offset;
        let offset_z = self.wrap(z) + self.z_offset;

        let floor_x = offset_x.floor();
        let floor_y = offset_y.floor();
        let floor_z = offset_z.floor();

        let dx = F::from_f64(offset_x - floor_x);
        let dy = F::from_f64(offset_y - floor_y);
        let dz = F::from_f64(offset_z - floor_z);

        self.sample_lerp(
            [floor_x as i32, floor_y as i32, floor_z as i32],
            [dx, dy, dz],
            dy,
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct NoiseLayer<F: Float>(PerlinNoise<F>, f64);

#[derive(Debug, Clone, Copy)]
pub struct NormalNoise<V: Version> {
    first: NoiseLayer<V::Float>,
    second: NoiseLayer<V::Float>,
    value_factor: V::Float,
    _version: PhantomData<V>,
}

impl<V: Version> NormalNoise<V> {
    const PERSISTENCE: f64 = 0.5;
    const LACUNARITY: f64 = 2.0;

    const PERLIN_DEVIATION: f64 = 0.2702247831245211;
    const TARGET_DEVIATION: f64 = 1.0 / 3.0;
    const PARITY_EXPECTED_DEVIATION: f64 = 0.2;

    const INPUT_FACTOR: f64 = 1.0181268882175227; // 337/331

    pub fn new(random: &mut JavaRandom, octave: i32, amplitude: f64) -> Self {
        let (perlin1, perlin2) = V::new_perlin_pair::<V::Float>(random, octave);

        let frequency = Self::get_frequency(octave);
        let base_amplitude = Self::parity_base_amplitude(amplitude);
        let normalization = Self::normalization(base_amplitude);
        let value_factor = V::Float::from_f64(base_amplitude * normalization);

        Self {
            first: NoiseLayer(perlin1, frequency),
            second: NoiseLayer(perlin2, frequency * Self::INPUT_FACTOR),
            value_factor,
            _version: PhantomData,
        }
    }

    pub fn unnormalized(mut self, amplitude: f64) -> Self {
        self.value_factor = V::Float::from_f64(Self::base_amplitude(amplitude));
        self
    }

    fn get_frequency(octave: i32) -> f64 {
        Self::LACUNARITY.powi(octave)
    }

    const fn base_amplitude(amplitude: f64) -> f64 {
        amplitude * Self::PERSISTENCE * Self::TARGET_DEVIATION / Self::PARITY_EXPECTED_DEVIATION
    }

    const fn normalization(target_amplitude: f64) -> f64 {
        let input_deviation: f64 = Self::PERLIN_DEVIATION * target_amplitude.abs() * SQRT_2;
        let target_deviation: f64 = target_amplitude * Self::TARGET_DEVIATION;
        target_deviation / input_deviation
    }

    const fn parity_base_amplitude(amplitude: f64) -> f64 {
        let target_amplitude = amplitude.abs();
        let new_normalization = Self::normalization(target_amplitude);
        let old_normalization = Self::base_amplitude(target_amplitude);
        old_normalization / new_normalization
    }

    pub fn get(&self, x: f64, y: f64, z: f64) -> V::Float {
        let freq1 = self.first.1;
        let layer1 = self.first.0.get(x * freq1, y * freq1, z * freq1);

        let freq2 = self.second.1;
        let layer2 = self.second.0.get(x * freq2, y * freq2, z * freq2);

        V::apply_amplitude(self.value_factor, layer1, layer2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::version::{MC17, MC18, MC263, Version};

    const SEED: i64 = 0xDEAD_BEEF;

    #[test]
    fn perlin_f64() {
        let perlin = PerlinNoise::<f64>::new(&mut JavaRandom::new(SEED));
        assert_eq!(
            perlin.get(0.0, 0.0, 0.0),
            0.11848800878838489,
            "perlin at 0.0 failed"
        );
        assert_eq!(
            perlin.get(-0.5, -0.5, -0.5),
            -0.08726343323690666,
            "perlin at -0.5 failed"
        );
        assert_eq!(
            perlin.get(-29_999_999.5, 64.25, 29_999_999.125),
            -0.34497749247741055,
            "perlin at world border failed"
        );
    }

    #[test]
    fn perlin_f32() {
        let perlin = PerlinNoise::<f32>::new(&mut JavaRandom::new(SEED));
        assert_eq!(
            perlin.get(0.0, 0.0, 0.0),
            0.118487954,
            "perlin at 0.0 failed"
        );
        assert_eq!(
            perlin.get(-0.5, -0.5, -0.5),
            -0.087263405,
            "perlin at -0.5 failed"
        );
        assert_eq!(
            perlin.get(-29_999_999.5, 64.25, 29_999_999.125),
            -0.3449775,
            "perlin at world border failed"
        );
    }

    #[derive(Debug, Clone, Copy)]
    struct NormalNoiseExpected<F: Float> {
        zero: F,
        negative_one: F,
        mixed: F,
    }

    fn test_normal_noise<V: Version>(expected: NormalNoiseExpected<V::Float>) {
        let noise: NormalNoise<V> =
            V::new_normal_noise(&mut JavaRandom::new(SEED), V::OCTAVE, V::AMPLITUDE);
        let scale = noise.second.1.recip();

        assert_eq!(
            noise.get(0.0, 0.0, 0.0),
            expected.zero,
            "normal noise at 0.0 failed"
        );
        assert_eq!(
            noise.get(-scale, -scale, -scale),
            expected.negative_one,
            "normal noise at -1.0 failed"
        );
        assert_eq!(
            noise.get(0.5 * scale, 0.25 * scale, 0.125 * scale),
            expected.mixed,
            "normal with (0.5, 0.25, 0.125) failed"
        );
    }

    #[test]
    fn normal_noise_17() {
        test_normal_noise::<MC17>(NormalNoiseExpected {
            zero: -0.13714612453645852,
            negative_one: 0.2836656908436528,
            mixed: 0.009630604775887711,
        });
    }

    #[test]
    fn normal_noise_18() {
        test_normal_noise::<MC18>(NormalNoiseExpected {
            zero: -0.29788683109182007,
            negative_one: -0.09486360661852414,
            mixed: 0.13893713975091695,
        });
    }

    #[test]
    fn normal_noise_26_3() {
        test_normal_noise::<MC263>(NormalNoiseExpected {
            zero: -0.2978868,
            negative_one: -0.094863616,
            mixed: 0.13893716,
        });
    }
}
