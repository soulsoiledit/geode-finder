use crate::math::{JavaRandom, Random};

#[derive(Clone, Copy)]
pub struct ImprovedNoise {
    x_offset: f64,
    y_offset: f64,
    z_offset: f64,
    permutation: [i8; 256],
}

pub trait Noise {
    fn get_value(&self, x: f64, y: f64, z: f64) -> f64;
}

impl ImprovedNoise {
    const GRADIENT: [[f64; 3]; 16] = [
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [1.0, -1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, 0.0, 1.0],
        [-1.0, 0.0, 1.0],
        [1.0, 0.0, -1.0],
        [-1.0, 0.0, -1.0],
        [0.0, 1.0, 1.0],
        [0.0, -1.0, 1.0],
        [0.0, 1.0, -1.0],
        [0.0, -1.0, -1.0],
        [1.0, 1.0, 0.0],
        [0.0, -1.0, 1.0],
        [-1.0, 1.0, 0.0],
        [0.0, -1.0, -1.0],
    ];

    pub fn new(random: &mut JavaRandom) -> Self {
        let x_offset = random.next_double() * 256.0;
        let y_offset = random.next_double() * 256.0;
        let z_offset = random.next_double() * 256.0;

        let mut permutation: [i8; 256] = std::array::from_fn(|i| i as i8);
        for i in 0..256 {
            let remaining = 256 - i;
            let offset = random.next_int(remaining as i32).cast_unsigned() as usize;
            permutation.swap(i, i + offset);
        }

        Self {
            x_offset,
            y_offset,
            z_offset,
            permutation,
        }
    }

    fn p(&self, i: i32) -> i32 {
        i32::from(self.permutation[(i & 0xFF) as usize])
    }

    fn grad_dot(&self, hash: i32, x: f64, y: f64, z: f64) -> f64 {
        let [gx, gy, gz] = Self::GRADIENT[(hash & 0xF) as usize];
        gx * x + gy * y + gz * z
    }

    fn smooth_step(&self, x: f64) -> f64 {
        x * x * x * (x * (x * 6.0 - 15.0) + 10.0)
    }

    #[allow(non_snake_case)]
    fn lerp(&self, ax: f64, x: f64, X: f64) -> f64 {
        x + ax * (X - x)
    }

    #[allow(non_snake_case)]
    fn lerp2(&self, ax: f64, ay: f64, xy: f64, Xy: f64, xY: f64, XY: f64) -> f64 {
        self.lerp(ay, self.lerp(ax, xy, Xy), self.lerp(ax, xY, XY))
    }

    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn lerp3(
        &self,
        ax: f64,
        ay: f64,
        az: f64,
        xyz: f64,
        Xyz: f64,
        xYz: f64,
        XYz: f64,
        xyZ: f64,
        XyZ: f64,
        xYZ: f64,
        XYZ: f64,
    ) -> f64 {
        self.lerp(
            az,
            self.lerp2(ax, ay, xyz, Xyz, xYz, XYz),
            self.lerp2(ax, ay, xyZ, XyZ, xYZ, XYZ),
        )
    }

    #[allow(non_snake_case)]
    fn sample_lerp(&self, x: i32, y: i32, z: i32, ax: f64, ay: f64, az: f64) -> f64 {
        let px = self.p(x);
        let pX = self.p(x + 1);

        let xy = self.p(px + y);
        let Xy = self.p(pX + y);
        let xY = self.p(px + y + 1);
        let XY = self.p(pX + y + 1);

        let xyz = self.grad_dot(self.p(xy + z), ax, ay, az);
        let Xyz = self.grad_dot(self.p(Xy + z), ax - 1.0, ay, az);
        let xYz = self.grad_dot(self.p(xY + z), ax, ay - 1.0, az);
        let XYz = self.grad_dot(self.p(XY + z), ax - 1.0, ay - 1.0, az);
        let xyZ = self.grad_dot(self.p(xy + z + 1), ax, ay, az - 1.0);
        let XyZ = self.grad_dot(self.p(Xy + z + 1), ax - 1.0, ay, az - 1.0);
        let xYZ = self.grad_dot(self.p(xY + z + 1), ax, ay - 1.0, az - 1.0);
        let XYZ = self.grad_dot(self.p(XY + z + 1), ax - 1.0, ay - 1.0, az - 1.0);

        let ax = self.smooth_step(ax);
        let ay = self.smooth_step(ay);
        let az = self.smooth_step(az);

        self.lerp3(ax, ay, az, xyz, Xyz, xYz, XYz, xyZ, XyZ, xYZ, XYZ)
    }
}

impl Noise for ImprovedNoise {
    fn get_value(&self, x: f64, y: f64, z: f64) -> f64 {
        let xo = x + self.x_offset;
        let yo = y + self.y_offset;
        let zo = z + self.z_offset;

        let x = xo.floor();
        let y = yo.floor();
        let z = zo.floor();

        let ax = xo - x;
        let ay = yo - y;
        let az = zo - z;

        self.sample_lerp(x as i32, y as i32, z as i32, ax, ay, az)
    }
}

#[derive(Clone, Copy)]
pub struct PerlinNoise {
    noise: ImprovedNoise,
}

impl PerlinNoise {
    const ROUND_OFF: f64 = (1 << 26) as f64;
    const FREQUENCY: f64 = 1.0 / (1 << 4) as f64;

    pub const fn new(noise: ImprovedNoise) -> Self {
        Self { noise }
    }

    fn wrap(&self, value: f64) -> f64 {
        value - f64::floor(value / Self::ROUND_OFF + 0.5) * Self::ROUND_OFF
    }
}

impl Noise for PerlinNoise {
    fn get_value(&self, x: f64, y: f64, z: f64) -> f64 {
        self.noise.get_value(
            self.wrap(x * Self::FREQUENCY),
            self.wrap(y * Self::FREQUENCY),
            self.wrap(z * Self::FREQUENCY),
        )
    }
}

#[derive(Clone, Copy)]
pub struct NormalNoise {
    pub first: PerlinNoise,
    pub second: PerlinNoise,
}

impl NormalNoise {
    const INPUT_FACTOR: f64 = 1.0181268882175227;
    // floating point precision prevents us from simplifying this
    const AMPLITUDE: f64 = 1.0 / 6.0 * 5.0;

    pub const fn new(first: PerlinNoise, second: PerlinNoise) -> Self {
        Self { first, second }
    }
}

impl Noise for NormalNoise {
    fn get_value(&self, x: f64, y: f64, z: f64) -> f64 {
        let x2 = x * Self::INPUT_FACTOR;
        let y2 = y * Self::INPUT_FACTOR;
        let z2 = z * Self::INPUT_FACTOR;

        (self.first.get_value(x, y, z) + self.second.get_value(x2, y2, z2)) * Self::AMPLITUDE
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::version::{MC17, MC18, Version};

    const SEED: i64 = 0;
    const ZERO: f64 = 0.0;
    const MIN: f64 = -30_000_000.0;
    const MAX: f64 = 30_000_000.0;

    #[derive(Debug, PartialEq)]
    struct TestNoiseResult {
        zero: f64,
        min: f64,
        max: f64,
        mix: f64,
    }

    fn test_noise(noise: impl Noise) -> TestNoiseResult {
        TestNoiseResult {
            zero: noise.get_value(ZERO, ZERO, ZERO),
            min: noise.get_value(MIN, MIN, MIN),
            max: noise.get_value(MAX, MAX, MAX),
            mix: noise.get_value(MIN, ZERO, MAX),
        }
    }

    #[test]
    fn test_improved_noise() {
        assert_eq!(
            test_noise(ImprovedNoise::new(&mut JavaRandom::new(SEED))),
            TestNoiseResult {
                zero: -0.09566354243549174,
                min: -0.2971040982148004,
                max: -0.2971040982148004,
                mix: 0.2263779027387616,
            }
        )
    }

    #[test]
    fn test_perlin_noise_17() {
        let mut random = JavaRandom::new(SEED);
        assert_eq!(
            test_noise(PerlinNoise::new(MC17::new_noise(&mut random))),
            TestNoiseResult {
                zero: -0.12439944493997021,
                min: 0.09218348114911419,
                max: -0.31362102310331874,
                mix: -0.04021242727494904,
            }
        )
    }

    #[test]
    fn test_perlin_noise() {
        let mut random = JavaRandom::new(SEED);
        assert_eq!(
            test_noise(PerlinNoise::new(MC18::new_noise(&mut random))),
            TestNoiseResult {
                zero: 0.22153284884252,
                min: 0.11128716610713107,
                max: 0.014484740544551677,
                mix: 0.32216788680541897,
            }
        )
    }

    #[test]
    fn test_normal_noise() {
        assert_eq!(
            test_noise(MC18::new_normal_noise(SEED)),
            TestNoiseResult {
                zero: 0.6284375015501312,
                min: 0.29515151085881075,
                max: -0.15898889860095572,
                mix: 0.37658086774180666,
            }
        )
    }
}
