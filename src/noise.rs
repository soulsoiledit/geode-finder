use std::num::Wrapping;

use crate::random::JavaRandom;
use crate::random::Random;

const NOISE_GRADIENTS: [[f64; 3]; 16] = [
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

pub struct PerlinNoiseSampler {
    origin_x: f64,
    origin_y: f64,
    origin_z: f64,
    permutations: [i8; 256],
}

pub struct OctavePerlinNoiseSampler {
    perlin_samplers: PerlinNoiseSampler,
}

pub struct DoublePerlinNoiseSampler {
    first_sampler: OctavePerlinNoiseSampler,
    second_sampler: OctavePerlinNoiseSampler,
}

impl PerlinNoiseSampler {
    fn new(random: &mut JavaRandom) -> Self {
        let origin_x = random.next_double() * 256.0;
        let origin_y = random.next_double() * 256.0;
        let origin_z = random.next_double() * 256.0;

        let mut permutations = [0i8; 256];

        for index in 0..256 {
            permutations[index] = Wrapping(index as i8).0;
        }

        for index in 0..256 {
            let random_index = random.next_int(256 - index as i32);
            let old = permutations[index];
            permutations[index] = permutations[index + random_index as usize];
            permutations[index + random_index as usize] = old;
        }

        PerlinNoiseSampler {
            origin_x,
            origin_y,
            origin_z,
            permutations,
        }
    }

    fn map(&self, input: i32) -> i32 {
        (Wrapping(self.permutations[(input & 0xFF) as usize]) & (Wrapping(0xFFu8 as i8))).0 as i32
    }

    fn grad(&self, hash: i32, x: f64, y: f64, z: f64) -> f64 {
        let gradient = NOISE_GRADIENTS[(hash & 0xF) as usize];
        gradient[0] * x + gradient[1] * y + gradient[2] * z
    }

    fn perlin_fade(&self, val: f64) -> f64 {
        val * val * val * (val * (val * 6.0 - 15.0) + 10.0)
    }

    fn lerp(&self, delta: f64, x0: f64, x1: f64) -> f64 {
        x0 + delta * (x1 - x0)
    }

    fn lerp2(&self, dx: f64, dy: f64, x0: f64, x1: f64, x2: f64, x3: f64) -> f64 {
        self.lerp(dy, self.lerp(dx, x0, x1), self.lerp(dx, x2, x3))
    }

    fn lerp3(
        &self,
        dx: f64,
        dy: f64,
        dz: f64,
        x0: f64,
        x1: f64,
        x2: f64,
        x3: f64,
        x4: f64,
        x5: f64,
        x6: f64,
        x7: f64,
    ) -> f64 {
        return self.lerp(
            dz,
            self.lerp2(dx, dy, x0, x1, x2, x3),
            self.lerp2(dx, dy, x4, x5, x6, x7),
        );
    }

    fn sample(&self, x: f64, y: f64, z: f64) -> f64 {
        let d = x + self.origin_x;
        let e = y + self.origin_y;
        let f = z + self.origin_z;
        let i = d.floor();
        let j = e.floor();
        let k = f.floor();
        let g = d - i;
        let h = e - j;
        let l = f - k;

        self.sample_base(i as i32, j as i32, k as i32, g, h, l, h)
    }

    fn sample_base(
        &self,
        secx: i32,
        secy: i32,
        secz: i32,
        locx: f64,
        locy: f64,
        locz: f64,
        fade_locy: f64,
    ) -> f64 {
        let i = self.map(secx);
        let j = self.map(secx + 1);
        let k = self.map(i + secy);
        let l = self.map(i + secy + 1);
        let m = self.map(j + secy);
        let n = self.map(j + secy + 1);

        let d = self.grad(self.map(k + secz), locx, locy, locz);
        let e = self.grad(self.map(m + secz), locx - 1.0, locy, locz);
        let f = self.grad(self.map(l + secz), locx, locy - 1.0, locz);
        let g = self.grad(self.map(n + secz), locx - 1.0, locy - 1.0, locz);

        let h = self.grad(self.map(k + secz + 1), locx, locy, locz - 1.0);
        let o = self.grad(self.map(m + secz + 1), locx - 1.0, locy, locz - 1.0);
        let p = self.grad(self.map(l + secz + 1), locx, locy - 1.0, locz - 1.0);
        let q = self.grad(self.map(n + secz + 1), locx - 1.0, locy - 1.0, locz - 1.0);

        let r = self.perlin_fade(locx);
        let s = self.perlin_fade(fade_locy);
        let t = self.perlin_fade(locz);

        self.lerp3(r, s, t, d, e, f, g, h, o, p, q)
    }
}

impl OctavePerlinNoiseSampler {
    fn new(random: &mut JavaRandom, is_17: bool) -> Self {
        let octave_samplers: PerlinNoiseSampler;
        if is_17 {
            random.skip(262 * 4);
            octave_samplers = PerlinNoiseSampler::new(random);
        } else {
            let mut random2 = random.next_split();
            octave_samplers = PerlinNoiseSampler::new(&mut random2);
        }

        OctavePerlinNoiseSampler {
            perlin_samplers: octave_samplers,
        }
    }

    fn maintain_precision(&self, value: f64) -> f64 {
        value - f64::floor(value / 3.3554432e7 + 0.5) * 3.3554432e7
    }

    fn sample(&self, x: f64, y: f64, z: f64) -> f64 {
        self.perlin_samplers.sample(
            self.maintain_precision(x * 0.0625),
            self.maintain_precision(y * 0.0625),
            self.maintain_precision(z * 0.0625),
        )
    }
}

impl DoublePerlinNoiseSampler {
    pub fn new(mut random: JavaRandom, is_17: bool) -> Self {
        let first_sampler = OctavePerlinNoiseSampler::new(&mut random, is_17);
        let second_sampler = OctavePerlinNoiseSampler::new(&mut random, is_17);

        DoublePerlinNoiseSampler {
            first_sampler,
            second_sampler,
        }
    }

    pub fn sample(&self, x: f64, y: f64, z: f64) -> f64 {
        let d = x * 1.0181268882175227;
        let e = y * 1.0181268882175227;
        let f = z * 1.0181268882175227;

        (self.first_sampler.sample(x, y, z) + self.second_sampler.sample(d, e, f))
            * ((1.0 / 6.0) * 5.0)
    }
}
