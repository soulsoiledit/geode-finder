pub fn inv_sqrt(x: f64) -> f64 {
    x.sqrt().recip()
}

#[derive(Debug, Clone, Copy)]
pub struct Block {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl Block {
    pub const fn new(x: i32, y: i32, z: i32) -> Self {
        Self { x, y, z }
    }

    pub const fn add(&self, dx: i32, dy: i32, dz: i32) -> Self {
        Self {
            x: self.x + dx,
            y: self.y + dy,
            z: self.z + dz,
        }
    }
}

pub trait Random {
    const FLOAT_MULTIPLIER: f32 = 1.0 / (1 << 24) as f32;
    const DOUBLE_MULTIPLIER: f64 = 1.0 / (1u64 << 53) as f64;

    fn new(seed: i64) -> Self;
    fn set_seed(&mut self, seed: i64);
    fn next_bits(&mut self, bits: u32) -> i32;

    fn next_int(&mut self, bound: i32) -> i32 {
        if bound.cast_unsigned().is_power_of_two() {
            return i64::from(bound)
                .wrapping_mul(i64::from(self.next_bits(31)))
                .wrapping_shr(31) as i32;
        }

        loop {
            let sample = self.next_bits(31);
            let modulo = sample % bound;
            if (sample.wrapping_sub(modulo).wrapping_add(bound - 1)) >= 0 {
                return modulo;
            }
        }
    }

    fn next_between(&mut self, mini: i32, maxi: i32) -> i32 {
        self.next_int(maxi - mini + 1) + mini
    }

    fn next_long(&mut self) -> i64 {
        let upper = i64::from(self.next_bits(32));
        let lower = i64::from(self.next_bits(32));
        upper.wrapping_shl(32).wrapping_add(lower)
    }

    fn next_float(&mut self) -> f32 {
        self.next_bits(24) as f32 * Self::FLOAT_MULTIPLIER
    }

    fn next_double(&mut self) -> f64 {
        let upper = i64::from(self.next_bits(26));
        let lower = i64::from(self.next_bits(27));
        upper.wrapping_shl(27).wrapping_add(lower) as f64 * Self::DOUBLE_MULTIPLIER
    }

    fn skip(&mut self, steps: i32) {
        for _ in 0..steps {
            self.next_bits(32);
        }
    }
}

#[derive(Clone, Copy)]
pub struct Xoroshiro128PlusPlusRandom {
    seed_hi: i64,
    seed_lo: i64,
}

impl Xoroshiro128PlusPlusRandom {
    const GOLDEN_RATIO: i64 = 0x9E37_79B9_7F4A_7C15_u64.cast_signed();
    const SILVER_RATIO: i64 = 0x6A09_E667_F3BC_C909_u64.cast_signed();

    const STAFFORD_1: u64 = 0xBF58_476D_1CE4_E5B9;
    const STAFFORD_2: u64 = 0x94D0_49BB_1331_11EB;

    // splitmix64
    const fn mix_stafford_13(&self, z: i64) -> i64 {
        let mut z = z.cast_unsigned();
        z = (z ^ z.wrapping_shr(30)).wrapping_mul(Self::STAFFORD_1);
        z = (z ^ z.wrapping_shr(27)).wrapping_mul(Self::STAFFORD_2);
        (z ^ z.wrapping_shr(31)).cast_signed()
    }
}

impl Random for Xoroshiro128PlusPlusRandom {
    fn new(seed: i64) -> Self {
        let mut state = Self {
            seed_lo: 0,
            seed_hi: 0,
        };
        state.set_seed(seed);
        state
    }

    fn set_seed(&mut self, seed: i64) {
        let seed_lo = seed ^ Self::SILVER_RATIO;
        let seed_hi = seed_lo.wrapping_add(Self::GOLDEN_RATIO);
        self.seed_lo = self.mix_stafford_13(seed_lo);
        self.seed_hi = self.mix_stafford_13(seed_hi);
    }

    fn next_bits(&mut self, bits: u32) -> i32 {
        let lo = self.seed_lo;
        let hi = self.seed_hi;
        let xor = hi ^ lo;
        self.seed_lo = lo.rotate_left(49) ^ xor ^ (xor.wrapping_shl(21));
        self.seed_hi = xor.rotate_left(28);
        let result = lo.wrapping_add(hi).rotate_left(17).wrapping_add(lo);
        result.cast_unsigned().wrapping_shr(64 - bits) as i32
    }
}

// LegacyRandomSource (LCG)
#[derive(Clone, Copy)]
pub struct JavaRandom {
    seed: i64,
}

impl JavaRandom {
    const MOD_BITS: u32 = 48;
    const MODULUS: i64 = 0xFFFF_FFFF_FFFF;
    const MULTIPLIER: i64 = 0x5DEECE66D;
    const INCREMENT: i64 = 11;
    // Java "octave_-4".hashCode()
    const HASH: i64 = 440898198;

    pub fn fork_from_hash(&mut self) -> Self {
        Self::new(Self::HASH ^ self.next_long())
    }
}

impl Random for JavaRandom {
    fn new(seed: i64) -> Self {
        let mut state = Self { seed: 0 };
        state.set_seed(seed);
        state
    }

    fn set_seed(&mut self, seed: i64) {
        self.seed = seed ^ Self::MULTIPLIER & Self::MODULUS;
    }

    fn next_bits(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(Self::MULTIPLIER) + Self::INCREMENT) & Self::MODULUS;
        self.seed.wrapping_shr(Self::MOD_BITS - bits) as i32
    }
}

#[derive(Clone, Copy)]
pub struct UniformInt {
    pub min: i32,
    pub max: i32,
}

impl UniformInt {
    pub const fn new(min: i32, max: i32) -> Self {
        Self { min, max }
    }

    pub fn sample(&self, random: &mut impl Random) -> i32 {
        random.next_between(self.min, self.max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEED: i64 = 0;

    #[derive(Debug, PartialEq)]
    struct TestRandomResult {
        bits: i32,
        int: i32,
        between: i32,
        long: i64,
        float: f32,
        double: f64,
        skip: i32,
    }

    fn test_random(random: &mut impl Random) -> TestRandomResult {
        TestRandomResult {
            bits: random.next_bits(32),
            int: random.next_int(256),
            between: random.next_between(16, 64),
            long: random.next_long(),
            float: random.next_float(),
            double: random.next_double(),
            skip: {
                random.skip(1);
                random.next_int(256)
            },
        }
    }

    #[test]
    fn test_xoroshiro_random() {
        assert_eq!(
            test_random(&mut Xoroshiro128PlusPlusRandom::new(SEED)),
            TestRandomResult {
                bits: 707568776,
                int: 204,
                between: 47,
                long: 2160572956399813066,
                float: 0.7566797,
                double: 0.7723285845227084,
                skip: 156,
            }
        );
    }

    #[test]
    fn test_java_random() {
        let mut random = JavaRandom::new(SEED);
        assert_eq!(
            test_random(&mut random),
            TestRandomResult {
                bits: -1155484576,
                int: 212,
                between: 41,
                long: -7261648964369397258,
                float: 0.30905056,
                double: 0.5504370051176339,
                skip: 200,
            }
        );
        assert_eq!(random.fork_from_hash().next_int(256), 2);
    }
}
