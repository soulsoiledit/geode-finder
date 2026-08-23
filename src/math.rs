use std::ops::RangeInclusive;

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
    const FLOAT_MULTIPLIER: f32 = 1.0 / (1u32 << 24) as f32;
    const DOUBLE_MULTIPLIER: f64 = 1.0 / (1u64 << 53) as f64;

    fn new(seed: i64) -> Self;
    fn set_seed(&mut self, seed: i64);
    fn next_bits(&mut self, bits: u32) -> i32;
    // skips progressing the seed
    fn next_seed_bits(seed: i64, bits: u32) -> i32;

    fn next_int(&mut self, bound: i32) -> i32 {
        assert!(bound > 0, "bound should always be positive: found {bound}");

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

    fn next_between(&mut self, range: RangeInclusive<i32>) -> i32 {
        let (min, max) = (range.start(), range.end());
        self.next_int(max - min + 1) + min
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
    const GOLDEN_RATIO: i64 = 0x9E37_79B9_7F4A_7C15_u64.cast_signed(); // -7046029254386353131
    const SILVER_RATIO: i64 = 0x6A09_E667_F3BC_C909_u64.cast_signed(); //  7640891576956012809

    const STAFFORD_MIX13_1: u64 = 0xBF58_476D_1CE4_E5B9;
    const STAFFORD_MIX13_2: u64 = 0x94D0_49BB_1331_11EB;

    const fn splitmix64(z: i64) -> i64 {
        let mut z = z.cast_unsigned();
        z = (z ^ z.wrapping_shr(30)).wrapping_mul(Self::STAFFORD_MIX13_1);
        z = (z ^ z.wrapping_shr(27)).wrapping_mul(Self::STAFFORD_MIX13_2);
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

        self.seed_lo = Self::splitmix64(seed_lo);
        self.seed_hi = Self::splitmix64(seed_hi);

        // Should never happen
        if (self.seed_lo | self.seed_hi) == 0 {
            self.seed_lo = Self::GOLDEN_RATIO;
            self.seed_hi = Self::SILVER_RATIO;
        }
    }

    fn next_bits(&mut self, bits: u32) -> i32 {
        let lo = self.seed_lo;
        let mut hi = self.seed_hi;
        let result = lo.wrapping_add(hi).rotate_left(17).wrapping_add(lo);

        hi ^= lo;
        self.seed_lo = lo.rotate_left(49) ^ hi ^ (hi.wrapping_shl(21));
        self.seed_hi = hi.rotate_left(28);

        result.cast_unsigned().wrapping_shr(64 - bits) as i32
    }

    fn next_seed_bits(seed: i64, bits: u32) -> i32 {
        let mut lo = seed ^ Self::SILVER_RATIO;
        let mut hi = lo.wrapping_add(Self::GOLDEN_RATIO);

        lo = Self::splitmix64(lo);
        hi = Self::splitmix64(hi);

        if (lo | hi) == 0 {
            lo = Self::GOLDEN_RATIO;
            hi = Self::SILVER_RATIO;
        }

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
    const MOD_MASK: i64 = 0xFFFF_FFFF_FFFF;
    const MULTIPLIER: i64 = 0x5DEECE66D;
    const INCREMENT: i64 = 11;

    pub fn fork_from_hash(&mut self) -> Self {
        // Java "octave_-4".hashCode()
        const HASH: i64 = 440898198;
        Self::new(HASH ^ self.next_long())
    }
}

impl Random for JavaRandom {
    fn new(seed: i64) -> Self {
        let mut state = Self { seed: 0 };
        state.set_seed(seed);
        state
    }

    fn set_seed(&mut self, seed: i64) {
        self.seed = (seed ^ Self::MULTIPLIER) & Self::MOD_MASK;
    }

    fn next_bits(&mut self, bits: u32) -> i32 {
        self.seed = (self
            .seed
            .wrapping_mul(Self::MULTIPLIER)
            .wrapping_add(Self::INCREMENT))
            & Self::MOD_MASK;
        self.seed.wrapping_shr(Self::MOD_BITS.wrapping_sub(bits)) as i32
    }

    fn next_seed_bits(seed: i64, bits: u32) -> i32 {
        let mut next_seed = (seed ^ Self::MULTIPLIER) & Self::MOD_MASK;
        next_seed = (next_seed.wrapping_mul(Self::MULTIPLIER) + Self::INCREMENT) & Self::MOD_MASK;
        next_seed.wrapping_shr(Self::MOD_BITS.wrapping_sub(bits)) as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEED: i64 = 0xDEAD_BEEF;

    #[test]
    fn random() {
        let mut random = JavaRandom::new(SEED);
        assert_eq!(
            random.next_int(256),
            33,
            "next_int failed for power-of-2 bound"
        );
        assert_eq!(
            random.next_int(192),
            152,
            "next_int failed for non power-of-2 bound"
        );
        assert_eq!(random.next_between(192..=256), 222, "next_between failed");
        assert_eq!(random.next_long(), 1683040361022731026, "next_long failed");
        assert_eq!(random.next_float(), 0.8564069, "next_float failed");
        assert_eq!(
            random.next_double(),
            0.544608645520025,
            "next_double failed"
        );
        random.skip(256);
        assert_eq!(random.next_int(256), 103, "skip failed");
    }

    #[test]
    fn java_random() {
        let mut random = JavaRandom::new(SEED);

        let bits = (0..256).map(|_| random.next_bits(32)).last();
        assert_eq!(bits, Some(-1183353915), "random progression failed");

        random.set_seed(SEED);
        for _ in 0..256 {
            let seed = random.next_long();
            assert_eq!(
                {
                    random.set_seed(seed);
                    random.next_bits(32)
                },
                JavaRandom::next_seed_bits(seed, 32),
                "next_bits and next_seed_bits mismatched"
            );
        }

        random.set_seed(SEED);
        assert_eq!(
            random.fork_from_hash().next_long(),
            -3702631657396612959,
            "forked next_int failed"
        );
    }

    #[test]
    fn xoroshiro_random() {
        let mut random = Xoroshiro128PlusPlusRandom::new(SEED);

        let bits = (0..256).map(|_| random.next_bits(32)).last();
        assert_eq!(bits, Some(1573955217), "random progression failed");

        random.set_seed(SEED);
        for _ in 0..256 {
            let seed = random.next_long();
            assert_eq!(
                {
                    random.set_seed(seed);
                    random.next_bits(32)
                },
                Xoroshiro128PlusPlusRandom::next_seed_bits(seed, 32),
                "next_bits and next_seed_bits mismatched"
            );
        }
    }
}
