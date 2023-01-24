use java_utils::HashCode;

pub struct Xoroshiro128PlusPlus {
    upper_bits: i64,
    lower_bits: i64,
}

pub struct JavaRandom {
    seed: i64,
}

pub trait Random {
    fn set_seed(&mut self, seed: i64);

    fn next_bits(&mut self, bits: u32) -> i32;

    fn next_long(&mut self) -> i64 {
        let i = self.next_bits(32);
        let j = self.next_bits(32);
        let l = (i as i64).wrapping_shl(32);
        l.wrapping_add(j as i64)
    }

    fn next_int(&mut self, max: i32) -> i32 {
        if (max & -max) == max {
            return (max as i64)
                .wrapping_mul(self.next_bits(31) as i64)
                .wrapping_shr(31) as i32;
        }

        let mut i = self.next_bits(31);
        let mut j = i % max;

        while (i - j + (max - 1)) < 0 {
            i = self.next_bits(31);
            j = i % max;
        }

        j
    }

    fn next_between(&mut self, min: i32, max: i32) -> i32 {
        self.next_int(max - min + 1) + min
    }

    fn next_float(&mut self) -> f32 {
        self.next_bits(24) as f32 * 5.9604645e-8f32
    }

    fn next_double(&mut self) -> f64 {
        let i = self.next_bits(26) as i64;
        let j = self.next_bits(27) as i64;
        i.wrapping_shl(27).wrapping_add(j) as f64 * 1.110223e-16f32 as f64
    }
}

impl Xoroshiro128PlusPlus {
    pub fn with_seed(seed: i64) -> Self {
        let mut random = Xoroshiro128PlusPlus {
            lower_bits: 0,
            upper_bits: 0,
        };

        random.set_seed(seed);

        random
    }
}

impl Random for Xoroshiro128PlusPlus {
    fn set_seed(&mut self, seed: i64) {
        self.lower_bits = seed ^ 7640891576956012809;
        self.upper_bits = self.lower_bits.wrapping_sub(7046029254386353131);

        self.lower_bits = (self.lower_bits ^ (self.lower_bits as u64 >> 30) as i64)
            .wrapping_mul(-4658895280553007687);
        self.lower_bits = (self.lower_bits ^ (self.lower_bits as u64 >> 27) as i64)
            .wrapping_mul(-7723592293110705685);
        self.lower_bits ^= (self.lower_bits as u64 >> 31) as i64;

        self.upper_bits = (self.upper_bits ^ (self.upper_bits as u64 >> 30) as i64)
            .wrapping_mul(-4658895280553007687);
        self.upper_bits = (self.upper_bits ^ (self.upper_bits as u64 >> 27) as i64)
            .wrapping_mul(-7723592293110705685);
        self.upper_bits ^= (self.upper_bits as u64 >> 31) as i64;
    }

    fn next_bits(&mut self, bits: u32) -> i32 {
        let xor = self.upper_bits ^ self.lower_bits;
        let value = {
            self.lower_bits
                .wrapping_add(self.upper_bits)
                .rotate_left(17)
                .wrapping_add(self.lower_bits)
        };

        self.lower_bits = self.lower_bits.rotate_left(49) ^ xor ^ (xor << 21);
        self.upper_bits = xor.rotate_left(28);

        (value as u64).wrapping_shr(64 - bits) as i32
    }
}

impl JavaRandom {
    pub fn with_seed(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ 0x5DEECE66D) & 0xFFFFFFFFFFFF,
        }
    }

    pub fn next_split(&mut self) -> Self {
        JavaRandom::with_seed("octave_-4".hash_code() as i64 ^ self.next_long())
    }

    pub fn skip(&mut self, next: i32) {
        for _ in 0..next {
            self.next_bits(32);
        }
    }
}

impl Random for JavaRandom {
    fn set_seed(&mut self, seed: i64) {
        self.seed = seed ^ 0x5DEECE66D & 0xFFFFFFFFFFFF
    }

    fn next_bits(&mut self, bits: u32) -> i32 {
        self.seed = self.seed.wrapping_mul(25214903917) + 11 & 0xFFFFFFFFFFFF;
        return self.seed.wrapping_shr(48 - bits) as i32;
    }
}
