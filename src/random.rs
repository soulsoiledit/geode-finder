#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xoroshiro_random() {
        let mut random = Xoroshiro128PlusPlusRandom::with_seed(1);
        assert_eq!(random.next_bits(32), -240669518);
        assert_eq!(random.next_long(), 6451672560482867317);
        assert_eq!(random.next_int(256), 12);
        assert_eq!(random.next_double(), 0.4388219140570252);
        assert_eq!(random.next_float(), 0.883254647);
        assert_eq!(random.next_between(16, 64), 42);
    }

    #[test]
    fn java_random() {
        let mut random = JavaRandom::with_seed(1);
        assert_eq!(random.next_bits(32), -1155869325);
        assert_eq!(random.next_long(), 1853403699951111791);
        assert_eq!(random.next_int(256), 104);
        assert_eq!(random.next_double(), 0.20771484130971707);
        assert_eq!(random.next_float(), 0.332717001);
        assert_eq!(random.next_between(16, 64), 17);
    }

    #[test]
    fn java_random_split() {
        let mut random = JavaRandom::with_seed(1);
        let mut split_random = random.next_split();
        assert_eq!(split_random.next_bits(32), 1217617485);
        assert_eq!(split_random.next_long(), -5858143487791806086);
        assert_eq!(split_random.next_int(256), 164);
        assert_eq!(split_random.next_double(), 0.04430274619703056);
        assert_eq!(split_random.next_float(), 0.801519394);
        assert_eq!(split_random.next_between(16, 64), 18);
    }
}

pub struct Xoroshiro128PlusPlusRandom {
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
        let i = self.next_bits(32) as i64;
        let j = self.next_bits(32) as i64;
        i.wrapping_shl(32).wrapping_add(j)
    }

    fn next_int(&mut self, max: i32) -> i32 {
        if (max & -max) == max {
            return (max as i64)
                .wrapping_mul(self.next_bits(31) as i64)
                .wrapping_shr(31) as i32;
        }

        let mut i = self.next_bits(31);
        let mut j = i % max;

        while (i.wrapping_sub(j).wrapping_add(max - 1)) < 0 {
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

impl Xoroshiro128PlusPlusRandom {
    pub fn with_seed(seed: i64) -> Self {
        let mut random = Xoroshiro128PlusPlusRandom {
            lower_bits: 0,
            upper_bits: 0,
        };

        random.set_seed(seed);

        random
    }
}

impl Random for Xoroshiro128PlusPlusRandom {
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
        let value = self
            .lower_bits
            .wrapping_add(self.upper_bits)
            .rotate_left(17)
            .wrapping_add(self.lower_bits);

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
        // 440898198 is Java hash code for "octave_-4"
        JavaRandom::with_seed(440898198 ^ self.next_long())
    }

    pub fn skip(&mut self, steps: i32) {
        for _ in 0..steps {
            self.next_bits(32);
        }
    }
}

impl Random for JavaRandom {
    fn set_seed(&mut self, seed: i64) {
        self.seed = seed ^ 0x5DEECE66D & 0xFFFFFFFFFFFF
    }

    fn next_bits(&mut self, bits: u32) -> i32 {
        self.seed = self.seed.wrapping_mul(0x5DEECE66D) + 11 & 0xFFFFFFFFFFFF;
        self.seed.wrapping_shr(48 - bits) as i32
    }
}
