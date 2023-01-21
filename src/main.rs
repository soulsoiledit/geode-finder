use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use java_utils::HashCode;
use std::{collections::HashMap, num::Wrapping};

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

struct Xoroshiro128PlusPlus {
    upper_bits: i64,
    lower_bits: i64,
}

struct JavaRandom {
    seed: i64,
}

trait Random {
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

struct PerlinNoiseSampler {
    origin_x: f64,
    origin_y: f64,
    origin_z: f64,
    permutations: [i8; 256],
}

struct OctavePerlinNoiseSampler {
    perlin_samplers: PerlinNoiseSampler,
}

struct DoublePerlinNoiseSampler {
    first_sampler: OctavePerlinNoiseSampler,
    second_sampler: OctavePerlinNoiseSampler,
}

struct GeodeGenerator {
    random: Box<dyn Random>,
    world_seed: i64,
    a: i64,
    b: i64,
    is_17: bool,
    noise_source: DoublePerlinNoiseSampler,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short = 'v', long, default_value_t = String::from("1.18"))]
    game_version: String,

    #[arg(short = 's', long, default_value_t = 011377)]
    seed: i64,

    #[arg(short = 'r', long, default_value_t = 10000)]
    search_radius: i64,

    #[arg(short = 'g', long, default_value_t = 26)]
    geode_threshold: i32,

    #[arg(short = 'b', long, default_value_t = 1000)]
    budding_threshold: i32,
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

impl Xoroshiro128PlusPlus {
    fn with_seed(seed: i64) -> Self {
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
    fn with_seed(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ 0x5DEECE66D) & 0xFFFFFFFFFFFF,
        }
    }

    fn next_split(&mut self) -> Self {
        JavaRandom::with_seed("octave_-4".hash_code() as i64 ^ self.next_long())
    }

    fn skip(&mut self, next: i32) {
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
    fn new(mut random: JavaRandom, is_17: bool) -> Self {
        let first_sampler = OctavePerlinNoiseSampler::new(&mut random, is_17);
        let second_sampler = OctavePerlinNoiseSampler::new(&mut random, is_17);

        DoublePerlinNoiseSampler {
            first_sampler,
            second_sampler,
        }
    }

    fn sample(&self, x: f64, y: f64, z: f64) -> f64 {
        let d = x * 1.0181268882175227;
        let e = y * 1.0181268882175227;
        let f = z * 1.0181268882175227;

        (self.first_sampler.sample(x, y, z) + self.second_sampler.sample(d, e, f))
            * ((1.0 / 6.0) * 5.0)
    }
}

impl GeodeGenerator {
    fn new(world_seed: i64, is_17: bool) -> Self {
        let mut random: Box<dyn Random>;

        if is_17 {
            random = Box::new(JavaRandom::with_seed(world_seed));
        } else {
            random = Box::new(Xoroshiro128PlusPlus::with_seed(world_seed));
        }

        let a = random.next_long() | 1;
        let b = random.next_long() | 1;

        GeodeGenerator {
            random,
            world_seed,
            a,
            b,
            is_17,
            noise_source: DoublePerlinNoiseSampler::new(JavaRandom::with_seed(world_seed), is_17),
        }
    }

    fn set_decorator_seed(&mut self, chunk_x: i64, chunk_z: i64, salt: i64) -> () {
        self.random.set_seed(
            (chunk_x
                .wrapping_shl(4)
                .wrapping_mul(self.a)
                .wrapping_add(chunk_z.wrapping_shl(4).wrapping_mul(self.b))
                ^ self.world_seed)
                .wrapping_add(salt),
        );
    }

    fn check_chunk(&mut self) -> bool {
        if self.is_17 {
            self.random.next_float() < (1.0 / 53.0)
        } else {
            self.random.next_float() < (1.0 / 24.0)
        }
    }

    fn fast_inv_sqrt(&mut self, x: f64) -> f64 {
        let y = f64::from_bits(6910469410427058090 - (x.to_bits() >> 1));
        y * (1.5 - 0.5 * x * y * y)
    }

    fn generate(&mut self, chunk_x: i64, chunk_z: i64) -> i32 {
        let origin = BlockPos {
            x: self.random.next_int(16) + 16 * chunk_x as i32,
            z: self.random.next_int(16) + 16 * chunk_z as i32,
            y: if self.is_17 {
                self.random.next_between(6, 46)
            } else {
                self.random.next_between(-58, 30)
            },
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
            let dx = if n % 2 == 0 { o } else { 0 };
            let dz = if n == 1 || n == 2 { o } else { 0 };

            block_list2.push(origin.add(dx, 7, dz));
            block_list2.push(origin.add(dx, 5, dz));
            block_list2.push(origin.add(dx, 1, dz));
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

                        if bl2 {
                            budding_count += 1;
                        }

                        if !bl2 || self.random.next_float() < 0.35 {
                            continue;
                        }
                    }
                }
            }
        }
        budding_count
    }
}

fn main() {
    let args = Args::parse();
    let seed = args.seed;
    let search_radius = args.search_radius;
    let geode_threshold = args.geode_threshold;
    let budding_threshold = args.budding_threshold;
    let is_17 = args.game_version.as_str() == "1.17";
    let salt: i64 = if is_17 { 20000 } else { 20002 };

    if is_17 {
        println!("Running geode search for version 1.17...");
    } else {
        println!("Running geode search for versions 1.18+...");
    }

    let mut finder = GeodeGenerator::new(seed, is_17);

    let progress_bar = ProgressBar::new(search_radius as u64 * 2 + 1).with_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed}] [{bar:.green/white}] ({eta_precise})",
        )
        .unwrap()
        .progress_chars("ùwú"),
    );

    let mut possible_locations: Vec<(i64, i64)> = vec![];
    for i in -search_radius..search_radius {
        let mut id = 0;
        let mut geode_count = [0; 15];
        let mut total_geodes = 0;

        for j in -search_radius..search_radius {
            total_geodes -= geode_count[id];
            geode_count[id] = 0;

            for x in 0..15 {
                finder.set_decorator_seed(i + x, j, salt);
                if finder.check_chunk() {
                    geode_count[id] += 1;
                }
            }

            total_geodes += geode_count[id];

            if total_geodes >= geode_threshold {
                possible_locations.push((i, j));
            }

            id = (id + 1) % 15;
        }

        progress_bar.inc(1);
    }
    progress_bar.finish();

    println!("{} possible locations found", possible_locations.len());

    let mut budding_amethyst_counts = HashMap::new();
    for loc in possible_locations {
        let minx = loc.0;
        let miny = loc.1 + 1;
        let mut budding = 0;
        for i in (minx)..(minx + 15) {
            for j in (miny - 15)..(miny) {
                if budding_amethyst_counts.contains_key(&(i, j)) {
                    budding += budding_amethyst_counts[&(i, j)];
                } else {
                    finder.set_decorator_seed(i, j, salt);
                    if finder.check_chunk() {
                        let budding_geode = finder.generate(i, j);
                        budding += budding_geode;
                        budding_amethyst_counts.insert((i, j), budding_geode);
                    }
                }
            }
        }

        if budding >= budding_threshold {
            println!(
                "Geode cluster with {} budding amethyst centered at {} {}",
                budding,
                (loc.0 + 7) * 16,
                (loc.1 - 7) * 16
            );
        }
    }
}
