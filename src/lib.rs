//! Rust implementation of a 64-bit Mersenne Twister
//! https://en.wikipedia.org/wiki/Mersenne_Twister
//! https://dl.acm.org/doi/pdf/10.1145/369534.369540
//!
//! Ported from C version sourced from here:
//! http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
//! and included intact in README_C.txt

const NN: usize = 312;
const MM: usize = 156;
const UPPER_MASK: u64 = 0xFFFFFFFF80000000; // Most significant 33 bits
const LOWER_MASK: u64 = 0x7FFFFFFF; // Least significant 31 bits
const MAG01: [u64; 2] = [0, 0xB5026F5AA96619E9];

pub struct MersenneTwister64 {
    v: [u64; NN],
    i: usize,
}

impl MersenneTwister64 {
    pub const fn new(seed: u64) -> Self {
        let mut v = [0; NN];
        v[0] = seed;
        let mut i = 1;
        while i < NN {
            //for i in 1..NN {
            v[i] = (v[i - 1] ^ (v[i - 1] >> 62)).wrapping_mul(6364136223846793005) + i as u64;
            i += 1;
        }
        MersenneTwister64 { v, i: NN }
    }

    pub fn new_init_by_array(a: &[u64]) -> Self {
        let mut mt = MersenneTwister64::new(19650218);
        let mut i = 1;
        let mut j = 0;
        for _ in 0..std::cmp::max(NN, a.len()) {
            mt.v[i] = (mt.v[i]
                ^ ((mt.v[i - 1] ^ (mt.v[i - 1].wrapping_shr(62)))
                    .wrapping_mul(3935559000370003845)))
            .wrapping_add(a[j])
            .wrapping_add(j as u64);
            j = (j + 1) % a.len();
            i += 1;
            if i >= NN {
                mt.v[0] = mt.v[NN - 1];
                i = 1;
            }
        }
        for _ in 0..NN - 1 {
            mt.v[i] = (mt.v[i]
                ^ ((mt.v[i - 1] ^ (mt.v[i - 1].wrapping_shr(62)))
                    .wrapping_mul(2862933555777941757)))
            .wrapping_sub(i as u64); // non linear
            i += 1;
            if i >= NN {
                mt.v[0] = mt.v[NN - 1];
                i = 1;
            }
        }
        mt.v[0] = 1 << 63; // MSB is 1; assuring non-zero initial array
        mt
    }

    /// return a random integer from [0, 2^64-1]
    pub fn genrand(&mut self) -> u64 {
        if self.i >= NN {
            // Generate next NN words
            for i in 0..(NN - MM) {
                let x = (self.v[i] & UPPER_MASK) | (self.v[i + 1] & LOWER_MASK);
                self.v[i] = self.v[i + MM] ^ (x >> 1) ^ MAG01[(x & 1) as usize];
            }
            for i in (NN - MM)..NN {
                let x = (self.v[i] & UPPER_MASK) | (self.v[(i + 1) % NN] & LOWER_MASK);
                self.v[i] = self.v[i + MM - NN] ^ (x >> 1) ^ MAG01[(x & 1) as usize];
            }
            self.i = 0;
        }

        let mut x = self.v[self.i];
        self.i += 1;

        x ^= x.wrapping_shr(29) & 0x5555555555555555;
        x ^= x.wrapping_shl(17) & 0x71D67FFFEDA60000;
        x ^= x.wrapping_shl(37) & 0xFFF7EEE000000000;
        x ^ x.wrapping_shr(43)
    }

    /// return an array of random integers
    pub fn genrand_array<const N: usize>(&mut self) -> Box<[u64; N]> {
        // start with Vec to avoid stack allocation at any point
        let mut a = vec![0; N];
        for x in a.iter_mut() {
            *x = self.genrand();
        }
        a.try_into().expect("Failed to create Box<array> from Vec")
    }

    /// return a real random number from the interval [0,1]
    pub fn genrand_real1(&mut self) -> f64 {
        (self.genrand() >> 11) as f64 * (1.0 / 9007199254740991.0)
    }

    /// return a real random number from the interval [0,1[
    pub fn genrand_real2(&mut self) -> f64 {
        (self.genrand() >> 11) as f64 * (1.0 / 9007199254740992.0)
    }

    /// return a real random number from the interval ]0,1[
    pub fn genrand_real3(&mut self) -> f64 {
        ((self.genrand() >> 12) as f64 + 0.5) * (1.0 / 4503599627370496.0)
    }
}

/// convert u64 to a float in the range [0,1[
pub const fn u64_to_f64(x: u64) -> f64 {
    // shift x because of mantissa bits
    (x >> 11) as f64 * (1.0 / 9007199254740992.0)
}

pub const fn genrand_matrix<const NROWS: usize, const NCOLS: usize>(
    seed: u64,
) -> [[u64; NCOLS]; NROWS] {
    let mut mt = MersenneTwister64::new(seed);
    let mut matrix = [[0u64; NCOLS]; NROWS];
    let mut row = 0;
    while row < NROWS {
        let mut col = 0;
        while col < NCOLS {
            if mt.i >= NN {
                // Generate next NN words
                let mut i = 0;
                while i < NN - MM {
                    let x = (mt.v[i] & UPPER_MASK) | (mt.v[i + 1] & LOWER_MASK);
                    mt.v[i] = mt.v[i + MM] ^ (x >> 1) ^ MAG01[(x & 1) as usize];
                    i += 1;
                }
                let mut i = NN - MM;
                while i < NN {
                    let x = (mt.v[i] & UPPER_MASK) | (mt.v[(i + 1) % NN] & LOWER_MASK);
                    mt.v[i] = mt.v[i + MM - NN] ^ (x >> 1) ^ MAG01[(x & 1) as usize];
                    i += 1;
                }
                mt.i = 0;
            }

            let mut x = mt.v[mt.i];
            mt.i += 1;

            x ^= x.wrapping_shr(29) & 0x5555555555555555;
            x ^= x.wrapping_shl(17) & 0x71D67FFFEDA60000;
            x ^= x.wrapping_shl(37) & 0xFFF7EEE000000000;
            matrix[row][col] = x ^ x.wrapping_shr(43);
            col += 1;
        }
        row += 1;
    }
    matrix
}
