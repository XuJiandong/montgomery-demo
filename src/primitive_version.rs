// based on https://research.nccgroup.com/2021/06/09/optimizing-pairing-based-cryptography-montgomery-arithmetic-in-rust/

// Extended Euclidean algorithm
// https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
pub fn egcd(a: i64, b: i64) -> (i64, i64, i64) {
    assert!(a < b);
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = egcd(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

pub struct Mont {
    pub np1: u32,
    pub r: u64,
    pub n: u32,
    pub bits: usize,
    pub init: bool,
}

impl Mont {
    pub fn new(n: u32) -> Self {
        let bits = n.to_be_bytes().len() * 8;
        let r = 2u64.pow(bits as u32);
        Mont {
            r,
            n,
            np1: 0,
            bits,
            init: false,
        }
    }

    pub fn precompute(&mut self) {
        let r = self.r as i64;
        let n = self.n as i64;
        let (gcd, np, rp) = egcd(n, r);
        assert_eq!(gcd, 1);

        let rp1 = rp + self.n as i64;
        assert!(rp1 >= 0);

        let np1 = self.r as i64 - np;
        assert!(np1 >= 0);
        self.np1 = np1 as u32; // can be truncated

        let r = self.r * rp1 as u64;
        let n = self.n as u64 * np1 as u64;
        assert_eq!(n + 1, r);

        self.init = true;
    }
    // m = T*Np1 mod R
    // U = (T + m * N) / R
    // The overall process delivers T · R−1 mod N without an expensive division operation!
    pub fn reduce(&self, t: u64) -> u32 {
        assert!(self.init);
        let t0 = t as u32 as u64; // low part of `t`, same as `% self.r`, avoid overflow
        let m = (t0 * self.np1 as u64) as u32 as u64; // `% self.r`
        let u = (t + m * (self.n as u64)) >> self.bits; // `/ self.r`
        if u >= self.n as u64 {
            (u - self.n as u64) as u32
        } else {
            u as u32
        }
    }
    pub fn to_mont(&self, x: u32) -> u32 {
        assert!(self.init);
        // divide n, need a lot of cycles
        let res = ((x as u64) * self.r) % (self.n as u64);
        res as u32
    }
    pub fn multi(&self, x: u32, y: u32) -> u32 {
        let xy = x as u64 * y as u64;
        self.reduce(xy)
    }
    pub fn pow(&self, x: u32, y: u32) -> u32 {
        let mut base = x;
        let mut res: u32 = 0;
        let mut first_time: bool = true;

        for index in 0..self.bits {
            // at most self.bits (32 here) multiplications
            if ((y >> index) & 1) == 1 {
                if first_time {
                    // note:
                    // `res = self.multi(1, base)` is not same
                    // as res = base;
                    res = base;
                    first_time = false;
                } else {
                    res = self.multi(res, base);
                }
            }
            base = self.multi(base, base); // at most self.bits (32 here) multiplications
        }
        res
    }
}
