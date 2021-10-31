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

struct Mont {
    pub rp1: u32,
    pub np1: u32,
    pub r: u64,
    pub n: u32,
    pub bits: usize,
}

impl Mont {
    pub fn new(n: u32) -> Self {
        let bits = n.to_be_bytes().len() * 8;
        let r = 2u64.pow(bits as u32);
        Mont {
            r,
            n,
            rp1: 0,
            np1: 0,
            bits,
        }
    }

    pub fn precompute(&mut self) {
        let r = self.r as i64;
        let n = self.n as i64;
        let (gcd, np, rp) = egcd(n, r);
        assert_eq!(gcd, 1);
        self.rp1 = if rp < 0 { -rp as u32 } else { rp as u32 };
        self.np1 = if np < 0 { -np as u32 } else { np as u32 };
    }
    // m = T*Np1 mod R
    // U = (T + m * N) / R
    // The overall process delivers T · R−1 mod N without an expensive division operation!
    pub fn reduce(&self, t: u64) -> u32 {
        let m = (t * (self.np1 as u64)) as u32;
        let u = (t + (m as u64) * (self.n as u64)) >> self.bits;
        u as u32
    }
    pub fn to_mont(&self, x: u32) -> u32 {
        // divide n, need a lot of cycles
        let res = ((x as u64) * self.r) % (self.n as u64);
        res as u32
    }
    pub fn multi(&self, x: u32, y: u32) -> u32 {
        let xy = x as u64 * y as u64;
        self.reduce(xy)
    }
}

#[test]
pub fn test() {
    let mut mont = Mont::new(17);
    mont.precompute();

    let x: u32 = 10;
    let x2 = mont.to_mont(x);
    let x3 = mont.reduce(x2 as u64);
    assert_eq!(x, x3);

    let x: u32 = 100;
    let y: u32 = 200;
    // into montgomery form
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    // do multiplication operation
    let xy2 = mont.multi(x2, y2);
    // into normal form
    let xy = mont.reduce(xy2 as u64);
    // the result should be same
    assert_eq!(xy, (x * y) % mont.n);
}

pub fn main() {}
