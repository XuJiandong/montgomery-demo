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
    pub init: bool,
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
        self.rp1 = rp1 as u32;
        let np1 = self.r as i64 - np;
        assert!(np1 >= 0);
        self.np1 = np1 as u32;

        self.init = true;
    }
    // m = T*Np1 mod R
    // U = (T + m * N) / R
    // The overall process delivers T · R−1 mod N without an expensive division operation!
    pub fn reduce(&self, t: u64) -> u32 {
        assert!(self.init);
        let t0 = t as u32 as u64; // low part of `t`, same as `mode self.r`, avoid overflow
        let m = (t0 * self.np1 as u64) as u32 as u64; // `mod self.r`
        let u = (t + m * (self.n as u64)) >> self.bits; // `/self.r`
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
        let mut res: u32 = 1;
        for index in 0..31 {
            // at most 32 multiplications
            if ((y >> index) & 1) == 1 {
                res = self.multi(res, base);
            }
            base = self.multi(base, base); // at most 32 multiplications
        }
        res
    }
}

fn test_n(n: u32) {
    let mut mont = Mont::new(n);
    mont.precompute();

    let x: u32 = 10;
    let x2 = mont.to_mont(x);
    let x3 = mont.reduce(x2 as u64);
    assert_eq!(x, x3);

    let x: u32 = 10;
    let y: u32 = 20;
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

fn test_xy(x: u32, y: u32) {
    let mut mont = Mont::new(1000001);
    mont.precompute();

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

#[test]
pub fn test_n_loops() {
    for n in 19..10001 {
        if n % 2 == 1 {
            test_n(n);
        }
    }
}

#[test]
pub fn test_xy_loops() {
    for x in 10000..20000 {
        test_xy(x, x + 20);
    }
}

#[test]
pub fn test_multiple() {
    let mut mont = Mont::new(17);
    mont.precompute();

    let x = 1;
    let y = 2;
    let z = 3;
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    let z2 = mont.to_mont(z);

    let res = mont.multi(mont.multi(x2, y2), z2);

    let res2 = mont.reduce(res as u64);
    assert_eq!(x * y * z % mont.n, res2);
}

#[test]
pub fn test_pow() {
    let mut mont = Mont::new(17);
    mont.precompute();
    for base in 2..5 {
        for n in 5..10 {
            let x = mont.to_mont(base);

            let v = mont.pow(x, n);
            let v = mont.reduce(v as u64);

            let expected = base.pow(n) % mont.n;
            assert_eq!(expected, v);
        }
    }
}

pub fn main() {}
