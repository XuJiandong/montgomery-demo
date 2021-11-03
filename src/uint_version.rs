#![allow(dead_code)]

use std::{
    cmp::Ordering,
    ops::{Add, Div, Mul, Rem, Sub},
};

use uint::construct_uint;

construct_uint! {
    pub struct U256(4);
}

construct_uint! {
    pub struct U512(8);
}

impl From<U256> for U512 {
    fn from(u: U256) -> Self {
        let U256(ref arr) = u;
        let mut val = [0; 8];
        val[0] = arr[0];
        val[1] = arr[1];
        val[2] = arr[2];
        val[3] = arr[3];
        Self(val)
    }
}

impl From<U512> for U256 {
    // use it with care
    fn from(u: U512) -> Self {
        let U512(ref arr) = u;
        let mut val = [0; 4];
        val[0] = arr[0];
        val[1] = arr[1];
        val[2] = arr[2];
        val[3] = arr[3];
        Self(val)
    }
}

#[derive(Copy, Clone, Eq, Debug)]
pub struct I512 {
    signed: bool,
    num: U512,
}

impl I512 {
    pub fn new(signed: bool, num: U512) -> Self {
        I512 { signed, num }
    }
}

impl From<U512> for I512 {
    fn from(n: U512) -> Self {
        Self::new(true, n)
    }
}

impl From<U256> for I512 {
    fn from(n: U256) -> Self {
        Self::new(true, n.into())
    }
}

impl From<u32> for I512 {
    fn from(n: u32) -> Self {
        Self::new(true, U512::from(n))
    }
}

impl Add for I512 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        match (self.signed, rhs.signed) {
            (true, true) => I512::new(true, self.num + rhs.num),
            (true, false) => {
                if self.num > rhs.num {
                    I512::new(true, self.num - rhs.num)
                } else {
                    I512::new(false, rhs.num - self.num)
                }
            }
            (false, true) => {
                if self.num > rhs.num {
                    I512::new(false, self.num - rhs.num)
                } else {
                    I512::new(true, rhs.num - self.num)
                }
            }
            (false, false) => I512::new(false, self.num + rhs.num),
        }
    }
}

impl Sub for I512 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        match (self.signed, rhs.signed) {
            (true, true) => {
                if self.num > rhs.num {
                    I512::new(true, self.num - rhs.num)
                } else {
                    I512::new(false, rhs.num - self.num)
                }
            }
            (true, false) => I512::new(true, self.num + rhs.num),
            (false, true) => I512::new(false, self.num + rhs.num),
            (false, false) => {
                if self.num > rhs.num {
                    I512::new(true, self.num - rhs.num)
                } else {
                    I512::new(false, rhs.num - self.num)
                }
            }
        }
    }
}

impl Mul for I512 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let signed = match (self.signed, rhs.signed) {
            (true, true) => true,
            (true, false) => false,
            (false, true) => false,
            (false, false) => false,
        };
        Self {
            signed,
            num: self.num * rhs.num,
        }
    }
}

impl Div for I512 {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let signed = match (self.signed, rhs.signed) {
            (true, true) => true,
            (true, false) => false,
            (false, true) => false,
            (false, false) => false,
        };
        Self {
            signed,
            num: self.num / rhs.num,
        }
    }
}

impl Rem for I512 {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self {
        assert!(rhs.num > 0.into());
        assert!(rhs.signed);
        if self.signed {
            I512::new(true, self.num % rhs.num)
        } else {
            let res = I512::new(false, self.num % rhs.num);
            res + I512::new(true, rhs.num)
        }
    }
}

impl Ord for I512 {
    fn cmp(&self, rhs: &Self) -> Ordering {
        match (self.signed, rhs.signed) {
            (true, true) => self.num.cmp(&rhs.num),
            (true, false) => Ordering::Greater,
            (false, true) => Ordering::Less,
            (false, false) => self.num.cmp(&rhs.num).reverse(),
        }
    }
}

impl PartialOrd for I512 {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}

impl PartialEq for I512 {
    fn eq(&self, rhs: &Self) -> bool {
        self.signed == rhs.signed && self.num == rhs.num
    }
}

// Extended Euclidean algorithm
// https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
pub fn egcd(a: I512, b: I512) -> (I512, I512, I512) {
    assert!(a < b);
    let zero = I512::from(0);
    let one = I512::from(1);
    if a == zero {
        (b, zero, one)
    } else {
        let (g, x, y) = egcd(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

struct Mont {
    pub np1: U256,
    pub r: U512,
    pub n: U256,
    pub bits: usize,
    pub init: bool,
}

impl Mont {
    pub fn new(n: U256) -> Self {
        let bits = 256;
        let base = U512::from(2);
        let r = base.pow(bits.into());
        Mont {
            r,
            n,
            np1: U256::from(0),
            bits,
            init: false,
        }
    }

    pub fn precompute(&mut self) {
        let one = I512::from(1);
        let zero = I512::from(0);

        let r: I512 = self.r.into();
        let n: I512 = self.n.into();
        let (gcd, np, rp) = egcd(n, r);
        assert_eq!(gcd, one);

        let rp1 = rp + self.n.into();
        assert!(rp1 >= zero);

        let np1 = I512::from(self.r) - np;
        assert!(np1 >= zero);
        self.np1 = np1.num.into(); // can be truncated

        let r = I512::from(self.r) * rp1;
        let n = I512::from(self.n) * np1;
        assert_eq!(n + one, r);

        self.init = true;
    }
    // m = T*Np1 mod R
    // U = (T + m * N) / R
    // The overall process delivers T · R−1 mod N without an expensive division operation!
    pub fn reduce(&self, t: U512) -> U256 {
        assert!(self.init);
        let t0: U512 = U256::from(t).into(); // low part of `t`, same as `% self.r`, avoid overflow
        let m: U512 = U256::from(t0 * U512::from(self.np1)).into(); // `% self.r`
        let u = (t + m * U512::from(self.n)) >> self.bits; // `/ self.r`
        if u >= U512::from(self.n) {
            U256::from(u - self.n)
        } else {
            U256::from(u)
        }
    }

    pub fn to_mont(&self, x: U256) -> U256 {
        assert!(self.init);
        let x2: U512 = x.into();
        let res = (x2 * self.r) % self.n;
        U256::from(res)
    }
    pub fn multi(&self, x: U256, y: U256) -> U256 {
        let xy = U512::from(x) * U512::from(y);
        self.reduce(xy)
    }

    pub fn pow(&self, x: U256, y: U256) -> U256 {
        let mut base = x;
        let one = U256::from(1);
        let mut res = U256::from(0);
        let mut first_time: bool = true;

        for index in 0..31 {
            // at most 32 multiplications
            if ((y >> index) & one) == one {
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
            base = self.multi(base, base); // at most 32 multiplications
        }
        res
    }
}

#[test]
pub fn test() {
    let a = U512::from(100);
    let b: U256 = a.into();
    let b2: U512 = b.into();
    assert_eq!(a, b2);
}

#[test]
pub fn test_i512() {
    let a = I512::from(2);
    let b = I512::from(10);
    let sum = a + b;
    assert_eq!(sum, I512::from(12));

    let sub = a - b;
    assert_eq!(sub, I512::new(false, U512::from(8)));

    let mul = a * b;
    assert_eq!(mul, I512::from(20));

    let div = b / a;
    assert_eq!(div, I512::from(5));

    let rem = b % a;
    assert_eq!(rem, I512::from(0));

    assert!(b > a);
    assert!(a < b);
    assert!(a != b);

    let x = I512::new(false, 11.into());
    let y = I512::from(7);
    let rem = x % y;
    assert_eq!(rem, I512::from(3));
}

#[test]
pub fn test_egcd() {
    let a = I512::from(11);
    let b = I512::from(17);
    let (gcd, x, y) = egcd(a, b);
    assert_eq!(gcd, I512::from(1));
    assert_eq!(a * x + b * y, I512::from(1));
}

fn test_n(n: U256) {
    let mut mont = Mont::new(n);
    mont.precompute();

    let x = U256::from(10);
    let x2 = mont.to_mont(x);
    let x3 = mont.reduce(x2.into());
    assert_eq!(x, x3);

    let x = U256::from(10);
    let y = U256::from(20);
    // into montgomery form
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    // do multiplication operation
    let xy2 = mont.multi(x2, y2);
    // into normal form
    let xy = mont.reduce(xy2.into());
    // the result should be same
    assert_eq!(xy, (x * y) % mont.n);
}

fn test_xy(x: U256, y: U256) {
    let mut mont = Mont::new(U256::from(1000001));
    mont.precompute();

    // into montgomery form
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    // do multiplication operation
    let xy2 = mont.multi(x2, y2);
    // into normal form
    let xy = mont.reduce(xy2.into());
    // the result should be same
    assert_eq!(xy, (x * y) % mont.n);
}

#[test]
pub fn test_n_loops() {
    for n in 19..10001 {
        if n % 2 == 1 {
            test_n(U256::from(n));
        }
    }
}

#[test]
pub fn test_xy_loops() {
    for x in 10000..20000 {
        test_xy(U256::from(x), U256::from(x + 20));
    }
}

#[test]
pub fn test_multiple() {
    let mut mont = Mont::new(U256::from(17));
    mont.precompute();

    let x = U256::from(1);
    let y = U256::from(2);
    let z = U256::from(3);
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    let z2 = mont.to_mont(z);

    let res = mont.multi(mont.multi(x2, y2), z2);

    let res2 = mont.reduce(res.into());
    assert_eq!(x * y * z % mont.n, res2);
}

#[test]
pub fn test_pow() {
    let mut mont = Mont::new(U256::from(17));
    mont.precompute();
    for base in 2..5 {
        for n in 5..10 {
            let base = U256::from(base);
            let n = U256::from(n);
            let x = mont.to_mont(base);

            let v = mont.pow(x, n);
            let v = mont.reduce(v.into());

            let expected = base.pow(n) % mont.n;
            assert_eq!(expected, v);
        }
    }
}

#[test]
pub fn test_pow2() {
    let mut mont = Mont::new(U256::from(33));
    mont.precompute();
    let base = U256::from(2);
    let x = mont.to_mont(base);
    let v = mont.pow(x, U256::from(7));
    let v = mont.reduce(v.into());
    assert_eq!(v, U256::from(29));
}

#[test]
pub fn test_pow3() {
    let mut mont = Mont::new(U256::from(33));
    mont.precompute();
    let x = U256::from(2);
    let y = U256::from(2);
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    let z2 = mont.multi(x2, y2);
    let z = mont.reduce(z2.into());
    assert_eq!(z, U256::from(4));

    let p2 = mont.pow(x2, U256::from(2));
    let p = mont.reduce(p2.into());
    assert_eq!(p, U256::from(4));
}
