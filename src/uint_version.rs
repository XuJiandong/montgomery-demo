#![allow(dead_code)]
use std::{
    cell::RefCell,
    ops::{Add, Mul, Sub},
    rc::Rc,
};

use crate::signed_integer::SignedInteger;

use uint::construct_uint;

#[macro_export]
macro_rules! U {
    ($n:expr) => {
        U256::from($n)
    };
}

construct_uint! {
    pub struct U256(4);
}

construct_uint! {
    pub struct U512(8);
}

type I512 = SignedInteger<U512>;

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

#[derive(Clone)]
pub struct Mont {
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
            np1: U!(0),
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
        let t0: U512 = U!(t).into(); // low part of `t`, same as `% self.r`, avoid overflow
        let m: U512 = U!(t0 * U512::from(self.np1)).into(); // `% self.r`
        let u = (t + m * U512::from(self.n)) >> self.bits; // `/ self.r`
        if u >= U512::from(self.n) {
            U!(u - self.n)
        } else {
            U!(u)
        }
    }

    pub fn to_mont(&self, x: U256) -> U256 {
        assert!(self.init);
        let x2: U512 = x.into();
        let res = (x2 * self.r) % self.n;
        U!(res)
    }
    pub fn multi(&self, x: U256, y: U256) -> U256 {
        let xy = U512::from(x) * U512::from(y);
        self.reduce(xy)
    }

    pub fn pow(&self, x: U256, y: U256) -> U256 {
        let mut base = x;
        let one = U!(1);
        let mut res = U!(0);
        let mut first_time: bool = true;

        for index in 0..self.bits {
            // at most self.bits(256 here) multiplications
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
            base = self.multi(base, base); // at most self.bits(256 here) multiplications
        }
        res
    }
}

#[derive(Clone)]
pub struct MontForm {
    mont: Rc<RefCell<Mont>>,
    num: U256,
}

impl MontForm {
    // it's not possible to implement `From<U256>` for `MontForm` because it requires extra `mont`
    pub fn new(num: U256, mont: Mont) -> Self {
        MontForm {
            mont: Rc::new(RefCell::new(mont)),
            num,
        }
    }
    pub fn derive(&self, num: U256) -> Self {
        MontForm {
            mont: self.mont.clone(),
            num,
        }
    }
    pub fn pow(&self, pow: U256) -> Self {
        let num = self.mont.borrow().pow(self.num, pow);
        self.derive(num)
    }
}

impl From<MontForm> for U256 {
    fn from(m: MontForm) -> Self {
        m.mont.borrow().reduce(m.num.into())
    }
}

impl Add for MontForm {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let sum = (self.num + rhs.num) % self.mont.borrow().n;
        self.derive(sum)
    }
}

impl<'a> Add<&'a MontForm> for MontForm {
    type Output = Self;
    fn add(self, rhs: &'a Self) -> Self::Output {
        let sum = (self.num + rhs.num) % self.mont.borrow().n;
        self.derive(sum)
    }
}

impl<'a> Add<&'a MontForm> for &'a MontForm {
    type Output = MontForm;
    fn add(self, rhs: Self) -> Self::Output {
        let sum = (self.num + rhs.num) % self.mont.borrow().n;
        self.derive(sum)
    }
}

impl Sub for MontForm {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let sub = (self.num - rhs.num) % self.mont.borrow().n;
        self.derive(sub)
    }
}

impl<'a> Sub<&'a MontForm> for MontForm {
    type Output = Self;
    fn sub(self, rhs: &'a Self) -> Self::Output {
        let sub = (self.num - rhs.num) % self.mont.borrow().n;
        self.derive(sub)
    }
}

impl<'a> Sub<&'a MontForm> for &'a MontForm {
    type Output = MontForm;
    fn sub(self, rhs: Self) -> Self::Output {
        let sub = (self.num - rhs.num) % self.mont.borrow().n;
        self.derive(sub)
    }
}

impl Mul for MontForm {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul = self.mont.borrow().multi(self.num, rhs.num);
        self.derive(mul)
    }
}

impl<'a> Mul<&'a MontForm> for MontForm {
    type Output = Self;
    fn mul(self, rhs: &'a Self) -> Self::Output {
        let mul = self.mont.borrow().multi(self.num, rhs.num);
        self.derive(mul)
    }
}

impl<'a> Mul<&'a MontForm> for &'a MontForm {
    type Output = MontForm;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul = self.mont.borrow().multi(self.num, rhs.num);
        self.derive(mul)
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

    let x = U!(10);
    let x2 = mont.to_mont(x);
    let x3 = mont.reduce(x2.into());
    assert_eq!(x, x3);

    let x = U!(10);
    let y = U!(20);
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
    let mut mont = Mont::new(U!(1000001));
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
            test_n(U!(n));
        }
    }
}

#[test]
pub fn test_xy_loops() {
    for x in 10000..20000 {
        test_xy(U!(x), U!(x + 20));
    }
}

#[test]
pub fn test_multiple() {
    let mut mont = Mont::new(U!(17));
    mont.precompute();

    let x = U!(1);
    let y = U!(2);
    let z = U!(3);
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    let z2 = mont.to_mont(z);

    let res = mont.multi(mont.multi(x2, y2), z2);

    let res2 = mont.reduce(res.into());
    assert_eq!(x * y * z % mont.n, res2);
}

#[test]
pub fn test_pow() {
    let mut mont = Mont::new(U!(17));
    mont.precompute();
    for base in 2..5 {
        for n in 5..10 {
            let base = U!(base);
            let n = U!(n);
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
    let mut mont = Mont::new(U!(33));
    mont.precompute();
    let base = U!(2);
    let x = mont.to_mont(base);
    let v = mont.pow(x, U!(7));
    let v = mont.reduce(v.into());
    assert_eq!(v, U!(29));
}

#[test]
pub fn test_pow3() {
    let mut mont = Mont::new(U!(33));
    mont.precompute();
    let x = U!(2);
    let y = U!(2);
    let x2 = mont.to_mont(x);
    let y2 = mont.to_mont(y);
    let z2 = mont.multi(x2, y2);
    let z = mont.reduce(z2.into());
    assert_eq!(z, U!(4));

    let p2 = mont.pow(x2, U!(2));
    let p = mont.reduce(p2.into());
    assert_eq!(p, U!(4));
}

pub fn test_ops() {
    let mut mont = Mont::new(U!(33));
    mont.precompute();

    let x = &MontForm::new(U!(2), mont);
    let y = &x.derive(U!(3));

    let sum = x + y;
    assert_eq!(U!(5), sum.into());

    let sub = y - x;
    assert_eq!(U!(1), sub.into());

    let mul = x * y;
    assert_eq!(U!(6), mul.into());

    let pow = x.pow(U!(3));
    assert_eq!(U!(8), pow.into());
}
