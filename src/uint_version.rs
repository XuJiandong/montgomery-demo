#![allow(dead_code)]
extern crate alloc;

use crate::signed_integer::SignedInteger;
use alloc::rc::Rc;
use core::cell::RefCell;
use core::ops::{Add, Mul, Sub};

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

pub type I512 = SignedInteger<U512>;

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
    pub mont: Rc<RefCell<Mont>>,
    pub num: U256,
}

impl MontForm {
    // it's not possible to implement `From<U256>` for `MontForm` because it requires extra `mont`
    pub fn from_u256(num: U256, mont: Mont) -> Self {
        let num = mont.to_mont(num);
        MontForm {
            mont: Rc::new(RefCell::new(mont)),
            num,
        }
    }
    pub fn derive_from_u256(&self, num: U256) -> Self {
        let num = self.mont.borrow().to_mont(num);
        MontForm {
            mont: self.mont.clone(),
            num,
        }
    }

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
        let sum = self.num + rhs.num;
        let sum = if sum > self.mont.borrow().n {
            sum - self.mont.borrow().n
        } else {
            sum
        };
        self.derive(sum)
    }
}

impl<'a> Add<&'a MontForm> for MontForm {
    type Output = Self;
    fn add(self, rhs: &'a Self) -> Self::Output {
        let sum = self.num + rhs.num;
        let sum = if sum > self.mont.borrow().n {
            sum - self.mont.borrow().n
        } else {
            sum
        };
        self.derive(sum)
    }
}

impl<'a> Add<&'a MontForm> for &'a MontForm {
    type Output = MontForm;
    fn add(self, rhs: Self) -> Self::Output {
        let sum = self.num + rhs.num;
        let sum = if sum > self.mont.borrow().n {
            sum - self.mont.borrow().n
        } else {
            sum
        };
        self.derive(sum)
    }
}

impl Sub for MontForm {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let sub = if self.num >= rhs.num {
            self.num - rhs.num
        } else {
            self.num + self.mont.borrow().n - rhs.num
        };
        self.derive(sub)
    }
}

impl<'a> Sub<&'a MontForm> for MontForm {
    type Output = Self;
    fn sub(self, rhs: &'a Self) -> Self::Output {
        let sub = if self.num >= rhs.num {
            self.num - rhs.num
        } else {
            self.num + self.mont.borrow().n - rhs.num
        };
        self.derive(sub)
    }
}

impl<'a> Sub<&'a MontForm> for &'a MontForm {
    type Output = MontForm;
    fn sub(self, rhs: Self) -> Self::Output {
        let sub = if self.num >= rhs.num {
            self.num - rhs.num
        } else {
            self.num + self.mont.borrow().n - rhs.num
        };
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
