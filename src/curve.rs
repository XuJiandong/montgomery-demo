use super::uint_version::{Mont, MontForm, U256};
use super::U;
use lazy_static::lazy_static;

use core::ops::{Add, Mul};
use core::str::FromStr;

lazy_static! {
    static ref MONT: Mont = {
        // https://github.com/XuJiandong/taproot-playground/blob/d5f48c5d5797395ce3f2e209cca29b01695a3d48/py/reference.py#L14
        let n = U256::from_str("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F").unwrap();
        let mut mont = Mont::new(n);
        mont.precompute();
        mont
    };
}

#[derive(Clone, Copy, Eq, Debug)]
pub struct Point {
    infinite: bool,
    x: U256,
    y: U256,
}

impl PartialEq for Point {
    fn eq(&self, rhs: &Self) -> bool {
        self.infinite == rhs.infinite && self.x == rhs.y && self.y == rhs.y
    }
}

impl Point {
    pub fn new(x: U256, y: U256) -> Self {
        Point {
            infinite: false,
            x,
            y,
        }
    }
    pub fn new_infinite() -> Self {
        Point {
            infinite: true,
            x: U!(0),
            y: U!(0),
        }
    }
    pub fn x(&self) -> U256 {
        assert!(!self.infinite);
        self.x
    }
    pub fn y(&self) -> U256 {
        assert!(!self.infinite);
        self.y
    }
}

impl Add for Point {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        if self.infinite {
            return rhs;
        }
        if rhs.infinite {
            return self;
        }
        if self.x() == rhs.x() && self.y() != rhs.y() {
            return Point::new_infinite();
        }

        let mont = MONT.clone();
        let two = MontForm::from_u256(U!(2), mont);
        let three = &two.derive_from_u256(U!(3));
        let x_p1 = &two.derive_from_u256(self.x());
        let y_p1 = &two.derive_from_u256(self.y());
        let x_p2 = &two.derive_from_u256(rhs.x());
        let y_p2 = &two.derive_from_u256(rhs.y());

        let lam = if self == rhs {
            //  lam = 3 * x(P1) * x(P1) * pow(2 * y(P1), p - 2, p)
            let pow = (two * y_p1).pow(MONT.n - 2);
            three * x_p1 * x_p1 * pow
        } else {
            // lam = (y(P2) - y(P1)) * pow(x(P2) - x(P1), p - 2, p)
            let base = x_p2 - x_p1;
            (y_p2 - y_p1) * base.pow(MONT.n - 2)
        };
        // x3 = (lam * lam - x(P1) - x(P2)) % p
        // y3 = (lam * (x(P1) - x3) - y(P1)) % p
        let x3 = &lam * &lam - x_p1 - x_p2;
        let y3 = lam * (x_p1 - &x3) - y_p1;

        Point::new(U!(x3), U!(y3))
    }
}

// https://github.com/XuJiandong/taproot-playground/blob/d5f48c5d5797395ce3f2e209cca29b01695a3d48/py/reference.py#L54
impl Mul<U256> for Point {
    type Output = Self;
    fn mul(self, n: U256) -> Self::Output {
        let mut result: Self = self.clone();
        let mut first_time = true;
        let mut p = self.clone();
        for index in 0..256 {
            if ((n >> index) & U!(1)) == U!(1) {
                if first_time {
                    result = p;
                    first_time = false;
                } else {
                    result = result + p;
                }
            }
            p = p + p;
        }
        result
    }
}
