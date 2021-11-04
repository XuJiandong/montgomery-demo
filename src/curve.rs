use super::uint_version::{Mont, MontForm, U256, U512};
use lazy_static::lazy_static;
use std::ops::{Add, Mul};
use std::str::FromStr;

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
            x: U256::from(0),
            y: U256::from(0),
        }
    }
    pub fn from_str(x: &str, y: &str) -> Self {
        Self::new(U256::from_str(x).unwrap(), U256::from_str(y).unwrap())
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
        let two = MontForm::new(U256::from(2), mont);
        let three = two.derive(U256::from(3));
        let x_p1 = two.derive(self.x());
        let y_p1 = two.derive(self.y());
        let x_p2 = two.derive(rhs.x());
        let y_p2 = two.derive(rhs.y());

        let lam = if self == rhs {
            //  lam = 3 * x(P1) * x(P1) * pow(2 * y(P1), p - 2, p)
            let pow = (two * y_p1.clone()).pow(MONT.n - 2);
            three * x_p1.clone() * x_p1.clone() * pow
        } else {
            // lam = (y(P2) - y(P1)) * pow(x(P2) - x(P1), p - 2, p)
            let base = x_p2.clone() - x_p1.clone();
            (y_p2 - y_p1.clone()) * base.pow(MONT.n - 2)
        };
        // x3 = (lam * lam - x(P1) - x(P2)) % p
        // y3 = (lam * (x(P1) - x3) - y(P1)) % p
        let x3 = lam.clone() * lam.clone() - x_p1.clone() - x_p2.clone();
        let y3 = lam * (x_p1 - x3.clone()) - y_p1.clone();

        Point::new(U256::from(x3), U256::from(y3))
    }
}

impl Mul for Point {
    type Output = Self;
    fn mul(self, _rhs: Self) -> Self::Output {
        unimplemented!();
    }
}
