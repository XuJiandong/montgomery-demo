// use std::println;
mod primitive_version_test {
    use rvv_algo_demo::primitive_version::*;
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
    pub fn test_xy_loops() {
        for x in 10000..20000 {
            test_xy(x, x + 20);
        }
    }

    pub fn test_n_loops() {
        for n in 19..10001 {
            if n % 2 == 1 {
                test_n(n);
            }
        }
    }

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
    pub fn test_addsub() {
        let mut mont = Mont::new(17);
        mont.precompute();

        let x = 15;
        let y = 7;
        let x2 = mont.to_mont(x);
        let y2 = mont.to_mont(y);
        let sum2 = x2 + y2;
        let sum = mont.reduce(sum2 as u64);
        assert_eq!(sum, (x + y) % mont.n);

        let sub2 = x2 - y2;
        let sub = mont.reduce(sub2 as u64);
        assert_eq!(sub, (x - y) % mont.n);
    }
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
    pub fn test_pow2() {
        let mut mont = Mont::new(33);
        mont.precompute();
        let base: u32 = 2;
        let x = mont.to_mont(base);
        // println!("x = {}", x);
        let v = mont.pow(x, 7);
        let v = mont.reduce(v as u64);
        assert_eq!(v, 29);
    }

    pub fn test_pow3() {
        let mut mont = Mont::new(33);
        mont.precompute();
        let x: u32 = 2;
        let y: u32 = 2;
        let x2 = mont.to_mont(x);
        let y2 = mont.to_mont(y);
        let z2 = mont.multi(x2, y2);
        let z = mont.reduce(z2 as u64);
        assert_eq!(z, 4);

        let p2 = mont.pow(x2, 2);
        let p = mont.reduce(p2 as u64);
        assert_eq!(p, 4);
    }
}

mod rsa_test {
    use rvv_algo_demo::rsa::*;
    pub fn test_rsa() {
        let rsa = Rsa::new();
        let plain_text: u32 = 2;

        let cipher_text = rsa.encrypt(plain_text);
        // println!("cipher_text = {}", cipher_text);
        let plain_text2 = rsa.decrypt(cipher_text);

        assert_eq!(plain_text, plain_text2);
    }
}

mod uint_version_test {
    use rvv_algo_demo::uint_version::*;
    use rvv_algo_demo::U;

    pub fn test() {
        let a = U512::from(100);
        let b: U256 = a.into();
        let b2: U512 = b.into();
        assert_eq!(a, b2);
    }

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

    pub fn test_n_loops() {
        for n in 19..10001 {
            if n % 2 == 1 {
                test_n(U!(n));
            }
        }
    }

    pub fn test_xy_loops() {
        for x in 10000..20000 {
            test_xy(U!(x), U!(x + 20));
        }
    }

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

    pub fn test_pow2() {
        let mut mont = Mont::new(U!(33));
        mont.precompute();
        let base = U!(2);
        let x = mont.to_mont(base);
        let v = mont.pow(x, U!(7));
        let v = mont.reduce(v.into());
        assert_eq!(v, U!(29));
    }

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

        let x = &MontForm::from_u256(U!(2), mont);
        let y = &x.derive_from_u256(U!(3));

        let sum = x + y;
        let res: U256 = sum.into();
        assert_eq!(U!(5), res);

        let sub = y - x;
        assert_eq!(U!(1), sub.into());

        let sub2 = x - y;
        assert_eq!(U!(32), sub2.into());

        let mul = x * y;
        assert_eq!(U!(6), mul.into());

        let pow = x.pow(U!(3));
        assert_eq!(U!(8), pow.into());
    }
    /*
    macro_rules! field_impl {
        ($name:ident, $modulus:expr, $rsquared:expr, $rcubed:expr, $one:expr, $inv:expr) => {

    field_impl!(
        Fr,
        [
            0x2833e84879b9709143e1f593f0000001,
            0x30644e72e131a029b85045b68181585d
        ],
        [
            0x53fe3ab1e35c59e31bb8e645ae216da7,
            0x0216d0b17f4e44a58c49833d53bb8085
        ],
        [
            0x2a489cbe1cfbb6b85e94d8e1b4bf0040,
            0x0cf8594b7fcc657c893cc664a19fcfed
        ],
        [
            0x36fc76959f60cd29ac96341c4ffffffb,
            0x0e0a77c19a07df2f666ea36f7879462e
        ],
        0x6586864b4c6911b3c2e1f593efffffff
    );
     */
    fn from_u128pair(n: &[u128; 2]) -> U256 {
        let mut buf = [0u8; 32];
        buf[0..16].copy_from_slice(&n[0].to_le_bytes());
        buf[16..32].copy_from_slice(&n[1].to_le_bytes());
        U256::from_little_endian(&buf)
    }

    #[test]
    pub fn test_number() {
        let n = from_u128pair(&[
            0x2833e84879b9709143e1f593f0000001u128,
            0x30644e72e131a029b85045b68181585du128,
        ]);
        let mut mont = Mont::new(n);
        mont.precompute();
        // println!("np1 = 0x{:x}", mont.np1);
        let np1 = 0x6586864b4c6911b3c2e1f593efffffffu128;
        assert_eq!(np1, mont.np1.low_u128());

        let r: U512 = mont.r.into();
        let r_square = (r % mont.n) * (r % mont.n);
        let r_square2 = r_square % mont.n;
        // println!("r_square = 0x{:x}", r_square2);
        let original_r_square = from_u128pair(&[
            0x53fe3ab1e35c59e31bb8e645ae216da7,
            0x0216d0b17f4e44a58c49833d53bb8085,
        ]);
        assert_eq!(original_r_square, r_square2.into());

        let r_cubed = r_square2 * (r % mont.n) % mont.n;
        let original_r_cubed = from_u128pair(&[
            0x2a489cbe1cfbb6b85e94d8e1b4bf0040,
            0x0cf8594b7fcc657c893cc664a19fcfed,
        ]);
        assert_eq!(original_r_cubed, r_cubed.into());
    }
    #[test]
    pub fn test_number2() {
        // inv = 0x9ede7d651eca6ac987d20782e4866389
        // modulo = 0x97816a916871ca8d3c208c16d87cfd47, 0x30644e72e131a029b85045b68181585d
        let n = from_u128pair(&[
            0x97816a916871ca8d3c208c16d87cfd47u128, 0x30644e72e131a029b85045b68181585du128
        ]);
        let mut mont = Mont::new(n);
        mont.precompute();
        // println!("np1 = 0x{:x}", mont.np1);
        let np1 = 0x9ede7d651eca6ac987d20782e4866389;
        assert_eq!(np1, mont.np1.low_u128());
        println!("mont.np1 = 0x{:x},0x{:x}", mont.np1.low_u128(), (mont.np1 >> 128).low_u128());
    }
    #[test]
    pub fn test_number3() {
        // inv = 0x6586864b4c6911b3c2e1f593efffffff
        // modulo = 0x2833e84879b9709143e1f593f0000001, 0x30644e72e131a029b85045b68181585d
        let n = from_u128pair(&[
            0x2833e84879b9709143e1f593f0000001, 0x30644e72e131a029b85045b68181585d
        ]);
        let mut mont = Mont::new(n);
        mont.precompute();
        // println!("np1 = 0x{:x}", mont.np1);
        let np1 = 0x6586864b4c6911b3c2e1f593efffffff;
        assert_eq!(np1, mont.np1.low_u128());
        println!("mont.np1 = 0x{:x},0x{:x}", mont.np1.low_u128(), (mont.np1 >> 128).low_u128());
    }
}

#[test]
pub fn test() {
    main()
}

pub fn main() {
    // #[test] is not supported in ckb-std
    primitive_version_test::test_multiple();
    primitive_version_test::test_addsub();
    primitive_version_test::test_pow();
    primitive_version_test::test_pow2();
    primitive_version_test::test_pow3();
    primitive_version_test::test_n_loops();
    primitive_version_test::test_xy_loops();

    rsa_test::test_rsa();

    uint_version_test::test();
    uint_version_test::test_i512();
    uint_version_test::test_egcd();
    uint_version_test::test_n_loops();
    uint_version_test::test_xy_loops();
    uint_version_test::test_multiple();
    uint_version_test::test_pow();
    uint_version_test::test_pow2();
    uint_version_test::test_pow3();
    uint_version_test::test_ops();
}
