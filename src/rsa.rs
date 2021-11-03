use crate::Mont;

// https://www.cs.utexas.edu/~mitra/honors/soln.html

// Choose p = 3 and q = 11
// Compute n = p * q = 3 * 11 = 33
// Compute φ(n) = (p - 1) * (q - 1) = 2 * 10 = 20
// Choose e such that 1 < e < φ(n) and e and φ (n) are co-prime. Let e = 7
// Compute a value for d such that (d * e) % φ(n) = 1. One solution is d = 3 [(3 * 7) % 20 = 1]
// Public key is (e, n) => (7, 33)
// Private key is (d, n) => (3, 33)

// The encryption of m = 2 is c = 2**7 % 33 = 29 (m ** e mod n)
// The decryption of c = 29 is m = 29**3 % 33 = 2 ( c ** d mod n)

pub struct Rsa {
    pub p: u32,
    pub q: u32,
    pub n: u32,
    pub e: u32,
    pub d: u32,
}

impl Rsa {
    pub fn new() -> Self {
        Rsa {
            p: 3,
            q: 11,
            n: 33,
            e: 7,
            d: 3,
        }
    }

    pub fn encrypt(&self, plain_text: u32) -> u32 {
        assert!(plain_text < self.n);
        let mut mont = Mont::new(self.n);
        mont.precompute();
        let plain_text2 = mont.to_mont(plain_text);
        let cipher_text2 = mont.pow(plain_text2, self.e);
        let cipher_text = mont.reduce(cipher_text2 as u64);
        cipher_text
    }
    pub fn decrypt(&self, cipher_text: u32) -> u32 {
        assert!(cipher_text < self.n);
        let mut mont = Mont::new(self.n);
        mont.precompute();
        let cipher_text2 = mont.to_mont(cipher_text);
        let plain_text2 = mont.pow(cipher_text2, self.d);
        let plain_text = mont.reduce(plain_text2 as u64);
        plain_text
    }
}

#[test]
pub fn test_rsa() {
    let rsa = Rsa::new();
    let plain_text: u32 = 2;

    let cipher_text = rsa.encrypt(plain_text);
    println!("cipher_text = {}", cipher_text);
    let plain_text2 = rsa.decrypt(cipher_text);

    assert_eq!(plain_text, plain_text2);
}
