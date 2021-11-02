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

#[test]
pub fn test() {
    let a = U512::from(100);
    let b: U256 = a.into();
    let b2: U512 = b.into();
    assert_eq!(a, b2);
}
