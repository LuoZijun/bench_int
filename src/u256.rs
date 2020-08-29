use uint::{construct_uint, uint_full_mul_reg};


construct_uint! {
    pub struct U256(4);
}

construct_uint! {
    pub struct U512(8);
}

// construct_uint! {
//     pub struct U1024(16);
// }

// construct_uint! {
//     pub struct U2048(32);
// }

// construct_uint! {
//     pub struct U4096(64);
// }

// construct_uint! {
//     pub struct U8192(128);
// }

impl U256 {
    #[inline(always)]
    pub fn full_mul(self, other: U256) -> U512 {
        U512(uint_full_mul_reg!(U256, 4, self, other))
    }
}