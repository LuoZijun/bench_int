
/// Calculate a - b - borrow, returning the result and modifying
/// the borrow value.
#[inline(always)]
pub fn sbb(a: u64, b: u64, borrow: &mut u64) -> u64 {
    let tmp = (1u128 << 64) + u128::from(a) - u128::from(b) - u128::from(*borrow);

    *borrow = if tmp >> 64 == 0 { 1 } else { 0 };

    tmp as u64
}

/// Calculate a + b + carry, returning the sum and modifying the
/// carry value.
#[inline(always)]
pub fn adc(a: u64, b: u64, carry: &mut u64) -> u64 {
    let tmp = u128::from(a) + u128::from(b) + u128::from(*carry);

    *carry = (tmp >> 64) as u64;

    tmp as u64
}

/// Calculate a + (b * c) + carry, returning the least significant digit
/// and setting carry to the most significant digit.
#[inline(always)]
pub fn mac_with_carry(a: u64, b: u64, c: u64, carry: &mut u64) -> u64 {
    let tmp = (u128::from(a)) + u128::from(b) * u128::from(c) + u128::from(*carry);

    *carry = (tmp >> 64) as u64;

    tmp as u64
}

#[cfg(target_arch = "x86_64")]
#[link(name = "ff-derive-crypto", kind = "static")]
extern "C" {
    fn mod_mul_4w(a: &[u64; 4], b: &[u64; 4], res: &mut [u64; 4]);
}

#[cfg(target_arch = "x86_64")]
#[inline]
pub fn mod_mul_4w_assign(a: &mut [u64; 4], b: &[u64; 4]) {
    let mut res = [0; 4];
    unsafe {
        mod_mul_4w(&*a, b, &mut res);
    }
    let _ = std::mem::replace(a, res);
}

// use crate::BaseFromRO;
// use fff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};
// use sha2ni::digest::generic_array::{typenum::U48, GenericArray};
// use std::io::{Cursor, Read};


// #[PrimeFieldModulus = "52435875175126190479447740508185965837690552500527637822603658699938581184513"]
// #[PrimeFieldGenerator = "7"]
#[derive(Copy, Clone, PartialEq, Eq, Default)]
pub struct Fr(FrRepr);

/// This is the modulus m of the prime field
pub const MODULUS: FrRepr = FrRepr([
    18446744069414584321u64,
    6034159408538082302u64,
    3691218898639771653u64,
    8353516859464449352u64,
]);

/// The number of bits needed to represent the modulus.
pub const MODULUS_BITS: u32 = 255u32;

/// The number of bits that must be shaved from the beginning of
/// the representation when randomly sampling.
pub const REPR_SHAVE_BITS: u32 = 1u32;

/// 2^{limbs*64} mod m
pub const R: FrRepr = FrRepr([
    8589934590u64,
    6378425256633387010u64,
    11064306276430008309u64,
    1739710354780652911u64,
]);

/// 2^{limbs*64*2} mod m
pub const R2: FrRepr = FrRepr([
    14526898881837571181u64,
    3129137299524312099u64,
    419701826671360399u64,
    524908885293268753u64,
]);

/// -(m^{-1} mod m) mod m
pub const INV: u64 = 18446744069414584319u64;

/// Multiplicative generator of `MODULUS` - 1 order, also quadratic
/// nonresidue.
pub const GENERATOR: FrRepr = FrRepr([
    64424509425u64,
    1721329240476523535u64,
    18418692815241631664u64,
    3824455624000121028u64,
]);

/// 2^s * t = MODULUS - 1 with t odd
pub const S: u32 = 32u32;

/// 2^s root of unity computed by GENERATOR^t
pub const ROOT_OF_UNITY: FrRepr = FrRepr([
    13381757501831005802u64,
    6564924994866501612u64,
    789602057691799140u64,
    6625830629041353339u64,
]);


#[derive(Copy, Clone, PartialEq, Eq, Default)]
pub struct FrRepr(pub [u64; 4usize]);

#[allow(dead_code)]
impl FrRepr {
    #[inline(always)]
    fn is_odd(&self) -> bool {
        self.0[0] & 1 == 1
    }

    #[inline(always)]
    fn is_even(&self) -> bool {
        !self.is_odd()
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0.iter().all(|&e| e == 0)
    }

    #[inline(always)]
    fn shr(&mut self, mut n: u32) {
        if n as usize >= 64 * 4usize {
            *self = Self::from(0);
            return;
        }
        while n >= 64 {
            let mut t = 0;
            for i in self.0.iter_mut().rev() {
                std::mem::swap(&mut t, i);
            }
            n -= 64;
        }
        if n > 0 {
            let mut t = 0;
            for i in self.0.iter_mut().rev() {
                let t2 = *i << (64 - n);
                *i >>= n;
                *i |= t;
                t = t2;
            }
        }
    }

    #[inline(always)]
    fn div2(&mut self) {
        let mut t = 0;
        for i in self.0.iter_mut().rev() {
            let t2 = *i << 63;
            *i >>= 1;
            *i |= t;
            t = t2;
        }
    }

    #[inline(always)]
    fn mul2(&mut self) {
        let mut last = 0;
        for i in &mut self.0 {
            let tmp = *i >> 63;
            *i <<= 1;
            *i |= last;
            last = tmp;
        }
    }

    #[inline(always)]
    fn shl(&mut self, mut n: u32) {
        if n as usize >= 64 * 4usize {
            *self = Self::from(0);
            return;
        }
        while n >= 64 {
            let mut t = 0;
            for i in &mut self.0 {
                std::mem::swap(&mut t, i);
            }
            n -= 64;
        }
        if n > 0 {
            let mut t = 0;
            for i in &mut self.0 {
                let t2 = *i >> (64 - n);
                *i <<= n;
                *i |= t;
                t = t2;
            }
        }
    }

    #[inline(always)]
    fn num_bits(&self) -> u32 {
        let mut ret = (4usize as u32) * 64;
        for i in self.0.iter().rev() {
            let leading = i.leading_zeros();
            ret -= leading;
            if leading != 64 {
                break;
            }
        }
        ret
    }

    #[inline(always)]
    pub fn add_nocarry(&mut self, other: &FrRepr) {
        let mut carry = 0;
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a = adc(*a, *b, &mut carry);
        }
    }

    #[inline(always)]
    pub fn sub_noborrow(&mut self, other: &FrRepr) {
        let mut borrow = 0;
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a = sbb(*a, *b, &mut borrow);
        }
    }
}

impl std::fmt::Debug for FrRepr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "0x")?;
        for i in self.0.iter().rev() {
            write!(f, "{:016x}", *i)?;
        }
        Ok(())
    }
}

impl std::fmt::Display for FrRepr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "0x")?;
        for i in self.0.iter().rev() {
            write!(f, "{:016x}", *i)?;
        }
        Ok(())
    }
}

impl AsRef<[u64]> for FrRepr {
    #[inline(always)]
    fn as_ref(&self) -> &[u64] {
        &self.0
    }
}

impl AsMut<[u64]> for FrRepr {
    #[inline(always)]
    fn as_mut(&mut self) -> &mut [u64] {
        &mut self.0
    }
}
impl From<u64> for FrRepr {
    #[inline(always)]
    fn from(val: u64) -> FrRepr {
        let mut repr = Self::default();
        repr.0[0] = val;
        repr
    }
}
impl Ord for FrRepr {
    #[inline(always)]
    fn cmp(&self, other: &FrRepr) -> std::cmp::Ordering {
        for (a, b) in self.0.iter().rev().zip(other.0.iter().rev()) {
            if a < b {
                return std::cmp::Ordering::Less;
            } else if a > b {
                return std::cmp::Ordering::Greater;
            }
        }
        std::cmp::Ordering::Equal
    }
}
impl PartialOrd for FrRepr {
    #[inline(always)]
    fn partial_cmp(&self, other: &FrRepr) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[allow(dead_code)]
impl Fr {
    pub const NUM_BITS: u32 = MODULUS_BITS;
    pub const CAPACITY: u32 = Self::NUM_BITS - 1;
    pub const S: u32        = S;

    // type Repr = FrRepr;
    pub fn from_repr(r: FrRepr) -> Result<Fr, ()> {
        let mut r = Fr(r);
        if r.is_valid() {
            r.mul_assign(&Fr(R2));
            Ok(r)
        } else {
            // Err(PrimeFieldDecodingError::NotInField(format!("{}", r.0)))
            println!("PrimeFieldDecodingError::NotInField: {:?}", format!("{}", r.0));
            Err(())
        }
    }

    pub fn into_repr(&self) -> FrRepr {
        let mut r = *self;
        r.mont_reduce(
            (self.0).0[0usize],
            (self.0).0[1usize],
            (self.0).0[2usize],
            (self.0).0[3usize],
            0,
            0,
            0,
            0,
        );
        r.0
    }

    fn char() -> FrRepr {
        MODULUS
    }

    fn multiplicative_generator() -> Self {
        Fr(GENERATOR)
    }
    
    fn root_of_unity() -> Self {
        Fr(ROOT_OF_UNITY)
    }
    
    // /// Computes a uniformly random element using rejection sampling.
    // fn random<R: ::rand_core::RngCore>(rng: &mut R) -> Self {
    //     loop {
    //         let mut tmp = {
    //             let mut repr = [0u64; 4usize];
    //             for i in 0..4usize {
    //                 repr[i] = rng.next_u64();
    //             }
    //             Fr(FrRepr(repr))
    //         };
    //         tmp.0.as_mut()[3usize] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;
    //         if tmp.is_valid() {
    //             return tmp;
    //         }
    //     }
    // }

    #[inline]
    pub fn zero() -> Self {
        Fr(FrRepr::from(0))
    }

    #[inline]
    pub fn one() -> Self {
        Fr(R)
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    pub fn add_assign_amd64(&mut self, other: &Fr) {
        use std::arch::x86_64::*;
        
        unsafe {
            let mut carry = _addcarry_u64(0, (self.0).0[0], (other.0).0[0], &mut (self.0).0[0]);
            carry = _addcarry_u64(carry, (self.0).0[1], (other.0).0[1], &mut (self.0).0[1]);
            carry = _addcarry_u64(carry, (self.0).0[2], (other.0).0[2], &mut (self.0).0[2]);
            _addcarry_u64(carry, (self.0).0[3], (other.0).0[3], &mut (self.0).0[3]);

            #[allow(deprecated)]
            let mut s_sub: [u64; 4] = std::mem::uninitialized();

            carry = _subborrow_u64(0, (self.0).0[0], MODULUS.0[0], &mut s_sub[0]);
            carry = _subborrow_u64(carry, (self.0).0[1], MODULUS.0[1], &mut s_sub[1]);
            carry = _subborrow_u64(carry, (self.0).0[2], MODULUS.0[2], &mut s_sub[2]);
            carry = _subborrow_u64(carry, (self.0).0[3], MODULUS.0[3], &mut s_sub[3]);
            if carry == 0 {
                (self.0).0[0] = s_sub[0];
                (self.0).0[1] = s_sub[1];
                (self.0).0[2] = s_sub[2];
                (self.0).0[3] = s_sub[3];
            }
        }
    }

    #[inline]
    pub fn add_assign_soft(&mut self, other: &Fr) {
        self.0.add_nocarry(&other.0);
        self.reduce();
    }

    #[inline]
    pub fn add_assign(&mut self, other: &Fr) {
        #[cfg(target_arch = "x86_64")]
        self.add_assign_amd64(other);

        #[cfg(not(target_arch = "x86_64"))]
        self.add_assign_soft(other);
    }

    #[inline]
    fn double(&mut self) {
        self.0.mul2();
        self.reduce();
    }

    #[inline]
    pub fn sub_assign(&mut self, other: &Fr) {
        if other.0 > self.0 {
            self.0.add_nocarry(&MODULUS);
        }
        self.0.sub_noborrow(&other.0);
    }

    #[inline]
    fn negate(&mut self) {
        if !self.is_zero() {
            let mut tmp = MODULUS;
            tmp.sub_noborrow(&self.0);
            self.0 = tmp;
        }
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let one = FrRepr::from(1);
        let mut u = self.0;
        let mut v = MODULUS;
        let mut b = Fr(R2);
        let mut c = Self::zero();
        while u != one && v != one {
            while u.is_even() {
                u.div2();
                if b.0.is_even() {
                    b.0.div2();
                } else {
                    b.0.add_nocarry(&MODULUS);
                    b.0.div2();
                }
            }
            while v.is_even() {
                v.div2();
                if c.0.is_even() {
                    c.0.div2();
                } else {
                    c.0.add_nocarry(&MODULUS);
                    c.0.div2();
                }
            }
            if v < u {
                u.sub_noborrow(&v);
                b.sub_assign(&c);
            } else {
                v.sub_noborrow(&u);
                c.sub_assign(&b);
            }
        }
        if u == one {
            Some(b)
        } else {
            Some(c)
        }
    }

    // #[inline(always)]
    // fn frobenius_map(&mut self, _: usize) { }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    pub fn mul_assign_amd64_adx_bmi2(&mut self, other: &Fr) {
        mod_mul_4w_assign(&mut (self.0).0, &(other.0).0);
    }

    #[inline]
    pub fn mul_assign_soft(&mut self, other: &Fr) {
        let mut carry = 0;
        let r0 = mac_with_carry(0, (self.0).0[0usize], (other.0).0[0usize], &mut carry);
        let r1 = mac_with_carry(0, (self.0).0[0usize], (other.0).0[1usize], &mut carry);
        let r2 = mac_with_carry(0, (self.0).0[0usize], (other.0).0[2usize], &mut carry);
        let r3 = mac_with_carry(0, (self.0).0[0usize], (other.0).0[3usize], &mut carry);
        let r4 = carry;

        let mut carry = 0;
        let r1 = mac_with_carry(r1, (self.0).0[1usize], (other.0).0[0usize], &mut carry);
        let r2 = mac_with_carry(r2, (self.0).0[1usize], (other.0).0[1usize], &mut carry);
        let r3 = mac_with_carry(r3, (self.0).0[1usize], (other.0).0[2usize], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[1usize], (other.0).0[3usize], &mut carry);
        let r5 = carry;

        let mut carry = 0;
        let r2 = mac_with_carry(r2, (self.0).0[2usize], (other.0).0[0usize], &mut carry);
        let r3 = mac_with_carry(r3, (self.0).0[2usize], (other.0).0[1usize], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[2usize], (other.0).0[2usize], &mut carry);
        let r5 = mac_with_carry(r5, (self.0).0[2usize], (other.0).0[3usize], &mut carry);
        let r6 = carry;

        let mut carry = 0;
        let r3 = mac_with_carry(r3, (self.0).0[3usize], (other.0).0[0usize], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[3usize], (other.0).0[1usize], &mut carry);
        let r5 = mac_with_carry(r5, (self.0).0[3usize], (other.0).0[2usize], &mut carry);
        let r6 = mac_with_carry(r6, (self.0).0[3usize], (other.0).0[3usize], &mut carry);
        let r7 = carry;
        self.mont_reduce(r0, r1, r2, r3, r4, r5, r6, r7);
    }

    #[inline]
    pub fn mul_assign(&mut self, other: &Fr) {
        if is_x86_feature_detected!("adx") && is_x86_feature_detected!("bmi2") {
            self.mul_assign_amd64_adx_bmi2(other);
        } else {
            self.mul_assign_soft(other);
        }
    }

    #[inline]
    pub fn square(&mut self) {
        let mut carry = 0;
        let r1 = mac_with_carry(0, (self.0).0[0usize], (self.0).0[1usize], &mut carry);
        let r2 = mac_with_carry(0, (self.0).0[0usize], (self.0).0[2usize], &mut carry);
        let r3 = mac_with_carry(0, (self.0).0[0usize], (self.0).0[3usize], &mut carry);
        let r4 = carry;
        let mut carry = 0;
        let r3 = mac_with_carry(r3, (self.0).0[1usize], (self.0).0[2usize], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[1usize], (self.0).0[3usize], &mut carry);
        let r5 = carry;
        let mut carry = 0;
        let r5 = mac_with_carry(r5, (self.0).0[2usize], (self.0).0[3usize], &mut carry);
        let r6 = carry;
        let r7 = r6 >> 63;
        let r6 = (r6 << 1) | (r5 >> 63);
        let r5 = (r5 << 1) | (r4 >> 63);
        let r4 = (r4 << 1) | (r3 >> 63);
        let r3 = (r3 << 1) | (r2 >> 63);
        let r2 = (r2 << 1) | (r1 >> 63);
        let r1 = r1 << 1;
        let mut carry = 0;
        let r0 = mac_with_carry(0, (self.0).0[0usize], (self.0).0[0usize], &mut carry);
        let r1 = adc(r1, 0, &mut carry);
        let r2 = mac_with_carry(r2, (self.0).0[1usize], (self.0).0[1usize], &mut carry);
        let r3 = adc(r3, 0, &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[2usize], (self.0).0[2usize], &mut carry);
        let r5 = adc(r5, 0, &mut carry);
        let r6 = mac_with_carry(r6, (self.0).0[3usize], (self.0).0[3usize], &mut carry);
        let r7 = adc(r7, 0, &mut carry);
        self.mont_reduce(r0, r1, r2, r3, r4, r5, r6, r7);
    }

    /// Determines if the element is really in the field. This is only used
    /// internally.
    #[inline(always)]
    pub fn is_valid(&self) -> bool {
        self.0 < MODULUS
    }
    /// Subtracts the modulus from this element if this element is not in the
    /// field. Only used interally.
    #[inline(always)]
    pub fn reduce(&mut self) {
        if !self.is_valid() {
            self.0.sub_noborrow(&MODULUS);
        }
    }

    #[inline(always)]
    pub fn mont_reduce(&mut self, r0: u64, mut r1: u64, mut r2: u64, mut r3: u64, mut r4: u64, mut r5: u64, mut r6: u64, mut r7: u64) {
        let k = r0.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r0, k, MODULUS.0[0], &mut carry);
        r1 = mac_with_carry(r1, k, MODULUS.0[1usize], &mut carry);
        r2 = mac_with_carry(r2, k, MODULUS.0[2usize], &mut carry);
        r3 = mac_with_carry(r3, k, MODULUS.0[3usize], &mut carry);
        r4 = adc(r4, 0, &mut carry);
        let carry2 = carry;
        let k = r1.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r1, k, MODULUS.0[0], &mut carry);
        r2 = mac_with_carry(r2, k, MODULUS.0[1usize], &mut carry);
        r3 = mac_with_carry(r3, k, MODULUS.0[2usize], &mut carry);
        r4 = mac_with_carry(r4, k, MODULUS.0[3usize], &mut carry);
        r5 = adc(r5, carry2, &mut carry);
        let carry2 = carry;
        let k = r2.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r2, k, MODULUS.0[0], &mut carry);
        r3 = mac_with_carry(r3, k, MODULUS.0[1usize], &mut carry);
        r4 = mac_with_carry(r4, k, MODULUS.0[2usize], &mut carry);
        r5 = mac_with_carry(r5, k, MODULUS.0[3usize], &mut carry);
        r6 = adc(r6, carry2, &mut carry);
        let carry2 = carry;
        let k = r3.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r3, k, MODULUS.0[0], &mut carry);
        r4 = mac_with_carry(r4, k, MODULUS.0[1usize], &mut carry);
        r5 = mac_with_carry(r5, k, MODULUS.0[2usize], &mut carry);
        r6 = mac_with_carry(r6, k, MODULUS.0[3usize], &mut carry);
        r7 = adc(r7, carry2, &mut carry);
        (self.0).0[0usize] = r4;
        (self.0).0[1usize] = r5;
        (self.0).0[2usize] = r6;
        (self.0).0[3usize] = r7;
        self.reduce();
    }

    // fn legendre(&self) -> ::fff::LegendreSymbol {
    //     let s = self.pow([
    //         9223372034707292160u64,
    //         12240451741123816959u64,
    //         1845609449319885826u64,
    //         4176758429732224676u64,
    //     ]);
    //     if s == Self::zero() {
    //         ::fff::LegendreSymbol::Zero
    //     } else if s == Self::one() {
    //         ::fff::LegendreSymbol::QuadraticResidue
    //     } else {
    //         ::fff::LegendreSymbol::QuadraticNonResidue
    //     }
    // }
    // fn sqrt(&self) -> Option<Self> {
    //     match self.legendre() {
    //         ::fff::LegendreSymbol::Zero => Some(*self),
    //         ::fff::LegendreSymbol::QuadraticNonResidue => None,
    //         ::fff::LegendreSymbol::QuadraticResidue => {
    //             let mut c = Fr(ROOT_OF_UNITY);
    //             let mut r = self.pow([
    //                 9223141137265459200u64,
    //                 347036667491570177u64,
    //                 10722717374829358084u64,
    //                 972477353u64,
    //             ]);
    //             let mut t = self.pow([
    //                 18446282274530918399u64,
    //                 694073334983140354u64,
    //                 2998690675949164552u64,
    //                 1944954707u64,
    //             ]);
    //             let mut m = S;
    //             while t != Self::one() {
    //                 let mut i = 1;
    //                 {
    //                     let mut t2i = t;
    //                     t2i.square();
    //                     loop {
    //                         if t2i == Self::one() {
    //                             break;
    //                         }
    //                         t2i.square();
    //                         i += 1;
    //                     }
    //                 }
    //                 for _ in 0..(m - i - 1) {
    //                     c.square();
    //                 }
    //                 r.mul_assign(&c);
    //                 c.square();
    //                 t.mul_assign(&c);
    //                 m = i;
    //             }
    //             Some(r)
    //         }
    //     }
    // }
}


/// Elements are ordered lexicographically.
impl Ord for Fr {
    #[inline(always)]
    fn cmp(&self, other: &Fr) -> std::cmp::Ordering {
        self.into_repr().cmp(&other.into_repr())
    }
}
impl PartialOrd for Fr {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fr) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl std::fmt::Display for Fr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Fr({})", self.into_repr())
    }
}

impl From<Fr> for FrRepr {
    fn from(e: Fr) -> FrRepr {
        e.into_repr()
    }
}


// impl BaseFromRO for Fr {
//     type Length = U48;
//     fn from_okm(okm: &GenericArray<u8, U48>) -> Fr {
//         const F_2_192: Fr = Fr(FrRepr([
//             0x59476ebc41b4528fu64,
//             0xc5a30cb243fcc152u64,
//             0x2b34e63940ccbd72u64,
//             0x1e179025ca247088u64,
//         ]));
//         let mut repr = FrRepr::default();
//         repr.read_be(Cursor::new([0; 8]).chain(Cursor::new(&okm[..24])))
//             .unwrap();
//         let mut elm = Fr::from_repr(repr).unwrap();
//         elm.mul_assign(&F_2_192);
//         repr.read_be(Cursor::new([0; 8]).chain(Cursor::new(&okm[24..])))
//             .unwrap();
//         elm.add_assign(&Fr::from_repr(repr).unwrap());
//         elm
//     }
// }
