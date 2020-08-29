#![feature(test)]
extern crate test;

#[allow(unused_imports)]
#[macro_use]
extern crate uint;
// extern crate num_bigint;
// extern crate rug;
extern crate bls12_381;
extern crate fff;


// use rug::Integer;
// use rug::integer::Order;
// use num_bigint::BigUint;
// use uint::{construct_uint, uint_full_mul_reg};
use bls12_381::Scalar;

pub mod u256;

use self::u256::U256;

// construct_uint! {
//     pub struct U256(4);
// }

// construct_uint! {
//     pub struct U512(8);
// }

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

// impl U256 {
//     #[inline(always)]
//     pub fn full_mul(self, other: U256) -> U512 {
//         U512(uint_full_mul_reg!(U256, 4, self, other))
//     }
// }
#[test]
fn test_mul() {
    let c = U256([1u64, 1, 1, 0]);
    let d = U256([1u64, 1, 0, 0]);
    let out1 = c * d;

    let c = Scalar::from_raw([1u64, 1, 1, 0]);
    let d = Scalar::from_raw([1u64, 1, 0, 0]);
    let out2 = c * d;

    use fff::Fr;
    use fff::FrRepr;
    
    let c = Fr::from_repr(FrRepr([1u64, 1, 1, 0])).unwrap();
    let d = Fr::from_repr(FrRepr([1u64, 1, 0, 0])).unwrap();
    let mut out3 = c.clone();
    out3.mul_assign_soft(&d);

    // let mut out_a = [0u64; 4];
    // let mut octets = [0u8; 32];
    // out1.to_little_endian(&mut octets);
    // out_a[0] = u64::from_le_bytes([octets[ 0], octets[ 1], octets[ 2], octets[ 3], octets[ 4], octets[ 5], octets[ 6], octets[ 7] ]);
    // out_a[1] = u64::from_le_bytes([octets[ 8], octets[ 9], octets[10], octets[11], octets[12], octets[13], octets[14], octets[15] ]);
    // out_a[2] = u64::from_le_bytes([octets[16], octets[17], octets[18], octets[19], octets[20], octets[21], octets[22], octets[23] ]);
    // out_a[3] = u64::from_le_bytes([octets[24], octets[25], octets[26], octets[27], octets[28], octets[29], octets[30], octets[31] ]);
    // println!("{:?}", &out_a);


    let mut out_b = [0u64; 4];
    let mut octets = out2.to_bytes();
    // println!("{:?}", &octets);
    out_b[0] = u64::from_le_bytes([octets[ 0], octets[ 1], octets[ 2], octets[ 3], octets[ 4], octets[ 5], octets[ 6], octets[ 7] ]);
    out_b[1] = u64::from_le_bytes([octets[ 8], octets[ 9], octets[10], octets[11], octets[12], octets[13], octets[14], octets[15] ]);
    out_b[2] = u64::from_le_bytes([octets[16], octets[17], octets[18], octets[19], octets[20], octets[21], octets[22], octets[23] ]);
    out_b[3] = u64::from_le_bytes([octets[24], octets[25], octets[26], octets[27], octets[28], octets[29], octets[30], octets[31] ]);
    println!("{:?}", &out_b);


    // println!("{:?}", &out3.0.0);
    println!("{:?}", &out3.into_repr().0);
    
    let a: [u64; 4] = out1.0;
    let b: [u64; 4] = out_b;
    let c: [u64; 4] = out3.into_repr().0;
    
    assert_eq!(a, b);
    assert_eq!(a, c);
    assert_eq!(b, c);
}

#[bench]
fn bench_fixedsize_mul(b: &mut test::Bencher) {
    let c = test::black_box(U256([1u64, 1, 1, 0]));
    let d = test::black_box(U256([1u64, 1, 0, 0]));
    
    b.iter(|| c * d)
}

#[bench]
fn bench_fixedsize_add(b: &mut test::Bencher) {
    let c = U256([255u64, 255, 255, 255]);
    let d = U256([255u64, 255, 255, 255]);

    b.iter(|| c + d)
}

#[bench]
fn bench_fixedsize_sub(b: &mut test::Bencher) {
    let c = U256([255u64, 255, 255, 255]);
    let d = U256([255u64, 255, 255, 255]);

    b.iter(|| c - d)
}

// #[bench]
// fn bench_fixedsize_div(b: &mut test::Bencher) {
//     let c = U256([255u64, 255, 255, 255]);
//     let d = U256([255u64, 255, 255, 255]);

//     b.iter(|| c / d)
// }

// #[bench]
// fn bench_fixedsize_rem(b: &mut test::Bencher) {
//     let c = U256([255u64, 255, 255, 255]);
//     let d = U256([255u64, 255, 255, 255]);

//     b.iter(|| c % d)
// }



//    ----------------  BigInt ---------------
// #[bench]
// fn bench_bigint_mul(b: &mut test::Bencher) {
//     let c = BigUint::new(vec![1u32, 0, 1, 0, 1, 0, 0, 0]);
//     let d = BigUint::new(vec![1u32, 0, 1, 0, 0, 0, 0, 0]);

//     b.iter(|| c.clone() * &d)
// }

// #[bench]
// fn bench_bigint_div(b: &mut test::Bencher) {
//     let c = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);
//     let d = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);

//     b.iter(|| c.clone() / &d)
// }

// #[bench]
// fn bench_bigint_rem(b: &mut test::Bencher) {
//     let c = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);
//     let d = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);

//     b.iter(|| c.clone() % &d)
// }

// #[bench]
// fn bench_bigint_add(b: &mut test::Bencher) {
//     let c = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);
//     let d = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);

//     b.iter(|| c.clone() + &d)
// }

// #[bench]
// fn bench_bigint_sub(b: &mut test::Bencher) {
//     let c = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);
//     let d = BigUint::new(vec![255u32, 0, 255, 0, 255, 0, 255, 0]);

//     b.iter(|| c.clone() - &d)
// }


// -------- GNU GMP -----------------
// #[bench]
// fn bench_gmp_mul(b: &mut test::Bencher) {
//     let mut c = Integer::new();
//     let mut d = Integer::new();

//     c.assign_digits(&[1u64, 1, 1, 0], Order::Lsf);
//     d.assign_digits(&[1u64, 1, 0, 0], Order::Lsf);

//     b.iter(|| c.clone() * &d)
// }

// #[bench]
// fn bench_gmp_div(b: &mut test::Bencher) {
//     let mut c = Integer::new();
//     let mut d = Integer::new();

//     c.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);
//     d.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);

//     b.iter(|| c.clone() / &d)
// }

// #[bench]
// fn bench_gmp_rem(b: &mut test::Bencher) {
//     let mut c = Integer::new();
//     let mut d = Integer::new();

//     c.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);
//     d.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);

//     b.iter(|| c.clone() % &d)
// }


// #[bench]
// fn bench_gmp_add(b: &mut test::Bencher) {
//     let mut c = Integer::new();
//     let mut d = Integer::new();

//     c.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);
//     d.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);

//     b.iter(|| c.clone() + &d)
// }

// #[bench]
// fn bench_gmp_sub(b: &mut test::Bencher) {
//     let mut c = Integer::new();
//     let mut d = Integer::new();

//     c.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);
//     d.assign_digits(&[255u64, 255, 255, 255], Order::Lsf);

//     b.iter(|| c.clone() - &d)
// }

// ---- BLS12-381 -------
#[bench]
fn bench_bls12_381_mul(b: &mut test::Bencher) {
    let c = Scalar::from_raw([1u64, 1, 1, 0]);
    let d = Scalar::from_raw([1u64, 1, 0, 0]);

    b.iter(|| c * d)
}

// #[bench]
// fn bench_bls12_381_div(b: &mut test::Bencher) {
//     let c = Scalar::from_raw([255u64, 255, 255, 255]);
//     let d = Scalar::from_raw([255u64, 255, 255, 255]);

//     b.iter(|| c / d)
// }

// #[bench]
// fn bench_bls12_381_rem(b: &mut test::Bencher) {
//     let c = Scalar::from_raw([255u64, 255, 255, 255]);
//     let d = Scalar::from_raw([255u64, 255, 255, 255]);

//     b.iter(|| c % d)
// }


#[bench]
fn bench_bls12_381_add(b: &mut test::Bencher) {
    let c = Scalar::from_raw([255u64, 255, 255, 255]);
    let d = Scalar::from_raw([255u64, 255, 255, 255]);

    b.iter(|| c + d)
}

#[bench]
fn bench_bls12_381_sub(b: &mut test::Bencher) {
    let c = Scalar::from_raw([255u64, 255, 255, 255]);
    let d = Scalar::from_raw([255u64, 255, 255, 255]);

    b.iter(|| c - d)
}



// --------- FILECOIN ff & Pairing -------------
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
#[bench]
fn bench_fff_mul_amd64_adx_bmi2(b: &mut test::Bencher) {
    use fff::Fr;
    use fff::FrRepr;
    
    let c = Fr::from_repr(FrRepr([1u64, 1, 1, 0])).unwrap();
    let d = Fr::from_repr(FrRepr([1u64, 1, 0, 0])).unwrap();

    b.iter(|| {
        let mut x = c.clone();
        x.mul_assign_amd64_adx_bmi2(&d);
        x
    })
}

#[bench]
fn bench_fff_mul_soft(b: &mut test::Bencher) {
    use fff::Fr;
    use fff::FrRepr;
    
    let c = Fr::from_repr(FrRepr([1u64, 1, 1, 0])).unwrap();
    let d = Fr::from_repr(FrRepr([1u64, 1, 0, 0])).unwrap();
    
    b.iter(|| {
        let mut x = c.clone();
        x.mul_assign_soft(&d);
        x
    })
}

#[cfg(target_arch = "x86_64")]
#[bench]
fn bench_fff_add_amd64(b: &mut test::Bencher) {
    use fff::Fr;
    use fff::FrRepr;
    
    let c = Fr::from_repr(FrRepr([255u64, 255, 255, 255])).unwrap();
    let d = Fr::from_repr(FrRepr([255u64, 255, 255, 255])).unwrap();

    b.iter(|| {
        let mut x = c.clone();
        x.add_assign_amd64(&d);
        x
    })
}

#[bench]
fn bench_fff_add_soft(b: &mut test::Bencher) {
    use fff::Fr;
    use fff::FrRepr;
    
    let c = Fr::from_repr(FrRepr([255u64, 255, 255, 255])).unwrap();
    let d = Fr::from_repr(FrRepr([255u64, 255, 255, 255])).unwrap();
    
    b.iter(|| {
        let mut x = c.clone();
        x.add_assign_soft(&d);
        x
    })
}

#[bench]
fn bench_fff_sub(b: &mut test::Bencher) {
    use fff::Fr;
    use fff::FrRepr;
    
    let c = Fr::from_repr(FrRepr([255u64, 255, 255, 255])).unwrap();
    let d = Fr::from_repr(FrRepr([255u64, 255, 255, 255])).unwrap();
    
    b.iter(|| {
        let mut x = c.clone();
        x.sub_assign(&d);
        x
    })
}
