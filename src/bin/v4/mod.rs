extern crate nalgebra as na;
use lib::mul_vec;
use lib::primes::get_primes;
use lib::timers::Timers;
use na::DMatrix;
use num_bigint::*;
use num_integer::*;
use num_traits::ToPrimitive;
use num_traits::*;
use std::collections::HashMap;
use std::rc::Rc;
use std::time::Instant;

#[cfg(test)]
mod test;

mod linalg;
mod poly;

#[derive(Debug)]
struct Relation {
    num: BigInt,
    factors: Vec<i64>,
}

fn small_factor(m: &BigInt, fb: &Vec<u64>) -> Vec<i64> {
    let mut factors = vec![];
    if m.is_zero() {
        return factors;
    } else if m.is_negative() {
        factors.push(-1);
    }
    let mut mm = m.abs();
    for x in fb.iter() {
        loop {
            let (mmx, rem) = mm.div_rem(&BigInt::from(*x as i64));
            if rem.is_zero() {
                mm = mmx;
                factors.push(*x as i64);
            } else {
                break;
            }
        }
    }
    factors
}

pub fn main_fun(
    num_: &BigInt,
    timers: Rc<Timers>,
    _use_simd: bool,
    _use_asm: bool,
    _la_lanes: u8,
    _scan_lanes: u8,
    _l1_cache_bits: u8,
    _optimize: bool,
) {
    let num = num_.clone();
    let (smoothness_bound, area) = match num.bits() {
        122.. => panic!("Number too large"),
        102..=121 => (28_000u64, 50_000_000usize),
        82..=101 => (5_500u64, 16_000_000usize),
        62..=81 => (2_000u64, 1_000_000usize),
        42..=61 => (1_200u64, 80_000usize),
        30..=41 => (380u64, 2_000usize),
        _ => panic!("Number too small"),
    };

    let primes = get_primes(smoothness_bound);

    // x = ceil(sqrt(num))
    let sqrt_num = if num.sqrt().pow(2u32) < num {
        num.sqrt() + 1
    } else {
        num.sqrt()
    };

    let sieve_start = Instant::now();

    let mut sieve = Vec::with_capacity(area);
    sieve.resize(area, 0f32);
    for p in primes.iter() {
        let sqrtp = (sqrt_num.clone() % p).to_i64().unwrap();
        let logp = (*p as f64).log2() as f32;
        let roots = poly::get_roots(&num, *p);
        line_sieve(area, &mut sieve, *p as isize, sqrtp as isize, logp, roots);
    }
    timers.record("line sieve", sieve_start);

    let scan_start = Instant::now();
    let mut smooth_points = vec![];
    let recalculate_every = 65536;
    let fudge = 5.0;
    for j in (0..area).step_by(recalculate_every) {
        let start =
            BigInt::from(j + sqrt_num.clone()) * BigInt::from(j + sqrt_num.clone()) - num.clone();
        let cutoff = start.to_f64().unwrap().log2() as f32 - fudge;
        for i in j..((j + recalculate_every).min(area)) {
            if sieve[i] > cutoff {
                smooth_points.push(i)
            }
        }
    }
    timers.record("scan", scan_start);
    let smooth_check_start = Instant::now();

    let mut relations = vec![];

    for point in smooth_points.iter() {
        let smooth_check = BigInt::from(point + sqrt_num.clone())
            * BigInt::from(point + sqrt_num.clone())
            - num.clone();
        let smooth_factors = small_factor(&smooth_check, &primes);
        let check = mul_vec(&smooth_factors);
        if check != smooth_check {
            continue;
        }
        // println!(
        //     "Smooth point {:?} factors {:?}",
        //     smooth_check, smooth_factors
        // );

        relations.push(Relation {
            num: point + sqrt_num.clone(),
            factors: smooth_factors,
        });
    }
    timers.record("smooth check", smooth_check_start);

    println!("Relations {}", relations.len());
    println!("Factor base size: {}", primes.len());

    let matrix_build_start = Instant::now();

    // build a binary matrix
    let cols = primes.len();

    let mut col_map = HashMap::new();
    // keep track of column counts to remove empty columns
    let mut col_counts: HashMap<usize, u32> = HashMap::new();

    // index
    for pi in 0..primes.len() {
        col_map.insert(primes[pi], pi);
    }
    let mut mat = DMatrix::<u8>::zeros(relations.len(), cols);
    for r in 0..relations.len() {
        for p in relations[r].factors.iter() {
            let pi = *col_map.get(&(*p as u64)).unwrap();
            mat[(r, pi)] ^= 1;
            col_counts.insert(pi, col_counts.get(&pi).unwrap_or(&0u32) + 1);
        }
    }
    let mut del_cols = vec![];
    for ci in (0..cols).rev() {
        if *col_counts.get(&ci).unwrap_or(&0u32) == 0 {
            //println!("Remove column {}", ci);
            del_cols.push(ci);
            mat = mat.remove_column(ci);
        }
    }
    timers.record("matrix build", matrix_build_start);

    println!(
        "Factor base size after removing empty columns: {}",
        mat.ncols()
    );

    // if mat.ncols() > mat.nrows() {
    //     panic!("Not enough relations; try increasing sieve area");
    // }

    let left_kernel_start = Instant::now();
    println!("Solving {}x{} matrix", mat.ncols(), mat.nrows());
    let kernel = linalg::left_kernel_gf2(&mat).unwrap();
    timers.record("left kernel", left_kernel_start);

    let sqrt_start = Instant::now();

    for ki in 0..kernel.len() {
        let result_vector = kernel.get(ki).unwrap();
        //println!("Result = {:?}", result_vector);
        let check = result_vector
            .map(|x| x as u32)
            .tr_mul(&mat.map(|x| x as u32))
            .iter()
            .all(|x| x % 2 == 0);
        if !check {
            panic!("Null vector was invalid");
        }
        println!("Null vector is valid");

        println!("Computing squares");

        let mut a = BigInt::one();
        let mut b = BigInt::one();

        let mut all_factors = HashMap::new();

        for i in 0..mat.nrows() {
            if result_vector[i] == 0 {
                continue;
            }

            // println!(
            //     "Mul {} and {}",
            //     relations[i].num.clone() * relations[i].num.clone() - num.clone(),
            //     relations[i].num.clone()
            // );
            for f in relations[i].factors.iter() {
                all_factors.insert(f, all_factors.get(&f).unwrap_or(&0u32) + 1);
            }
            b = b * (relations[i].num.clone());
        }

        for (f, e) in all_factors.iter() {
            assert!(e % 2 == 0);
            a = a * BigInt::from_i64(**f)
                .unwrap()
                .modpow(&BigInt::from_u32(e / 2).unwrap(), &num)
                % num.clone();
        }
        a = a % num.clone();

        println!("a is ~2^{}, b is 2^{}", a.bits(), b.bits());
        println!("a = {}", a);
        println!("b = {}", b);

        // if a_r.clone() * a_r.clone() != a {
        //     panic!("a is not square!");
        // }
        // if b_r.clone() * b_r.clone() != b {
        //     panic!("b is not square!");
        // }
        // if !(a - b).rem(&num).is_zero() {
        //     panic!("number does not divide a^2 - b^2!")
        // }

        let g1 = (a.clone() + b.clone()).gcd(&num.clone());
        let g2 = (a - b).gcd(&num.clone());

        let mut factor_found = false;
        println!("g1 = {}, g2 = {}", g1, g2);

        if !g1.is_one() && g1 != num && (num.clone() % g1.clone()).is_zero() {
            println!("Factor found! {}", g1);
            factor_found = true;
        }
        if !g2.is_one() && g2 != num && (num.clone() % g2.clone()).is_zero() {
            println!("Factor found! {}", g2);
            factor_found = true;
        }
        if factor_found {
            break;
        } else {
            println!("No factor found :(");
        }
    }
    timers.record("sqrt", sqrt_start);
}

#[inline(always)]
fn line_sieve(
    area: usize,
    sieve: &mut Vec<f32>,
    p: isize,
    sqrtp: isize,
    logp: f32,
    roots: Vec<u64>,
) {
    if roots.len() == 2 {
        let i0 = (roots[0] as isize - sqrtp) % p;
        let i1 = (roots[1] as isize - sqrtp) % p;
        let (i0, i1) = (i0.min(i1), i0.max(i1));
        unsafe {
            line_sieve_double(sieve.as_mut_ptr(), i0, i1, logp, p, area as isize);
        }
        return;
    }
    for r in roots.iter() {
        let mut idx_i64 = (*r as i64 - sqrtp as i64) % p as i64;
        while idx_i64 < 0 {
            idx_i64 += p as i64;
        }
        let mut idx = idx_i64 as usize;
        while idx < area {
            sieve[idx] += logp;
            idx += p as usize;
        }
    }
}

#[inline(always)]
pub unsafe fn line_sieve_double(
    sieve: *mut f32,
    mut i0: isize,
    mut i1: isize,
    logp: f32,
    p_i: isize,
    block_size: isize,
) {
    let sieve0 = sieve.offset(i0 - i1);
    let sieve1 = sieve;
    let init_i1 = i1;
    while i1 < block_size {
        *sieve0.offset(i1) += logp;
        *sieve1.offset(i1) += logp;
        i1 += p_i;
    }
    i0 += i1 - init_i1;
    while i0 < block_size {
        *sieve1.offset(i0) += logp;
        i0 += p_i;
    }
}
