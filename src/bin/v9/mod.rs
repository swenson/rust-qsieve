extern crate nalgebra as na;
use crate::v9::smooth::{jacobi_bn, TinyPrimes};
use lib::mul_vec_rug;
use lib::primes::get_primes;
use lib::timers::Timers;
use na::DMatrix;
use num_bigint::BigInt;
use num_traits::*;
use rug::{Complete, Integer as Int};
use std::collections::{HashMap, HashSet};
use std::rc::Rc;
use std::time::Instant;

mod linalg;
mod line_sieve;
mod poly;
mod smooth;

#[derive(Debug)]
struct Relation {
    num: Int,
    factors: Vec<i64>,
}

pub fn main_fun(
    num_: &BigInt,
    timers: Rc<Timers>,
    _use_simd: bool,
    _use_asm: bool,
    _la_lanes: u8,
    _scan_lanes: u8,
    l1_cache_bits: u8,
    _optimize: bool,
) {
    let num = Int::parse(num_.to_str_radix(10)).unwrap().complete();
    let (smoothness_bound, area) = match num_.bits() {
        162.. => panic!("Number too large"),
        142..=161 => (68_000u64, 2_000_000_000usize),
        122..=141 => (38_000u64, 250_000_000usize),
        102..=121 => (28_000u64, 50_000_000usize),
        82..=101 => (5_500u64, 16_000_000usize),
        62..=81 => (2_000u64, 1_000_000usize),
        42..=61 => (1_200u64, 80_000usize),
        30..=41 => (380u64, 2_000usize),
        _ => panic!("Number too small"),
    };
    let area = area.next_multiple_of(1 << l1_cache_bits);

    let primes = Rc::new(get_primes(smoothness_bound));

    // x = ceil(sqrt(num))
    let sqrt_num: Int = num.clone().sqrt() + 1;

    let sieve_start = Instant::now();
    let mut sieve = vec![0u8; 1 << l1_cache_bits];
    timers.record("line sieve - alloc", sieve_start);
    let sieve_fun = match l1_cache_bits {
        14 => line_sieve::split_line_sieve::<14>,
        15 => line_sieve::split_line_sieve::<15>,
        16 => line_sieve::split_line_sieve::<16>,
        17 => line_sieve::split_line_sieve::<17>,
        18 => line_sieve::split_line_sieve::<18>,
        19 => line_sieve::split_line_sieve::<19>,
        20 => line_sieve::split_line_sieve::<20>,
        21 => line_sieve::split_line_sieve::<21>,
        22 => line_sieve::split_line_sieve::<22>,
        _ => panic!("L1 cache bits wrong size: {}", l1_cache_bits),
    };
    let smooth_points = sieve_fun(&mut sieve, primes.clone(), &num, area, &timers);
    let smooth_check_start = Instant::now();

    let mut relations = vec![];

    let tiny_primes = TinyPrimes::generate(&num);
    let primes_except_tiny: Vec<u64> = get_primes(smoothness_bound)
        .into_iter()
        .filter(|p| jacobi_bn(&num, *p) == 1)
        .filter(|p| p > tiny_primes.primes.last().unwrap_or(&0))
        .collect();
    let prime_set: HashSet<u64> = primes.iter().map(|x| *x).collect();

    for point in smooth_points.iter() {
        let smooth_check =
            Int::from(point + sqrt_num.clone()) * Int::from(point + sqrt_num.clone()) - num.clone();
        let (smooth_factors, _, _) = smooth::small_factor(
            &smooth_check,
            primes_except_tiny.as_slice(),
            &tiny_primes,
            &prime_set,
        );
        let check = mul_vec_rug(&smooth_factors);
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

        let mut a = Int::one();
        let mut b = Int::one();

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
            a = (a * Int::from(**f).pow_mod(&Int::from(e / 2), &num).unwrap()) % num.clone();
        }
        a = a % num.clone();

        println!(
            "a is ~2^{}, b is 2^{}",
            a.to_f64().abs().log2(),
            b.to_f64().abs().log2()
        );
        // println!("a = {}", a);
        // println!("b = {}", b);

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
pub unsafe fn line_sieve_double(
    sieve: *mut u8,
    mut i0: isize,
    mut i1: isize,
    logp: u8,
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
    // while i0 < block_size {
    //     *sieve.offset(i0) += logp;
    //     i0 += p_i;
    // }
    // while i1 < block_size {
    //     *sieve.offset(i1) += logp;
    //     i1 += p_i;
    // }
}
