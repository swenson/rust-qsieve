use crate::latest::{multiple_poly::QuadraticPoly, poly::jacobi_bn};
use lib::{primes::get_primes, timers::Timers};
use memoize::lazy_static::lazy_static;
use num_integer::Integer;
use num_traits::{FromPrimitive, One, Signed, Zero};
use rug::{Complete, Integer as Int};
use std::simd::{LaneCount, Simd, SupportedLaneCount};
use std::{collections::HashSet, rc::Rc, simd::SimdPartialOrd, time::Instant};

#[derive(Debug, Clone)]
pub struct Relation {
    pub poly: Rc<QuadraticPoly>,
    pub num: Int,
    pub factors: Vec<i64>,
}

#[derive(Debug, Clone)]
pub struct TinyPrimes {
    pub prod: u64,
    pub big_prod: Int,
    pub primes: Vec<u64>,
}

impl TinyPrimes {
    pub fn generate(modulo: &Int) -> TinyPrimes {
        let mut prod: u64 = 1;
        let mut primes = vec![];
        for p in get_primes(5000) {
            if p == 2 {
                continue;
            }
            if jacobi_bn(modulo, p) != 1 {
                continue;
            }
            match prod.checked_mul(p) {
                Some(newprod) => {
                    prod = newprod;
                    primes.push(p);
                }
                None => {
                    break;
                }
            }
        }
        TinyPrimes {
            prod,
            big_prod: Int::from(prod),
            primes,
        }
    }
}

fn big_gcd(a: &Int, mut v: u64) -> u64 {
    let mut u = (a % v).complete().to_u64().unwrap();
    //u >>= u.trailing_zeros();
    while u & 1 == 0 {
        u >>= 1;
    }
    // we know both params are odd and non-zero
    // while v != 0 {
    //     (u, v) = (v, u % v);
    // }
    // u
    // binary GCD is slightly faster
    loop {
        if u > v {
            (u, v) = (v, u);
        }
        v -= u;
        if v == 0 {
            return u;
        }
        while v & 1 == 0 {
            v >>= 1;
        }
        //v >>= v.trailing_zeros();
    }
}

pub fn tiny_factor(num: &Int, tiny_primes: &TinyPrimes) -> (Int, Vec<i64>) {
    let mut m: Int = num.clone();
    let twos = m.find_one(0).unwrap_or(0);
    let mut result = vec![];
    result.resize(twos as usize, 2i64);
    m >>= twos;
    loop {
        let mut g = big_gcd(&m, tiny_primes.prod);
        if g == 1 {
            break;
        }
        m /= g;
        for p in tiny_primes.primes.iter() {
            while g % p == 0 {
                g /= p;
                result.push(*p as i64);
            }
            if g == 1 {
                break;
            }
        }
    }
    (m, result)
}

const POLLARD_B: u64 = 50;
lazy_static! {
    static ref POLLARD_PRIMES: (Vec<u64>, Vec<u32>) = pollard_make_powers(POLLARD_B);
}

fn pollard_make_powers(b: u64) -> (Vec<u64>, Vec<u32>) {
    let primes = get_primes(POLLARD_B);
    let mut result_p = vec![];
    let mut result_e = vec![];
    let log_b = (b as f64).ln();
    for p in primes {
        let e = (log_b / (p as f64).ln()).floor() as u32;
        assert!(e >= 1);
        result_p.push(p);
        result_e.push(e);
    }
    (result_p, result_e)
}

pub fn scan(
    poly: Rc<QuadraticPoly>,
    center: &Int,
    sieve: &[u8],
    add_to_point: usize,
    large_prime_bound: u64,
    log_scale: f64,
    timers: &Timers,
) -> Vec<usize> {
    let mut smooth_points = vec![];
    let cutoff = calculate_cutoff(poly, center, add_to_point, large_prime_bound, log_scale);
    // always push the 0 point to check since it will have the smallest size
    //smooth_points.push(0);
    let scan_start = Instant::now();
    assert_eq!(0, sieve.len() % 64);
    for i in 0..sieve.len() {
        if sieve[i] >= cutoff {
            smooth_points.push(i + add_to_point);
        }
    }
    timers.record("scan", scan_start);
    smooth_points
}

fn calculate_cutoff(
    poly: Rc<QuadraticPoly>,
    center: &Int,
    add_to_point: usize,
    large_prime_bound: u64,
    log_scale: f64,
) -> u8 {
    let lg2_large_prime_bound = (large_prime_bound as f64).log2();
    // always push the 0 point to check since it will have the smallest size
    //smooth_points.push(0);
    let fudge = -1.0;
    ((poly
        .eval(&(center + add_to_point).complete())
        .abs()
        .to_f64()
        .log2()
        - lg2_large_prime_bound
        - fudge)
        / log_scale) as u8
}

pub fn scan_simd_specific<const LANES: usize>(
    poly: Rc<QuadraticPoly>,
    center: &Int,
    sieve: &[u8],
    add_to_point: usize,
    large_prime_bound: u64,
    log_scale: f64,
    timers: &Timers,
) -> Vec<usize>
where
    LaneCount<LANES>: SupportedLaneCount,
{
    let mut smooth_points = vec![];
    let cutoff = calculate_cutoff(poly, center, add_to_point, large_prime_bound, log_scale);
    let scan_start = Instant::now();

    assert_eq!(0, sieve.len() % 64);
    let sieved = sieve.as_ptr() as *const Simd<u8, LANES>;
    let mask = Simd::<u8, LANES>::splat(cutoff);

    for i in 0..sieve.len() >> LANES.trailing_zeros() {
        let x = unsafe { (*(sieved.offset(i as isize))).simd_ge(mask) };
        if x.any() {
            for k in 0..LANES {
                unsafe {
                    if x.test_unchecked(k) {
                        smooth_points.push((i << LANES.trailing_zeros()) + k + add_to_point);
                    }
                }
            }
        }
    }
    let name = format!("scan simd u8x{}", LANES);
    timers.record(name.as_str(), scan_start);
    smooth_points
}

//#[cfg(target_arch = "x86_64")]
pub fn scan_asm(
    poly: Rc<QuadraticPoly>,
    center: &Int,
    sieve: &[u8],
    add_to_point: usize,
    large_prime_bound: u64,
    log_scale: f64,
    timers: &Timers,
) -> Vec<usize> {
    scan_simd_specific::<64>(
        poly,
        center,
        sieve,
        add_to_point,
        large_prime_bound,
        log_scale,
        timers,
    )
}

// #[cfg(target_arch = "aarch64")]
// pub fn scan_asm(
//     poly: Rc<QuadraticPoly>,
//     center: &Int,
//     sieve: &[u8],
//     add_to_point: usize,
//     large_prime_bound: u64,
//     log_scale: f64,
//     timers: &Timers,
// ) -> Vec<usize> {
//     use std::arch::asm;
//     let mut smooth_points = vec![];
//     let cutoff = calculate_cutoff(poly, center, add_to_point, large_prime_bound, log_scale);
//     assert_eq!(0, sieve.len() % 64);
//     let mut success: u64 = 0;
//     let mut ptr: *const u8 = &sieve[0] as *const u8;
//     let scan_start = Instant::now();
//
//     unsafe {
//         asm!(
//             "mov {ptr}, {sieve}",
//             "add {max}, {ptr}, {len}",
//             "ld1r {vcutoff}.8b, {cutoff}",
//
//             "2:",
//             "ldr {temp}, {ptr}",
//             "mov {vtemp}, {temp}",
//             "cmhs.8b {vtemp2}, {vtemp}, {vcutoff}",
//             "mov {temp}, {vtemp2}",
//             "cmp {temp}, 0",
//             "b.eq 3",
//             "add {ptr}, {ptr}, 8",
//             "cmp {ptr}, {max}",
//             "b.ne 2",
//             "b 4",
//             "3:",
//             "mov {success}, 1",
//             "4:",
//             success = out(reg) success,
//             ptr = inout(reg) ptr,
//             max = out(reg) _,
//             temp = out(reg) _,
//             vtemp = out(vreg) _,
//             vtemp2 = out(vreg) _,
//             vcutoff = out(vreg) _,
//             cutoff = in(reg) cutoff,
//             sieve = in(reg) &sieve[0] as *const u8,
//             len = in(reg) sieve.len(),
//         );
//     }
//
//     if success != 0 {
//         let offset = ptr as usize - (&sieve[0] as *const u8) as usize;
//         for i in offset..offset + 8 {
//             if sieve[i] >= cutoff {
//                 smooth_points.push(unsafe { i.unchecked_add(add_to_point) });
//             }
//         }
//         // TODO: repeat scan from offset
//     }
//     timers.record("scan asm", scan_start);
//     smooth_points
// }

pub fn smooth_check(
    smooth_points: &[usize],
    poly: Rc<QuadraticPoly>,
    center: &Int,
    offset: isize,
    primes: &Vec<u64>,
    tiny_primes: &TinyPrimes,
    prime_set: &HashSet<u64>,
    large_prime_bound: u64,
    timers: &Timers,
    a_factors: Rc<Vec<u32>>,
) -> (Box<Vec<Relation>>, Box<HashSet<u64>>) {
    let mut large_primes = Box::new(HashSet::new());
    let mut relations = Box::new(vec![]);
    let small_prime_bound = *match primes.last() {
        Some(p) => p,
        None => {
            return (relations, large_primes);
        }
    };

    let a = poly.a.as_ref();
    let asq = a.clone().sqrt().to_i64().unwrap_or(-2);
    let b = poly.b.as_ref();
    let b2 = &(b >> 1u32).complete();

    let smooth_start = Instant::now();

    for ppoint in smooth_points.iter() {
        let point = *ppoint as isize + offset;
        let x = (center + point).complete();
        let smooth_check = poly.eval(&x);
        let (mut smooth_factors, smooth, left) =
            small_factor(&smooth_check, &primes, tiny_primes, &prime_set);
        // weird edge case, skip these
        if smooth_factors.contains(&asq) {
            continue;
        }
        if !a.is_one() {
            for a_factor in a_factors.iter() {
                smooth_factors.push(*a_factor as i64);
            }
        }
        if !smooth {
            if left.to_u64().unwrap_or(large_prime_bound + 1) > large_prime_bound {
                continue;
            }
            let lp = left.to_u64().unwrap();
            smooth_factors.push(lp as i64);
            large_primes.insert(lp);
        }
        let count_set = smooth_factors
            .iter()
            .filter(|q| **q != -1)
            .filter(|q| **q != asq)
            .filter(|q| **q as u64 >= small_prime_bound)
            .map(|q| {
                large_primes.insert(*q as u64);
                *q as u64
            })
            .collect::<HashSet<u64>>();
        if count_set.len() > 1 {
            continue;
        }
        relations.push(Relation {
            poly: poly.clone(),
            num: a * x + b2,
            factors: smooth_factors,
        });
    }
    timers.record("smooth check", smooth_start);
    return (relations, large_primes);
}

fn small_factor(
    m: &Int,
    primes: &[u64],
    tiny_primes: &TinyPrimes,
    prime_set: &HashSet<u64>,
) -> (Vec<i64>, bool, Int) {
    let mut factors = vec![];
    if m.is_zero() {
        return (factors, true, Int::one());
    } else if m.is_negative() {
        factors.push(-1);
    }
    let (mut mm, mut new_factors) = tiny_factor(&m.abs(), tiny_primes);
    if !new_factors.is_empty() {
        factors.append(&mut new_factors);
    }
    match mm.to_u64() {
        Some(mut mmu) => {
            let largest_prime = *primes.last().unwrap();
            let pollard_factors = vec![]; //pollard_p1(mmu).unwrap_or(vec![]);
            let mut left = 1u64;
            for f in pollard_factors.iter() {
                if *f < largest_prime {
                    factors.push(*f as i64);
                    mmu /= f;
                } else {
                    left *= *f;
                }
            }
            let new_factors = small_factor64(mmu, primes, prime_set);
            for f in new_factors.iter() {
                if *f < largest_prime {
                    factors.push(*f as i64);
                } else {
                    left *= *f;
                }
            }
            return (factors, left == 1, Int::from_u64(left).unwrap());
        }
        _ => {}
    }
    for xi in 0..primes.len() {
        let x = primes[xi];
        loop {
            if !(&mm % x).complete().is_zero() {
                break;
            }
            mm /= x;
            factors.push(x as i64);
            match mm.to_u64() {
                Some(n) => {
                    if prime_set.contains(&n) {
                        factors.push(n as i64);
                        return (factors, true, Int::one());
                    } else {
                        let largest_prime = *primes.last().unwrap();
                        let mut m = n;
                        let mut left = 1u64;
                        let pollard_factors = vec![]; //pollard_p1(m).unwrap_or(vec![]);
                        for f in pollard_factors.iter() {
                            if *f < largest_prime {
                                factors.push(*f as i64);
                                m /= f;
                            } else {
                                m /= f;
                                left *= *f;
                            }
                        }
                        let new_factors = small_factor64(m, &primes[xi + 1..], prime_set);
                        for f in new_factors.iter() {
                            if f < primes.last().unwrap() {
                                factors.push(*f as i64);
                            } else {
                                left *= *f;
                            }
                        }
                        return (factors, left == 1, Int::from(left));
                    }
                }
                _ => {}
            }
        }
        if mm.is_one() {
            return (factors, true, Int::one());
        }
    }
    (factors, false, mm)
}

fn small_factor64(num: u64, primes: &[u64], prime_set: &HashSet<u64>) -> Vec<u64> {
    let mut factors = vec![];
    if num <= 1 {
        return factors;
    }
    if primes.len() == 0 {
        return factors;
    }
    let mut mm = num;
    let maxprime = *primes.last().unwrap();
    for xi in 0..primes.len() {
        let x = primes[xi];
        loop {
            let (mmx, mmd) = mm.div_rem(&x);
            if mmd != 0 {
                break;
            }
            mm = mmx;
            factors.push(x);
            if mm <= maxprime {
                if prime_set.contains(&mm) {
                    factors.push(mm);
                    return factors;
                }
            }
        }
        if mm == 1 {
            return factors;
        }
    }
    if mm != 1 {
        factors.push(mm);
    }
    factors
}

#[allow(dead_code)]
fn mul_vec(v: &[i64]) -> Int {
    if v.len() == 0 {
        return Int::one();
    } else if v.len() == 1 {
        return Int::from_i64(v[0]).unwrap();
    } else {
        let a = mul_vec(&v[..v.len() / 2]);
        let b = mul_vec(&v[v.len() / 2..]);
        return a * b;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn small_factor_test() {
        let primes = get_primes(1_000_000);
        let prime_set = primes.iter().map(|x| *x).collect();
        assert_eq!(
            small_factor(
                &Int::from(4),
                &primes,
                &TinyPrimes::generate(&Int::from(4)),
                &prime_set
            )
            .0,
            vec![2, 2]
        );
        assert_eq!(
            small_factor(
                &Int::from(-15),
                &primes,
                &TinyPrimes::generate(&Int::from(15)),
                &prime_set
            )
            .0,
            vec![-1, 3, 5]
        );
        assert_eq!(
            small_factor(
                &Int::from(1023),
                &primes,
                &TinyPrimes::generate(&Int::from(1023)),
                &prime_set
            )
            .0,
            vec![3, 11, 31]
        );
        assert_eq!(
            small_factor(
                &Int::from(7708410801051940u64),
                &primes,
                &TinyPrimes::generate(&Int::from(7708410801051940u64)),
                &prime_set
            )
            .0,
            vec![2, 2, 5, 7]
        );
        assert_eq!(
            small_factor(
                &Int::from(13301488142305634u64),
                &primes,
                &TinyPrimes::generate(&Int::from(13301488142305634u64)),
                &prime_set
            )
            .0,
            vec![2, 46271, 338383, 424769]
        );
        assert_eq!(
            small_factor(
                &Int::from(9582621056685920u64),
                &primes,
                &TinyPrimes::generate(&Int::from(9582621056685920u64)),
                &prime_set
            )
            .0,
            vec![2, 2, 2, 2, 2, 5, 17, 19, 19, 5443]
        );
    }
}
