use lib::primes::get_primes;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use rug::{Complete, Integer as Int};
use std::collections::HashSet;

pub fn jacobi_bn(aa: &Int, mm: u64) -> i8 {
    return jacobi((aa % Int::from(mm)).to_u64().unwrap(), mm);
}

pub fn jacobi(aa: u64, mm: u64) -> i8 {
    Int::from(aa).jacobi(&Int::from(mm)) as i8
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

pub fn small_factor(
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
        Some(mmu) => {
            let largest_prime = *primes.last().unwrap();
            let mut left = 1u64;
            let new_factors = small_factor64(mmu, primes, prime_set);
            for f in new_factors.iter() {
                if *f < largest_prime {
                    factors.push(*f as i64);
                } else {
                    left *= *f;
                }
            }
            return (factors, left == 1, Int::from(left));
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
                        let m = n;
                        let mut left = 1u64;
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
