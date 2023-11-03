use lib::nums::{modpow32, mulmod32};
use memoize::lazy_static::lazy_static;
use num_integer::Integer;
use rug::{ops::RemRounding, Complete, Integer as Int};

use crate::latest::primes::chinese_remainder_theorem;

pub fn jacobi_bn(aa: &Int, mm: u64) -> i8 {
    return jacobi((aa % Int::from(mm)).to_u64().unwrap(), mm);
}

pub fn jacobi32(aa: u32, mm: u32) -> i8 {
    // algorithm 2.3.5 in Prime Numbers: A Computational Perspective
    let mut m = mm;
    let mut a = aa % m;
    let mut t = 1;
    while a != 0 {
        while a & 1 == 0 {
            a >>= 1;
            if m & 7 == 3 || m & 7 == 5 {
                t = -t;
            }
        }
        (a, m) = (m, a);
        if a & 3 == 3 && m & 3 == 3 {
            t = -t;
        }
        a %= m;
    }
    if m == 1 {
        return t;
    }
    return 0;
}

pub fn jacobi(aa: u64, mm: u64) -> i8 {
    Int::from(aa).jacobi(&Int::from(mm)) as i8
    // algorithm 2.3.5 in Prime Numbers: A Computational Perspective
    // let mut m = mm;
    // let mut a = aa % m;
    // let mut t = 1;
    // while a != 0 {
    //     while a & 1 == 0 {
    //         a >>= 1;
    //         if m & 7 == 3 || m & 7 == 5 {
    //             t = -t;
    //         }
    //     }
    //     (a, m) = (m, a);
    //     if a & 3 == 3 && m & 3 == 3 {
    //         t = -t;
    //     }
    //     a %= m;
    // }
    // if m == 1 {
    //     return t;
    // }
    // return 0;
}

pub fn jacobi_nomod(aa: u64, mm: u64) -> i8 {
    // algorithm 2.3.5 in Prime Numbers: A Computational Perspective
    let mut m = mm;
    let mut a = aa;
    let mut t = 1;
    while a != 0 {
        while a & 1 == 0 {
            a >>= 1;
            if m & 7 == 3 || m & 7 == 5 {
                t = -t;
            }
        }
        (a, m) = (m, a);
        if a & 3 == 3 && m & 3 == 3 {
            t = -t;
        }
        a %= m;
    }
    if m == 1 {
        return t;
    }
    return 0;
}

#[allow(dead_code)]
pub fn is_quadratic_residue64_nomod(num: u64, p: u64) -> bool {
    jacobi_nomod(num, p) == 1
}

#[derive(Clone, Debug, PartialEq)]
pub struct DoubleRoot {
    pub r0: u64,
    pub r1: u64,
    pub prime: u32,
    pub logp: u8,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Roots {
    Zero(),
    One {
        r: u64,
        prime: u32,
        logp: u8,
    },
    Two {
        r0: u64,
        r1: u64,
        prime: u32,
        logp: u8,
    },
}

pub fn roots_mod_p2(a: &Int, b: &Int, c: &Int, p: i64) -> Roots {
    let ainv_opt = mod_inverse(a, p);
    if ainv_opt.is_none() {
        return Roots::Zero();
    }
    let rootsp = roots_mod_p(a, b, c, p, ainv_opt.unwrap(), None);

    // hensel lift
    // Algorithm 2.3.11 in P:ACP
    match rootsp {
        Roots::Zero() => Roots::Zero(),
        Roots::One { r, logp, .. } => Roots::One {
            r: hensel_lift(a, b, c, p as u64, r as u64),
            prime: p as u32,
            logp,
        },
        Roots::Two { r0, r1, logp, .. } => {
            let q1 = hensel_lift(a, b, c, p as u64, r0 as u64);
            let q2 = hensel_lift(a, b, c, p as u64, r1 as u64);
            Roots::Two {
                r0: q1.min(q2),
                r1: q1.max(q2),
                prime: p as u32,
                logp,
            }
        }
    }
}

fn hensel_lift(a: &Int, b: &Int, c: &Int, p: u64, r0: u64) -> u64 {
    // Algorithm 2.3.11 in P:ACP
    //let p2 = BigInt::from_u64(p).unwrap() * BigInt::from_u64(p).unwrap();
    let r = r0 as i64;
    let x = ((a * (r * r)).complete() + b * r + c) / p;
    let arb = ((a * (2 * r) + b).complete() % p)
        .to_i64()
        .unwrap()
        .rem_euclid((p) as i64);
    let z = arb.extended_gcd(&(p as i64)).x;
    let y = ((-x * z) % p).to_i64().unwrap().rem_euclid((p) as i64) as u64;
    let r1 = r0 + y * p;
    r1 as u64
}

pub fn mod_inverse64(x: i64, y: i64) -> Option<i64> {
    // // Algorithm 9.4.3 in P:AC, modified knowing y is prime
    // x %= y;
    // if x == 0 {
    //     return None;
    // }
    // let (mut a, mut b, mut h) = (1, 0, x);
    // let (mut v1, mut v2, mut v3) = (y, 1 - x, y);
    // let (mut t1, mut t2, mut t3, mut skip_first_half) = if x & 1 == 0 {
    //     (1, 0, x, false)
    // } else {
    //     (0, -1, -y, true)
    // };

    // loop {
    //     loop {
    //         // halve t3
    //         if skip_first_half {
    //             skip_first_half = false;
    //         } else {
    //             if t1 & 1 == 1 || t2 & 1 == 1 {
    //                 t1 += y;
    //                 t2 -= x;
    //             }
    //             t1 >>= 1;
    //             t2 >>= 1;
    //             t3 >>= 1;
    //         }
    //         // check even
    //         if t3 & 1 != 0 {
    //             break;
    //         }
    //     }
    //     // reset max
    //     if t3 > 0 {
    //         (a, b, h) = (t1, t2, t3);
    //     } else {
    //         (v1, v2, v3) = (y - t1, -x - t2, -t3);
    //     }
    //     // subtract
    //     (t1, t2, t3) = (a - v1, b - v2, h - v3);
    //     if t1 < 0 {
    //         t1 += y;
    //         t2 -= x;
    //     }
    //     if t3 == 0 {
    //         return Some(a);
    //     }
    // }

    if x == 0 {
        return None;
    }
    let mut u = x;
    let mut v = y;
    let mut x1 = 1;
    let mut x2 = 0i64;
    while u != 1 {
        let q = v / u;
        let r = v % u; //v - q * u;
        let xx = x2 - (q * x1);
        v = u;
        u = r;
        x2 = x1;
        x1 = xx;
    }
    Some(x1.rem_euclid(y))

    // let p = y;
    // let mut z = x;
    // if z == 0 {
    //     return None;
    // }
    // let mut a = 1;
    // while z != 1 {
    //     let q = p / z;
    //     z = p % z;
    //     a = (-q * a).rem_euclid(p);
    // }
    // Some(a)

    // known good
    // let p = y;
    // let mut z = x;
    // if z == 0 {
    //     return None;
    // }
    // let mut a = 1;
    // while z != 1 {
    //     let q = -(p / z);
    //     z = p + q * z;
    //     a = (q * a).rem_euclid(p);
    // }
    // Some(a)
}

pub fn mod_inverse(aa: &Int, b: i64) -> Option<i64> {
    match aa.to_i64() {
        Some(a) => mod_inverse64(a % b, b),
        None => mod_inverse64((aa % b).complete().to_i64().unwrap(), b),
    }
}

pub fn roots_mod_p(
    a: &Int,
    b: &Int,
    c: &Int,
    p: i64,
    a_inv: i64,
    n_sqrt_opt: Option<i64>,
) -> Roots {
    if p == 2 {
        let r0 = !c.get_bit(0);
        let r1 = !(a.get_bit(0) ^ b.get_bit(0) ^ c.get_bit(0));
        if r0 && r1 {
            return Roots::Two {
                r0: 0,
                r1: 1,
                prime: 2,
                logp: 2,
            };
        } else if r0 || r1 {
            return Roots::One {
                r: if r0 { 0 } else { 1 },
                prime: 2,
                logp: 2,
            };
        } else {
            return Roots::Zero();
        }
    }
    //r = (-b +- sqrt(b^2 - 4 a c)) / 2 a
    // r = -(b/2) +- sqrt((b/2)^2 - a c) / a
    assert!(!b.get_bit(0));
    let half_b = b >> 1u32;
    let bp = (half_b.complete() % p).to_i64().unwrap();
    let n_sqrt = match n_sqrt_opt {
        Some(n) => n,
        None => {
            //let a2_inv = a2_inv_opt.unwrap();
            let cp = (c % p).complete().to_i64().unwrap();
            let d = bp * bp - (a * cp).complete();
            let dp = d.rem_euc(p);
            let sqroots = square_roots_modp(dp.to_u32().unwrap(), p as u32);
            if sqroots.len() == 0 {
                return Roots::Zero();
            } else if sqroots[0] == 0 {
                if cp != 0 {
                    return Roots::Zero();
                }
            }
            sqroots[0] as i64
        }
    };
    //let a2: Int = (a << 1u32).complete();
    //let a2_inv_opt = mod_inverse(&a2, p);
    // if a2_inv_opt.is_none() {
    //     return Roots::Zero();
    // }
    let mb = p - bp;
    let r0 = ((mb - n_sqrt) * a_inv).rem_euclid(p);
    let r1 = ((mb + n_sqrt) * a_inv).rem_euclid(p);
    let logp = (64u32 - p.leading_zeros()) as u8;
    if r0 == r1 {
        Roots::One {
            r: r0 as u64,
            prime: p as u32,
            logp,
        }
    } else {
        Roots::Two {
            r0: r0.min(r1) as u64,
            r1: r0.max(r1) as u64,
            prime: p as u32,
            logp,
        }
    }
}

lazy_static! {
    static ref INT_ONE: Int = Int::from(1);
}

#[allow(dead_code)]
pub fn square_roots_mod_n(a: &Int, n_factors: &[u32]) -> Vec<Int> {
    let mut pieces: Vec<Int> = vec![];
    let mut a_s: Vec<Vec<u64>> = vec![];
    let mut n_s: Vec<u64> = vec![];
    for f in n_factors.iter() {
        let modf = (a % *f).complete();
        let roots = square_roots_modp(modf.to_u32().unwrap(), *f);

        let mut new_a_s = vec![];

        if a_s.len() == 0 {
            new_a_s.push(vec![roots[0]]);
            new_a_s.push(vec![roots[1]]);
        } else {
            for i in 0..a_s.len() {
                let mut v = a_s[i].clone();
                let mut w = a_s[i].clone();
                v.push(roots[0]);
                w.push(roots[1]);
                new_a_s.push(v);
                new_a_s.push(w);
            }
        }
        a_s = new_a_s;
        n_s.push(*f as u64);
    }

    for i in 0..a_s.len() {
        let soln = chinese_remainder_theorem(a_s[i].as_ref(), n_s.as_ref());
        pieces.push(soln);
    }
    pieces
}

#[allow(dead_code)]
pub fn square_roots_mod_n2(a: &Int, n_factors: &[u32]) -> Vec<Int> {
    let neg_a = (-a).complete();
    let mut pieces: Vec<Int> = vec![];
    let mut a_s: Vec<Vec<u64>> = vec![];
    let mut n_s: Vec<u64> = vec![];
    for f in n_factors.iter() {
        let modf = (a % *f).complete();
        let roots = square_roots_modp(modf.to_u32().unwrap(), *f);
        let roots2: Vec<u64> = roots
            .iter()
            .map(|r| hensel_lift(&INT_ONE, &Int::ZERO, &neg_a, *f as u64, *r))
            .collect();

        let mut new_a_s = vec![];

        if a_s.len() == 0 {
            new_a_s.push(vec![roots2[0]]);
            new_a_s.push(vec![roots2[1]]);
        } else {
            for i in 0..a_s.len() {
                let mut v = a_s[i].clone();
                let mut w = a_s[i].clone();
                v.push(roots2[0]);
                w.push(roots2[1]);
                new_a_s.push(v);
                new_a_s.push(w);
            }
        }
        a_s = new_a_s;
        n_s.push((*f as u64) * (*f as u64));
    }

    for i in 0..a_s.len() {
        let soln = chinese_remainder_theorem(a_s[i].as_ref(), n_s.as_ref());
        pieces.push(soln);
    }
    pieces
}

// get_roots returns the roots of x^2 - num modulo p.
// double roots are only included once.
#[allow(non_snake_case)]
pub fn square_roots_modp(a: u32, p: u32) -> Vec<u64> {
    if a == 0 {
        return vec![0];
    }
    // this never happens if the code is correct elsewhere
    // if !is_quadratic_residue64_nomod(a as u64, p as u64) {
    //     return vec![];
    // }

    // algorithm 2.3.8 in Prime Numbers: A Computational Perspective
    let p8 = p & 7;
    match p8 {
        3 | 7 => {
            let r = modpow32(a, (p + 1) >> 2, p);
            let p_minus_r = p - r;
            if r == p_minus_r {
                return vec![r as u64];
            }
            return vec![r.min(p_minus_r) as u64, r.max(p_minus_r) as u64];
        }
        5 => {
            let mut x = modpow32(a, (p + 3) >> 3, p);
            // TODO: do this correctly for nums >= 32 bits
            let c = mulmod32(x, x, p);
            if c != a {
                x = mulmod32(x, modpow32(2, (p - 1) >> 2, p), p);
            }
            let p_minus_x = p - x;
            if x == p_minus_x {
                return vec![x as u64];
            }
            return vec![x.min(p_minus_x) as u64, x.max(p_minus_x) as u64];
        }
        _ => {} // fall through to hard case
    }

    // hard case: p == 1 (mod 8)
    let mut d = 0;
    for dd in 2..p {
        if jacobi32(dd, p) == -1 {
            d = dd;
            break;
        }
    }
    let p1 = p - 1;
    let s = p1.trailing_zeros();
    let t = p1 >> s;
    let A = modpow32(a, t, p);
    let D = modpow32(d, t, p);
    let mut m = 0u32;
    let mut admp = A;
    for i in 0..s {
        if modpow32(admp, 1u32 << (s - 1 - i), p) == p1 {
            m += 1 << i;
            admp = mulmod32(admp, modpow32(D, 1u32 << i, p), p);
        }
    }
    let r = mulmod32(modpow32(a, (t + 1) >> 1, p), modpow32(D, m >> 1, p), p);
    let p_minus_r = p - r;
    if r == p_minus_r {
        return vec![r as u64];
    }
    return vec![r.min(p_minus_r) as u64, r.max(p_minus_r) as u64];
}

#[cfg(test)]
mod test {
    use lib::nums::modpow64;
    #[allow(unused_imports)]
    use lib::primes::get_primes;

    use super::*;
    #[allow(unused_imports)]
    use num_traits::Signed;
    #[allow(unused_imports)]
    use num_traits::{One, Zero};
    #[allow(unused_imports)]
    use rand::prelude::*;
    #[allow(unused_imports)]
    use rug::rand::RandState;

    #[allow(dead_code)]
    fn slow_roots_modp(num: &Int, p: u64) -> Vec<u64> {
        let np = (num % p).complete().to_i64().unwrap();
        let mut roots: Vec<u64> = vec![];
        for a in 0..p as i64 {
            if (a * a - np) % p as i64 == 0 {
                roots.push(a as u64);
            }
        }
        roots
    }

    #[test]
    fn square_roots_mod_n_test() {
        let mut results = square_roots_mod_n(&Int::from(123), vec![17, 19].as_ref());
        results.sort();
        assert_eq!(
            vec![
                Int::from(117),
                Int::from(155),
                Int::from(168),
                Int::from(206)
            ],
            results
        );
    }

    #[test]
    fn square_roots_mod_n2_test() {
        let mut results = square_roots_mod_n2(&Int::from(123), vec![17, 19].as_ref());
        results.sort();
        assert_eq!(
            vec![
                Int::from(14057),
                Int::from(26692),
                Int::from(77637),
                Int::from(90272)
            ],
            results
        );
    }

    #[test]
    fn test_modpow64() {
        assert_eq!(8, modpow64(2, 3, 11));
        let a = 123;
        let b = 32093;
        let c = 704838212083171;
        let aa = Int::from(a);
        let bb = Int::from(b);
        let cc = Int::from(c);
        assert_eq!(
            aa.pow_mod(&bb, &cc).unwrap().to_u64().unwrap(),
            modpow64(a, b, c)
        );
        assert_eq!(modpow64(2, 0, 3), 1);
    }

    #[test]
    fn get_roots_test_equals_slow() {
        let primes = get_primes(1000);
        let mut rng = RandState::new();
        for i in 3..100u32 {
            let num: &Int = &(((Int::one() << i).random_below(&mut rng) << 1) + 1);
            for p in primes.iter() {
                if num.to_u64().unwrap_or(0) == *p {
                    continue;
                }
                if !is_quadratic_residue64_nomod((num % p).complete().to_u64().unwrap(), *p) {
                    continue;
                }
                println!("Checking roots of {} mod {}", num, p);
                assert_eq!(
                    square_roots_modp((num % p).complete().to_u32().unwrap(), *p as u32),
                    slow_roots_modp(&num, *p)
                );
            }
        }
    }

    #[test]
    fn get_roots_test() {
        assert_eq!(square_roots_modp(123, 2), vec![1]);
        assert_eq!(square_roots_modp(123, 7), vec![2, 5]);
    }

    #[test]
    fn roots_mod_p_test() {
        assert_eq!(
            roots_mod_p(
                &Int::from(289),
                &Int::from(84),
                &Int::from(-47),
                47,
                mod_inverse64(289, 47).unwrap(),
                None
            ),
            Roots::Two {
                r0: 0,
                r1: 35,
                prime: 47,
                logp: 6,
            }
        );
        assert_eq!(
            roots_mod_p(&Int::one(), &Int::zero(), &Int::from(-15347), 47, 1, None),
            Roots::Two {
                r0: 5,
                r1: 42,
                prime: 47,
                logp: 6,
            }
        );

        assert_eq!(
            roots_mod_p(&Int::from(4), &Int::zero(), &Int::from(-15347), 2, 0, None),
            Roots::Zero()
        );

        assert_eq!(
            roots_mod_p(
                &Int::from(289),
                &Int::zero(),
                &Int::from(-15347),
                17,
                0,
                None
            ),
            Roots::Zero()
        );
    }

    #[test]
    fn roots_mod_p2_test() {
        assert_eq!(
            roots_mod_p2(&Int::from(289), &Int::from(84), &Int::from(-47), 47),
            Roots::Two {
                r0: 129,
                r1: 658,
                prime: 47,
                logp: 6,
            }
        );
        assert_eq!(
            roots_mod_p2(&Int::one(), &Int::zero(), &Int::from(-15347), 47),
            Roots::Two {
                r0: 230,
                r1: 1979,
                prime: 47,
                logp: 6,
            }
        );
        assert_eq!(
            roots_mod_p2(&Int::from(4), &Int::zero(), &Int::from(-15347), 47),
            Roots::Two {
                r0: 115,
                r1: 2094,
                prime: 47,
                logp: 6,
            }
        );
        assert_eq!(
            roots_mod_p2(&Int::from(4), &Int::zero(), &Int::from(-15347), 2),
            Roots::Zero()
        );
        assert_eq!(
            roots_mod_p2(&Int::from(289), &Int::zero(), &Int::from(-15347), 2),
            Roots::One {
                r: 1,
                prime: 2,
                logp: 2,
            }
        );
    }

    #[test]
    fn inverse_test() {
        assert_eq!(mod_inverse64(3, 13), Some(9));
        assert_eq!(mod_inverse64(7, 13), Some(2));
    }
}
