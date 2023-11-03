use num_bigint::BigInt;
use num_traits::Pow;
use rug::Integer as Int;

pub mod nums;
pub mod opt;
pub mod params;
pub mod primes;
pub mod timers;

pub fn mul_vec_rug(v: &Vec<i64>) -> Int {
    let mut n = Int::from(1);
    for a in v.iter() {
        n *= *a;
    }
    n
}

pub fn mul_vec(v: &Vec<i64>) -> BigInt {
    let mut n = BigInt::from(1);
    for a in v.iter() {
        n *= *a;
    }
    n
}

#[allow(dead_code)]
pub fn ell(f: f64, nf: f64) -> f64 {
    let ln = nf.ln();
    let lln = ln.ln();
    return (ln.pow(f) * lln.pow(1.0 - f)).exp();
}
