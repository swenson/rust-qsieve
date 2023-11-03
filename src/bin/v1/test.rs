use lib::primes::get_primes;

use super::*;

#[test]
fn small_factor_test() {
    let primes = get_primes(1024);
    assert_eq!(small_factor(&BigInt::from(4), &primes), vec![2, 2]);
    assert_eq!(small_factor(&BigInt::from(-15), &primes), vec![-1, 3, 5]);
    assert_eq!(small_factor(&BigInt::from(1023), &primes), vec![3, 11, 31]);
}
