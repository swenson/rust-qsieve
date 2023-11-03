use super::*;

#[test]
fn small_factor_test() {
    let primes = get_primes(1024);
    assert_eq!(small_factor(&Int::from(4), &primes), vec![2, 2]);
    assert_eq!(small_factor(&Int::from(-15), &primes), vec![-1, 3, 5]);
    assert_eq!(small_factor(&Int::from(1023), &primes), vec![3, 11, 31]);
}
