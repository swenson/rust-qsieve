use num_bigint::*;
use num_traits::ToPrimitive;

// get_roots returns the roots of x^2 - num modulo p.
// double roots are only included once.
pub fn get_roots(num: &BigInt, p: u64) -> Vec<u64> {
    let np = (num % p).to_i64().unwrap();
    let mut roots: Vec<u64> = vec![];
    for a in 0..p as i64 {
        if (a * a - np) % p as i64 == 0 {
            roots.push(a as u64);
        }
    }
    roots
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn get_roots_test() {
        assert_eq!(get_roots(&BigInt::from(123), 2), vec![1]);
        assert_eq!(get_roots(&BigInt::from(123), 7), vec![2, 5]);
    }
}
