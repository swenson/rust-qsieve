use std::collections::HashMap;

use num_traits::Pow;
use rug::{Complete, Integer};

use crate::latest::poly::mod_inverse;

pub fn chinese_remainder_theorem(nums: &[u64], mods: &[u64]) -> Integer {
    let mut modulus = Integer::from(1);
    for m in mods.iter() {
        modulus = modulus * m;
    }
    let mut crt_inverses = HashMap::new();
    for p in mods.iter() {
        crt_inverses.insert(
            *p,
            mod_inverse(&(&modulus / p).complete(), *p as i64).unwrap(),
        );
    }
    let mut s = Integer::from(0);
    for (x, n) in nums.iter().zip(mods.iter()) {
        let n_inv = crt_inverses.get(n).unwrap();
        s += *x * (&modulus / *n).complete() * n_inv;
    }
    s % modulus
}

pub fn sqrt_ceil(num: Integer) -> Integer {
    if num.clone().sqrt().pow(2u32) < num {
        num.sqrt() + 1
    } else {
        num.sqrt()
    }
}

mod test {
    #[allow(unused_imports)]
    use super::*;

    #[test]
    fn chinese_remainder_theorem_test() {
        assert_eq!(
            Integer::from(39),
            chinese_remainder_theorem(vec![0, 3, 4].as_ref(), vec![3, 4, 5].as_ref())
        );
    }
}
