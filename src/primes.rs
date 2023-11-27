use std::collections::HashSet;

use bitvector::BitVector;
use memoize::lazy_static::lazy_static;

use crate::nums::{modpow64, mulmod64};

const PRIME_SIEVE_SIZE: usize = 65536;

pub fn get_n_primes(n: u64) -> Vec<u64> {
    get_all_primes().take(n as usize).collect()
}

pub fn get_primes(bound: u64) -> Vec<u64> {
    let mut primes: Vec<u64> = vec![];
    if bound < 2 {
        return primes;
    }
    let maxi = (bound - 3) >> 1;
    let mut mem = BitVector::ones(maxi as usize + 1);
    primes.push(2);

    let mut next = 0;
    loop {
        if next > maxi {
            break;
        }
        if mem.contains(next as usize) {
            let q = (next << 1) + 3;
            primes.push(q);
            let mut i = next;
            while i <= maxi {
                mem.remove(i as usize);
                i += q;
            }
        }
        next += 1;
    }

    primes
}

#[derive(Clone, Debug)]
pub struct AllPrimes {
    primes: Vec<u64>,
    current_idx: usize,
    block_offset: usize,
    block: [u64; PRIME_SIEVE_SIZE],
}

impl Iterator for AllPrimes {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        while self.current_idx >= self.primes.len() {
            self.block.fill(0);
            let offset = self.block_offset * 64;
            for p in self.primes.iter() {
                let mut idx = offset.next_multiple_of(*p as usize) - offset;
                while idx < PRIME_SIEVE_SIZE * 64 {
                    self.block[idx >> 6] |= 1<<(idx & 0x3f);
                    idx += *p as usize;
                }
            }
            for idx in 0..PRIME_SIEVE_SIZE {
                if self.block[idx].count_zeros() == 0 {
                    continue;
                }
                for i in 0..64 {
                    if self.block[idx] & (1 << i) == 0 {
                        let q = offset + idx * 64 + i;
                        if q == 1 {
                            continue;
                        }
                        self.primes.push(q as u64);
                    }
                }
            }
            self.block_offset += PRIME_SIEVE_SIZE;
        }
        self.current_idx += 1;
        Some(self.primes[self.current_idx - 1])
    }
}

pub fn get_all_primes() -> AllPrimes {
    AllPrimes {
        primes: get_primes(PRIME_SIEVE_SIZE as u64 * 64),
        current_idx: 0,
        block_offset: 0,
        block: [0u64; PRIME_SIEVE_SIZE],
    }
}

const TRIAL_DIVIDE_LIMIT: u64 = 100;
lazy_static! {
    static ref TRIAL_DIVIDE_PRIMES: Vec<u64> =
        get_primes(TRIAL_DIVIDE_LIMIT).into_iter().skip(1).collect();
    static ref TRIAL_DIVIDE_PRIME_SET: HashSet<u64> = get_primes(100_000).into_iter().collect();
}

const MILLER_RABIN_POWERS: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

pub fn is_prime64(n: u64) -> bool {
    if TRIAL_DIVIDE_PRIME_SET.contains(&n) {
        return true;
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    for a in MILLER_RABIN_POWERS.iter() {
        let mut x = modpow64(*a, d, n);
        let mut y = 0;
        for _ in 0..s {
            y = mulmod64(x, x, n);
            if y == 1 && x != 1 && x != n - 1 {
                return false;
            }
            x = y;
        }
        if y != 1 {
            return false;
        }
    }
    true
}

mod test {
    #[allow(unused_imports)]
    use super::*;

    #[test]
    fn get_n_primes_test() {
        assert_eq!(get_n_primes(2), vec![2, 3]);
        assert_eq!(get_n_primes(3), vec![2, 3, 5]);
        assert_eq!(get_n_primes(4), vec![2, 3, 5, 7]);
        assert_eq!(get_n_primes(5), vec![2, 3, 5, 7, 11]);
        assert_eq!(
            get_n_primes(25),
            vec![
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97
            ]
        );
        assert_eq!(
            get_n_primes(172),
            vec![
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
                173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
                263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
                359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449,
                457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563,
                569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653,
                659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
                769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877,
                881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991,
                997, 1009, 1013, 1019, 1021
            ]
        );
    }

    #[test]
    fn get_primes_test() {
        assert_eq!(get_primes(3), vec![2, 3]);
        assert_eq!(get_primes(5), vec![2, 3, 5]);
        assert_eq!(get_primes(11), vec![2, 3, 5, 7, 11]);
        assert_eq!(get_primes(12), vec![2, 3, 5, 7, 11]);
        assert_eq!(
            get_primes(100),
            vec![
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97
            ]
        );
        assert_eq!(
            get_primes(1024),
            vec![
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
                173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
                263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
                359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449,
                457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563,
                569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653,
                659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
                769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877,
                881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991,
                997, 1009, 1013, 1019, 1021
            ]
        );
    }
}
