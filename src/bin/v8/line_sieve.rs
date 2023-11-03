use crate::v8::{line_sieve_double, poly};
use lib::timers::Timers;
use rug::Integer as Int;
use std::rc::Rc;
use std::time::Instant;

// according to msieve, we want to split the factor base into tiny primes,
// whose updates are sieved into the block, and medium primes, whose updates
// are put into hash tables to be added when we sieve that block later.
// we start with the medium primes, then do per block.
pub fn split_line_sieve<const L1_SHIFT: u32>(
    sieve: &mut [u8],
    primes: Rc<Vec<u64>>,
    num: &Int,
    area: usize,
    timers: &Timers,
) -> Vec<usize> {
    let fudge = 6;
    let recalculate_every = 65536;
    let mut smooth_points = vec![];
    let l1_mask = (1 << L1_SHIFT) - 1;
    let l1_cache_size = 1 << L1_SHIFT;
    let roots_start = Instant::now();
    let mut roots = vec![];
    let mut medium_prime_offset = 0;
    let sqrt_num: Int = num.clone().sqrt() + 1;
    for p in primes.iter() {
        let r = poly::get_roots(num, *p);
        if r.len() != 2 {
            continue;
        }
        if *p > 3 * l1_cache_size && medium_prime_offset == 0 {
            medium_prime_offset = roots.len();
        }
        let r0 = r[0] as i64;
        let r1 = r[1] as i64;
        let sqrtp = (sqrt_num.clone() % p).to_i64().unwrap();
        let i0 = (r0 - sqrtp).rem_euclid(*p as i64) as u32;
        let i1 = (r1 - sqrtp).rem_euclid(*p as i64) as u32;
        roots.push((*p as u32, i0, i1));
    }
    if medium_prime_offset == 0 {
        medium_prime_offset = roots.len();
    }
    let num_medium_primes = roots.len() - medium_prime_offset;

    timers.record("roots", roots_start);

    let medium_start = Instant::now();
    let num_buckets = (area >> L1_SHIFT) + 1;
    let mut hash_tables = Vec::with_capacity(num_buckets);
    for _ in 0..num_buckets {
        hash_tables.push(Vec::with_capacity(
            num_medium_primes * 2 / l1_cache_size as usize + 4,
        ));
    }
    for (p, r0, r1) in roots.iter().skip(medium_prime_offset) {
        let p = *p as usize;
        let logp = p.ilog2() as u8 + 1;
        // sieve the medium primes into the buckets
        let r0 = *r0 as usize;
        let r1 = *r1 as usize;
        let mut idx0 = r0.min(r1);
        let mut idx1 = r0.max(r1);

        while idx1 < area {
            unsafe {
                hash_tables
                    .get_unchecked_mut(idx0 >> L1_SHIFT)
                    .push((idx0 & l1_mask, logp));
                hash_tables
                    .get_unchecked_mut(idx1 >> L1_SHIFT)
                    .push((idx1 & l1_mask, logp));
            }
            idx0 += p;
            idx1 += p;
        }
    }
    timers.record("line sieve - medium calculate", medium_start);

    for (block_i, start) in (0..area as isize)
        .step_by(l1_cache_size as usize)
        .enumerate()
    {
        let tiny_sieve_start = Instant::now();
        sieve.fill(0);
        for (p, s0, s1) in roots.iter().take(medium_prime_offset) {
            let p = *p as isize;
            let logp = p.ilog2() as u8 + 1;
            let mut i0 = *s0 as isize;
            i0 += (-start).rem_euclid(p);
            if i0 >= p {
                i0 -= p;
            }
            let mut i1 = *s1 as isize;
            i1 += (-start).rem_euclid(p);
            if i1 >= p {
                i1 -= p;
            }
            unsafe {
                line_sieve_double(
                    sieve.as_mut_ptr(),
                    i0.min(i1),
                    i0.max(i1),
                    logp,
                    p,
                    l1_cache_size as isize,
                );
            }
        }
        timers.record("line sieve - small", tiny_sieve_start);
        // now empty the bucket
        let medium_sieve_start = Instant::now();
        let s = unsafe { sieve.as_mut_ptr().offset(start) };
        for (idx, logp) in hash_tables[block_i].iter() {
            unsafe {
                //*sieve.get_unchecked_mut(start as usize + idx) += *logp;
                *s.offset(*idx as isize) += *logp;
            }
        }
        timers.record("line sieve - medium", medium_sieve_start);

        let scan_start = Instant::now();
        for j in (0..l1_cache_size as usize).step_by(recalculate_every) {
            let s = Int::from(j as isize + start + sqrt_num.clone())
                * Int::from(j as isize + start + sqrt_num.clone())
                - num.clone();
            let cutoff = s.to_f64().log2().ceil() as u8 - fudge;
            for i in j..((j + recalculate_every).min(l1_cache_size as usize)) {
                if sieve[i] > cutoff {
                    smooth_points.push(i + start as usize)
                }
            }
        }
        timers.record("scan", scan_start);
    }
    println!("Found {} potential points", smooth_points.len());
    smooth_points
}
