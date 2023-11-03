use crate::latest::multiple_poly::QuadraticPoly;
use crate::latest::poly;
use lib::timers::Timers;
use rug::{Complete, Integer};
use std::{collections::HashMap, rc::Rc, time::Instant};

fn compute_root_map(
    a: Rc<Integer>,
    b: Rc<Integer>,
    c: Rc<Integer>,
    primes: Rc<Vec<u64>>,
    a_inverses: Rc<HashMap<u32, u32>>,
    sqrt_ns: Rc<HashMap<u32, u32>>,
) -> Rc<(poly::Roots, Vec<poly::DoubleRoot>)> {
    let mut root_map = Vec::with_capacity(primes.len() / 2);
    let mut root2 = poly::Roots::Zero();
    let skip_first = if primes[0] == 2 {
        root2 = poly::roots_mod_p(
            &a, &b, &c, 2, 0, /* inverse doesn't matter for 2 */
            None,
        );
        1
    } else {
        0
    };
    for pi in skip_first..primes.len() {
        let p = primes[pi];
        let ainv = a_inverses.get(&(p as u32));
        if ainv.is_none() {
            continue;
        }
        let sqn = sqrt_ns.get(&(p as u32));
        if sqn.is_none() {
            continue;
        }
        let roots = poly::roots_mod_p(
            &a,
            &b,
            &c,
            p as i64,
            *ainv.unwrap() as i64,
            Some(*sqn.unwrap() as i64),
        );
        match roots {
            poly::Roots::Zero() => {}
            poly::Roots::One { .. } => {
                panic!("Expected more than 1 root");
            }
            poly::Roots::Two {
                r0,
                r1,
                prime,
                logp,
            } => {
                root_map.push(poly::DoubleRoot {
                    r0,
                    r1,
                    prime,
                    logp,
                });
            }
        }
    }
    Rc::new((root2, root_map))
}

fn two_fill(sieve: &mut [u8], roots: &poly::Roots, offset: isize, center: &Integer) {
    // special init of 2
    // let mut sieve2 =
    //     unsafe { &mut *slice_from_raw_parts_mut(sieve.as_mut_ptr() as *mut u16, sieve.len() / 2) };
    match roots {
        poly::Roots::Zero() => {
            // otherwise init with 0s if no roots
            // sieve2.fill(0);
            sieve.fill(0);
        }
        poly::Roots::One { r, .. } => {
            let mask = if (center.get_bit(0) as isize + offset - *r as isize) & 1 == 1 {
                1
            } else {
                0
            };
            let mut idx = 0;
            if mask == 0 {
                // sieve2.fill(2);
                for _ in 0..sieve.len() / 2 {
                    unsafe {
                        *sieve.get_unchecked_mut(idx) = 2;
                        *sieve.get_unchecked_mut(idx + 1) = 0;
                    }
                    idx += 2;
                }
            } else {
                // sieve2.fill(0x20);
                for _ in 0..sieve.len() / 2 {
                    unsafe {
                        *sieve.get_unchecked_mut(idx) = 0;
                        *sieve.get_unchecked_mut(idx + 1) = 2;
                    }
                    idx += 2;
                }
            }
        }
        poly::Roots::Two { .. } => {
            // sieve2.fill(0x22);
            sieve.fill(2);
        }
    }
}

// according to msieve, we want to split the factor base into tiny primes,
// whose updates are sieved into the block, and medium primes, whose updates
// are put into hash tables to be added when we sieve that block later.
// we start with the medium primes, then do per block.
pub fn split_line_sieve<const L1_SHIFT: u32>(
    center: &Integer,
    offset: isize,
    area: usize,
    primes: Rc<Vec<u64>>,
    a: Rc<Integer>,
    b: Rc<Integer>,
    c: Rc<Integer>,
    sieve: &mut [u8],
    timers: &Timers,
    poly: Rc<QuadraticPoly>,
    large_prime_bound: u64,
    log_scale: f64,
    a_inverses: Rc<HashMap<u32, u32>>,
    sqrt_ns: Rc<HashMap<u32, u32>>,
    scan_fun: fn(Rc<QuadraticPoly>, &Integer, &[u8], usize, u64, f64, &Timers) -> Vec<usize>,
) -> Vec<usize> {
    let l1_mask: u32 = (1 << L1_SHIFT) - 1;
    let l1_cache_size: u32 = 1 << L1_SHIFT;
    let fb_start = Instant::now();
    let rm = compute_root_map(
        a.clone(),
        b.clone(),
        c.clone(),
        primes.clone(),
        a_inverses,
        sqrt_ns,
    );
    let root2 = &rm.0;
    let root_map = &rm.1;
    timers.record("root map", fb_start);
    let mut sieve_start = Instant::now();
    let neg_offset_center = (-offset - center).complete();
    let mut med_prime_start = 0;
    for pi in 0..root_map.len() {
        let rm = &root_map[pi];
        let prime = rm.prime;
        if prime >= l1_cache_size / 2 {
            med_prime_start = pi;
            break;
        }
    }
    if med_prime_start == 0 {
        med_prime_start = root_map.len();
    }

    // compute offsets
    let mut offsets_r0 = Vec::with_capacity(root_map.len());
    offsets_r0.resize(root_map.len(), 0usize);
    let mut offsets_r1 = Vec::with_capacity(root_map.len());
    offsets_r1.resize(root_map.len(), 0usize);
    // TODO: these only need to be computed once
    // let mut block_size_mods = Vec::with_capacity(root_map.len());
    // block_size_mods.resize(root_map.len(), 0usize);
    // let mut block_size_mods_origs = Vec::with_capacity(root_map.len());
    // block_size_mods_origs.resize(root_map.len(), 0usize);

    for pi in 0..root_map.len() {
        let rm = &root_map[pi];
        let r0 = rm.r0;
        let r1 = rm.r1;
        let p = rm.prime;

        let sp = (&neg_offset_center).mod_u(p) as i64;
        // adjust for offset
        let mut idx0 = (r0 as i64 + sp) as usize;
        if idx0 >= p as usize {
            idx0 -= p as usize;
        }
        let mut idx1 = (r1 as i64 + sp) as usize;
        if idx1 >= p as usize {
            idx1 -= p as usize;
        }
        offsets_r0[pi] = idx0.min(idx1);
        offsets_r1[pi] = idx0.max(idx1);
    }

    let num_buckets = area.div_ceil(l1_cache_size as usize);
    let mut hash_tables = Vec::with_capacity(num_buckets);
    for _ in 0..num_buckets {
        hash_tables.push(vec![]);
    }

    timers.record("line sieve - offsets", sieve_start);
    sieve_start = Instant::now();

    // sieve the medium primes
    for pi in med_prime_start..root_map.len() {
        let rm = &root_map[pi];
        let p = rm.prime;
        let logp = rm.logp;
        let mut idx0 = offsets_r0[pi].min(offsets_r1[pi]);
        let mut idx1 = offsets_r0[pi].max(offsets_r1[pi]);
        let p_u = p as usize;
        while idx1 < area {
            // for _ in 0..(area.saturating_sub(idx1)).div_ceil(p_u) {
            unsafe {
                hash_tables
                    .get_unchecked_mut(idx0 >> L1_SHIFT)
                    .push((idx0 as u32 & l1_mask, logp));
                hash_tables
                    .get_unchecked_mut(idx1 >> L1_SHIFT)
                    .push((idx1 as u32 & l1_mask, logp));
            }
            idx0 += p_u;
            idx1 += p_u;
        }
        // guaranteed to only need at most one more update
        if idx0 < area {
            hash_tables[idx0 >> L1_SHIFT].push((idx0 as u32 & l1_mask, logp));
        }
        // // while idx0 < area {
        // // for unknown reasons, this is faster than using a while loop, or one loop
        // for _ in 0..(area.saturating_sub(idx0)).div_ceil(p_u) {
        //     unsafe {
        //         hash_tables
        //             .get_unchecked_mut(idx0 >> L1_SHIFT)
        //             .push((idx0 as u32 & L1_MASK, logp))
        //     };
        //     // println!(
        //     //     "actual bucket size = {}",
        //     //     hash_tables[idx0 >> L1_SHIFT].len()
        //     // );
        //     idx0 += p_u;
        // }
        // for _ in 0..(area.saturating_sub(idx1)).div_ceil(p_u) {
        //     unsafe {
        //         hash_tables
        //             .get_unchecked_mut(idx1 >> L1_SHIFT)
        //             .push((idx1 as u32 & L1_MASK, logp))
        //     };
        //     idx1 += p_u;
        // }
    }
    let mut smooth_points = vec![];
    timers.record("line sieve - med prime updates", sieve_start);

    // now sieve by blocks
    for block_start in (0..area).step_by(l1_cache_size as usize) {
        sieve_start = Instant::now();
        let end = (block_start + l1_cache_size as usize).min(area);
        if end <= block_start {
            break;
        }
        let block_size = end - block_start;
        // we can reset the sieve and sieve the 2s at the same time
        two_fill(sieve, root2, offset, center);
        timers.record("fill 2s", sieve_start);
        sieve_start = Instant::now();

        // do the tiny primes
        let ptr = &mut sieve[0] as *mut u8;
        for pi in 0..med_prime_start {
            let rm = &root_map[pi];
            let p_i = rm.prime as isize;
            let logp = root_map[pi].logp;
            let i0 = offsets_r0[pi] as isize;
            let i1 = offsets_r1[pi] as isize;
            let (new_i0, new_i1) = unsafe {
                tiny_prime_double_sieve::<L1_SHIFT>(ptr, i0, i1, logp, p_i, block_size as isize)
            };
            // let new_i0 = unsafe { tiny_prime_sieve(ptr, i0, logp, p_i, block_size as isize) };
            // let new_i1 = unsafe { tiny_prime_sieve(ptr, i1, logp, p_i, block_size as isize) };
            // indices come back ready for the next block
            offsets_r0[pi] = new_i0.min(new_i1) as usize;
            offsets_r1[pi] = new_i0.max(new_i1) as usize;
        }
        timers.record("line sieve - small", sieve_start);
        sieve_start = Instant::now();
        // now the medium prime buckets
        for x in hash_tables[block_start >> L1_SHIFT].iter() {
            unsafe {
                *sieve.get_unchecked_mut(x.0 as usize) =
                    (*sieve.get_unchecked_mut(x.0 as usize)) + x.1
            };
        }
        timers.record("line sieve - medium", sieve_start);
        //println!("scanning {:?} block {}", poly, block_start);
        smooth_points.append(&mut scan_fun(
            poly.clone(),
            &(center + (offset + block_start as isize)).complete(),
            sieve,
            block_start,
            large_prime_bound,
            log_scale,
            timers,
        ));
    }

    smooth_points
}

#[inline(always)]
#[allow(dead_code)]
pub unsafe fn tiny_prime_sieve<const L1_SHIFT: u32>(
    ptr: *mut u8,
    mut i: isize,
    logp: u8,
    p_i: isize,
    block_size: isize,
) -> isize {
    let l1_mask = (1 << L1_SHIFT) - 1;
    while i < block_size {
        *ptr.offset(i) += logp;
        i += p_i;
    }
    i & l1_mask as isize
}

#[inline(always)]
pub unsafe fn tiny_prime_double_sieve<const L1_SHIFT: u32>(
    ptr: *mut u8,
    mut i0: isize,
    mut i1: isize,
    logp: u8,
    p_i: isize,
    block_size: isize,
) -> (isize, isize) {
    let l1_mask = (1 << L1_SHIFT) - 1;
    let ptr0 = ptr.offset(i0 - i1);
    let init_i1 = i1;
    while i1 < block_size {
        *ptr0.offset(i1) += logp;
        *ptr.offset(i1) += logp;
        i1 += p_i;
    }
    i0 += i1 - init_i1;
    while i0 < block_size {
        *ptr.offset(i0) += logp;
        i0 += p_i;
    }
    (i0 & l1_mask as isize, i1 & l1_mask as isize)
    // let mut ptr0 = ptr.offset(i0);
    // let mut ptr1 = ptr.offset(i1);
    // let end = ptr.offset(block_size);
    // while ptr1.lt(&end) {
    //     *ptr0 += logp;
    //     *ptr1 += logp;
    //     ptr0 = ptr0.offset(p_i);
    //     ptr1 = ptr1.offset(p_i);
    // }
    // while ptr0.lt(&end) {
    //     *ptr0 += logp;
    //     ptr0 = ptr0.offset(p_i);
    // }
    // (
    //     ptr0.sub_ptr(ptr) as isize & l1_mask as isize,
    //     ptr1.sub_ptr(ptr) as isize & l1_mask as isize,
    // )
}

// #[cfg(not(target_arch = "aarch64"))]
// #[inline(always)]
// unsafe fn tiny_prime_sieve(
//     ptr: *mut u8,
//     mut i0: isize,
//     mut i1: isize,
//     logp: u8,
//     p_i: isize,
//     block_size: isize,
// ) {
//     // while i1 < block_size {
//     //     *(ptr.offset(i0)) = (*(ptr.offset(i0))).unchecked_add(logp);
//     //     *(ptr.offset(i1)) = (*(ptr.offset(i1))).unchecked_add(logp);
//     //     i0 = i0.unchecked_add(p_i);
//     //     i1 = i1.unchecked_add(p_i);
//     // }
//     // while i0 < block_size {
//     //     *(ptr.offset(i0)) = (*(ptr.offset(i0))).unchecked_add(logp);
//     //     i0 = i0.unchecked_add(p_i);
//     // }
//     let logp8 = logp as u64;

//     // println!(
//     //     "Starting assembly routine prime {} size {} i0 {} i1 {}",
//     //     p_i, block_size, i0, i1,
//     // );
//     use std::arch::asm;
//     unsafe {
//         asm!(
//             "2:
//             ; while i1 < block_size {{
//             cmp {i1}, {block_size}
//             b.gt 3f

//             5:
//             ; *(ptr.offset(i0)) = (*(ptr.offset(i0))).unchecked_add(logp);
//             ; *(ptr.offset(i1)) = (*(ptr.offset(i1))).unchecked_add(logp);
//             ldrb {t0:w}, [{i0}, {ptr}]
//             add {t0:w}, {t0:w}, {logp:w}
//             strb {t0:w}, [{i0}, {ptr}]

//             ldrb {t1:w}, [{i1}, {ptr}]
//             add {t1:w}, {t1:w}, {logp:w}
//             strb {t1:w}, [{i1}, {ptr}]

//             ; i0 = i0.unchecked_add(p_i);
//             ; i1 = i1.unchecked_add(p_i);
//             add {i1}, {i1}, {p_i}
//             add {i0}, {i0}, {p_i}
//             cmp {i1}, {block_size}
//             b.le 5b

//             3:
//             cmp {i0}, {block_size}
//             b.gt 4f
//             6:
//             ldrb {t0:w}, [{i0}, {ptr}]
//             add {t0:w}, {t0:w}, {logp:w}
//             strb {t0:w}, [{i0}, {ptr}]
//             add {i0}, {i0}, {p_i}
//             cmp {i0}, {block_size}
//             b.le 6b
//             4:",
//             i0 = inout(reg) i0,
//             i1 = inout(reg) i1,
//             ptr = in(reg) ptr,
//             t0 = out(reg) _,
//             t1 = out(reg) _,
//             logp = in(reg) logp8,
//             p_i = in(reg) p_i,
//             block_size = in(reg) block_size,
//             options(nostack),
//         );
//     }
//     // println!("Ending assembly routine");
// }

// #[allow(dead_code)]
// pub fn simple_line_sieve(
//     center: &Integer,
//     offset: isize,
//     area: usize,
//     primes: Rc<Vec<u64>>,
//     a: Rc<Integer>,
//     b: Rc<Integer>,
//     c: Rc<Integer>,
//     sieve: &mut [u8],
//     timers: &mut Timers,
//     poly: Rc<crate::multiple_poly::QuadraticPoly>,
//     large_prime_bound: u64,
//     log_scale: f64,
// ) -> Vec<usize> {
//     let fb_start = Instant::now();
//     let rm = compute_root_map(a.clone(), b.clone(), c.clone(), primes.clone());
//     let root2 = rm.0;
//     let root_map = rm.1;
//     timers.record("root map", fb_start);

//     let sieve_start = Instant::now();
//     let mut skip_first = 0;
//     // we can reset the sieve and sieve the 2s at the same time
//     if *primes.first().unwrap_or(&0) == 2 {
//         skip_first = 1;
//         two_fill(sieve, &root_map[0], offset, center);
//     } else {
//         sieve.fill(0);
//     }
//     timers.record("fill 2s", sieve_start);

//     let neg_offset_center = (-offset - center).complete();

//     for pi in skip_first..root_map.len() {
//         match root_map[pi] {
//             poly::Roots::Zero() => {}
//             poly::Roots::One { r, prime: p, logp } => {
//                 let p_i64 = p as i64;
//                 let p_u = p as usize;
//                 let mut sp = (&neg_offset_center % p).complete().to_i64().unwrap();
//                 if sp < 0 {
//                     sp += p_i64;
//                 }
//                 // adjust for offset
//                 let mut idx = (r as i64 + sp) as usize;
//                 if idx >= p as usize {
//                     idx -= p as usize;
//                 }
//                 for _ in 0..(area.saturating_sub(idx)) / p_u {
//                     unsafe {
//                         *sieve.get_unchecked_mut(idx) += logp;
//                     }
//                     idx += p_u;
//                 }
//             }
//             poly::Roots::Two {
//                 r0: r1,
//                 r1: r2,
//                 prime: p,
//                 logp,
//             } => {
//                 let p_i64 = p as i64;
//                 let p_u = p as usize;
//                 let mut sp = (&neg_offset_center % p).complete().to_i64().unwrap();
//                 if sp < 0 {
//                     sp += p_i64;
//                 }
//                 // adjust for offset
//                 let mut idx1 = (r1 as i64 + sp) as usize;
//                 if idx1 >= p as usize {
//                     idx1 -= p as usize;
//                 }
//                 let mut idx2 = (r2 as i64 + sp) as usize;
//                 if idx2 >= p as usize {
//                     idx2 -= p as usize;
//                 }
//                 // let mut idx1 = (*r1 as i64 + sp as i64).rem_euclid(p_i64) as usize;
//                 // let mut idx2 = (*r2 as i64 + sp as i64).rem_euclid(p_i64) as usize;
//                 //let mut iters = 0;
//                 // for _ in 0..(area.saturating_sub(idx1.max(idx2))) / p_u {
//                 //     unsafe {
//                 //         *sieve.get_unchecked_mut(idx1) += logp;
//                 //         *sieve.get_unchecked_mut(idx2) += logp;
//                 //     }
//                 //     idx1 += p_u;
//                 //     idx2 += p_u;
//                 //     iters += 1;
//                 // }

//                 unsafe {
//                     let mut i1 = idx1.min(idx2) as isize;
//                     let mut i2 = idx1.max(idx2) as isize;
//                     let p_i = p_u as isize;
//                     let ptr = &mut sieve[0] as *mut u8;
//                     while i2 < area as isize {
//                         *(ptr.offset(i1)) = (*(ptr.offset(i1))).unchecked_add(logp);
//                         *(ptr.offset(i2)) = (*(ptr.offset(i2))).unchecked_add(logp);
//                         i1 = i1.unchecked_add(p_i);
//                         i2 = i2.unchecked_add(p_i);
//                         //iters += 1;
//                     }
//                 }
//                 // unsafe {
//                 //     let mut idx = idx1 as isize;
//                 //     let ptr = &mut sieve[0] as *mut u8;
//                 //     let p_i = p_u as isize;
//                 //     while idx < area as isize {
//                 //         idx += p_i;
//                 //     }
//                 // }
//                 //                println!("Done {} with {:?}", p, poly);
//             }
//         }
//     }

//     timers.record("line sieve", sieve_start);
//     scan(
//         poly,
//         &(center + offset).complete(),
//         sieve,
//         0,
//         large_prime_bound,
//         log_scale,
//         timers,
//     )
// }
