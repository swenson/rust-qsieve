use lib::opt::particle_swarm_optimize;
use lib::params::get_sieve_params;
use lib::primes::{get_n_primes, get_primes};
use num_traits::ToPrimitive;
use num_traits::*;
use rug::{Complete, Integer as Int, Integer};
use std::collections::{HashMap, HashSet};
use std::rc::Rc;
use std::time::Instant;

mod binary_matrix;
mod binary_matrix_simd;
mod dmatrix;
mod line_sieve;
mod multiple_poly;
mod poly;
mod primes;
mod smooth;

use binary_matrix::*;
use smooth::*;

use crate::latest::binary_matrix_simd::BinaryMatrixSimd;
use crate::latest::line_sieve::split_line_sieve;
use crate::latest::multiple_poly::generate_polys;
use crate::latest::poly::jacobi_bn;
use crate::latest::primes::sqrt_ceil;
use lib::timers::Timers;

pub fn main_fun(
    num_: &num_bigint::BigInt,
    timers: Rc<Timers>,
    use_simd: bool,
    use_asm: bool,
    la_lanes: u8,
    scan_lanes: u8,
    l1_cache_bits: u8,
    optimize: bool,
) {
    let num = Int::parse(num_.to_str_radix(10)).unwrap().complete();

    let params = get_sieve_params(num_.bits() as u32);
    println!("Params {:?}", params);
    let fb = get_n_primes(params.fb_size as u64);
    let smoothness_bound = *fb.last().unwrap();
    let large_mult = params.large_mult;
    let per_poly_area = params.sieve_size;
    let a_factors: usize = match num_.bits() {
        202.. => 6,
        78.. => 3,
        _ => 1,
    };
    if optimize {
        let threads = (std::thread::available_parallelism().unwrap().get() / 2).min(8);
        println!("Launching PSO with {} threads", threads);
        particle_swarm_optimize(
            run_factor,
            vec![1000.0, 65536.0, 1.0],
            vec![
                500_000.0.min(smoothness_bound as f64 * 10.0),
                10_000_000.0.min(per_poly_area as f64 * 10.0),
                1000.0.min(large_mult as f64 * 100.0),
            ],
            vec![
                smoothness_bound as f64,
                per_poly_area as f64,
                large_mult as f64,
            ],
            12,
            30,
            threads,
            &num,
            a_factors,
            use_simd,
            use_asm,
            la_lanes,
            scan_lanes,
            l1_cache_bits,
        );
    } else {
        run_factor(
            &num,
            &vec![
                smoothness_bound as f64,
                per_poly_area as f64,
                large_mult as f64,
            ],
            a_factors,
            use_simd,
            use_asm,
            la_lanes,
            scan_lanes,
            l1_cache_bits,
            timers,
        );
    }
}
struct Relations {
    raw: Vec<Relation>,
    filtered: Vec<usize>, // indices
    rel_set: HashSet<Integer>,
    n_set: HashSet<Integer>,
}

impl Relations {
    fn estimated_rows(&self, large_primes: &HashMap<u64, usize>) -> usize {
        let mut count = 0;
        for ri in self.filtered.iter() {
            let r: &Relation = &self.raw[*ri];
            if *large_primes
                .get(&(r.factors.last().unwrap_or(&0).abs() as u64))
                .unwrap_or(&0)
                == 1
            {
                // skip
            } else {
                count += 1;
            }
        }
        count
    }
}

impl Relations {
    fn add(&mut self, other: &mut Vec<Relation>) {
        let last = self.raw.len();
        self.raw.append(other);

        for i in last..self.raw.len() {
            let rel = &self.raw[i];
            let x = rel.num.clone();
            let n = rel.poly.eval(&x).abs();
            let xabs = rel.num.clone().abs();

            if self.rel_set.contains(&xabs) || self.n_set.contains(&n) {
            } else {
                self.rel_set.insert(xabs);
                self.n_set.insert(n.clone());
                self.filtered.push(i);
            }
        }
    }
}

fn run_factor(
    num: &Int,
    v: &Vec<f64>,
    a_factors: usize,
    use_simd: bool,
    use_asm: bool,
    la_lanes: u8,
    scan_lanes: u8,
    l1_cache_bits: u8,
    timers: Rc<Timers>,
) -> f64 {
    println!("Running with parameters {:?}, simd {}", v, use_simd);
    let smoothness_bound = v[0].round() as u64;
    let per_poly_area = v[1].round() as u64;
    let large_mult = v[2];
    let start = Instant::now();

    println!(
        "per poly area = 2^{}",
        per_poly_area.to_f64().unwrap().log2()
    );

    let large_prime_bound = ((smoothness_bound as f64) * large_mult).round() as u64;
    println!("Large prime bound {}", large_prime_bound);

    let primes: Rc<Vec<u64>> = Rc::new(
        get_primes(smoothness_bound)
            .into_iter()
            .filter(|p| jacobi_bn(&num, *p) == 1)
            .collect(),
    );
    let tiny_primes = TinyPrimes::generate(num);
    let primes_except_tiny: Vec<u64> = get_primes(smoothness_bound)
        .into_iter()
        .filter(|p| jacobi_bn(&num, *p) == 1)
        .filter(|p| p > tiny_primes.primes.last().unwrap_or(&0))
        .collect();
    let prime_set: HashSet<u64> = primes.iter().map(|x| *x).collect();

    let log_scale = 1.0;

    let mut rels = Relations {
        raw: vec![],
        filtered: vec![],
        rel_set: HashSet::new(),
        n_set: HashSet::new(),
    };
    let mut relations = vec![];
    let max_area = 8 * 1024 * 1024;
    let mut large_primes: HashMap<u64, usize> = HashMap::new();
    let sqrt_num = &sqrt_ceil(num.clone());
    let split = per_poly_area.div_ceil(max_area);
    println!("Splitting into {} areas", split);
    // round up to a multiple of 64
    let leftover = (per_poly_area / split) % 64;
    let ppas = ((per_poly_area / split) + (64 - leftover) % 64) as usize;
    println!("ppas = {}", ppas);
    let l1_cache_size = 1 << l1_cache_bits;

    let sieve = unsafe { &mut Box::new_uninit_slice(l1_cache_size as usize).assume_init() };
    println!(
        "Splitting into {} blocks per polynomial",
        ppas.div_ceil(l1_cache_size as usize)
    );

    let scanner = if use_simd {
        match scan_lanes {
            1 => scan_simd_specific::<1>,
            2 => scan_simd_specific::<2>,
            4 => scan_simd_specific::<4>,
            8 => scan_simd_specific::<8>,
            16 => scan_simd_specific::<16>,
            32 => scan_simd_specific::<32>,
            64 => scan_simd_specific::<64>,
            _ => panic!("Invalid scan_lanes: {}", scan_lanes),
        }
    } else if use_asm {
        scan_asm
    } else {
        scan
    };

    let sieve_fun = match l1_cache_bits {
        10 => split_line_sieve::<10>,
        11 => split_line_sieve::<11>,
        12 => split_line_sieve::<12>,
        13 => split_line_sieve::<13>,
        14 => split_line_sieve::<14>,
        15 => split_line_sieve::<15>,
        16 => split_line_sieve::<16>,
        17 => split_line_sieve::<17>,
        18 => split_line_sieve::<18>,
        19 => split_line_sieve::<19>,
        20 => split_line_sieve::<20>,
        21 => split_line_sieve::<21>,
        22 => split_line_sieve::<22>,
        23 => split_line_sieve::<23>,
        24 => split_line_sieve::<24>,
        _ => panic!("Unsupported L1 cache bits {}", l1_cache_bits),
    };

    let mut poly_num = 0;
    let mut last_estimated_rows = 0;
    let mut last_poly_estimated_rows_change = 0;
    let mut last_estimated_cols = 0;
    for poly_batch in generate_polys(&num, primes.clone(), timers.clone(), a_factors) {
        for poly in poly_batch.polys.iter() {
            let area = ppas;
            let a = poly.a.clone();
            let b = poly.b.clone();
            let c = poly.c.clone();
            let center = if a.is_one() {
                sqrt_num.clone()
            } else {
                (sqrt_num - b.as_ref()).complete() / a.as_ref()
            };
            for o in 0..split {
                let offset = -(per_poly_area as isize) / 2 + (ppas as u64 * o) as isize;
                let smooth_points = sieve_fun(
                    &center,
                    offset,
                    area,
                    primes.clone(),
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    sieve,
                    &timers,
                    poly.clone(),
                    large_prime_bound,
                    log_scale,
                    poly_batch.a_inverses.clone(),
                    poly_batch.sqrt_ns.clone(),
                    scanner,
                );
                let (mut new_relations, new_large_primes) = smooth_check(
                    smooth_points.as_ref(),
                    poly.clone(),
                    &center,
                    offset,
                    &primes_except_tiny,
                    &tiny_primes,
                    &prime_set,
                    large_prime_bound,
                    &timers,
                    poly_batch.a_factors.clone(),
                );
                new_large_primes.iter().for_each(|x| {
                    large_primes.insert(*x, 1 + *(large_primes.get(x).unwrap_or(&0usize)));
                });
                rels.add(&mut new_relations);
            }
            poly_num += 1;
        }
        let estimated_cols = primes.len() + large_primes.values().filter(|x| **x > 1).count();
        let estimated_rows = rels.estimated_rows(&large_primes);
        // println!(
        //     "Poly {}, estimated rows {}, estimated cols {}",
        //     poly_num, estimated_rows, estimated_cols,
        // );
        if estimated_rows > estimated_cols {
            break;
        }
        if estimated_rows != last_estimated_rows {
            if last_poly_estimated_rows_change > 0 {
                if poly_num > 1000 && poly_num > last_poly_estimated_rows_change + 2 {
                    if estimated_rows - last_estimated_rows
                        < (estimated_cols - last_estimated_cols) / 2
                    {
                        println!("Making negative progress; aborting");
                        return f64::INFINITY;
                    }
                }
            }
            last_estimated_rows = estimated_rows;
            last_poly_estimated_rows_change = poly_num;
            last_estimated_cols = estimated_cols;
        } else if poly_num > last_poly_estimated_rows_change + 10000 {
            println!("No progress after 10000 polys; aborting");
            return f64::INFINITY;
        }
    }

    for r in 0..rels.filtered.len() {
        relations.push(rels.raw[rels.filtered[r]].clone());
    }
    println!("Scanned {} polynomials", poly_num);

    let dupe_time = Instant::now();
    println!("All relations (including large primes) {}", relations.len());
    let all_relations = relations.len();
    relations = relations
        .into_iter()
        .filter(|rel| {
            *large_primes
                .get(&(*rel.factors.last().unwrap_or(&0) as u64))
                .unwrap_or(&0)
                != 1
        })
        .collect();
    println!(
        "Useful relations from large primes {}",
        large_primes
            .len()
            .saturating_sub(all_relations - relations.len())
    );
    println!("Relations {}", relations.len());
    println!("Factor base size: {}", primes.len());
    if relations.len() > 100000 {
        println!("Too many relations; aborting");
        return f64::INFINITY;
    }
    timers.record("deduplication", dupe_time);

    let build_start = Instant::now();
    // build a binary matrix
    let cols = primes.len() + large_primes.len() + 1;

    let mut col_map: HashMap<i64, usize> = HashMap::new();
    for pi in 0..primes.len() {
        col_map.insert(primes[pi] as i64, pi);
    }
    let large_primes_vec: Vec<u64> = large_primes
        .iter()
        .filter(|k| *k.1 > 1)
        .map(|x| *x.0)
        .collect();
    for i in 0..large_primes_vec.len() {
        col_map.insert(large_primes_vec[i] as i64, i + primes.len());
    }
    col_map.insert(-1, cols - 1);

    let mut mat_cols: HashMap<usize, Box<Vec<usize>>> = HashMap::new();
    for r in 0..relations.len() {
        for p in relations[r].factors.iter() {
            if !col_map.contains_key(p) {
                continue;
                // println!("Unknown prime {}; aborting", p);
                // return f64::INFINITY;
            }
            let pi = col_map.get(p).unwrap();
            if !mat_cols.contains_key(pi) {
                mat_cols.insert(*pi, Box::new(vec![]));
            }
            let v = mat_cols.get_mut(pi).unwrap();
            if v.is_empty() {
                v.push(r);
            } else if *v.last().unwrap() == r {
                v.pop();
            } else {
                v.push(r);
            }
        }
    }
    let mut col_indices = mat_cols.keys().map(|x| *x).collect::<Vec<usize>>();
    col_indices.sort();
    let mut cci = 0;
    let mut mat_rows: HashMap<usize, Box<Vec<usize>>> = HashMap::new();
    for ci in 0..col_indices.len() {
        let cx = col_indices[ci];
        let rows = mat_cols.get(&cx).unwrap();
        if rows.is_empty() {
            continue;
        }
        for r in rows.iter() {
            if !mat_rows.contains_key(r) {
                mat_rows.insert(*r, Box::new(vec![]));
            }
            mat_rows.get_mut(r).unwrap().push(cci);
        }
        cci += 1;
    }
    let mut mat = if use_simd {
        match la_lanes {
            1 => BinaryMatrixSimd::<1>::new(),
            2 => BinaryMatrixSimd::<2>::new(),
            4 => BinaryMatrixSimd::<4>::new(),
            8 => BinaryMatrixSimd::<8>::new(),
            16 => BinaryMatrixSimd::<16>::new(),
            32 => BinaryMatrixSimd::<32>::new(),
            64 => BinaryMatrixSimd::<64>::new(),
            _ => panic!("Invalid la lanes {}", la_lanes),
        }
    } else {
        BinaryMatrix64::new()
        //Box::new(DMatrix::<u8>::zeros(0, 0))
    };
    // TODO: we could repeat this since this reduces the column weight
    let extra_rows = (relations.len() as isize - cci as isize - 10).max(0) as usize;
    println!("Removing {} extra rows", extra_rows);
    if extra_rows > 0 {
        let len = relations.len();
        relations = relations.into_iter().take(len - extra_rows).collect();

        // TODO: DRY code
        // TODO: could maybe use dancing links to speed this up
        mat_cols.clear();
        for r in 0..relations.len() {
            for p in relations[r].factors.iter() {
                if !col_map.contains_key(p) {
                    continue;
                    // println!("Unknown prime {}; aborting", p);
                    // return f64::INFINITY;
                }
                let pi = col_map.get(p).unwrap();
                if !mat_cols.contains_key(pi) {
                    mat_cols.insert(*pi, Box::new(vec![]));
                }
                let v = mat_cols.get_mut(pi).unwrap();
                if v.is_empty() {
                    v.push(r);
                } else if *v.last().unwrap() == r {
                    v.pop();
                } else {
                    v.push(r);
                }
            }
        }
        col_indices = mat_cols.keys().map(|x| *x).collect::<Vec<usize>>();
        col_indices.sort();
        cci = 0;
        mat_rows.clear();
        for ci in 0..col_indices.len() {
            let cx = col_indices[ci];
            let rows = mat_cols.get(&cx).unwrap();
            if rows.is_empty() {
                continue;
            }
            for r in rows.iter() {
                if !mat_rows.contains_key(r) {
                    mat_rows.insert(*r, Box::new(vec![]));
                }
                mat_rows.get_mut(r).unwrap().push(cci);
            }
            cci += 1;
        }
    }
    mat.insert_rows(relations.len(), cci);
    for r in 0..relations.len() {
        if !mat_rows.contains_key(&r) {
            println!("Unknown row {}; aborting", r);
            return f64::INFINITY;
        }
        for c in mat_rows.get(&r).unwrap().iter() {
            mat.set(r, *c, 1);
        }
    }
    timers.record("matrix build", build_start);

    if mat.nrows() > 100_000 || mat.nrows() < 10 {
        println!("Matrix too big or small {}; aborting", mat.nrows());
        return f64::INFINITY;
    }
    println!(
        "Factor base size after removing empty columns: {}",
        mat.ncols()
    );
    println!("Solving {}x{} matrix", mat.ncols(), mat.nrows());
    let la_start = Instant::now();
    let kernel = mat.left_kernel().unwrap();
    timers.record("left kernel", la_start);
    let sqrt_start = Instant::now();

    let mut factor_found = false;

    for ki in 0..kernel.len() {
        let result_vector = kernel.get(ki).unwrap();
        let check = mat.left_mul(result_vector).all_zero();
        if !check {
            panic!("Null vector was invalid");
        }
        println!("Null vector is valid");

        let mut a = Int::one();
        let mut b = Int::one();

        let mut all_factors = HashMap::new();

        for i in 0..mat.nrows() {
            if result_vector.get(i) == 0 {
                continue;
            }
            for f in relations[i].factors.iter() {
                all_factors.insert(f, all_factors.get(&f).unwrap_or(&0u32) + 1);
            }
            b = b * (&relations[i].num) % num;
        }

        for (f, e) in all_factors.iter() {
            assert!(e % 2 == 0);
            a = a * Int::from(**f).pow_mod(&Int::from(e / 2), &num).unwrap() % num;
        }

        if !(((&a * &a).complete() - (&b * &b).complete()) % num).is_zero() {
            println!("Something is very broken -- (a^2 - b^2) % n != 0");
            return f64::INFINITY;
        }

        let g1 = (a.clone() + b.clone()).gcd(&num);
        let g2 = (a - b).gcd(&num.clone());

        if !g1.is_one() && &g1 != num && (num % &g1).complete().is_zero() {
            println!("Factor found! {}", g1);
            factor_found = true;
        }
        if !g2.is_one() && &g2 != num && (num % &g2).complete().is_zero() {
            println!("Factor found! {}", &g2);
            factor_found = true;
        }
        if factor_found {
            break;
        } else {
            println!("No factor found :( trying another vector");
        }
    }
    timers.record("sqrt", sqrt_start);
    println!("");
    if factor_found {
        // let fitness = ((mat.ncols() * mat.nrows()) as f64).pow(3)
        //     + (full_area as f64) * (full_area as f64).ln();
        let fitness = start.elapsed().as_secs_f64();
        println!("Fitness = {}", fitness);
        return fitness;
    }
    println!("Fitness = {}", f64::INFINITY);
    return f64::INFINITY;
}
