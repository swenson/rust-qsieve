use crate::timers::Timers;
use crossbeam::thread;
use rug::Integer;
use std::rc::Rc;

pub fn particle_swarm_optimize(
    f: fn(&Integer, &Vec<f64>, usize, bool, bool, u8, u8, u8, Rc<Timers>) -> f64,
    min: Vec<f64>,
    max: Vec<f64>,
    init: Vec<f64>,
    num_particles: usize,
    max_iters: usize,
    max_threads: usize,
    num: &Integer,
    a_factors: usize,
    use_simd: bool,
    use_asm: bool,
    la_lanes: u8,
    scan_lanes: u8,
    l1_cache_bits: u8,
) -> (Vec<f64>, f64) {
    let w = 0.5;
    let phip = 1.0;
    let phig = 1.0;
    use rand::prelude::*;
    let mut rng = SmallRng::from_entropy();

    let mut particles = vec![];
    let mut velocities = vec![];

    let mut fitness = vec![];
    let mut best_f = vec![];
    let mut best_x: Vec<Vec<f64>> = vec![];
    let mut overall_best_f = f64::INFINITY;
    let mut overall_best_x = vec![];

    for i in 0..num_particles {
        let mut p = vec![];
        let mut v = vec![];
        for d in 0..min.len() {
            if i == 0 {
                p.push(init[d]);
            } else {
                p.push(rng.gen_range(min[d]..max[d]));
            }
            let x = (max[d] - min[d]).abs();
            v.push(rng.gen_range(-x..x));
        }
        velocities.push(v);
        particles.push(p.clone());
        fitness.push(f64::INFINITY);
        best_f.push(f64::INFINITY);
        best_x.push(p.clone());
    }

    thread::scope(|s| {
        let mut i = 0;
        while i < num_particles {
            let mut results = vec![];
            for j in 0..(max_threads.min(num_particles - i)) {
                let p = particles[i + j].clone();
                let x = s.spawn(move |_| {
                    f(
                        num,
                        &p,
                        a_factors,
                        use_simd,
                        use_asm,
                        la_lanes,
                        scan_lanes,
                        l1_cache_bits,
                        Rc::new(Timers::new()),
                    )
                });
                results.push(x);
            }
            for result in results {
                fitness[i] = result.join().unwrap_or(f64::INFINITY);
                i += 1;
            }
        }
    })
    .unwrap();
    for i in 0..num_particles {
        if fitness[i] < overall_best_f {
            overall_best_x = particles[i].clone();
            overall_best_f = fitness[i];
        }
        if fitness[i] < best_f[i] {
            best_x[i] = particles[i].clone();
            best_f[i] = fitness[i];
        }
    }
    //     fitness.push(f(&p));
    //     best_f.push(*fitness.last().unwrap());
    //     best_x.push(p.clone());
    //     if i > 0 {
    //         if *fitness.last().unwrap() < fitness[g] {
    //             g = i;
    //         }
    //     }
    // }

    for _ in 0..max_iters {
        for i in 0..num_particles {
            let v = &mut velocities[i];
            let l = particles[0].len();
            {
                let p = &particles[i];
                for d in 0..l {
                    let rp = rng.gen_range(0.0..1.0);
                    let rg = rng.gen_range(0.0..1.0);
                    v[d] = w * v[d]
                        + phip * rp * (best_x[i][d] - p[d])
                        + phig * rg * (overall_best_x[d] - p[d]);
                }
            }
            let p = &mut particles[i];
            for d in 0..l {
                p[d] += v[d];
                p[d] = p[d].max(min[d]).min(max[d]);
            }
        }

        thread::scope(|s| {
            let mut i = 0;
            while i < num_particles {
                let mut results = vec![];
                for j in 0..(max_threads.min(num_particles - i)) {
                    let p = particles[i + j].clone();
                    let x = s.spawn(move |_| {
                        f(
                            num,
                            &p,
                            a_factors,
                            use_simd,
                            use_asm,
                            la_lanes,
                            scan_lanes,
                            l1_cache_bits,
                            Rc::new(Timers::new()),
                        )
                    });
                    results.push(x);
                }
                for result in results {
                    fitness[i] = result.join().unwrap_or(f64::INFINITY);
                    i += 1;
                }
            }
        })
        .unwrap();

        for i in 0..num_particles {
            if fitness[i] < overall_best_f {
                overall_best_x = particles[i].clone();
                overall_best_f = fitness[i];
            }
            if fitness[i] < best_f[i] {
                best_x[i] = particles[i].clone();
                best_f[i] = fitness[i];
            }
        }

        // for i in 0..particles.len() {
        //     let v = &mut velocities[i];
        //     let l = particles[0].len();
        //     {
        //         let p = &particles[i];
        //         for d in 0..l {
        //             let rp = rng.gen_range(0.0..1.0);
        //             let rg = rng.gen_range(0.0..1.0);
        //             v[d] = w * v[d]
        //                 + phip * rp * (best_x[i][d] - p[d])
        //                 + phig * rg * (particles[g][d] - p[d]);
        //         }
        //     }
        //     let p = &mut particles[i];
        //     for d in 0..l {
        //         p[d] += v[d];
        //         p[d] = p[d].max(min[d]).min(max[d]);
        //     }
        //     fitness[i] = f(p);
        //     if fitness[i] < fitness[g] {
        //         g = i;
        //     }
        //     if fitness[i] < best_f[i] {
        //         best_x[i] = p.clone();
        //         best_f[i] = fitness[i];
        //     }
        // }
        println!("Best (overall) {:?} -> {}", overall_best_x, overall_best_f)
    }

    (overall_best_x, overall_best_f)
}
