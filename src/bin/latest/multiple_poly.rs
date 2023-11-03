use itertools::{Itertools, MultiPeek};
use num_traits::{One, Zero};
use rug::{Complete, Integer};
use std::{collections::HashMap, rc::Rc, time::Instant};

use crate::latest::poly::{
    jacobi_bn, mod_inverse, roots_mod_p2, square_roots_mod_n, square_roots_modp, Roots,
};
use lib::{primes::get_all_primes, timers::Timers};

#[derive(Clone, Debug, PartialEq)]
pub struct QuadraticPoly {
    pub a: Rc<Integer>,
    pub b: Rc<Integer>,
    pub c: Rc<Integer>,
}

#[derive(Debug)]
pub struct PolyBatch {
    pub a_factors: Rc<Vec<u32>>,
    pub a_inverses: Rc<HashMap<u32, u32>>,
    pub sqrt_ns: Rc<HashMap<u32, u32>>,
    pub polys: Vec<Rc<QuadraticPoly>>,
}

impl QuadraticPoly {
    pub fn new(a: &Integer, b: &Integer, c: &Integer) -> Rc<QuadraticPoly> {
        Rc::new(QuadraticPoly {
            a: Rc::new(a.clone()),
            b: Rc::new(b.clone()),
            c: Rc::new(c.clone()),
        })
    }

    #[allow(dead_code)]
    pub fn new_i64(a: i64, b: i64, c: i64) -> Rc<QuadraticPoly> {
        Rc::new(QuadraticPoly {
            a: Rc::new(Integer::from(a)),
            b: Rc::new(Integer::from(b)),
            c: Rc::new(Integer::from(c)),
        })
    }

    #[allow(dead_code)]
    pub fn sqrt_disc(&self) -> Rc<Integer> {
        Rc::new(
            (self.b.as_ref() * self.b.as_ref()
                - ((self.a.as_ref() * self.c.as_ref()).complete() << 2) as Integer)
                .sqrt(),
        )
    }

    pub fn eval(&self, x: &Integer) -> Integer {
        (self.a.as_ref() * x + self.b.as_ref()).complete() * x + self.c.as_ref()
    }
}

pub struct MultiplePolys {
    factor_base: Rc<Vec<u64>>,
    num: Integer,
    curr: usize,
    primes: MultiPeek<Box<dyn Iterator<Item = u64>>>,
    timers: Rc<Timers>,
    sqrt_ns: Rc<HashMap<u32, u32>>,
    num_a_factors: usize,
}

// #[cfg(test)]
// mod test {
//     use super::*;

//     #[test]
//     fn generate_polys_test() {
//         let g: Vec<Rc<QuadraticPoly>> = generate_polys(
//             &Integer::from(15347),
//             4,
//             Rc::new(vec![2, 3, 5, 7, 11]),
//             Rc::new(Timers::new()),
//         )
//         .into_iter()
//         .flat_map(|x| {
//             x.polys
//                 .iter()
//                 .map(|y| y.to_owned())
//                 .collect::<Vec<Rc<QuadraticPoly>>>()
//         })
//         .collect();
//         assert_eq!(
//             g[1..],
//             vec![
//                 QuadraticPoly::new_i64(289, 84, -47),
//                 QuadraticPoly::new_i64(529, 252, 1),
//                 QuadraticPoly::new_i64(841, 770, 158),
//             ]
//         );
//     }
// }

pub fn generate_polys(
    num: &Integer,
    factor_base: Rc<Vec<u64>>,
    timers: Rc<Timers>,
    num_a_factors: usize,
) -> MultiplePolys {
    let n = num.clone();
    let n2 = num.clone();
    let mut sqrt_ns = HashMap::new();
    for q in factor_base.iter() {
        sqrt_ns.insert(
            *q as u32,
            square_roots_modp((num % q).complete().to_u32().unwrap(), *q as u32)[0] as u32,
        );
    }
    let iter: Box<dyn Iterator<Item = u64>> = Box::new(
        get_all_primes()
            .filter(move |x| jacobi_bn(&n, *x) == 1)
            .filter(move |p| {
                let rmp2 = roots_mod_p2(&Integer::one(), &Integer::zero(), &n2, *p as i64);
                rmp2 != Roots::Zero()
            })
            .filter(|p| *p != 2),
    );
    MultiplePolys {
        factor_base,
        num: num.clone(),
        curr: 0,
        primes: iter.multipeek(),
        timers,
        sqrt_ns: Rc::new(sqrt_ns),
        num_a_factors,
    }
}

impl MultiplePolys {
    fn initial_poly(&mut self) -> Rc<PolyBatch> {
        let mut a_inverses = HashMap::new();
        for z in self.factor_base.iter() {
            a_inverses.insert(*z as u32, 1);
        }
        Rc::new(PolyBatch {
            a_factors: Rc::new(vec![]),
            a_inverses: Rc::new(a_inverses),
            sqrt_ns: self.sqrt_ns.clone(),
            polys: vec![QuadraticPoly::new(
                &Integer::one(),
                &Integer::zero(),
                &-self.num.clone(),
            )],
        })
    }

    fn compute_a_inverses(&mut self) -> (Integer, HashMap<u32, u32>, Vec<u32>) {
        let mut p_factors = vec![];
        p_factors.push(self.primes.next().unwrap() as u32);
        let mut a = Integer::from(p_factors[0]);
        for _ in 0..self.num_a_factors - 1 {
            let q = self.primes.peek().unwrap();
            p_factors.push(*q as u32);
            a = a * *q;
        }
        let mut a_inverses = HashMap::new();
        // precompute a_inverses
        for z in self.factor_base.iter() {
            if p_factors.contains(&(*z as u32)) {
                continue;
            }
            a_inverses.insert(*z as u32, mod_inverse(&a, *z as i64).unwrap() as u32);
        }

        (a, a_inverses, p_factors)
    }
}

impl Iterator for MultiplePolys {
    type Item = Rc<PolyBatch>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr == 0 {
            self.curr += 1;
            return Some(self.initial_poly());
        }
        let (a, a_inverses, p_factors) = self.compute_a_inverses();
        let start = Instant::now();
        // solve b^2 = num (mod a)
        let a_2: &Integer = &(&a >> 1u32).complete();
        let b_s: Vec<Integer> = square_roots_mod_n(&self.num, p_factors.as_ref())
            .into_iter()
            .filter(|b| *b < *a_2)
            .collect();

        // solve b^2 - a * c = num
        let polys: Vec<Rc<QuadraticPoly>> = b_s
            .iter()
            .map(|b| {
                //assert!(((b * b - &self.num).complete() % &a) == 0);
                let c = (b * b - &self.num).complete() / &a;
                QuadraticPoly::new(&a, &(b << 1u32).complete(), &c)
            })
            .collect();
        if self.curr == 1 {
            println!("Processing polynomials in blocks of {}", polys.len());
        }
        self.curr += polys.len();
        self.timers.record("poly generate", start);
        Some(Rc::new(PolyBatch {
            a_inverses: Rc::new(a_inverses),
            sqrt_ns: self.sqrt_ns.clone(),
            a_factors: Rc::new(p_factors),
            polys,
        }))
    }
}
