use std::{cell::RefCell, collections::HashMap, time::Instant};

pub struct Timers {
    timers: RefCell<HashMap<String, (u64, f64)>>,
}

impl Timers {
    pub fn new() -> Timers {
        Timers {
            timers: RefCell::new(HashMap::new()),
        }
    }

    pub fn record(&self, name: &str, start: Instant) {
        self.timers
            .borrow_mut()
            .entry(String::from(name))
            .and_modify(|f| {
                let (a, b) = *f;
                *f = (a + 1, b + start.elapsed().as_secs_f64())
            })
            .or_insert((1, start.elapsed().as_secs_f64()));
    }

    pub fn report(&self) {
        let mut keys: Vec<String> = self
            .timers
            .borrow()
            .keys()
            .map(|x| String::from(x))
            .collect();
        keys.sort();
        let total: f64 = self.timers.borrow().values().map(|f| f.1).sum();
        for k in keys.iter() {
            let v = *self.timers.borrow().get(k).unwrap();
            println!(
                "{:32} time: {:10.6} seconds ({:5.2}%)  count {:7}",
                k,
                v.1,
                100.0 * v.1 / total,
                v.0,
            );
        }
        println!("--> total: {:.6} seconds", total);
    }

    pub fn report_latex(&self) {
        let mut keys: Vec<String> = self
            .timers
            .borrow()
            .keys()
            .map(|x| String::from(x))
            .collect();
        keys.sort();
        let total: f64 = self.timers.borrow().values().map(|f| f.1).sum();
        println!("\\begin{{tabular}}{{lrr}}");
        println!("\\toprule");
        for k in keys.iter() {
            let v = *self.timers.borrow().get(k).unwrap();
            println!("{} & {:.3}s & {:5.2}\\% \\\\", k, v.1, 100.0 * v.1 / total,);
        }
        println!("\\midrule");
        println!("total & {:.3}s & \\\\", total);
        println!("\\bottomrule");
        println!("\\end{{tabular}}")
    }
}
