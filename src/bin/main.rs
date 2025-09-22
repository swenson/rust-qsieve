#![feature(test)]
#![feature(portable_simd)]
#![feature(bigint_helper_methods)]
#![feature(int_roundings)]
#![feature(ptr_sub_ptr)]

use clap::{Parser, Subcommand};
use lib::timers::Timers;
use num_bigint::*;
use num_traits::ToPrimitive;
use rug::rand::RandState;
use rug::Integer;
use std::time::SystemTime;
use std::{fmt, rc::Rc};

mod latest;
mod v1;
mod v2;
mod v3;
mod v4;
mod v5;
mod v6;
mod v7;
mod v8;
mod v9;

const DEFAULT_LA_LANES: u8 = 64;
#[cfg(target_arch = "aarch64")]
const DEFAULT_SCAN_LANES: u8 = 16;
#[cfg(target_arch = "x86_64")]
const DEFAULT_SCAN_LANES: u8 = 64;

// M2 Max performance core has a 128 KB data cache
#[cfg(target_arch = "aarch64")]
pub const DEFAULT_L1_CACHE_BITS: u8 = 17;

// Threadripper performs best with this set very high. TODO: figure out why
#[cfg(target_arch = "x86_64")]
pub const DEFAULT_L1_CACHE_BITS: u8 = 23;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(short, long, default_value_t = 100)]
    size: u8,
    #[arg(long)]
    use_simd: bool,
    #[arg(long)]
    use_asm: bool,
    #[arg(long, default_value_t = DEFAULT_LA_LANES)]
    la_lanes: u8,
    #[arg(long, default_value_t = DEFAULT_SCAN_LANES)]
    scan_lanes: u8,
    #[arg(long, default_value_t = DEFAULT_L1_CACHE_BITS)]
    l1_cache_bits: u8,
    #[arg(short, long)]
    optimize: bool,
    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Subcommand, Debug)]
enum Command {
    V1,
    V2,
    V3,
    V4,
    V5,
    V6,
    V7,
    V8,
    V9,
    Latest,
}

impl fmt::Display for Command {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let x = match self {
            Command::V1 => "v1",
            Command::V2 => "v2",
            Command::V3 => "v3",
            Command::V4 => "v4",
            Command::V5 => "v5",
            Command::V6 => "v6",
            Command::V7 => "v7",
            Command::V8 => "v8",
            Command::V9 => "v9",
            Command::Latest => "latest",
        };
        write!(f, "{}", x)
    }
}

fn random_prime(bits: usize, rand: &mut RandState) -> Integer {
    let mut p = Integer::new();
    while p.significant_bits() != bits as u32 {
        p = Integer::from(Integer::random_bits(bits as u32, rand)).next_prime();
    }
    p
}

fn random_semiprime(bits: usize, rand: &mut RandState) -> Integer {
    let mut n = Integer::new();
    while n.significant_bits() != bits as u32 {
        let a = random_prime(bits / 2, rand);
        let b = random_prime(bits / 2, rand);
        n = a * b;
    }
    n
}

fn main() {
    let args = Cli::parse();

    let cmd = args.command.unwrap_or(Command::Latest);
    let main_fun = match cmd {
        Command::V1 => v1::main_fun,
        Command::V2 => v2::main_fun,
        Command::V3 => v3::main_fun,
        Command::V4 => v4::main_fun,
        Command::V5 => v5::main_fun,
        Command::V6 => v6::main_fun,
        Command::V7 => v7::main_fun,
        Command::V8 => v8::main_fun,
        Command::V9 => v9::main_fun,
        Command::Latest => latest::main_fun,
    };

    let mut rand = RandState::new();
    rand.seed(&Integer::from(
        SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .unwrap()
            .as_secs(),
    ));
    //let num_str = random_semiprime(args.size as usize, &mut rand).to_string();
    let num_str = if args.size == 160 {
        "1004244384729623680896301239440627501113842567031".to_string()
    } else {
        random_semiprime(args.size as usize, &mut rand).to_string()
    };
    let num = BigInt::parse_bytes(num_str.as_bytes(), 10).unwrap();

    println!(
        "Factoring {} ({} bits) with {}, use_simd? {}, la lanes, scan lanes {}, L1 cache bits {}",
        num,
        num.to_f64().unwrap().log2(),
        cmd,
        args.la_lanes,
        args.scan_lanes,
        args.l1_cache_bits,
    );

    let timers = Rc::new(Timers::new());
    main_fun(
        &num,
        timers.clone(),
        args.use_simd,
        args.use_asm,
        args.la_lanes,
        args.scan_lanes,
        args.l1_cache_bits,
        args.optimize,
    );
    if !args.optimize {
        timers.report();
        timers.report_latex();
    }
}

#[cfg(test)]
mod bench {
    extern crate test;

    use crate::latest::main_fun as latest_main_fun;
    use crate::v2::main_fun as v2_main_fun;
    use crate::v3::main_fun as v3_main_fun;
    use crate::v4::main_fun as v4_main_fun;
    use crate::v5::main_fun as v5_main_fun;
    use crate::v6::main_fun as v6_main_fun;
    use crate::v7::main_fun as v7_main_fun;
    use crate::v8::main_fun as v8_main_fun;
    use crate::v9::main_fun as v9_main_fun;
    use crate::{DEFAULT_L1_CACHE_BITS, DEFAULT_LA_LANES, DEFAULT_SCAN_LANES};
    use lib::timers::Timers;
    use num_bigint::BigInt;
    use std::rc::Rc;
    use std::str::FromStr;
    use test::Bencher;

    const N80: &str = "1122813492023266650422531";
    const N100: &str = "921440771183728043661609546101";
    const N120: &str = "1240013853791532781287252774232659019";
    const N140: &str = "1133807386503023029911374595596154922160401";
    const N160: &str = "1202648794555690927303943531657894025303077839469";
    const N180: &str = "1422757630709378770473285749241082264027413063256050623";
    const N200: &str = "1008584493700239131569706984984577648342145878959825008194317";
    const N220: &str = "1378171389487685613300642617577088401175427867275667544964725680831";
    const N240: &str = "1261735615599667314838352547765439669337684592314874195620466176142898337";

    // M2 Max:
    // test bench::bench_v1_80                                      ... bench: 2,453,125,041 ns/iter (+/- 156,475,088)
    // #[bench]
    // fn bench_v1_80(b: &mut Bencher) {
    //     let num = BigInt::from_str(N80).unwrap();
    //     let timers = Rc::new(Timers::new());
    //     b.iter(|| {
    //         test::black_box(v1_main_fun(&num, timers.clone(), false, false, 1, 1, 17, false));
    //     });
    // }

    #[bench]
    fn bench_v2_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v2_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v3_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v3_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    // #[bench]
    // fn bench_v4_80(b: &mut Bencher) {
    //     let num = BigInt::from_str(N80).unwrap();
    //     let timers = Rc::new(Timers::new());
    //     b.iter(|| {
    //         test::black_box(v4_main_fun(&num, timers.clone(), false, false, 1, 1, DEFAULT_L1_CACHE_BITS, false));
    //     });
    // }

    #[bench]
    fn bench_v5_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v5_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v6_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v6_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v7_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v7_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v8_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v8_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v9_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v9_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_80(b: &mut Bencher) {
        let num = BigInt::from_str(N80).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    // From M2 Max:
    // test bench::bench_v1_100                                     ... bench: 3,851,756,054 ns/iter (+/- 182,516,831)
    // #[bench]
    // fn bench_v1_100(b: &mut Bencher) {
    //     let num = BigInt::from_str(N100).unwrap();
    //     let timers = Rc::new(Timers::new());
    //     b.iter(|| {
    //         test::black_box(v1_main_fun(&num, timers.clone(), false, false, 1, 1, 17, false));
    //     });
    // }

    #[bench]
    fn bench_v2_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v2_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v3_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v3_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v4_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v4_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v5_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v5_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v6_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v6_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v7_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v7_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v8_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v8_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v9_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v9_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_100(b: &mut Bencher) {
        let num = BigInt::from_str(N100).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                64,
                16,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v2_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v2_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v3_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v3_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v4_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v4_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v5_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v5_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v5_140(b: &mut Bencher) {
        let num = BigInt::from_str(N140).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v5_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v6_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v6_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v6_140(b: &mut Bencher) {
        let num = BigInt::from_str(N140).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v6_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v7_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v7_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v7_140(b: &mut Bencher) {
        let num = BigInt::from_str(N140).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v7_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v8_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v8_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v8_140(b: &mut Bencher) {
        let num = BigInt::from_str(N140).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v8_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v9_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v9_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v9_140(b: &mut Bencher) {
        let num = BigInt::from_str(N140).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v9_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_v9_160(b: &mut Bencher) {
        let num = BigInt::from_str(N160).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(v9_main_fun(
                &num,
                timers.clone(),
                false,
                false,
                1,
                1,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_120(b: &mut Bencher) {
        let num = BigInt::from_str(N120).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_140(b: &mut Bencher) {
        let num = BigInt::from_str(N140).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_160(b: &mut Bencher) {
        let num = BigInt::from_str(N160).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_180(b: &mut Bencher) {
        let num = BigInt::from_str(N180).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_200(b: &mut Bencher) {
        let num = BigInt::from_str(N200).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_220(b: &mut Bencher) {
        let num = BigInt::from_str(N220).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }

    #[bench]
    fn bench_latest_240(b: &mut Bencher) {
        let num = BigInt::from_str(N240).unwrap();
        let timers = Rc::new(Timers::new());
        b.iter(|| {
            test::black_box(latest_main_fun(
                &num,
                timers.clone(),
                true,
                true,
                DEFAULT_LA_LANES,
                DEFAULT_SCAN_LANES,
                DEFAULT_L1_CACHE_BITS,
                false,
            ));
        });
    }
}
