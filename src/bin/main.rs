#![feature(portable_simd)]
#![feature(bigint_helper_methods)]
#![feature(int_roundings)]
#![feature(new_uninit)]
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
