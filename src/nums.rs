#[inline(always)]
pub fn mulmod32(a: u32, b: u32, p: u32) -> u32 {
    ((a as u64) * (b as u64) % p as u64) as u32
}
#[inline(always)]
pub fn mulmod64(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128) * (b as u128) % p as u128) as u64
}

pub fn modpow64(a: u64, e: u64, p: u64) -> u64 {
    let mut result = 1;

    let mut pow = a;
    let mut ee = e;

    while ee != 0 {
        if ee & 1 == 1 {
            result = mulmod64(result, pow, p);
        }
        pow = mulmod64(pow, pow, p);
        ee >>= 1;
    }
    return result;
}

pub fn modpow32(a: u32, e: u32, p: u32) -> u32 {
    let mut result = 1;

    let mut pow = a;
    let mut ee = e;

    while ee != 0 {
        if ee & 1 == 1 {
            result = ((result as u64) * (pow as u64) % (p as u64)) as u32;
        }
        pow = ((pow as u64) * (pow as u64) % (p as u64)) as u32;
        ee >>= 1;
    }
    return result;
}

#[cfg(test)]
mod bench {
    extern crate test;

    use crate::nums::modpow32;
    use crate::nums::modpow64;
    use test::Bencher;

    #[bench]
    fn bench_modpow32(b: &mut Bencher) {
        b.iter(|| test::black_box(modpow32(1234088719, 4109293439, 3989641949)))
    }

    #[bench]
    fn bench_modpow64(b: &mut Bencher) {
        b.iter(|| {
            test::black_box(modpow64(
                12422686881634578097,
                2245559321401646113,
                13729269479923652819,
            ))
        })
    }
}
