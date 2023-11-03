// adapted from Jason Papadopoulos's msieve

#[derive(Copy, Clone, Debug)]
pub struct SieveParam {
    pub bits: u32,       /* size of integer this config info applies to */
    pub fb_size: u32,    /* number of factor base primes */
    pub large_mult: u32, /* the large prime multiplier */
    pub sieve_size: u32, /* the size of the sieve (actual sieve is 2x this) */
}

pub fn get_msieve_params(bits: u32) -> SieveParam {
    if bits <= MSIEVE_PREBUILT_PARAMS[0].bits {
        return MSIEVE_PREBUILT_PARAMS[0];
    }

    if bits >= MSIEVE_PREBUILT_PARAMS[MSIEVE_PREBUILT_PARAMS.len() - 1].bits {
        return MSIEVE_PREBUILT_PARAMS[MSIEVE_PREBUILT_PARAMS.len() - 1];
    }

    let mut i = 0;
    for j in 0..MSIEVE_PREBUILT_PARAMS.len() {
        if MSIEVE_PREBUILT_PARAMS[j].bits == bits {
            return MSIEVE_PREBUILT_PARAMS[j];
        }
        if MSIEVE_PREBUILT_PARAMS[j].bits > bits {
            i = j - 1;
            break;
        }
    }

    let low = &MSIEVE_PREBUILT_PARAMS[i];
    let high = &MSIEVE_PREBUILT_PARAMS[i + 1];

    SieveParam {
        bits,
        fb_size: weighted_sum(low.bits, bits, high.bits, low.fb_size, high.fb_size),
        large_mult: weighted_sum(low.bits, bits, high.bits, low.large_mult, high.large_mult),
        sieve_size: weighted_sum(low.bits, bits, high.bits, low.sieve_size, high.sieve_size),
    }
}

pub fn get_sieve_params(bits: u32) -> SieveParam {
    if bits <= OUR_PREBUILT_PARAMS[0].bits {
        return OUR_PREBUILT_PARAMS[0];
    }

    if bits >= OUR_PREBUILT_PARAMS[OUR_PREBUILT_PARAMS.len() - 1].bits {
        return OUR_PREBUILT_PARAMS[OUR_PREBUILT_PARAMS.len() - 1];
    }

    let mut i = 0;
    for j in 0..OUR_PREBUILT_PARAMS.len() {
        if OUR_PREBUILT_PARAMS[j].bits == bits {
            return OUR_PREBUILT_PARAMS[j];
        }
        if OUR_PREBUILT_PARAMS[j].bits > bits {
            i = j - 1;
            break;
        }
    }

    let low = &OUR_PREBUILT_PARAMS[i];
    let high = &OUR_PREBUILT_PARAMS[i + 1];

    SieveParam {
        bits,
        fb_size: weighted_sum(low.bits, bits, high.bits, low.fb_size, high.fb_size),
        large_mult: weighted_sum(low.bits, bits, high.bits, low.large_mult, high.large_mult),
        sieve_size: weighted_sum(low.bits, bits, high.bits, low.sieve_size, high.sieve_size),
    }
}

fn weighted_sum(low: u32, curr: u32, high: u32, a: u32, b: u32) -> u32 {
    let i = (curr - low) as f64;
    let j = (high - curr) as f64;
    let dist = (high - low) as f64;
    ((a as f64 * j + b as f64 * i) / dist + 0.5).round() as u32
}

const OUR_PREBUILT_PARAMS: [SieveParam; 11] = [
    SieveParam {
        bits: 40,
        fb_size: 185,
        large_mult: 1,
        sieve_size: 65536,
    },
    SieveParam {
        bits: 60,
        fb_size: 185,
        large_mult: 1,
        sieve_size: 65536,
    },
    SieveParam {
        bits: 80,
        fb_size: 185,
        large_mult: 1,
        sieve_size: 65536,
    },
    SieveParam {
        bits: 100,
        fb_size: 590,
        large_mult: 4,
        sieve_size: 2 * 65536,
    },
    SieveParam {
        bits: 120,
        fb_size: 717,
        large_mult: 74,
        sieve_size: 315_180,
    },
    SieveParam {
        bits: 140,
        fb_size: 1961,
        large_mult: 17,
        sieve_size: 636_045,
    },
    SieveParam {
        bits: 160,
        fb_size: 1900,
        large_mult: 213,
        sieve_size: 351_793,
    },
    SieveParam {
        bits: 180,
        fb_size: 5197,
        large_mult: 50,
        sieve_size: 881_126,
        // fb_size: 6036,
        // large_mult: 59,
        // sieve_size: 615_526,
    },
    SieveParam {
        bits: 200,
        fb_size: 10819,
        large_mult: 338,
        sieve_size: 6_236_467,
    },
    SieveParam {
        bits: 220,
        fb_size: 14859,
        large_mult: 403,
        sieve_size: 4_963_353,
        // fb_size: 12997,
        // large_mult: 187,
        // sieve_size: 1 << 20,
    },
    SieveParam {
        bits: 240,
        // fb_size: 27668,
        // large_mult: 939,
        // sieve_size: 2794129,
        fb_size: 29650,
        large_mult: 151,
        sieve_size: 8_175_424,
        //fb_size: 20000,
        //large_mult: 600,
        // sieve_size: 8194304, // 156.937474333
        //sieve_size: 8291456, // 206.637348
        //sieve_size: 8340032, // 139.151366
        //sieve_size: 8364320, // 224.538690
        // sieve_size: 8376464, // 144.869937
        //sieve_size: 8382536, // 167.955816
        //sieve_size: 8385572, // 114.574321
        //sieve_size: 8387090, // 167.691650
        //sieve_size: 8387849, // 191.777508
        //sieve_size: 8388229, // 139.421019
        //sieve_size: 8388419, // 177.529210
        //sieve_size: 8388514, // 130.616710
        //sieve_size: 8388561, // 86.800616
        // sieve_size: 8388585, // broken
        // sieve_size: 8388573, // 243.056202
        //sieve_size: 8388579, // 182.563095
        //sieve_size: 8388582, // 156.192051
        //sieve_size: 8388583, // broken
        // TODO: figure out why 1 << 24 fails exactly
        // has something to do with the max_area calculation almost surely
    },
];

const MSIEVE_PREBUILT_PARAMS: [SieveParam; 21] = [
    SieveParam {
        bits: 64,
        fb_size: 100,
        large_mult: 40,
        sieve_size: 1 * 65536,
    },
    SieveParam {
        bits: 128,
        fb_size: 450,
        large_mult: 40,
        sieve_size: 1 * 65536,
    },
    SieveParam {
        bits: 183,
        fb_size: 2000,
        large_mult: 40,
        sieve_size: 1 * 65536,
    },
    SieveParam {
        bits: 200,
        fb_size: 3000,
        large_mult: 50,
        sieve_size: 1 * 65536,
    },
    SieveParam {
        bits: 212,
        fb_size: 5400,
        large_mult: 50,
        sieve_size: 3 * 65536,
    },
    SieveParam {
        bits: 233,
        fb_size: 10000,
        large_mult: 100,
        sieve_size: 3 * 65536,
    },
    SieveParam {
        bits: 249,
        fb_size: 27000,
        large_mult: 100,
        sieve_size: 3 * 65536,
    },
    SieveParam {
        bits: 266,
        fb_size: 50000,
        large_mult: 100,
        sieve_size: 3 * 65536,
    },
    SieveParam {
        bits: 283,
        fb_size: 55000,
        large_mult: 80,
        sieve_size: 3 * 65536,
    },
    SieveParam {
        bits: 298,
        fb_size: 60000,
        large_mult: 80,
        sieve_size: 9 * 65536,
    },
    SieveParam {
        bits: 315,
        fb_size: 80000,
        large_mult: 150,
        sieve_size: 9 * 65536,
    },
    SieveParam {
        bits: 332,
        fb_size: 100000,
        large_mult: 150,
        sieve_size: 9 * 65536,
    },
    SieveParam {
        bits: 348,
        fb_size: 140000,
        large_mult: 150,
        sieve_size: 9 * 65536,
    },
    SieveParam {
        bits: 363,
        fb_size: 210000,
        large_mult: 150,
        sieve_size: 13 * 65536,
    },
    SieveParam {
        bits: 379,
        fb_size: 300000,
        large_mult: 150,
        sieve_size: 17 * 65536,
    },
    SieveParam {
        bits: 395,
        fb_size: 400000,
        large_mult: 150,
        sieve_size: 21 * 65536,
    },
    SieveParam {
        bits: 415,
        fb_size: 500000,
        large_mult: 150,
        sieve_size: 25 * 65536,
    }, /* beyond this point you're crazy */
    SieveParam {
        bits: 440,
        fb_size: 700000,
        large_mult: 150,
        sieve_size: 33 * 65536,
    },
    SieveParam {
        bits: 465,
        fb_size: 900000,
        large_mult: 150,
        sieve_size: 50 * 65536,
    },
    SieveParam {
        bits: 490,
        fb_size: 1100000,
        large_mult: 150,
        sieve_size: 75 * 65536,
    },
    SieveParam {
        bits: 512,
        fb_size: 1300000,
        large_mult: 150,
        sieve_size: 100 * 65536,
    },
];
