// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Youtao Lu@Kim Lab, 2016-2020

use std::env;
use std::process::exit;
use std::mem::take;
use std::path::Path;
use std::io::Write;
use rust_htslib::{bam, bam::Read, bam::Record};
use rand::prelude::*;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use chrono::Local;
use getopts::Options;
use env_logger::{self, Builder};
use log::{error, warn, info, debug, trace, LevelFilter};


type RecordSet = Vec<Record>; 
static VERSION: &str = "0.1.0";

struct Params {
    infile: String,
    outfile: String,
    num: usize,
    seed: u64,
    level: String,
}

fn init_logger(level: &str) {
    Builder::new()
    .format(|buf, record| {
        writeln!(
            buf,
            "[{} {}] {}",
            Local::now().format("%Y-%m-%d %H:%M:%S%.3f %z"),
            record.level(),
            record.args()
        )
    })
    .filter(None, 
        match level {
            "error" => LevelFilter::Error, 
            "warn"  => LevelFilter::Warn,
            "info"  => LevelFilter::Info,
            "debug" => LevelFilter::Debug,
            "trace" => LevelFilter::Trace,
            _ => panic!("unknown debug level!"),
            })
    .init();
}

fn usage(prog: &str, opts: Options) {
    let s = format!("\
Summary:
Random sample --num reads (SE) or read pairs (PE) from BAM or SAM

Usage:
{} --infile input.[bam|sam] --outfile output.bam [--num 5000] [--seed 43] [--help] [--version] [--debug error|warn|info|debug|trace]",
prog);
    println!("{}", opts.usage(&s));
}

fn parse_args(args: &Vec<String>, mut opts: Options) -> Params {
    opts.optopt("i", "infile", "input BAM/SAM, queryname sorted", "FILE");
    opts.optopt("o", "outfile", "output BAM", "FILE");
    opts.optopt("n", "num", "number of reads (read pairs if PE) to downsample (default: 5000)", "INTEGER");
    opts.optopt("s", "seed", "seed (default: None)", "INTEGER");
    opts.optopt("", "level", "level of debugging info, choose from 'error', 'warn', 'info', 'debug', 'trace'", "");
    opts.optflag("h", "help", "print usage");
    opts.optflag("v", "version", "print version");

    let m = opts.parse(&args[1..]).expect("failed to parse arguments!");
    if m.opt_present("h") {
        usage(&args[0], opts);
        exit(0);
    }
    if m.opt_present("v") {
        println!("v{}", VERSION);
        exit(0);
    }
    let infile = match m.opt_str("infile") {
        Some(f) => match Path::new(&f).exists() {
            true => match &*(f
                .split('.')
                .last()
                .expect("Faied to find the file extension!")
                .to_lowercase())
            {
                "sam" | "bam" => f,
                _ => panic!("{} does not seem to be a SAM or BAM!", f),
            },
            false => panic!("{} does not exist!", f),
        },
        None => panic!("--infile is empty!"),
    };
    let outfile = m.opt_str("outfile").expect("invalid --outfile");
    let num = m.opt_get_default("num", 5000).expect("invalid --num");
    let seed = m.opt_get::<u64>("seed").expect("invalid --seed, must be integer");
    let seed = match seed {
        Some(x) => x,
        None => Local::now().timestamp_millis() as u64,
    };
    let level = m.opt_get_default("level", String::from("info")).expect("invalid --level, choose from 'info', 'warn', 'error', 'debug', 'trace'");
    Params {
        infile, 
        outfile,
        num,
        seed,
        level,
    }
}

fn check_header(header: &bam::Header) {
    let header = header.to_hashmap();
    let so = match header.get("HD") {
        Some(a) => a,
        None => { 
            error!("'@HD' not found in header!");
            panic!()
        },
    };

    let so = match so[0].get("SO") {
        Some(a) => a,
        None => { 
            error!("'SO' not found in '@HD'!");
            panic!();
        },
    };

    if so != "queryname" {
        error!("Not sorted by queryname! Please run 'samtools sort -n -o output.bam input.bam' first!");
        panic!();
    }

}

fn main() {
    let args: Vec<String> = env::args().collect();
    let params = parse_args(&args, Options::new());

    let infile = params.infile;
    let outfile = params.outfile;
    let num = params.num;
    let seed = params.seed;
    let level = params.level;
    init_logger(&level);
    info!("{{ infile = {}, outfile = {}, num = {}, seed = {}, level = {} }}", infile, outfile, num, seed, level);

    let mut infh = match bam::Reader::from_path(&infile) {
            Ok(f) => f,
            Err(e) => {
                error!("failed to read {}: {}", &infile, e);
                panic!()
            },
    };

    let header = bam::Header::from_template(infh.header());
    let mut outfh = match bam::Writer::from_path(&outfile, &header, bam::Format::Bam) {
        Ok(f) => f,
        Err(e) => {
            error!("failed to write {}: {}", &outfile, e);
            panic!();
        },
    };
    check_header(&header);

    let mut k = 0;
    let mut v = Vec::<RecordSet>::new();
    let mut rs: RecordSet = RecordSet::new();
    let mut rng = Pcg64::seed_from_u64(seed);
    let mut rid: Option<String>;
    let mut rid_prev: Option<String> = None;

    info!("Iteration starts.");

    for rec in infh.records() {
        match rec {
            Ok(r) => { 
                rid = Some(String::from_utf8(r.qname().to_vec()).expect("invalid qname!"));
                if rid_prev.is_none() || rid_prev.as_ref().expect("invalid qname!").eq(rid.as_ref().expect("invalid qname!")) { 
                    // first record or current record has same qname as previous one; cache it
                    rid_prev = rid.take();
                    rs.push(r);
                    continue;
                } else { // current record is a new template; process the cached; cache it
                    if k < num {
                        v.push(take(&mut rs));
                    } else {
                        let f: f64 = rng.gen();
                        let i = (f * (k as f64)) as usize;
                        if i < num {
                            v[i] = take(&mut rs);
                        }
                    }
                    rid_prev = rid.take();
                    rs.clear();
                    rs.push(r);
                    k += 1;
                    if k % 1_000_000 == 0 {
                        info!("{} reads (read pairs) processed...", k);
                    }
                }
            },
            Err(e) => { 
                error!("empty record: {}", e); 
                panic!();    
            }
        }
    }
    // last record; process the cached
    if k < num {
        v.push(take(&mut rs));
        warn!("--num exceeds the input read counts! output all.");
    } else {
        let f: f64 = rng.gen();
        let i = (f * (k as f64)) as usize;
        if i < num {
            v[i] = take(&mut rs);
        }
    }
    rs.clear();
    for rs in &v {
        for r in rs {
            outfh.write(&r).unwrap();
        }
    }
    info!("All done.");
}
