use std::env;
use std::io;
use std::io::prelude::*;
use std::process::exit;
use std::fs::read_to_string;

#[macro_use]
extern crate serde_derive; // all
extern crate serde; // all

mod test_1;

const USAGE: &str = "
Usage: simulation <test_number> [ options ]
       simulation --help

Description of current simulation area

Alternatively, you can run individual benchmarks by running `simulation foo`,
where `foo` is the test number such as `test_1` of a benchmark. Each benchmark
has its own options and modes, so try `simulation foo --help`.

Benchmarks:

  - test_1 : Test 1 description
  - all    : Run all tests
";

fn usage() -> ! {
    let _ = writeln!(&mut io::stderr(), "{}", USAGE);
    exit(1);
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        usage();
    }

    let bench_name = &args[1];
    match &bench_name[..] {
        "test_1" => test_1::main(&args[1..]),
        "all" => {
            test_1::main(&args[1..]);
        }
        _ => usage(),
    }
}


fn read_xy_pairs(filename: &str) -> (Vec<f64>, Vec<f64>) {
    let data = read_to_string(filename);
    let xy = match data {
        Ok(content) => content,
        Err(error) => {
            panic!(
                "Could not open or find file: {}. File name is {:?}",
                error, filename
            );
        }
    };

    let xy_pairs: Vec<&str> = xy.trim().split("\n").collect();
    let mut x: Vec<f64> = Vec::new();
    let mut y: Vec<f64> = Vec::new();

    for pair in xy_pairs {
        let p: Vec<&str> = pair.trim().split(", ").collect();

        x.push(p[0].parse().unwrap());
        y.push(p[1].parse().unwrap());
    }

    (x, y)
}
