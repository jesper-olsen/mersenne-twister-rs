use mersenne_twister_rs::*;
use rand::Rng;

fn main() {
    //const RMAT: [[u64; 4]; 5] = genrand_matrix::<5, 4>(42);
    const RMAT: [[u64; 4]; 5] = genrand_matrix::<5, 4>(42);

    println!("const 5x4 u64 matrix");
    for row in &RMAT {
        for e in row {
            print!("{:20} ", *e);
        }
        println!()
    }

    println!("\nconst 5x4 f64 matrix");
    for row in &RMAT {
        for e in row {
            print!("{:2.2} ", u64_to_f64(*e));
        }
        println!()
    }

    let mut mt =
        //MersenneTwister64::new(43);
        MersenneTwister64::new_init_by_array(&[0x12345, 0x23456, 0x34567, 0x45678]);

    println!("\n500 u64");
    for _ in 0..500 {
        //println!("{:>20}", mt.genrand());
        println!("{:>20}", mt.random::<u64>());
    }

    println!("\nconst 500 u64 array");
    for x in mt.genrand_array::<500>().iter() {
        println!("{x:>20}");
    }

    println!("\n1000 f64");
    for _ in 0..1000 {
        println!("{:>10.8}", mt.genrand_real2());
    }
}
