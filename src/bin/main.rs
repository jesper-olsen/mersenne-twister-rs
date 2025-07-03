use mersenne_twister_rs::*;

fn main() {
    const RMAT: [[u64; 4]; 5] = genrand_matrix::<5, 4>(42);

    for row in &RMAT {
        for e in row {
            print!("{:20} ", *e);
        }
        println!()
    }

    let mut mt =
        //MersenneTwister64::new(43);
        MersenneTwister64::new_init_by_array(&[0x12345, 0x23456, 0x34567, 0x45678]);

    for _ in 0..500 {
        println!("{:>20}", mt.genrand());
    }

    for x in mt.genrand_array::<500>().iter() {
        println!("{x:>20}");
    }

    for _ in 0..1000 {
        println!("{:>10.8}", mt.genrand_real2());
    }
}
