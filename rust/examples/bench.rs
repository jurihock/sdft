use sdft::{SDFT,Window};

use num::Zero;
use std::time::Instant;

#[allow(non_camel_case_types)]
type c64 = num::complex::Complex<f64>;

fn main() {
    println!("RUST;\tSDFT;\tISDFT");

    let samplerate = 44100.0;
    let dftsize = 1000;
    let window = Window::Hann;
    let latency = 1.0;

    let t0 = Instant::now();
    let mut sdft = SDFT::<f64, f64>::new(
        dftsize,
        window,
        latency);
    let e0 = t0.elapsed().as_micros();

    println!("0;\t{e0};\t{e0}");

    let n = 1 * samplerate as usize;
    let m = sdft.size();

    let mut x = vec![f64::zero(); n];
    let mut y = vec![c64::zero(); n * m];

    let runs = 10;

    for run in 1 .. runs + 1 {
        x.fill(f64::zero());
        y.fill(c64::zero());

        let t1 = Instant::now();
        sdft.sdft(&x, &mut y);
        let e1 = t1.elapsed().as_micros();

        let t2 = Instant::now();
        sdft.isdft(&y, &mut x);
        let e2 = t2.elapsed().as_micros();

        println!("{run};\t{e1};\t{e2}");
    }
}
