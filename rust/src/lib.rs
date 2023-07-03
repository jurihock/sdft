#![allow(dead_code)]

use num::{Float,One,Zero};
use num::complex::Complex;
use std::marker::PhantomData;

/// Cast between primitive types such as `f32` and `f64`.
pub trait CastFrom<T> {
    /// Casts the specified source value
    /// to the destination primitive type
    /// by using the `as` expression.
    fn cast(value: T) -> Self;
}

/// Implements the [`CastFrom`] trait for the specified
/// source `x` and destination `y` types.
macro_rules! impl_cast_from_x_for_y {
    ($x:ty, $y:ty) => {
        impl CastFrom<$x> for $y {
            #[inline]
            fn cast(value: $x) -> $y {
                value as $y
            }
        }
    }
}

// Since only f32 and f64 types are intended,
// the cast trait must only be implemented
// for these two types:
impl_cast_from_x_for_y!(f32, f32);
impl_cast_from_x_for_y!(f32, f64);
impl_cast_from_x_for_y!(f64, f32);
impl_cast_from_x_for_y!(f64, f64);

const SDFT_CONVOLUTION_KERNEL_SIZE: usize = 2;

#[derive(Clone,Copy)]
pub enum Window {
    Boxcar,
    Hann,
    Hamming,
    Blackman,
}

pub struct Analysis<T, F>
    where T: Float + CastFrom<F>,
          F: Float + CastFrom<T> + CastFrom<f64> {
    window: Window,

    weight: F,
    roi: (usize, usize),
    twiddles: Vec<Complex<F>>,

    cursor: usize,
    maxcursor: usize,
    input: Vec<T>,

    accoutput: Vec<Complex<F>>,
    auxoutput: Vec<Complex<F>>,
    fiddles: Vec<Complex<F>>,
}

pub struct Synthesis<T, F>
    where T: Float + CastFrom<F>,
          F: Float + CastFrom<T> + CastFrom<f64> {
    latency: f64,

    weight: F,
    roi: (usize, usize),
    twiddles: Vec<Complex<F>>,

    dontcare: PhantomData<T>,
}

pub struct SDFT<T, F>
    where T: Float + CastFrom<F>,
          F: Float + CastFrom<T> + CastFrom<f64> {
    kernelsize: usize,
    dftsize: usize,
    analysis: Analysis<T, F>,
    synthesis: Synthesis<T, F>,
}

impl<T, F> SDFT<T, F>
    where T: Float + CastFrom<F>,
          F: Float + CastFrom<T> + CastFrom<f64> {

    pub fn new(dftsize: usize,
               window: Window,
               latency: f64
    ) -> Self {
        let kernelsize = SDFT_CONVOLUTION_KERNEL_SIZE;
        let pi = std::f64::consts::PI; // acos(-1)
        let omega = -2.0 * pi / ((dftsize as f64) * 2.0);
        let weight = 2.0 / (1.0 - (omega * (dftsize as f64) * latency).cos());

        SDFT {
            kernelsize: kernelsize,
            dftsize: dftsize,
            analysis: Analysis::<T, F> {
                window,

                weight: F::cast(1.0 / ((dftsize as f64) * 2.0)),
                roi: (0, dftsize),

                twiddles: (0 .. dftsize)
                    .map(|i| Complex::<F>::from_polar(
                        F::one(),
                        F::cast(omega * (i as f64))))
                    .collect(),

                cursor: 0,
                maxcursor: dftsize * 2 - 1,
                input: vec![T::zero(); dftsize * 2],

                accoutput: vec![Complex::<F>::zero(); dftsize],
                auxoutput: vec![Complex::<F>::zero(); dftsize + kernelsize * 2],
                fiddles: vec![Complex::<F>::one(); dftsize],
            },
            synthesis: Synthesis::<T, F> {
                latency,

                weight: F::cast(2.0),
                roi: (0, dftsize),

                twiddles: (0 .. dftsize)
                    .map(|i| Complex::<F>::from_polar(
                        F::cast(weight),
                        F::cast(omega * (i as f64) * (dftsize as f64) * latency)))
                    .collect(),

                dontcare: PhantomData::<T>,
            }
        }
    }

    pub fn size(&self) -> usize { self.dftsize }

    pub fn sdft_scalar(&mut self, sample: &T, dft: &mut [Complex::<F>]) {
        assert_eq!(dft.len(), self.dftsize);

        let newsample = *sample;
        let oldsample = self.analysis.input[self.analysis.cursor];
        self.analysis.input[self.analysis.cursor] = newsample;

        let delta = F::cast(newsample - oldsample);

        if self.analysis.cursor >= self.analysis.maxcursor {
            self.analysis.cursor = 0;

            let mut i = self.analysis.roi.0;
            let mut j = i + self.kernelsize;

            while i < self.analysis.roi.1 {
                self.analysis.accoutput[i] = self.analysis.accoutput[i] + self.analysis.fiddles[i] * delta;
                self.analysis.fiddles[i]   = Complex::<F>::one();
                self.analysis.auxoutput[j] = self.analysis.accoutput[i];

                i += 1;
                j += 1;
            }
        }
        else {
            self.analysis.cursor += 1;

            let mut i = self.analysis.roi.0;
            let mut j = i + self.kernelsize;

            while i < self.analysis.roi.1 {
                self.analysis.accoutput[i] = self.analysis.accoutput[i] + self.analysis.fiddles[i] * delta;
                self.analysis.fiddles[i]   = self.analysis.fiddles[i] * self.analysis.twiddles[i];
                self.analysis.auxoutput[j] = self.analysis.accoutput[i] * self.analysis.fiddles[i].conj();

                i += 1;
                j += 1;
            }
        }

        let auxoffset = (self.kernelsize, self.kernelsize + (self.dftsize - 1));

        for i in 1 .. self.kernelsize + 1  {
            self.analysis.auxoutput[auxoffset.0 - i] = self.analysis.auxoutput[auxoffset.0 + i].conj();
            self.analysis.auxoutput[auxoffset.1 + i] = self.analysis.auxoutput[auxoffset.1 - i].conj();
        }

        self.convolve(&self.analysis.auxoutput, dft);
    }

    pub fn isdft_scalar(&mut self, dft: &[Complex::<F>], sample: &mut T) {
        assert_eq!(dft.len(), self.dftsize);

        let mut result = F::zero();

        if self.synthesis.latency == 1.0 {
            for i in self.synthesis.roi.0 .. self.synthesis.roi.1 {
                let twiddle = F::cast(if i % 2 != 0 { -1.0 } else { 1.0 });
                result = result + dft[i].re * twiddle;
            }
        }
        else {
            for i in self.synthesis.roi.0 .. self.synthesis.roi.1 {
                result = result + (dft[i] * self.synthesis.twiddles[i]).re;
            }
        }

        result = result * self.synthesis.weight;

        *sample = T::cast(result);
    }

    /// Estimate the DFT matrix for the given sample array.
    #[inline]
    pub fn sdft_vector(&mut self, samples: &[T], dfts: &mut [Complex::<F>]) {
        assert_eq!(dfts.len(), samples.len() * self.dftsize);
        for i in 0 .. samples.len() {
            let j = i * self.dftsize .. (i + 1) * self.dftsize;
            self.sdft_scalar(&samples[i], &mut dfts[j]);
        }
    }

    /// Synthesize the sample array from the given DFT matrix.
    #[inline]
    pub fn isdft_vector(&mut self, dfts: &[Complex::<F>], samples: &mut [T]) {
        assert_eq!(dfts.len(), samples.len() * self.dftsize);
        for i in 0 .. samples.len() {
            let j = i * self.dftsize .. (i + 1) * self.dftsize;
            self.isdft_scalar(&dfts[j], &mut samples[i]);
        }
    }

    /// Estimate the DFT matrix for the given sample array.
    /// This is a shortcut for the function [`sdft_vector`].
    #[inline]
    pub fn sdft(&mut self, samples: &[T], dfts: &mut [Complex::<F>]) {
        self.sdft_vector(samples, dfts);
    }

    /// Synthesize the sample array from the given DFT matrix.
    /// This is a shortcut for the function [`isdft_vector`].
    #[inline]
    pub fn isdft(&mut self, dfts: &[Complex::<F>], samples: &mut [T]) {
        self.isdft_vector(dfts, samples);
    }

    #[inline]
    fn convolve(&self, input: &[Complex::<F>], output: &mut [Complex::<F>]) {
        let roi = self.analysis.roi;
        let window = self.analysis.window;
        let weight = self.analysis.weight;
        let kernelsize = self.kernelsize;

        let l2 = kernelsize - 2;
        let l1 = kernelsize - 1;
        let m  = kernelsize;
        let r1 = kernelsize + 1;
        let r2 = kernelsize + 2;

        for i in roi.0 .. roi.1 {
            match window {
                Window::Hann => {
                    let a = input[i + m] + input[i + m];
                    let b = input[i + l1] + input[i + r1];
                    output[i] = (a - b) * weight * F::cast(0.25);
                },
                Window::Hamming => {
                    let a = input[i + m] * F::cast(0.54);
                    let b = (input[i + l1] + input[i + r1]) * F::cast(0.23);
                    output[i] = (a - b) * weight;
                },
                Window::Blackman => {
                    let a = input[i + m] * F::cast(0.42);
                    let b = (input[i + l1] + input[i + r1]) * F::cast(0.25);
                    let c = (input[i + l2] + input[i + r2]) * F::cast(0.04);
                    output[i] = (a - b + c) * weight;
                },
                _ => {
                    output[i] = input[i + m] * weight;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn todo() {
        // let result = add(2, 2);
        // assert_eq!(result, 4);
    }
}
