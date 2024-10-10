use crate::error::{Error, Result};
use crate::dotprod::DotProd;
use crate::filter::iir::iirfilt::IirFilt;
use crate::filter::iir::design::{IirdesFilterType, IirdesBandType, IirdesFormat};
use std::collections::VecDeque;
use num_complex::{ComplexFloat, Complex32};

#[derive(Clone, Debug)]
pub struct IirDecim<T, Coeff = T> {
    decimation_factor: usize,
    iirfilt: IirFilt<T, Coeff>,
}


impl<T, Coeff> IirDecim<T, Coeff>
where
    T: Copy + Default + ComplexFloat<Real = f32> + std::ops::Mul<Coeff, Output = T> + From<Coeff>,
    Coeff: Copy + Default + ComplexFloat<Real = f32> + std::ops::Mul<T, Output = T> + Into<Complex32>,
    VecDeque<T>: DotProd<Coeff, Output = T>,
    f32: Into<Coeff>,
{
    pub fn new(decimation_factor: usize, b: &[Coeff], a: &[Coeff]) -> Result<Self> {
        if decimation_factor < 2 {
            return Err(Error::Config("decimation factor must be greater than 1".into()));
        }

        let iirfilt = IirFilt::new(b, a)?;

        Ok(Self { decimation_factor, iirfilt })
    }

    pub fn new_default(decimation_factor: usize, order: usize) -> Result<Self> {
        Self::new_prototype(
            decimation_factor,
            IirdesFilterType::Butter,
            IirdesBandType::Lowpass,
            IirdesFormat::SecondOrderSections,
            order,
            0.5 / decimation_factor as f32,
            0.0,
            0.1,
            60.0,
        )
    }

    pub fn new_prototype(
        decimation_factor: usize,
        ftype: IirdesFilterType,
        btype: IirdesBandType,
        format: IirdesFormat,
        order: usize,
        fc: f32,
        f0: f32,
        ap: f32,
        as_: f32,
    ) -> Result<Self> {
        if decimation_factor < 2 {
            return Err(Error::Config("interp factor must be greater than 1".into()));
        }

        let iirfilt = IirFilt::new_prototype(ftype, btype, format, order, fc, f0, ap, as_)?;

        Ok(Self { decimation_factor, iirfilt })
    }

    pub fn reset(&mut self) -> () {
        self.iirfilt.reset();
    }

    pub fn execute(&mut self, x: &[T]) -> T {
        let mut y = T::default();
        for (i, &xi) in x.iter().enumerate() {
            let v = self.iirfilt.execute(xi);
            if i == 0 {
                y = v;
            }
        }
        y
    }

    pub fn execute_block(&mut self, x: &[T], y: &mut [T]) -> () {
        for (i, xi) in x.chunks(self.decimation_factor).enumerate() {
            y[i] = self.execute(xi);
        }
    }

    pub fn groupdelay(&self, fc: f32) -> Result<f32> {
        self.iirfilt.groupdelay(fc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test_macro::autotest_annotate;
    use crate::random::randnf;

    #[test]
    #[autotest_annotate(autotest_iirdecim_copy)]
    fn test_iirdecim_copy() {
        // create base object
        let mut q0 = IirDecim::<Complex32, f32>::new_default(3, 7).unwrap();

        // run samples through filter
        let mut buf = [Complex32::default(); 3];
        for _ in 0..20 {
            for j in 0..3 {
                buf[j] = Complex32::new(randnf(), randnf());
            }
            let _y0 = q0.execute(&buf);
        }

        // copy object
        let mut q1 = q0.clone();

        // run samples through both filters in parallel
        for _ in 0..60 {
            for j in 0..3 {
                buf[j] = Complex32::new(randnf(), randnf());
            }
            let y0 = q0.execute(&buf);
            let y1 = q1.execute(&buf);

            assert_eq!(y0, y1);
        }

        // objects are automatically destroyed when they go out of scope
    }
}