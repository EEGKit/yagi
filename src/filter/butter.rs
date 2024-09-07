use crate::error::{Error, Result};
use num_complex::Complex32;
use std::f32::consts::PI;

// Compute analog zeros, poles, gain of low-pass Butterworth
// filter, grouping complex conjugates together. If filter
// order is odd, the single real pole (-1) is at the end of
// the array.  There are no zeros for the analog Butterworth
// filter.  The gain is unity.
pub fn butter_azpkf(n: usize) -> Result<(Vec<Complex32>, Vec<Complex32>, Complex32)> {
    if n == 0 {
        return Err(Error::Config("filter order must be greater than zero".to_string()));
    }

    let r = n % 2;
    let l = (n - r) / 2;
    
    let mut pa = Vec::with_capacity(n);
    
    for i in 0..l {
        let theta = (2.0 * (i as f32 + 1.0) + n as f32 - 1.0) * PI / (2.0 * n as f32);
        pa.push(Complex32::from_polar(1.0, theta));
        pa.push(Complex32::from_polar(1.0, -theta));
    }

    if r == 1 {
        pa.push(Complex32::new(-1.0, 0.0));
    }

    if pa.len() != n {
        return Err(Error::Internal("filter order mismatch".to_string()));
    }

    Ok((Vec::new(), pa, Complex32::new(1.0, 0.0)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_butter_azpkf() {
        // test odd/even filter orders
        for n in 1..=10 {
            let (za, pa, ka) = butter_azpkf(n).unwrap();

            // check lengths
            assert_eq!(za.len(), 0);
            assert_eq!(pa.len(), n);

            // check gain
            assert_relative_eq!(ka.re, 1.0);
            assert_relative_eq!(ka.im, 0.0);

            // check poles
            for i in 0..n {
                // poles should be on unit circle
                assert_relative_eq!(pa[i].norm(), 1.0, epsilon = 1e-6);

                // poles should be in left half of s-plane
                assert!(pa[i].re < 0.0);

                // check values for even orders
                if n % 2 == 0 {
                    // poles should be in conjugate pairs
                    if i % 2 == 0 {
                        assert_relative_eq!(pa[i].re, pa[i+1].re, epsilon = 1e-6);
                        assert_relative_eq!(pa[i].im, -pa[i+1].im, epsilon = 1e-6);
                    }
                } else {
                    // odd orders should have one real pole at -1
                    if i == n-1 {
                        assert_relative_eq!(pa[i].re, -1.0, epsilon = 1e-6);
                        assert_relative_eq!(pa[i].im, 0.0, epsilon = 1e-6);
                    } else if i % 2 == 0 {
                        // remaining poles should be in conjugate pairs
                        assert_relative_eq!(pa[i].re, pa[i+1].re, epsilon = 1e-6);
                        assert_relative_eq!(pa[i].im, -pa[i+1].im, epsilon = 1e-6);
                    }
                }
            }
        }
    }
}