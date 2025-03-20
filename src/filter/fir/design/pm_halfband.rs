use crate::error::Result;
use crate::filter::fir::design::estimate_req_filter_transition_bandwidth;
use crate::filter::fir::design::pm::{fir_design_pm, FirPmBandType, FirPmWeightType};
use crate::fft::{Fft, Direction};
use crate::optim::qs1dsearch::{Qs1dSearch, OptimDirection};

use num_complex::Complex32;

// Structured data type
struct FirdespmHalfband {
    // top-level filter design parameters
    m: usize,          // filter semi-length
    h_len: usize,      // filter length, 4*m+1
    ft: f32,           // desired transition band
    h: Vec<f32>,       // resulting filter coefficients

    // utility calculation
    nfft: usize,                // transform size for analysis
    buf_time: Vec<Complex32>,   // time buffer
    buf_freq: Vec<Complex32>,   // frequency buffer
    fft: Fft<f32>,              // transform object
    n: usize,                   // number of points to evaluate
}

impl FirdespmHalfband {
    fn new(m: usize, h_len: usize, nfft: usize, ft: f32) -> Result<Self> {
        let mut nfft = nfft;
        while nfft < 20 * m {
            nfft <<= 1;
        }
        let n = (nfft as f32 * (0.25 - 0.5 * ft)) as usize;

        Ok(Self {
            m,
            h_len,
            ft,
            h: vec![0.0; h_len],
            nfft,
            buf_time: vec![Complex32::new(0.0, 0.0); nfft],
            buf_freq: vec![Complex32::new(0.0, 0.0); nfft],
            fft: Fft::new(nfft, Direction::Forward),
            n,
        })
    }
}

fn firdespm_halfband_utility(gamma: f32, userdata: &mut Option<&mut dyn std::any::Any>) -> f32 {
    // design filter
    let userdata = userdata.as_mut().unwrap().downcast_mut::<FirdespmHalfband>().unwrap();
    let f0 = 0.25 - 0.5 * userdata.ft * gamma;
    let f1 = 0.25 + 0.5 * userdata.ft;
    let bands = [0.00, f0, f1, 0.50];
    let des = [1.0, 0.0];
    let weights = [1.0, 1.0]; // best with {1, 1}
    let wtype = [FirPmWeightType::Flat, FirPmWeightType::Flat]; // best with {flat, flat}
    let h = fir_design_pm(userdata.h_len, 2, &bands, &des, Some(&weights), Some(&wtype), FirPmBandType::Bandpass)
        .expect("firdespm_run failed");
    userdata.h = h;

    // compute utility; copy ideal non-zero coefficients and compute transform
    // force zeros for even coefficients
    for i in 0..userdata.m {
        userdata.h[2 * i] = 0.0;
        userdata.h[userdata.h_len - 2 * i - 1] = 0.0;
    }
    // copy coefficients to input buffer
    for i in 0..userdata.nfft {
        userdata.buf_time[i] = if i < userdata.h_len {
            Complex32::new(userdata.h[i], 0.0)
        } else {
            Complex32::new(0.0, 0.0)
        };
    }
    // compute transform
    userdata.fft.run(&mut userdata.buf_time, &mut userdata.buf_freq);

    // compute metric: power in stop-band
    let u: f32 = (0..userdata.n)
        .map(|i| {
            let idx = userdata.nfft / 2 - i;
            let u_test = userdata.buf_freq[idx].norm();
            u_test * u_test
        })
        .sum();

    // return utility in dB
    10.0 * (u / userdata.n as f32).log10()
}

/// Design halfband filter using Parks-McClellan algorithm given the
/// filter length and desired transition band
///
/// # Arguments
/// * `m` : filter semi-length
/// * `ft` : desired transition band
///
/// # Returns
/// 
/// A vec of filter coefficients
pub fn fir_design_pm_halfband_ft(m: usize, ft: f32) -> Result<Vec<f32>> {
    // create and initialize object
    let mut q = FirdespmHalfband::new(m, 4 * m + 1, 1200, ft)?;

    // create and run search
    {
        let mut optim = Qs1dSearch::new(
            |gamma, userdata| firdespm_halfband_utility(gamma, userdata),
            Some(&mut q),
            OptimDirection::Minimize,
        );
        optim.init_bounds(1.0, 0.9)?;
        for _ in 0..32 {
            optim.step()?;
        }
    }

    Ok(q.h)
}

/// Design halfband filter using Parks-McClellan algorithm given the
/// filter length and desired stop-band suppression
///
/// # Arguments
/// * `m` : filter semi-length
/// * `as_` : desired stop-band suppression
///
/// # Returns
/// 
/// A vec of filter coefficients
pub fn fir_design_pm_halfband_stopband_attenuation(m: usize, as_: f32) -> Result<Vec<f32>> {
    // estimate transition band given other parameters
    let ft = estimate_req_filter_transition_bandwidth(as_, 4 * m + 1)?;

    // return filter design
    fir_design_pm_halfband_ft(m, ft)
}

#[cfg(test)]
mod tests {
    use super::*;
    use test_macro::autotest_annotate;
    use crate::filter::fir::design::estimate_req_filter_stopband_attenuation;
    use crate::utility::test_helpers::{PsdRegion, validate_psd_signalf};

    // test halfband filter design by specifying filter semi-length and transition bandwidth
    fn testbench_firdespm_halfband_ft(m: usize, ft: f32) {
        let h_len = 4 * m + 1;
        let h = fir_design_pm_halfband_ft(m, ft).unwrap();

        // estimate stop band suppression
        let as_ = estimate_req_filter_stopband_attenuation(ft, h_len).unwrap();

        // verify resulting spectrum
        let f0 = 0.25 - 0.5 * ft;
        let f1 = 0.25 + 0.5 * ft;
        let regions = [
            PsdRegion { fmin: -0.5, fmax: -f1, pmin: 0.0,  pmax: -as_, test_lo: false, test_hi: true },
            PsdRegion { fmin: -f0,  fmax:  f0, pmin: -0.1, pmax:  0.1, test_lo: true,  test_hi: true },
            PsdRegion { fmin:  f1,  fmax: 0.5, pmin: 0.0,  pmax: -as_, test_lo: false, test_hi: true },
        ];

        assert!(validate_psd_signalf(&h, &regions).unwrap());
    }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m2_ft400)]
    fn test_firdespm_halfband_m2_ft400() { testbench_firdespm_halfband_ft( 3, 0.400); }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m4_ft400)]
    fn test_firdespm_halfband_m4_ft400() { testbench_firdespm_halfband_ft( 4, 0.400); }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m4_ft200)]
    fn test_firdespm_halfband_m4_ft200() { testbench_firdespm_halfband_ft( 4, 0.200); }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m10_ft200)]
    fn test_firdespm_halfband_m10_ft200() { testbench_firdespm_halfband_ft(10, 0.200); }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m12_ft100)]
    fn test_firdespm_halfband_m12_ft100() { testbench_firdespm_halfband_ft(12, 0.100); }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m20_ft050)]
    fn test_firdespm_halfband_m20_ft050() { testbench_firdespm_halfband_ft(20, 0.050); }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m40_ft050)]
    fn test_firdespm_halfband_m40_ft050() { testbench_firdespm_halfband_ft(40, 0.050); }

    #[test]
    #[autotest_annotate(autotest_firdespm_halfband_m80_ft010)]
    fn test_firdespm_halfband_m80_ft010() { testbench_firdespm_halfband_ft(80, 0.010); }
}