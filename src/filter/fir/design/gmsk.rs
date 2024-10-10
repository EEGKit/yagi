use crate::error::{Error, Result};
use crate::fft::{fft_run, Direction};
use crate::filter::fir::design;
use crate::math::qf;
use std::f32::consts::PI;

use num_complex::Complex32;

/// Design GMSK transmit filter
///
/// # Arguments
/// * `k`      : samples/symbol
/// * `m`      : symbol delay
/// * `beta`   : rolloff factor (0 < beta <= 1)
/// * `dt`     : fractional sample delay
///
/// # Returns
/// * `h`      : output coefficient buffer (length: 2*k*m+1)
pub fn fir_design_gmsktx(k: usize, m: usize, beta: f32, dt: f32) -> Result<Vec<f32>> {
    // validate input
    if k < 1 {
        return Err(Error::Config("k must be greater than 0".into()));
    }
    if m < 1 {
        return Err(Error::Config("m must be greater than 0".into()));
    }
    if beta < 0.0 || beta > 1.0 {
        return Err(Error::Config("beta must be in [0,1]".into()));
    }

    // derived values
    let h_len = 2 * k * m + 1;

    // compute filter coefficients
    let mut h = vec![0.0; h_len];
    let c0 = 1.0 / (2.0_f32.ln()).sqrt();
    for i in 0..h_len {
        let t = i as f32 / k as f32 - m as f32 + dt;

        h[i] = qf(2.0 * PI * beta * (t - 0.5) * c0) -
               qf(2.0 * PI * beta * (t + 0.5) * c0);
    }

    // normalize filter coefficients such that the filter's
    // integral is pi/2
    let e: f32 = h.iter().sum();
    for h_i in h.iter_mut() {
        *h_i *= PI / (2.0 * e) * k as f32;
    }

    Ok(h)
}

/// Design GMSK receive filter
///
/// # Arguments
/// * `k`      : samples/symbol
/// * `m`      : symbol delay
/// * `beta`   : rolloff factor (0 < beta <= 1)
/// * `dt`     : fractional sample delay
///
/// # Returns
/// * `h`      : output coefficient buffer (length: 2*k*m+1)
pub fn fir_design_gmskrx(k: usize, m: usize, beta: f32, _dt: f32) -> Result<Vec<f32>> {
    // validate input
    if k < 1 {
        return Err(Error::Config("k must be greater than 0".into()));
    }
    if m < 1 {
        return Err(Error::Config("m must be greater than 0".into()));
    }
    if beta < 0.0 || beta > 1.0 {
        return Err(Error::Config("beta must be in [0,1]".into()));
    }

    let bt = beta;

    // internal options
    let beta = bt;                // prototype filter cut-off
    let delta = 1e-3;             // filter design correction factor
    let prototype = design::FirdesFilterType::Kaiser;    // Nyquist prototype

    // derived values
    let h_len = 2 * k * m + 1;   // filter length

    // arrays
    let mut hr = vec![0.0; h_len];         // receive filter coefficients

    // design transmit filter
    let ht = fir_design_gmsktx(k, m, bt, 0.0)?;

    //
    // start of filter design procedure
    //

    // 'internal' arrays
    let mut h_tx = vec![Complex32::new(0.0, 0.0); h_len];      // impulse response of transmit filter
    let mut h_prime = vec![Complex32::new(0.0, 0.0); h_len];   // impulse response of 'prototype' filter
    let mut g_prime = vec![Complex32::new(0.0, 0.0); h_len];   // impulse response of 'gain' filter
    let mut h_hat = vec![Complex32::new(0.0, 0.0); h_len];     // impulse response of receive filter
    
    let mut h_freq_tx = vec![Complex32::new(0.0, 0.0); h_len];      // frequency response of transmit filter
    let mut h_freq_prime = vec![Complex32::new(0.0, 0.0); h_len];   // frequency response of 'prototype' filter
    let mut g_freq_prime = vec![Complex32::new(0.0, 0.0); h_len];   // frequency response of 'gain' filter
    let mut h_freq_hat = vec![Complex32::new(0.0, 0.0); h_len];     // frequency response of receive filter

    // create 'prototype' matched filter
    let h_primef = design::fir_design_prototype(prototype, k, m, beta, 0.0)?;

    // create 'gain' filter to improve stop-band rejection
    let fc = (0.7 + 0.1 * beta) / k as f32;
    let as_ = 60.0;
    let g_primef = design::kaiser::fir_design_kaiser(h_len, fc, as_, 0.0)?;

    // copy to fft input buffer, shifting appropriately
    for i in 0..h_len {
        h_prime[i] = Complex32::new(h_primef[(i + k * m) % h_len], 0.0);
        g_prime[i] = Complex32::new(g_primef[(i + k * m) % h_len], 0.0);
        h_tx[i] = Complex32::new(ht[(i + k * m) % h_len], 0.0);
    }

    // run ffts
    fft_run(&h_prime, &mut h_freq_prime, Direction::Forward);
    fft_run(&g_prime, &mut g_freq_prime, Direction::Forward);
    fft_run(&h_tx, &mut h_freq_tx, Direction::Forward);

    // find minimum of responses
    let mut h_freq_tx_min = f32::MAX;
    let mut h_freq_prime_min = f32::MAX;
    let mut g_freq_prime_min = f32::MAX;
    for i in 0..h_len {
        h_freq_tx_min = h_freq_tx_min.min(h_freq_tx[i].re);
        h_freq_prime_min = h_freq_prime_min.min(h_freq_prime[i].re);
        g_freq_prime_min = g_freq_prime_min.min(g_freq_prime[i].re);
    }

    // compute 'prototype' response, removing minima, and add correction factor
    for i in 0..h_len {
        // compute response necessary to yield prototype response (not exact, but close)
        h_freq_hat[i] = Complex32::new(
            (h_freq_prime[i].re - h_freq_prime_min + delta) / (h_freq_tx[i].re - h_freq_tx_min + delta),
            0.0
        );

        // include additional term to add stop-band suppression
        h_freq_hat[i] *= (g_freq_prime[i].re - g_freq_prime_min) / g_freq_prime[0].re;
    }

    // compute ifft and copy response
    fft_run(&h_freq_hat, &mut h_hat, Direction::Backward);
    for i in 0..h_len {
        hr[i] = h_hat[(i + k * m + 1) % h_len].re / (k * h_len) as f32;
    }

    // copy result, scaling by (samples/symbol)^2
    let mut h = vec![0.0; h_len];
    for i in 0..h_len {
        h[i] = hr[i] * (k * k) as f32;
    }

    Ok(h)
}