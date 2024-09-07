use num_complex::Complex;
use std::f32::consts::PI;

use crate::error::{Error, Result};

// Constants
const NCO_STATIC_LUT_WORDBITS: u32 = 32;
const NCO_STATIC_LUT_NBITS: u32 = 10;
const NCO_STATIC_LUT_SIZE: usize = 1 << NCO_STATIC_LUT_NBITS;

const NCO_PLL_BANDWIDTH_DEFAULT: f32 = 0.1;

// Main NCO structure
pub struct Nco {
    theta: u32,
    d_theta: u32,
    nco_sintab: Vec<f32>,
    alpha: f32,
    beta: f32,
}

impl Nco {
    // Create a new NCO
    pub fn new() -> Self {
        let mut nco = Nco { theta: 0, d_theta: 0, nco_sintab: vec![0.0; NCO_STATIC_LUT_SIZE], alpha: 0.0, beta: 0.0 };

        // Initialize sine table
        for i in 0..NCO_STATIC_LUT_SIZE {
            nco.nco_sintab[i] = (2.0 * PI * i as f32 / NCO_STATIC_LUT_SIZE as f32).sin();
        }

        // Set default PLL bandwidth
        nco.pll_set_bandwidth(NCO_PLL_BANDWIDTH_DEFAULT);

        // Reset object
        nco.reset();

        nco
    }

    // Reset internal state
    pub fn reset(&mut self) {
        self.theta = 0;
        self.d_theta = 0;
        self.pll_reset();
    }

    // Set frequency
    pub fn set_frequency(&mut self, dtheta: f32) {
        self.d_theta = Self::constrain(dtheta);
    }

    // Adjust frequency
    pub fn adjust_frequency(&mut self, df: f32) {
        self.d_theta = self.d_theta.wrapping_add(Self::constrain(df));
    }

    // Set phase
    pub fn set_phase(&mut self, phi: f32) {
        self.theta = Self::constrain(phi);
    }

    // Adjust phase
    pub fn adjust_phase(&mut self, dphi: f32) {
        self.theta = self.theta.wrapping_add(Self::constrain(dphi));
    }

    // Increment internal phase
    pub fn step(&mut self) {
        self.theta = self.theta.wrapping_add(self.d_theta);
    }

    // Get phase
    pub fn get_phase(&self) -> f32 {
        2.0 * PI * self.theta as f32 / ((1u64 << 32) as f32)
    }

    // Get frequency
    pub fn get_frequency(&self) -> f32 {
        let d_theta = 2.0 * PI * self.d_theta as f32 / (1u64 << 32) as f32;
        if d_theta > PI {
            d_theta - 2.0 * PI
        } else {
            d_theta
        }
    }

    // Compute sine of internal phase
    pub fn sin(&self) -> f32 {
        let index = self.static_index();
        self.nco_sintab[index]
    }

    // Compute cosine of internal phase
    pub fn cos(&self) -> f32 {
        let index = self.static_index();
        let index_pi2 = (index + NCO_STATIC_LUT_SIZE / 4) % NCO_STATIC_LUT_SIZE;
        self.nco_sintab[index_pi2]
    }

    // Compute sine and cosine of internal phase
    pub fn sin_cos(&self) -> (f32, f32) {
        let index = self.static_index();
        let index_pi2 = (index + NCO_STATIC_LUT_SIZE / 4) % NCO_STATIC_LUT_SIZE;
        (self.nco_sintab[index], self.nco_sintab[index_pi2])
    }

    // Compute complex exponential of internal phase
    pub fn cexp(&self) -> Complex<f32> {
        let (sin, cos) = self.sin_cos();
        Complex::new(cos, sin)
    }

    // PLL methods

    // Reset PLL state
    pub fn pll_reset(&mut self) {
        // In this implementation, we don't need to do anything here
    }

    // Set PLL bandwidth
    pub fn pll_set_bandwidth(&mut self, bw: f32) {
        if bw < 0.0 {
            panic!("Bandwidth must be positive");
        }
        self.alpha = bw;
        self.beta = bw.sqrt();
    }

    // Advance PLL phase
    pub fn pll_step(&mut self, dphi: f32) {
        self.adjust_frequency(dphi * self.alpha);
        self.adjust_phase(dphi * self.beta);
    }

    // Mixing methods

    // Mix up
    pub fn mix_up(&self, input: Complex<f32>) -> Complex<f32> {
        let (sin, cos) = self.sin_cos();
        input * Complex::new(cos, sin)
    }

    // Mix block up
    pub fn mix_block_up(&mut self, input: &[Complex<f32>], output: &mut [Complex<f32>]) -> Result<()> {
        if input.len() != output.len() {
            return Err(Error::Range("Input and output slices must have the same length".to_owned()));
        }
        for (x, y) in input.iter().zip(output.iter_mut()) {
            *y = self.mix_up(*x);
            self.step();
        }
        Ok(())
    }

    // Mix down
    pub fn mix_down(&self, input: Complex<f32>) -> Complex<f32> {
        let (sin, cos) = self.sin_cos();
        input * Complex::new(cos, sin).conj()
    }

    // Mix block down
    pub fn mix_block_down(&mut self, input: &[Complex<f32>], output: &mut [Complex<f32>]) -> Result<()> {
        if input.len() != output.len() {
            return Err(Error::Range("Input and output slices must have the same length".to_owned()));
        }
        for (x, y) in input.iter().zip(output.iter_mut()) {
            *y = self.mix_down(*x);
            self.step();
        }
        Ok(())
    }

    // Helper functions
    fn constrain(theta: f32) -> u32 {
        let mut theta = theta;
        while theta >= 2.0 * PI {
            theta -= 2.0 * PI;
        }
        while theta < 0.0 {
            theta += 2.0 * PI;
        }
        ((theta / (2.0 * PI)) * (u32::MAX as f32)) as u32
    }

    fn static_index(&self) -> usize {
        (((self.theta as usize) + (1 << (NCO_STATIC_LUT_WORDBITS - NCO_STATIC_LUT_NBITS - 1)))
            >> (NCO_STATIC_LUT_WORDBITS - NCO_STATIC_LUT_NBITS)) as usize
            & (NCO_STATIC_LUT_SIZE - 1)
    }
}

// Additional implementations...

#[cfg(test)]
mod tests {
    use crate::nco::Nco;
    use lazy_static::lazy_static;
    use num_complex::Complex;
    use std::f32::consts::PI;

    // Helper function to calculate phase/frequency error
    fn pll_error(a: f32, b: f32) -> f32 {
        let mut error = a - b;
        while error >= 2.0 * PI {
            error -= 2.0 * PI;
        }
        while error <= -2.0 * PI {
            error += 2.0 * PI;
        }
        error
    }

    // Test phase-locked loop
    fn nco_crcf_pll_test(phase_offset: f32, freq_offset: f32, pll_bandwidth: f32, num_iterations: usize, tol: f32) {
        // Create NCO objects
        let mut nco_tx = Nco::new();
        let mut nco_rx = Nco::new();

        // Initialize objects
        nco_tx.set_phase(phase_offset);
        nco_tx.set_frequency(freq_offset);
        nco_rx.pll_set_bandwidth(pll_bandwidth);

        // Run loop
        for _ in 0..num_iterations {
            // Received complex signal
            let r = nco_tx.cexp();
            let v = nco_rx.cexp();

            // Error estimation
            let phase_error = (r * v.conj()).arg();

            // Update PLL
            nco_rx.pll_step(phase_error);

            // Update NCO objects
            nco_tx.step();
            nco_rx.step();
        }

        // Ensure phase of oscillators is locked
        let phase_error = pll_error(nco_tx.get_phase(), nco_rx.get_phase());
        assert!((phase_error).abs() < tol, "Phase error: {}", phase_error);

        // Ensure frequency of oscillators is locked
        let freq_error = pll_error(nco_tx.get_frequency(), nco_rx.get_frequency());
        assert!((freq_error).abs() < tol, "Frequency error: {}", freq_error);

        println!(
            "nco[bw:{:.4},n={}], phase:{:.6},e={:.4e}, freq:{:.6},e={:.4e}",
            pll_bandwidth, num_iterations, phase_offset, phase_error, freq_offset, freq_error
        );
    }

    #[test]
    fn test_nco_crcf_pll_phase() {
        let bandwidths = [0.1, 0.01, 0.001, 0.0001];
        let tol = 1e-2;

        for &bw in &bandwidths {
            // Adjust number of steps according to loop bandwidth
            let num_steps = (32.0 / bw) as usize;

            // Test various phase offsets
            nco_crcf_pll_test(-PI / 1.1, 0.0, bw, num_steps, tol);
            nco_crcf_pll_test(-PI / 2.0, 0.0, bw, num_steps, tol);
            nco_crcf_pll_test(-PI / 4.0, 0.0, bw, num_steps, tol);
            nco_crcf_pll_test(-PI / 8.0, 0.0, bw, num_steps, tol);
            nco_crcf_pll_test(PI / 8.0, 0.0, bw, num_steps, tol);
            nco_crcf_pll_test(PI / 4.0, 0.0, bw, num_steps, tol);
            nco_crcf_pll_test(PI / 2.0, 0.0, bw, num_steps, tol);
            nco_crcf_pll_test(PI / 1.1, 0.0, bw, num_steps, tol);
        }
    }

    #[test]
    fn test_nco_crcf_pll_freq() {
        let bandwidths = [0.1, 0.05, 0.02, 0.01];
        let tol = 1e-2;

        for &bw in &bandwidths {
            // Adjust number of steps according to loop bandwidth
            let num_steps = (32.0 / bw) as usize;

            // Test various frequency offsets
            nco_crcf_pll_test(0.0, -0.8, bw, num_steps, tol);
            nco_crcf_pll_test(0.0, -0.4, bw, num_steps, tol);
            nco_crcf_pll_test(0.0, -0.2, bw, num_steps, tol);
            nco_crcf_pll_test(0.0, -0.1, bw, num_steps, tol);
            nco_crcf_pll_test(0.0, 0.1, bw, num_steps, tol);
            nco_crcf_pll_test(0.0, 0.2, bw, num_steps, tol);
            nco_crcf_pll_test(0.0, 0.4, bw, num_steps, tol);
            nco_crcf_pll_test(0.0, 0.8, bw, num_steps, tol);
        }
    }

    // Autotest helper function
    fn nco_crcf_phase_test(theta: f32, expected_cos: f32, expected_sin: f32, tol: f32) {
        // Create object
        let mut nco = Nco::new();

        // Set phase
        nco.set_phase(theta);

        // Compute cosine and sine outputs
        let c = nco.cos();
        let s = nco.sin();

        println!(
            "cos({:8.5}) = {:8.5} ({:8.5}) e:{:8.5}, sin({:8.5}) = {:8.5} ({:8.5}) e:{:8.5}",
            theta,
            expected_cos,
            c,
            expected_cos - c,
            theta,
            expected_sin,
            s,
            expected_sin - s
        );

        // Run tests
        assert!((c - expected_cos).abs() < tol, "Cosine error: expected {}, got {}", expected_cos, c);
        assert!((s - expected_sin).abs() < tol, "Sine error: expected {}, got {}", expected_sin, s);
    }

    #[test]
    fn test_nco_crcf_phase() {
        // Error tolerance (higher for NCO)
        let tol = 0.02;

        nco_crcf_phase_test(-6.283185307, 1.000000000, 0.000000000, tol);
        nco_crcf_phase_test(-6.195739393, 0.996179042, 0.087334510, tol);
        nco_crcf_phase_test(-5.951041106, 0.945345356, 0.326070787, tol);
        nco_crcf_phase_test(-5.131745978, 0.407173250, 0.913350943, tol);
        nco_crcf_phase_test(-4.748043551, 0.035647016, 0.999364443, tol);
        nco_crcf_phase_test(-3.041191113, -0.994963998, -0.100232943, tol);
        nco_crcf_phase_test(-1.947799864, -0.368136099, -0.929771914, tol);
        nco_crcf_phase_test(-1.143752030, 0.414182352, -0.910193924, tol);
        nco_crcf_phase_test(-1.029377689, 0.515352252, -0.856978446, tol);
        nco_crcf_phase_test(-0.174356887, 0.984838307, -0.173474811, tol);
        nco_crcf_phase_test(-0.114520496, 0.993449692, -0.114270338, tol);
        nco_crcf_phase_test(0.000000000, 1.000000000, 0.000000000, tol);
        nco_crcf_phase_test(1.436080000, 0.134309213, 0.990939471, tol);
        nco_crcf_phase_test(2.016119855, -0.430749878, 0.902471353, tol);
        nco_crcf_phase_test(2.996498473, -0.989492293, 0.144585621, tol);
        nco_crcf_phase_test(3.403689755, -0.965848729, -0.259106603, tol);
        nco_crcf_phase_test(3.591162483, -0.900634128, -0.434578148, tol);
        nco_crcf_phase_test(5.111428476, 0.388533479, -0.921434607, tol);
        nco_crcf_phase_test(5.727585681, 0.849584319, -0.527452828, tol);
        nco_crcf_phase_test(6.283185307, 1.000000000, -0.000000000, tol);
    }

    #[test]
    fn test_nco_basic() {
        let mut nco = Nco::new();

        let tol = 1e-4; // Error tolerance
        let f = 2.0 * PI / 64.0; // Frequency to test

        nco.set_phase(0.0);
        assert!((nco.cos() - 1.0).abs() < tol, "Cosine at phase 0 error");
        assert!(nco.sin().abs() < tol, "Sine at phase 0 error");

        let (s, c) = nco.sin_cos();
        assert!(s.abs() < tol, "Sine at phase 0 error (sin_cos)");
        assert!((c - 1.0).abs() < tol, "Cosine at phase 0 error (sin_cos)");

        nco.set_phase(PI / 2.0);
        assert!(nco.cos().abs() < tol, "Cosine at phase PI/2 error");
        assert!((nco.sin() - 1.0).abs() < tol, "Sine at phase PI/2 error");

        let (s, c) = nco.sin_cos();
        assert!((s - 1.0).abs() < tol, "Sine at phase PI/2 error (sin_cos)");
        assert!(c.abs() < tol, "Cosine at phase PI/2 error (sin_cos)");

        // Cycle through one full period in 64 steps
        nco.set_phase(0.0);
        nco.set_frequency(f);
        for i in 0..128 {
            let (s, c) = nco.sin_cos();
            assert!((s - (i as f32 * f).sin()).abs() < tol, "Sine error at step {}", i);
            assert!((c - (i as f32 * f).cos()).abs() < tol, "Cosine error at step {}", i);
            nco.step();
        }

        // Double frequency: cycle through one full period in 32 steps
        nco.set_phase(0.0);
        nco.set_frequency(2.0 * f);
        for i in 0..128 {
            let (s, c) = nco.sin_cos();
            assert!((s - (i as f32 * 2.0 * f).sin()).abs() < tol, "Sine error at step {} (double frequency)", i);
            assert!((c - (i as f32 * 2.0 * f).cos()).abs() < tol, "Cosine error at step {} (double frequency)", i);
            nco.step();
        }
    }

    #[test]
    fn test_nco_mixing() {
        // frequency, phase
        let f = 0.1;
        let phi = PI;

        // error tolerance (high for NCO)
        let tol = 0.05;

        // initialize nco object
        let mut nco = Nco::new();
        nco.set_frequency(f);
        nco.set_phase(phi);

        for _ in 0..64 {
            // generate sin/cos
            let (nco_q, nco_i) = nco.sin_cos();

            // mix back to zero phase
            let nco_cplx_in = Complex::new(nco_i, nco_q);
            let nco_cplx_out = nco.mix_down(nco_cplx_in);

            // assert mixer output is correct
            assert!((nco_cplx_out.re - 1.0).abs() < tol, "Real part mixing error");
            assert!(nco_cplx_out.im.abs() < tol, "Imaginary part mixing error");

            // step nco
            nco.step();
        }
    }

    #[test]
    fn test_nco_block_mixing() {
        // frequency, phase
        let f = 0.1;
        let phi = PI;

        // error tolerance (high for NCO)
        let tol = 0.05;

        // number of samples
        const NUM_SAMPLES: usize = 1024;

        // store samples
        let mut x = [Complex::new(0.0, 0.0); NUM_SAMPLES];
        let mut y = [Complex::new(0.0, 0.0); NUM_SAMPLES];

        // generate complex sin/cos
        for i in 0..NUM_SAMPLES {
            x[i] = Complex::new(0.0, f * i as f32 + phi).exp();
        }

        // initialize nco object
        let mut nco = Nco::new();
        nco.set_frequency(f);
        nco.set_phase(phi);

        // mix signal back to zero phase (in pieces)
        let mut i = 0;
        while i < NUM_SAMPLES {
            let n = std::cmp::min(7, NUM_SAMPLES - i);
            nco.mix_block_down(&x[i..i + n], &mut y[i..i + n]).unwrap();
            i += n;
        }

        // assert mixer output is correct
        for i in 0..NUM_SAMPLES {
            assert!((y[i].re - 1.0).abs() < tol, "Real part mixing error at index {}", i);
            assert!(y[i].im.abs() < tol, "Imaginary part mixing error at index {}", i);
        }
    }

    fn testbench_nco_crcf_mix(phase: f32, frequency: f32) {
        use rand::Rng;
        // options
        let buf_len = 1200;
        let tol = 1e-2;

        // create and initialize object
        let mut nco = Nco::new();
        nco.set_phase(phase);
        nco.set_frequency(frequency);

        // generate signal (pseudo-random)
        let mut rng = rand::thread_rng();
        let buf_0: Vec<Complex<f32>> = (0..buf_len).map(|_| Complex::new(0.0, 2.0 * PI * rng.gen::<f32>()).exp()).collect();

        // mix signal
        let mut buf_1 = vec![Complex::new(0.0, 0.0); buf_len];
        nco.mix_block_up(&buf_0, &mut buf_1).unwrap();

        // compare result to expected
        let mut theta = phase;
        for i in 0..buf_len {
            let v = buf_0[i] * Complex::new(0.0, theta).exp();
            assert!((buf_1[i].re - v.re).abs() < tol, "Real part mixing error at index {}", i);
            assert!((buf_1[i].im - v.im).abs() < tol, "Imaginary part mixing error at index {}", i);

            // update and constrain phase
            theta += frequency;
            while theta > PI {
                theta -= 2.0 * PI;
            }
            while theta < -PI {
                theta += 2.0 * PI;
            }
        }
    }

    // test NCO mixing
    #[test]
    fn test_nco_crcf_mix_nco_0() {
        testbench_nco_crcf_mix(0.000, 0.000);
    }

    #[test]
    fn test_nco_crcf_mix_nco_1() {
        testbench_nco_crcf_mix(1.234, 0.000);
    }

    #[test]
    fn test_nco_crcf_mix_nco_2() {
        testbench_nco_crcf_mix(-1.234, 0.000);
    }

    #[test]
    fn test_nco_crcf_mix_nco_3() {
        testbench_nco_crcf_mix(99.000, 0.000);
    }

    #[test]
    fn test_nco_crcf_mix_nco_4() {
        testbench_nco_crcf_mix(PI, 0.000);
    }

    #[test]
    fn test_nco_crcf_mix_nco_5() {
        testbench_nco_crcf_mix(0.000, PI);
    }

    #[test]
    fn test_nco_crcf_mix_nco_6() {
        testbench_nco_crcf_mix(0.000, -PI);
    }

    #[test]
    fn test_nco_crcf_mix_nco_7() {
        testbench_nco_crcf_mix(0.000, 0.123);
    }

    #[test]
    fn test_nco_crcf_mix_nco_8() {
        testbench_nco_crcf_mix(0.000, -0.123);
    }

    #[test]
    fn test_nco_crcf_mix_nco_9() {
        testbench_nco_crcf_mix(0.000, 1e-5);
    }

    // autotest helper function
    fn nco_crcf_frequency_test(phase: f32, frequency: f32, sincos: &[Complex<f32>], num_samples: usize, tol: f32) {
        // create object
        let mut nco = Nco::new();

        // set phase and frequency
        nco.set_phase(phase);
        nco.set_frequency(frequency);

        // run trials
        for i in 0..num_samples {
            // compute complex output
            let y_test = nco.cexp();

            // compare to expected output
            let y = sincos[i];

            // run tests
            assert!((y_test.re - y.re).abs() < tol, "Real part error at index {}: expected {}, got {}", i, y.re, y_test.re);
            assert!((y_test.im - y.im).abs() < tol, "Imaginary part error at index {}: expected {}, got {}", i, y.im, y_test.im);

            // step oscillator
            nco.step();
        }
    }

    #[test]
    fn test_nco_crcf_frequency() {
        // error tolerance (higher for NCO)
        let tol = 0.04;

        // test frequencies with irrational values
        nco_crcf_frequency_test(0.0, 1.0 / 2.0_f32.sqrt(), &NCO_SINCOS_FSQRT1_2, 256, tol); // 1/sqrt(2)
        nco_crcf_frequency_test(0.0, 1.0 / 3.0_f32.sqrt(), &NCO_SINCOS_FSQRT1_3, 256, tol); // 1/sqrt(3)
        nco_crcf_frequency_test(0.0, 1.0 / 5.0_f32.sqrt(), &NCO_SINCOS_FSQRT1_5, 256, tol); // 1/sqrt(5)
        nco_crcf_frequency_test(0.0, 1.0 / 7.0_f32.sqrt(), &NCO_SINCOS_FSQRT1_7, 256, tol);
        // 1/sqrt(7)
    }

    pub fn generate_sincos(frequency: f32, num_samples: usize) -> Vec<Complex<f32>> {
        (0..num_samples)
            .map(|i| {
                let phase = i as f32 * frequency;
                Complex::new(phase.cos(), phase.sin())
            })
            .collect()
    }

    lazy_static! {
        pub static ref NCO_SINCOS_FSQRT1_2: Vec<Complex<f32>> = generate_sincos(1.0 / 2.0_f32.sqrt(), 256);
        pub static ref NCO_SINCOS_FSQRT1_3: Vec<Complex<f32>> = generate_sincos(1.0 / 3.0_f32.sqrt(), 256);
        pub static ref NCO_SINCOS_FSQRT1_5: Vec<Complex<f32>> = generate_sincos(1.0 / 5.0_f32.sqrt(), 256);
        pub static ref NCO_SINCOS_FSQRT1_7: Vec<Complex<f32>> = generate_sincos(1.0 / 7.0_f32.sqrt(), 256);
    }
}
