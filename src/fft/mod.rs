extern crate rustfft;

use std::sync::Arc;
use num_complex::Complex;
pub enum Direction {
    Forward,
    Backward,
}

impl From<Direction> for rustfft::FftDirection {
    fn from(direction: Direction) -> Self {
        match direction {
            Direction::Forward => rustfft::FftDirection::Forward,
            Direction::Backward => rustfft::FftDirection::Inverse,
        }
    }
}

pub struct Fft<T> {
    fft: Arc<dyn rustfft::Fft<T>>,
}

impl<T: rustfft::FftNum + std::default::Default> Fft<T> {
    pub fn new(n: usize, direction: Direction) -> Self {
        let mut planner = rustfft::FftPlanner::new();
        let fft = planner.plan_fft(n, direction.into());
        Self { fft }
    }

    pub fn run(&self, input: &[Complex<T>], output: &mut [Complex<T>]) {
        output.copy_from_slice(input);
        self.fft.process(output);
    }

    pub fn shift(&self, input: &mut [Complex<T>], n: usize) {
        let n2 = if n % 2 == 0 { n / 2 } else { (n - 1) / 2 };
        for i in 0..n2 {
            let temp = input[i];
            input[i] = input[i + n2];
            input[i + n2] = temp;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test_macro::liquid_test_annotate;
    use approx::assert_relative_eq;

    #[test]
    #[liquid_test_annotate(autotest_fft_shift_4)]
    fn test_shift_4() {
        let mut input = vec![
            Complex::new(0, 0),
            Complex::new(1, 1),
            Complex::new(2, 2),
            Complex::new(3, 3),
        ];
        let expected = vec![
            Complex::new(2, 2),
            Complex::new(3, 3),
            Complex::new(0, 0),
            Complex::new(1, 1),
        ];
        let fft = Fft::new(4, Direction::Forward);
        fft.shift(&mut input, 4);
        assert_eq!(input, expected);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_shift_8)]
    fn test_shift_8() {
        let mut input = vec![
            Complex::new(0, 0),
            Complex::new(1, 1),
            Complex::new(2, 2),
            Complex::new(3, 3),
            Complex::new(4, 4),
            Complex::new(5, 5),
            Complex::new(6, 6),
            Complex::new(7, 7),
        ];
        let expected = vec![
            Complex::new(4, 4),
            Complex::new(5, 5),
            Complex::new(6, 6),
            Complex::new(7, 7),
            Complex::new(0, 0),
            Complex::new(1, 1),
            Complex::new(2, 2),
            Complex::new(3, 3),
        ];
        let fft = Fft::new(8, Direction::Forward);
        fft.shift(&mut input, 8);
        assert_eq!(input, expected);
    }

    fn fft_test_runner(x: &[Complex<f32>], test: &[Complex<f32>], n: usize) {
        let tol = 2e-4;

        let mut y = vec![Complex::<f32>::new(0.0, 0.0); n];
        let mut z = vec![Complex::<f32>::new(0.0, 0.0); n];

        // compute FFT
        let fft_forward = Fft::new(n, Direction::Forward);
        fft_forward.run(x, &mut y);

        // compute IFFT
        let fft_backward = Fft::new(n, Direction::Backward);
        fft_backward.run(&y, &mut z);

        // normalize inverse
        for z_i in z.iter_mut() {
            *z_i /= n as f32;
        }

        // validate results
        for i in 0..n {
            let fft_error = (y[i] - test[i]).norm();
            let ifft_error = (x[i] - z[i]).norm();
            assert_relative_eq!(fft_error, 0.0, epsilon = tol);
            assert_relative_eq!(ifft_error, 0.0, epsilon = tol);
        }
    }

    include!("test_data.rs");

    #[test]
    #[liquid_test_annotate(autotest_fft_2)]
    fn test_fft_2() {
        fft_test_runner(&FFT_TEST_X2, &FFT_TEST_Y2, 2);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_3)]
    fn test_fft_3() {
        fft_test_runner(&FFT_TEST_X3, &FFT_TEST_Y3, 3);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_4)]
    fn test_fft_4() {
        fft_test_runner(&FFT_TEST_X4, &FFT_TEST_Y4, 4);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_5)]
    fn test_fft_5() {
        fft_test_runner(&FFT_TEST_X5, &FFT_TEST_Y5, 5);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_6)]
    fn test_fft_6() {
        fft_test_runner(&FFT_TEST_X6, &FFT_TEST_Y6, 6);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_7)]
    fn test_fft_7() {
        fft_test_runner(&FFT_TEST_X7, &FFT_TEST_Y7, 7);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_8)]
    fn test_fft_8() {
        fft_test_runner(&FFT_TEST_X8, &FFT_TEST_Y8, 8);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_9)]
    fn test_fft_9() {
        fft_test_runner(&FFT_TEST_X9, &FFT_TEST_Y9, 9);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_10)]
    fn test_fft_10() {
        fft_test_runner(&FFT_TEST_X10, &FFT_TEST_Y10, 10);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_16)]
    fn test_fft_16() {
        fft_test_runner(&FFT_TEST_X16, &FFT_TEST_Y16, 16);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_17)]
    fn test_fft_17() {
        fft_test_runner(&FFT_TEST_X17, &FFT_TEST_Y17, 17);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_20)]
    fn test_fft_20() {
        fft_test_runner(&FFT_TEST_X20, &FFT_TEST_Y20, 20);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_21)]
    fn test_fft_21() {
        fft_test_runner(&FFT_TEST_X21, &FFT_TEST_Y21, 21);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_22)]
    fn test_fft_22() {
        fft_test_runner(&FFT_TEST_X22, &FFT_TEST_Y22, 22);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_24)]
    fn test_fft_24() {
        fft_test_runner(&FFT_TEST_X24, &FFT_TEST_Y24, 24);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_26)]
    fn test_fft_26() {
        fft_test_runner(&FFT_TEST_X26, &FFT_TEST_Y26, 26);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_30)]
    fn test_fft_30() {
        fft_test_runner(&FFT_TEST_X30, &FFT_TEST_Y30, 30);
    }

    #[test]
    #[liquid_test_annotate(autotest_fft_32)]
    fn test_fft_32() {
        fft_test_runner(&FFT_TEST_X32, &FFT_TEST_Y32, 32);
    }
}