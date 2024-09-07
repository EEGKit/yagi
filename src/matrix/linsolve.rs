use crate::matrix::{matrix_gjelim, FloatComplex};
use crate::error::Result;

/// Solve linear system of n equations: Ax = b
///
/// # Arguments
///
/// * `a` - System matrix [size: n x n]
/// * `n` - System size
/// * `b` - Equality vector [size: n x 1]
/// * `x` - Solution vector [size: n x 1]
/// * `opts` - Options (ignored for now)
///
/// # Returns
///
/// `Ok(())` if successful, `Err(...)` otherwise
pub fn matrix_linsolve<T>(a: &[T], n: usize, b: &[T], x: &mut [T], _opts: Option<&()>) -> Result<()>
where
    T: FloatComplex,
{
    // Compute augmented matrix M [size: n x n+1]
    let mut m = Vec::with_capacity(n * (n + 1));
    for r in 0..n {
        m.extend_from_slice(&a[r * n..(r + 1) * n]);
        m.push(b[r]);
    }

    // Run Gauss-Jordan elimination on M
    matrix_gjelim(&mut m, n, n + 1)?;

    // Copy result from right-most column of M
    for r in 0..n {
        x[r] = m[(n + 1) * (r + 1) - 1];
    }

    Ok(())
}
