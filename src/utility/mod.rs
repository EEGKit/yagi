// Utility module
// Current state:
// - bits ready to use (+autotests)
// - other parts of utility TBD

pub mod bits;

pub use bits::*;

#[cfg(test)]
pub(crate) mod test_helpers;