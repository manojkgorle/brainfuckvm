#![allow(dead_code, unused_imports)]

static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
struct ConsoleLogger;

pub mod channel;
pub mod evaluation_argument;
pub mod fields;
pub mod fri;
pub mod merkle;
pub mod multivariate_polynomial;
pub mod ntt;
pub mod permutation_argument;
pub mod stark;
pub mod tables;
pub mod univariate_polynomial;
pub mod vm;

fn main() {
    println!("Hello, world!");
}
