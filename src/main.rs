#![allow(dead_code, unused_imports)]

static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
struct ConsoleLogger;

pub mod fields;
pub mod multivariate_polynomial;
pub mod tables;
pub mod univariate_polynomial;
pub mod vm;
pub mod ntt;
pub mod evaluation_argument;
pub mod permutation_argument;
pub mod merkle;

fn main() {
    println!("Hello, world!");
}
