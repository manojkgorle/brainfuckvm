static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
struct ConsoleLogger;

pub mod fields;
pub mod multivariate_polynomial;
pub mod tables;
pub mod univariate_polynomial;
pub mod vm;

fn main() {
    println!("Hello, world!");
}
