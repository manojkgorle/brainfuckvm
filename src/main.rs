static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
struct ConsoleLogger;
pub mod fields;

pub mod vm;

pub mod univariate_polynomial;
pub mod multivariate_polynomial;

fn main() {
    println!("Hello, world!");
}
