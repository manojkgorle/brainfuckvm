#![allow(dead_code, unused_imports)]

static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
struct ConsoleLogger;
use chrono::Local;
use log::{Level, LevelFilter, Metadata, Record};

impl log::Log for ConsoleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Debug
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!(
                "{} [{}] {}:{} - {}",
                Local::now().format("%Y-%m-%dT%H:%M:%S"),
                record.level(),
                record.module_path().unwrap(),
                record.line().unwrap(),
                record.args()
            );
        }
    }

    fn flush(&self) {}
}
use pprof::protos::Message;
use pprof::ProfilerGuard;
use std::fs::File;
use std::io::Write;

pub mod channel;
pub mod evaluation_argument;
pub mod fields;
pub mod fri;
pub mod merkle;
pub mod ntt;
pub mod stark;
pub mod tables;
pub mod univariate_polynomial;
pub mod vm;

pub use crate::fields::{Field, FieldElement};
pub use crate::fri::FriDomain;
pub use crate::stark::{prove, verify_proof};
pub use crate::tables::derive_omicron;
pub use crate::tables::Table;
pub use crate::vm::VirtualMachine;
use rayon::ThreadPoolBuilder;
fn main() {
    env_logger::init();
    let guard = ProfilerGuard::new(100000).unwrap();
    // ThreadPoolBuilder::new()
    // .thread_name(|i| format!("par-iter-{}", i))
    // .build_global()
    // .unwrap();
    let field = Field(18446744069414584321);
    let vm = VirtualMachine::new(field);
    let generator = field.generator().pow((1 << 32) - 1);
    let order = 1 << 32;
    // let code = "++>+-[+--]++.".to_string();
    // let code = "++>+++++[<+>-]++++++++[<++++++>-]<.".to_string();
    let code =  "++++++++[>++++[>++>+++>+++>+<<<<-]>+>+>->>+[<]<-]>>.>---.+++++++..+++.>>.<-.<.+++.------.--------.>>+.>++.".to_string();
    let program = vm.compile(code);

    let (running_time, input_symbols, _output_symbols) =
        vm.run(&program, "1121231241223".to_string());
    log::info!("Running time: {}", running_time);
    let (processor_matrix, memory_matrix, instruction_matrix, input_matrix, output_matrix) =
        vm.simulate(&program, "1121231241223".to_string());
    assert_eq!(running_time as usize, processor_matrix.len());

    let offset = FieldElement::one(field);
    let expansion_f = 1;
    let num_queries = 1;

    let v: &[&[Vec<FieldElement>]] = &[
        &processor_matrix,
        &memory_matrix,
        &instruction_matrix,
        &input_matrix,
        &output_matrix,
    ];
    let (degree_bound, compressed_proof, tp, tm, tins, ti, to, fri_d) =
        prove(v, input_symbols, field, offset, expansion_f, num_queries);

    let maximum_random_int =
        ((degree_bound + 1) * expansion_f as u128 - expansion_f as u128) as u64;

    let _domain = FriDomain::new(
        offset,
        derive_omicron(generator, order, (degree_bound + 1) * expansion_f as u128),
        (degree_bound + 1) * expansion_f as u128,
    );
    verify_proof(
        num_queries as usize,
        maximum_random_int,
        expansion_f as usize,
        field,
        &fri_d,
        &compressed_proof,
        tp,
        tins,
        tm,
        ti,
        to,
        degree_bound as usize,
    );
    match guard.report().build() {
        Ok(report) => {
            let mut file = File::create("profile.pb").unwrap();
            let profile = report.pprof().unwrap();

            let mut content = Vec::new();
            profile.write_to_vec(&mut content).unwrap();
            file.write_all(&content).unwrap();
        }
        Err(_) => {}
    };
}
