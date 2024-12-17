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

use crate::fields::{Field, FieldElement};
use crate::fri::FriDomain;
use crate::stark::{prove, verify_proof};
use crate::tables::derive_omicron;
use crate::vm::VirtualMachine;

fn main() {
    let field = Field(18446744069414584321);
    let vm = VirtualMachine::new(field);
    let generator = field.generator().pow((1 << 32) - 1);
    let order = 1 << 32;
    let code = "++>+-[+--]++.".to_string();
    //let code = "++>+++++[<+>-]++++++++[<++++++>-]<.".to_string();
    let program = vm.compile(code);

    let (running_time, input_symbols, output_symbols) = vm.run(&program, "112".to_string());

    let (processor_matrix, memory_matrix, instruction_matrix, input_matrix, output_matrix) =
        vm.simulate(&program, "112".to_string());
    assert_eq!(running_time as usize, processor_matrix.len());

    let offset = FieldElement::one(field);
    let expansion_f = 1;
    let num_queries = 1;

    let v = vec![
        processor_matrix,
        memory_matrix,
        instruction_matrix,
        input_matrix,
        output_matrix,
    ];
    let (degree_bound, compressed_proof, Tp, Tm, Tins, Ti, To, fri_d) =
        prove(v, input_symbols, field, offset, expansion_f, num_queries);
    // println!("\ncompressed proof returned");
    // for i in 0..compressed_proof.len(){
    //     for j in 0..compressed_proof[i].len(){
    //         print!("{} ", compressed_proof[i][j]);
    //     }
    //     println!("\n{}: " ,i);
    // }
    let maximum_random_int =
        ((degree_bound + 1) * expansion_f as u128 - expansion_f as u128) as u64;

    let domain = FriDomain::new(
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
        Tp,
        Tins,
        Tm,
        Ti,
        To,
        degree_bound as usize,
    );
}
