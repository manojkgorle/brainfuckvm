#![allow(unused_variables)]
use io::IOTable;
use memory::MemoryTable;
use processor::ProcessorTable;
use instruction::InstructionTable;

use crate::merkle::*;
use crate::channel::*;
use crate::fri::*;
use crate::fields::Field;
use crate::fields::FieldElement;
use crate::tables::*;
use crate::univariate_polynomial::*;


//@todo boundary, transition and terminal constraints: in all tables: should we be adding them? does that ensure they are individually zero if the sum is zero? check once
//@todo Tipa, Tmpa, Tiea, Toea, Tpea, Tppai, Tppam, Tea, Tea' -> have to write equality amongst them, not written in terminal constraints 
pub struct Stark<'a> {
    pub running_time: i32,
    pub memory_length: usize,
    pub program: &'a [FieldElement],
    pub input_symbols: String,
    pub output_symbols: String,
    expansion_factor: u32,
    security_level: u32,
    num_collinearity_checks: u32,
}

// prove:
// prove parameter - matrices, inputs
// matrices -> processor, memory, instruction, i, o -> in this order
pub fn prove(matrices: Vec<Vec<Vec<FieldElement>>>, inputs: Vec<FieldElement>, field: Field, offset: FieldElement, expansion_f: usize){
    let generator = field.generator().pow((1<<32)-1);
    let order = 1<<32;

    let mut processor_table  = ProcessorTable::new(field, matrices[0].clone().len() as u128, generator, order, matrices[0].clone());
    let mut memory_table = MemoryTable::new(field, matrices[1].len() as u128, generator, order, matrices[1].clone());
    let mut instruction_table = InstructionTable::new(field, matrices[2].len() as u128, generator, order, matrices[2].clone());
    let mut input_table = IOTable::new(field, matrices[3].len() as u128, generator, order, matrices[3].clone());
    let mut output_table = IOTable::new(field, matrices[4].len() as u128, generator, order, matrices[4].clone());

    processor_table.pad();
    memory_table.pad();
    instruction_table.pad();
    input_table.pad();
    output_table.pad();

    let processor_interpol_columns = processor_table.table.clone().interpolate_columns(vec![0,1,2,3,4,5,6]);
    let memory_interpol_columns = memory_table.table.clone().interpolate_columns(vec![0,1,2]);
    let instruction_interpol_columns = instruction_table.table.clone().interpolate_columns(vec![0,1,2]);

    let initial_length = instruction_table.table.clone().height;
    //all codewords are evaluated on this expanded domain that has length expanded_length
    let expanded_length = initial_length*(expansion_f as u128);

    let domain = FriDomain::new(offset, derive_omicron(generator, order, expanded_length), expanded_length);
    
    let mut basecodewords: Vec<Vec<FieldElement>> = Vec::new();

    // basecodewords vector order:
    // processor: clk, ip, ci, ni, mp, mv, inv
    // memory: clk, mp, mv
    // instruction: ip, ci, ni
    // input and output tables are public, we dont commit to those, we only check their termnal extensions after extending

    for i in 0..processor_interpol_columns.clone().len(){
        basecodewords.push(domain.evaluate(processor_interpol_columns[i].clone()));
    }

    for i in 0..memory_interpol_columns.clone().len(){
        basecodewords.push(domain.evaluate(memory_interpol_columns[i].clone()));
    }

    for i in 0..instruction_interpol_columns.clone().len(){
        basecodewords.push(domain.evaluate(instruction_interpol_columns[i].clone()));
    }

    //we are zipping all the base codewords (for each index in order) using concatenation
    //@todo to_bytes function of field element is not working properly? check once

    let mut basecodeword: Vec<Vec<u8>> = Vec::new();

    // for i in 0..expanded_length as usize{
    //     let mut x: Vec<Vec<u8>> = vec![];
    //     for j in 0..basecodewords.len(){
    //         x.push(basecodewords[j][i].to_bytes().iter().map(y:));
    //     }
    // }

    //@todo could not find a function in channel for fiat shamir, ie sending data as string and then getting random element
    //@todo make extend columns function return Terminal value , eg. Tipa, for every table and store it, use it to compare

}

// commit this codeword in merkle tree -> send to verifier, and use the merkle root in fiat shamir
// get 11 challenges array from fiat shamir
// use extend column function on tables -> extends the base columns to extension columns
// interpolate extension columns of all matrices
// evaluate these polynomials on expanded evaluation domains to give extension codewords
// zip/concatenate the extension codewords to give one extension codeword
// commit this codeword in merkle tree -> send to verifier, and use the merkle root in fiat shamir
// get 2 challenges array from fiat shamir
// use generate AIR -> generate zerofier -> generate quotient: on all tables
// form combination polynomial from quotient polynomials and challenges array
// evaluate combination polynomials on expanded evaluation domains to get combination codeword
// perform fri :D, send commitments of fri functions (written in fri module) 
// lessgooo

//@todo IMP - we have interpolated columns of processor table already for commitment and fiat shamir, no need to do it again in AIR

// verifier
// verifier knows - 
// constraints (therefore AIR)
// zerofiers (instruction zerofiers //@todo discuss once)
// combination polynomial equation
// challenges of extension columns
// challenges of composition polynomial
// 
// prover sends to verifier - 
// height (whose correctness is indirectly verified through fri and degree bound)
// base codewords merkle root, extension codewords merkle root
// for each query (index) of verifier, prover sends respective evaluation and merkle authentication path of evaluation
// written in fri decommit_on_query
//
//verifier will perform IOTable computations like extension of columns, will then send those values to prover via channel


impl Stark<'_> {}

#[cfg(test)]
mod stark_test {
    use crate::fields::{Field, FieldElement};
    use crate::vm::VirtualMachine;
    #[test]
    fn test_proving() {
        let field = Field(18446744069414584321);
        let vm = VirtualMachine::new(field);
        let code = "++++".to_string();
        let program = vm.compile(code);
        let (running_time, input_symbols, output_symbols) = vm.run(&program, "".to_string());
        let (processor_matrix, memory_matrix, instruction_matrix, input_matrix, output_matrix) =
            vm.simulate(&program, "".to_string());
        assert_eq!(running_time as usize, processor_matrix.len());
    }
    #[test]
    fn helper_tests() {
        let x = FieldElement::new(318, Field::new(421));
        println!("{}", x);
        let y = x.to_bytes();
        for i in 0..y.len(){
            print!("{}, ", y[i]);
        }
    }
}

