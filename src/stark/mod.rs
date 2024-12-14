#![allow(unused_variables)]
use std::io::Read;

use instruction::Indices;
use instruction::InstructionTable;
use io::IOTable;
use memory::MemoryTable;
use processor::ProcessorTable;

use crate::channel;
use crate::channel::*;
use crate::fields::Field;
use crate::fields::FieldElement;
use crate::fri::*;
use crate::merkle::*;
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

pub enum ChallengeIndices {
    A,
    B,
    C,
    D,
    E,
    F,
    Alpha,
    Beta,
    Delta,
    Gamma,
    Eta,
}

// prove:
// prove parameter - matrices, inputs
// matrices -> processor, memory, instruction, i, o -> in this order
pub fn prove(
    matrices: Vec<Vec<Vec<FieldElement>>>,
    inputs: Vec<FieldElement>,
    field: Field,
    offset: FieldElement,
    expansion_f: usize,
) {
    let generator = field.generator().pow((1 << 32) - 1);
    let order = 1 << 32;

    let mut processor_table = ProcessorTable::new(
        field,
        matrices[0].clone().len() as u128,
        generator,
        order,
        matrices[0].clone(),
    );
    let mut memory_table = MemoryTable::new(
        field,
        matrices[1].len() as u128,
        generator,
        order,
        matrices[1].clone(),
    );
    let mut instruction_table = InstructionTable::new(
        field,
        matrices[2].len() as u128,
        generator,
        order,
        matrices[2].clone(),
    );
    let mut input_table = IOTable::new(
        field,
        matrices[3].len() as u128,
        generator,
        order,
        matrices[3].clone(),
    );
    let mut output_table = IOTable::new(
        field,
        matrices[4].len() as u128,
        generator,
        order,
        matrices[4].clone(),
    );

    //@todo instruction table height passed as parameter
    processor_table.pad();
    memory_table.pad();
    instruction_table.pad();
    input_table.pad();
    output_table.pad();

    let processor_interpol_columns = processor_table
        .table
        .clone()
        .interpolate_columns(vec![0, 1, 2, 3, 4, 5, 6]);
    let memory_interpol_columns = memory_table
        .table
        .clone()
        .interpolate_columns(vec![0, 1, 2]);
    let instruction_interpol_columns = instruction_table
        .table
        .clone()
        .interpolate_columns(vec![0, 1, 2]);

    let initial_length = instruction_table.table.clone().height;
    //all codewords are evaluated on this expanded domain that has length expanded_length
    let expanded_length = initial_length * (expansion_f as u128);

    let domain = FriDomain::new(
        offset,
        derive_omicron(generator, order, expanded_length),
        expanded_length,
    );

    let mut basecodewords: Vec<Vec<FieldElement>> = Vec::new();

    // basecodewords vector order:
    // processor: clk, ip, ci, ni, mp, mv, inv
    // memory: clk, mp, mv
    // instruction: ip, ci, ni
    // input and output tables are public, we dont commit to those, we only check their termnal extensions after extending

    for i in 0..processor_interpol_columns.clone().len() {
        basecodewords.push(domain.evaluate(processor_interpol_columns[i].clone()));
    }

    for i in 0..memory_interpol_columns.clone().len() {
        basecodewords.push(domain.evaluate(memory_interpol_columns[i].clone()));
    }

    for i in 0..instruction_interpol_columns.clone().len() {
        basecodewords.push(domain.evaluate(instruction_interpol_columns[i].clone()));
    }

    //we are zipping all the base codewords (for each index in order) using concatenation

    let mut basecodeword: Vec<Vec<u8>> = Vec::new();

    for i in 0..expanded_length as usize {
        let mut x: Vec<u8> = vec![];
        for j in 0..basecodewords.len() {
            x.extend(basecodewords[j][i].to_bytes().iter().map(|&x| x));
        }
        basecodeword.push(x);
    }

    //@todo make extend columns function return Terminal value , eg. Tipa, for every table and store it, use it to compare
    let mut data1 = vec![];

    for i in 0..basecodeword.len() {
        // difficulty in implementing -> let n = basecodewords[0].len();
        // so hardcoded the value to 32*13 = 416 -> where 13 => clk, ip, ci, ni, mp, mv, inv, clk, mp, mv, ip, ci, ni
        let array: &[u8] = &basecodeword[i].to_vec();

        data1.push(FieldElement::from_bytes(array));
    }

    // get 11 challenges array from fiat shamir
    let mut channel = Channel::new();
    let merkle1 = MerkleTree::new(&data1);
    channel.send(merkle1.inner.root().unwrap().to_vec());

    let mut challenges_extension = vec![];

    for i in 0..10 {
        let x = channel.receive_random_field_element(field);
        challenges_extension.push(x);
        channel.send(x.to_bytes());
    }
    challenges_extension.push(channel.receive_random_field_element(field));

    // use extend column function on tables -> extends the base columns to extension columns
    let Terminal_processor = processor_table.extend_columns(challenges_extension.clone());
    let Terminal_memory = memory_table.extend_column_ppa(1, challenges_extension.clone());
    let Terminal_instruction = instruction_table.extend_column(1, challenges_extension.clone());
    let Terminal_input = input_table
        .extend_column_ea(1, challenges_extension[ChallengeIndices::Gamma as usize])
        .clone();
    let Terminal_output = output_table
        .extend_column_ea(1, challenges_extension[ChallengeIndices::Delta as usize])
        .clone();

    //These contain polynomials for interpolation of extension columns
    let processor_interpol_columns_2 = processor_table
        .table
        .clone()
        .interpolate_columns(vec![7, 8, 9, 10]);
    let memory_interpol_columns_2 = memory_table.table.clone().interpolate_columns(vec![3]);
    let instruction_interpol_columns_2 = instruction_table
        .table
        .clone()
        .interpolate_columns(vec![3, 4]);

    let mut extension_codewords: Vec<Vec<FieldElement>> = Vec::new();

    // extensioncodewords vector order:
    // processor: ipa, mpa, iea, oea
    // memory: ppa
    // instruction: ppa, pea
    // input and output tables are public, we dont commit to those, we only check their termnal extensions after extending

    for i in 0..processor_interpol_columns_2.clone().len() {
        extension_codewords.push(domain.evaluate(processor_interpol_columns_2[i].clone()));
    }

    for i in 0..memory_interpol_columns_2.clone().len() {
        extension_codewords.push(domain.evaluate(memory_interpol_columns_2[i].clone()));
    }

    for i in 0..instruction_interpol_columns_2.clone().len() {
        extension_codewords.push(domain.evaluate(instruction_interpol_columns_2[i].clone()));
    }

    let mut extension_codeword: Vec<Vec<u8>> = Vec::new();

    for i in 0..expanded_length as usize {
        let mut x: Vec<u8> = vec![];
        for j in 0..extension_codewords.len() {
            x.extend(extension_codewords[j][i].to_bytes().iter().map(|&x| x));
        }
        extension_codeword.push(x);
    }

    let mut data2 = vec![];

    for i in 0..extension_codeword.len() {
        let array: &[u8] = &extension_codeword[i].to_vec();

        data2.push(FieldElement::from_bytes(array));
    }

    let merkle2 = MerkleTree::new(&data2);
    channel.send(merkle2.inner.root().unwrap().to_vec());

    let mut challenges_combination = vec![];
    let x = channel.receive_random_field_element(field);
    challenges_combination.push(x);
    channel.send(x.to_bytes());
    challenges_combination.push(channel.receive_random_field_element(field));

    let eval = FieldElement::zero(field);

    let processor_quotients = processor_table.generate_quotients(
        challenges_extension.clone(),
        Terminal_processor[0],
        Terminal_processor[1],
        Terminal_processor[2],
        Terminal_processor[3],
    );
    let memory_quotients =
        memory_table.generate_quotients(challenges_extension.clone(), Terminal_memory[0]);
    let instruction_quotients = instruction_table.generate_quotients(
        challenges_extension,
        Terminal_instruction[0],
        Terminal_instruction[1],
    );

    //for inter table arguments constraints
    assert_eq!(Terminal_processor[0], Terminal_instruction[0]); //Tipa = Tppa
    assert_eq!(Terminal_processor[1], Terminal_memory[0]); //Tmpa = Tppa
    assert_eq!(Terminal_processor[2], Terminal_input[0]); //Tipa = Tea input
    assert_eq!(Terminal_processor[3], Terminal_output[0]); //Tipa = Tea output
                                                           //let this be for now:- assert_eq!(Terminal_instruction[1], Tpea); //Tpea = program evaluation

    //form combination polynomial
    let combination = combination_polynomial(
        processor_quotients,
        memory_quotients,
        instruction_quotients,
        challenges_combination,
        instruction_table.table.height as usize,
        field,
    );
    let combination_codeword = domain.evaluate(combination.clone());

    let merkle_combination = MerkleTree::new(&combination_codeword);
    channel.send(merkle_combination.inner.root().unwrap().to_vec());

    let (fri_polys, fri_domains, fri_layers, fri_merkles) = fri_commit(
        combination.clone(),
        domain,
        combination_codeword,
        merkle_combination,
        &mut channel,
    );

    let no_of_queries = 5;
    decommit_fri(
        no_of_queries,
        expansion_f,
        1 << 64 - 1 << 32 + 1,
        vec![&data1, &data2],
        vec![&merkle1, &merkle2],
        &fri_layers,
        &fri_merkles,
        &mut channel,
    );

    //print channel proof, proofsize, time taken for running prover, space taken etc etc.
}

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
        for i in 0..y.len() {
            print!("{}, ", y[i]);
        }
    }
}
