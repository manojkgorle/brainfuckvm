#![allow(unused_variables)]
use std::io::Read;
use std::f64;

use instruction::Indices;
use io::IOTable;
use memory::MemoryTable;
use processor::ProcessorTable;
use instruction::InstructionTable;

use crate::channel;
use crate::merkle::*;
use crate::channel::*;
use crate::fri::*;
use crate::fields::Field;
use crate::fields::FieldElement;
use crate::tables::*;
use crate::univariate_polynomial::*;
static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
use log::{Level, LevelFilter, Metadata, Record};
use chrono::Local;
struct ConsoleLogger;

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
#[warn(non_snake_case)]
pub fn prove(matrices: Vec<Vec<Vec<FieldElement>>>, inputs: Vec<FieldElement>, field: Field, offset: FieldElement, expansion_f: usize){
 
    let generator = field.generator().pow((1<<32)-1);
    let order = 1<<32;
    log::info!("Generating tables");

    let mut processor_table  = ProcessorTable::new(field, matrices[0].clone().len() as u128, generator, order, matrices[0].clone());
    let mut memory_table = MemoryTable::new(field, matrices[1].len() as u128, generator, order, matrices[1].clone());
    let mut instruction_table = InstructionTable::new(field, matrices[2].len() as u128, generator, order, matrices[2].clone());
    let mut input_table = IOTable::new(field, matrices[3].len() as u128, generator, order, matrices[3].clone());
    let mut output_table = IOTable::new(field, matrices[4].len() as u128, generator, order, matrices[4].clone());
    log::info!("padding all tables");
    processor_table.pad();
    memory_table.pad();
    instruction_table.pad();
    input_table.pad();
    output_table.pad();
    log::info!("Interpolating all tables");

    let processor_interpol_columns = processor_table.table.clone().interpolate_columns(vec![0,1,2,3,4,5,6]);
    let memory_interpol_columns = memory_table.table.clone().interpolate_columns(vec![0,1,2]);
    let instruction_interpol_columns = instruction_table.table.clone().interpolate_columns(vec![0,1,2]);

    let initial_length = instruction_table.table.clone().height;
    //all codewords are evaluated on this expanded domain that has length expanded_length
    let expanded_length = initial_length*(expansion_f as u128);
    log::info!("Extending the domain");

    let domain = FriDomain::new(offset, derive_omicron(generator, order, expanded_length), expanded_length);
    
    let mut basecodewords: Vec<Vec<FieldElement>> = Vec::new();

    // basecodewords vector order:
    // processor: clk, ip, ci, ni, mp, mv, inv
    // memory: clk, mp, mv
    // instruction: ip, ci, ni
    // input and output tables are public, we dont commit to those, we only check their terminal extensions after extending
    log::info!("evaluating on the extended domain");

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

    let mut basecodeword: Vec<Vec<u8>> = Vec::new();
    log::info!("zipping all the codewords on the extended domain");
    for i in 0..expanded_length as usize{
        let mut x: Vec<u8> = vec![];
        for j in 0..basecodewords.len(){
            x.extend(basecodewords[j][i].to_bytes().iter().map(|&x| x));
        }
        basecodeword.push(x);
    }

    //@todo make extend columns function return Terminal value , eg. Tipa, for every table and store it, use it to compare
    let mut data1 = vec![];

    for i in 0..basecodeword.len(){ 
        // difficulty in implementing -> let n = basecodewords[0].len(); 
        // so hardcoded the value to 32*13 = 416 -> where 13 => clk, ip, ci, ni, mp, mv, inv, clk, mp, mv, ip, ci, ni
        let array: &[u8] = &basecodeword[i].to_vec();

        data1.push(FieldElement::from_bytes(array));
    }
    
    // get 11 challenges array from fiat shamir
    let mut channel = Channel::new();
    log::info!("commiting the base codewords");
    let merkle1 = MerkleTree::new(&data1);
    channel.send(merkle1.inner.root().unwrap().to_vec());

    let mut challenges_extension = vec![];

    for i in 0..10{
        let x = channel.receive_random_field_element(field);
        challenges_extension.push(x);
        // channel.send(x.to_bytes());
    }
    challenges_extension.push(channel.receive_random_field_element(field));

    // use extend column function on tables -> extends the base columns to extension columns
    log::info!("generating the extension coloumn using the fiat-shamir challenges");
    let Terminal_processor = processor_table.extend_columns(challenges_extension.clone());
    let Terminal_memory = memory_table.extend_column_ppa(1, challenges_extension.clone());
    let Terminal_instruction = instruction_table.extend_column(1, challenges_extension.clone());
    let Terminal_input = input_table.extend_column_ea(1, challenges_extension[ChallengeIndices::Gamma as usize]).clone();
    let Terminal_output = output_table.extend_column_ea(1, challenges_extension[ChallengeIndices::Delta as usize]).clone();

    //These contain polynomials for interpolation of extension columns 
    log::info!("interpolating the extension columns");
    let processor_interpol_columns_2 = processor_table.table.clone().interpolate_columns(vec![7, 8, 9, 10]);
    let memory_interpol_columns_2 = memory_table.table.clone().interpolate_columns(vec![3]);
    let instruction_interpol_columns_2 = instruction_table.table.clone().interpolate_columns(vec![3, 4]);

    let mut extension_codewords: Vec<Vec<FieldElement>> = Vec::new();

    // extensioncodewords vector order:
    // processor: ipa, mpa, iea, oea
    // memory: ppa
    // instruction: ppa, pea
    // input and output tables are public, we dont commit to those, we only check their terminal extensions after extending
    log::info!("evaluating the extendion columns on the extended domain");

    for i in 0..processor_interpol_columns_2.clone().len(){
        extension_codewords.push(domain.evaluate(processor_interpol_columns_2[i].clone()));
    }

    for i in 0..memory_interpol_columns_2.clone().len(){
        extension_codewords.push(domain.evaluate(memory_interpol_columns_2[i].clone()));
    }

    for i in 0..instruction_interpol_columns_2.clone().len(){
        extension_codewords.push(domain.evaluate(instruction_interpol_columns_2[i].clone()));
    }

    let mut extension_codeword: Vec<Vec<u8>> = Vec::new();
    log::info!("zipping all the extension codewords");
    for i in 0..expanded_length as usize{
        let mut x: Vec<u8> = vec![];
        for j in 0..extension_codewords.len(){
            x.extend(extension_codewords[j][i].to_bytes().iter().map(|&x| x));
        }
        extension_codeword.push(x);
    }

    let mut data2 = vec![];

    for i in 0..extension_codeword.len(){ 
        let array: &[u8] = &extension_codeword[i].to_vec();

        data2.push(FieldElement::from_bytes(array));
    }
    log::info!("commiting the extension codewords");
    let merkle2 = MerkleTree::new(&data2);
    channel.send(merkle2.inner.root().unwrap().to_vec());

    let mut challenges_combination = vec![];
    let x = channel.receive_random_field_element(field);
    challenges_combination.push(x);
    // channel.send(x.to_bytes());
    challenges_combination.push(channel.receive_random_field_element(field));

    let eval = FieldElement::zero(field);


    let processor_quotients = processor_table.generate_quotients(challenges_extension.clone(), Terminal_processor[0], Terminal_processor[1], Terminal_processor[2], Terminal_processor[3]);
    let memory_quotients = memory_table.generate_quotients(challenges_extension.clone(), Terminal_memory[0]);
    let instruction_quotients = instruction_table.generate_quotients(challenges_extension, Terminal_instruction[0], Terminal_instruction[1]);

    //for inter table arguments constraints
    assert_eq!(Terminal_processor[0], Terminal_instruction[0]); //Tipa = Tppa
    assert_eq!(Terminal_processor[1], Terminal_memory[0]); //Tmpa = Tppa
    assert_eq!(Terminal_processor[2], Terminal_input[0]); //Tipa = Tea input
    assert_eq!(Terminal_processor[3], Terminal_output[0]); //Tipa = Tea output
    //let this be for now:- assert_eq!(Terminal_instruction[1], Tpea); //Tpea = program evaluation

    //form combination polynomial
    let combination= combination_polynomial(processor_quotients, memory_quotients, instruction_quotients, challenges_combination, instruction_table.table.height as usize, field);
    let combination_codeword = domain.evaluate(combination.clone());
    
    let merkle_combination = MerkleTree::new(&combination_codeword);
    channel.send(merkle_combination.inner.root().unwrap().to_vec());

    let (fri_polys, fri_domains, fri_layers, fri_merkles) =  fri_commit(combination.clone(), domain, combination_codeword, merkle_combination, &mut channel);

    let no_of_queries = 5;
    decommit_fri(no_of_queries, expansion_f, 1<<64-1<<32+1, vec![&data1, &data2], vec![&merkle1, &merkle2], &fri_layers, &fri_merkles, &mut channel);

    //print channel proof, proofsize, time taken for running prover, space taken etc etc.

}

// use generate AIR -> generate zerofier -> generate quotient: on all tables
// form combination polynomial from quotient polynomials and challenges array
// evaluate combination polynomials on expanded evaluation domains to get combination codeword
// perform fri :D, send commitments of fri functions (written in fri module) 
// lessgooo

//@todo IMP - we have interpolated columns of processor table already for commitment and fiat shamir, no need to do it again in AIR

 // ii) Decommitment -> query phase
    // Decommitment involves verifier sending random elements from evaluation domain to prover. and prover responding with decommitments to the evaluations, which involve sending merkle paths along with evaluations.
    //
    //  The commitments made are: basecodewords, extension codewords, logd fri layers.
    //
    // with each successful query and valid decommitment, verifiers confidence in the proof increases.
    // decommit the basewords and the extension codewords and also the combination polynomial to satisfy the verifier.

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
//pub fn verify proof{collect all betas to check for fri layer}
//pub fn verify_queries{verify queries on the zipped value of base codewords and the extension codeowrds and also the terminal values }
//pub fn verify_frilayer{verify the consistency of all the fri_layers with the given betas}


//pub fn verify proof{collect all betas to check for fri layer}
pub fn verify_proof(
    num_of_queries: usize,
    maximum_random_int: u64,
    blow_up_factor: usize,//expansion_factor
    field: Field,
    fri_domains: &[Vec<FieldElement>],
    compressed_proof: &[Vec<u8>],
    terminal_processor:Vec<FieldElement>,
    terminal_instruction:Vec<FieldElement>,
    terminal_memory:Vec<FieldElement>,
    terminal_input:Vec<FieldElement>,
    terminal_output:Vec<FieldElement>,
height:FieldElement){
        let mut channel=Channel::new();
        // merkle root of the zipped base codewords
        let base_merkle_root =compressed_proof[0].clone();
        channel.send(base_merkle_root.clone());

        // get challenges for the extension columns
        let mut challenges_extension = vec![];

        for i in 0..10{
            let x = channel.receive_random_field_element(field);
            challenges_extension.push(x);
          
        }
        challenges_extension.push(channel.receive_random_field_element(field));

        let exten_merkle_root=compressed_proof[1].clone();
        channel.send(exten_merkle_root.clone());

        let combination_merkle_root=compressed_proof[2].clone();
        channel.send(combination_merkle_root);

        //commit to fri
        let mut fri_merkle_roots:Vec<Vec<u8>> = Vec::new();
        let mut betas:Vec<FieldElement>= Vec::new();
        // height should be power of 2
    // for fri_layer degree of the combination polynomial should be less then the height and the domain size will be the height*expansion_fact
    // fri.layer.len = 1+ log(height)/log2 
    let number = height.0;
    let base = 2.0;
    let log_base_2 = (number as f64).log2();
    let fri_layer_length:usize=(log_base_2+1 as f64) as usize;
    for i in 0..fri_layer_length-1 as usize{
        let beta = channel.receive_random_field_element(field);
        betas.push(beta);
        let fri_root=compressed_proof[3+i].clone();
        channel.send(fri_root.clone());
        fri_merkle_roots.push(fri_root); }
        let last_layer_free_term = 
        // last root will be 3+fri_layer_length-1
        //last term of the constant polynomial
        compressed_proof[2 as usize + fri_layer_length].clone();
    channel.send(last_layer_free_term.clone());
    // base_idx will be the point where the end of the compressed_proof indices for thr fri-layer_root commitment after this we have added the element and there authentication path we can see that in the utils.rs of the fri_layer decommit
    let mut base_idx = 3 as usize + fri_layer_length;
    for i in 0..num_of_queries{
        let idx = channel.receive_random_int(0, maximum_random_int, true) as usize;
        // verify_queries
        verify_queries(
            base_idx + i,
            idx,
            blow_up_factor,
            field,
            &fri_merkle_roots,
            &fri_domains,
            compressed_proof,
            &betas,
            &mut channel,
            terminal_processor.clone(),
            terminal_instruction.clone(),
            terminal_memory.clone(),
            terminal_input.clone(),
            terminal_output.clone(),
            height,
            fri_layer_length
        );
        // why 46 ??// here 46 is consistence acc to stark101 6 commitment of the f(x), f(gx), f(g^2x) for it's elem and the authentication path and other 41 for the fri layer 4 for all 10 layers and 1 for the last layer the constant term
        // in our case it will be 8 (for the base_x , base_gx,extenion_x, extension_gx) + 4*(fri_layer_length -1)+1 for the constant term 
        base_idx+= 8 +(4*(fri_layer_length-1));
}
}
//pub fn verify_queries{verify queries on the zipped value of base codewords and the extension codeowrds and also the terminal values }
pub fn verify_queries( base_idx: usize,
    idx: usize,
    blow_up_factor: usize,
    field: Field,
    fri_merkle_roots: &[Vec<u8>],
    fri_domains: &[Vec<FieldElement>],
    compressed_proof: &[Vec<u8>],
    betas: &[FieldElement],
    channel: &mut Channel,
    terminal_processor:Vec<FieldElement>,
    terminal_instruction:Vec<FieldElement>,
    terminal_memory:Vec<FieldElement>,
    terminal_input:Vec<FieldElement>,
    terminal_output:Vec<FieldElement>,
height:FieldElement,
fri_layer_length:usize
){
    // length of the eval_domain
    let len =height.0 as usize *blow_up_factor ;
    let base_merkle_root = compressed_proof[0].clone();
    // doubt here how this compressed proof is storing the leaf and the proof_bytes // leaf of the zipped base codewords
    // doubt solved
     let base_x=compressed_proof[base_idx].clone();
     channel.send(base_x.clone());
     //doubt here this is the proof_bytes
     let base_x_auth=compressed_proof[base_idx + 1].clone();
     channel.send(base_x_auth.clone());
     assert!(MerkleTree::validate(
        base_merkle_root.clone(),
        base_x_auth,
        idx,
        base_x,
        len as usize
    ));
    let base_gx = compressed_proof[base_idx + 2].clone();
    channel.send(base_gx.clone());
    let base_gx_auth = compressed_proof[base_idx + 3].clone();
    channel.send(base_gx_auth.clone());
    assert!(MerkleTree::validate(
        base_merkle_root.clone(),
        base_gx_auth,
        idx + blow_up_factor,
        base_gx,
        len as usize
    ));
    let exten_merkle_root =compressed_proof[1].clone();
    let exten_x=compressed_proof[base_idx+4].clone();
     channel.send(exten_x.clone());
     
     let exten_x_auth=compressed_proof[base_idx + 5].clone();
     channel.send(exten_x_auth.clone());
     assert!(MerkleTree::validate(
        exten_merkle_root.clone(),
        exten_x_auth,
        idx,// doubt this idx should change or not acc to me merkle tree is different so idx will remain same
        exten_x,
        len as usize
    ));
    let exten_gx = compressed_proof[base_idx + 6].clone();
    channel.send(exten_gx.clone());
    let exten_gx_auth = compressed_proof[base_idx + 7].clone();
    channel.send(exten_gx_auth.clone());
    assert!(MerkleTree::validate(
        exten_merkle_root.clone(),
        exten_gx_auth,
        idx + blow_up_factor,
        exten_gx,
        len as usize
    ));
    //for inter table arguments constraints
    assert_eq!(terminal_processor[0], terminal_instruction[0]); //Tipa = Tppa
    assert_eq!(terminal_processor[1], terminal_memory[0]); //Tmpa = Tppa
    assert_eq!(terminal_processor[2], terminal_input[0]); //Tipa = Tea input
    assert_eq!(terminal_processor[3], terminal_output[0]); //Tipa = Tea output
    //@todo
    //let this be for now:- assert_eq!(Terminal_instruction[1], Tpea); //Tpea = program evaluation
    let intial_length=blow_up_factor*height.0 as usize;

    verify_fri_layers(
        base_idx + 8,
        idx,
        field,
        fri_merkle_roots,
        fri_domains,
        compressed_proof,
        betas,
        channel,
        len,
        fri_layer_length


    );
}
//pub fn verify_frilayer{verify the consistency of all the fri_layers with the given betas}
pub fn verify_fri_layers(
    base_idx: usize,
    idx: usize,
    field: Field,
    fri_merkle_roots: &[Vec<u8>],
    fri_domains: &[Vec<FieldElement>],
    compressed_proof: &[Vec<u8>],
    betas: &[FieldElement],
    channel: &mut Channel,
    intial_length:usize,
    fri_layer_length:usize)

{
    let mut lengths:Vec<usize> =vec![0_usize;fri_layer_length-1 as usize];
    for i in 0..fri_layer_length-1 as usize{
        // let len_value =intial_length/2_usize.pow(i as u32);
        // length.push(len_value);
        lengths[i]=intial_length/2_usize.pow(i as u32);
     }
    for i in 0..fri_layer_length-1 as usize{
        let length = lengths[i];
        let elem_idx=idx%length;
        let elem = compressed_proof[base_idx + 4 * i].clone();
        channel.send(elem.clone());
        let elem_proof = compressed_proof[base_idx + 4 * i + 1].clone();
        channel.send(elem_proof.clone());
        let merkle_root = if i == 0 {
            compressed_proof[2].clone()
        } else {
            fri_merkle_roots[i - 1].clone()
        };
        if i != 0 {
            // checking fri polynomials consistency.
            let prev_elem =
                FieldElement::from_bytes(&compressed_proof[base_idx + 4 * (i - 1)].clone());
            let prev_sibling =
                FieldElement::from_bytes(&compressed_proof[base_idx + 4 * (i - 1) + 2].clone());
            let two = FieldElement::new(2, field);
            let computed_elem = (prev_elem + prev_sibling) / two
                + (betas[i - 1] * (prev_elem - prev_sibling)
                    / (two * fri_domains[i - 1][idx % lengths[i-1]]));
            assert!(computed_elem.0 == FieldElement::from_bytes(&elem).0);
        }
        assert!(MerkleTree::validate(
            merkle_root.clone(),
            elem_proof,
            elem_idx,
            elem.clone(),
            length,
        ));
        let sibling_idx = (idx + length / 2) % length;
        let sibling = compressed_proof[base_idx + 4 * i + 2].clone();
        channel.send(sibling.clone());
        let sibling_proof = compressed_proof[base_idx + 4 * i + 3].clone();
        channel.send(sibling_proof.clone());

        assert!(MerkleTree::validate(
            merkle_root,
            sibling_proof,
            sibling_idx,
            sibling.clone(),
            length,
        ));

    }

}



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

