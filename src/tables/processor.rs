
use crate::fields::{FieldElement, Field};
use super::Table;
use super::{roundup_npow2, derive_omicron};
use crate::univariate_polynomial::*;

pub struct ProcessorTable {
    table: Table
}

pub enum Indices {
    // Named indices for registers
    Cycle,
    InstructionPointer,
    CurrentInstruction,
    NextInstruction,
    MemoryPointer,
    MemoryValue,
    MemoryValueInvers,
    // Named indices for extension columns
    InstructionPermutaion,
    MemoryPermuation,
    InputEvaluation,
    OutputEvaluation,
}
#[allow(non_camel_case_types)]
pub enum IndicesPoly {
    //i0 = [, 
    //i1 = ],
    // i2 = <,
    // i3 = >,
    // i4 = +,
    // i5 = -,
    // i6 = ",",
    // i7 = ".",
    Boundary,
    Transition_i0,
    Transition_i1,
    Transition_i2,
    Transition_i3,
    Transition_i4,
    Transition_i5,
    Transition_i6,
    Transition_i7,
    Transition_all,
    //clk*-clk-1
    //inv.(1-inv.mv)
    //ci.(ipa.(a.ip+b.ci+c.ni-alpha)-ipa*)+(ipa*-ipa).deselector
    //mpa.(d.clk+e.mp+f.mv-beta)-mpa*
    Transition_iea, 
    Transition_oea,
    Terminal,
    Difference,
}
pub enum ChallengeIndices{
    A, B, C, D, E, F, Alpha, Beta, Delta, Gamma, Eta
}

impl ProcessorTable {
    pub fn new(field: Field, length:u128,  generator: FieldElement, order: u128) -> Self {
        let base_width = 7;
        let full_width = base_width + 4;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let matrix = vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        let table = Table::new(field, base_width, full_width, length,  height, omicron, generator, order, matrix);
        Self { table: table }
    }

    // Note: Before padding initiate the matrix in table.
    // Add padding rows to convert the matrix derived from trace to a matrix of length of a power of 2
    pub fn pad(&mut self) {
        let matrix = &mut self.table.matrix;
        let field = self.table.field;
        let zero = FieldElement::zero(field);
        while matrix.len() & (matrix.len() -1 ) != 0 {
            let last_row = matrix.last().unwrap();
            let mut new_row = vec![FieldElement::zero(self.table.field); 7];
            // Increment cycle by one for every new padded row.
            new_row[Indices::Cycle as usize] = last_row[Indices::Cycle as usize] + FieldElement::one(self.table.field);
            // keep the instruction pointer same as the last row in every padded row.
            new_row[Indices::InstructionPointer as usize] = last_row[Indices::InstructionPointer as usize];
            // keep the current instruction and next instruction as zero.
            new_row[Indices::CurrentInstruction as usize] = zero;
            new_row[Indices::NextInstruction as usize] = zero;
            // keep the memory pointer and memory value as last.
            new_row[Indices::MemoryPointer as usize] = last_row[Indices::MemoryPointer as usize];
            new_row[Indices::MemoryValue as usize] = last_row[Indices::MemoryValue as usize];
            new_row[Indices::MemoryValueInvers as usize] = last_row[Indices::MemoryValueInvers as usize];
            matrix.push(new_row);
        }
    }


    // define a selector polynomial for a specific instruction.
    //todo for a set of instructions.
   pub fn deselector_polynomial(
        instruction: char, 
        indeterminate: FieldElement, 
        field: Field,
    ) -> FieldElement {
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, field) };
        let mut acc = FieldElement::new(1, field); // Start with the multiplicative identity (1)
        
        for c in "[]<>,.+-".chars() {
            if c != instruction {
                acc = acc * (indeterminate - f(c));
            }
        }
          acc
    }
    // define a selector polynomial for a valid set
   pub fn universal_deselector(
        indeterminate: FieldElement, 
        field: Field,
        char:Vec<char>
    ) -> Vec<(char, FieldElement)> {
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, field) };
        let mut deselectors = Vec::new();
    
        for target_char in "[]<>,.+-".chars() {
            let mut acc = FieldElement::new(1, field); // Start with the multiplicative identity (1)
    
            for c in char.iter() {
                if *c != target_char {
                    acc = acc * (indeterminate - f(*c));
                }
            }
    
            deselectors.push((target_char, acc));
        }
    
        deselectors
    }
 // this is after padding and extension
//  pub fn generate_AIR(&self,challenges:Vec<FieldElement>)->Vec<Polynomial>{
//     let interpolated = self.table.clone().interpolate_columns(vec![Indices::Cycle as u128, Indices::InstructionPointer as u128, Indices::CurrentInstruction as u128, Indices::NextInstruction as u128, Indices::MemoryPointer as u128, Indices::MemoryValue as u128, Indices::MemoryValueInvers as u128, Indices::InstructionPermutaion as u128, Indices::MemoryPermuation as u128, Indices::InputEvaluation as u128, Indices::OutputEvaluation as u128]);
//     let clk=interpolated[Indices::Cycle as usize].clone();
//     let ip=interpolated[Indices::InstructionPointer as usize].clone();
//     let ci=interpolated[Indices::CurrentInstruction as usize].clone();
//     let ni=interpolated[Indices::NextInstruction as usize].clone();
//     let mp=interpolated[Indices::MemoryPointer as usize].clone();
//     let mv=interpolated[Indices::MemoryValue as usize].clone();
//     let inv_mv=interpolated[Indices::MemoryValueInvers as usize].clone();
//     let ipa=interpolated[Indices::InstructionPermutaion as usize].clone();
//     let mpa=interpolated[Indices::MemoryPermuation as usize].clone();
//     let iea=interpolated[Indices::InputEvaluation as usize].clone();
//     let oea=interpolated[Indices::OutputEvaluation as usize].clone();
//     let mut next_poly:Vec<Polynomial>=Vec::new();

  

    

    

//  }











    

}
    

// @todo test processor table padding