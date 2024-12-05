
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
        Self { table }
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

    //the matrix taken here is padded
    pub fn extend_columns(&mut self, challenges: Vec<FieldElement>){
        //@todo taking init 1 for now, change to random secret initial value which we check by difference constraint of Tmpa = Tppa
        let mut ipa = FieldElement::one(self.table.field);
        let mut mpa = FieldElement::one(self.table.field);
        let mut iea = FieldElement::zero(self.table.field);
        let mut oea = FieldElement::zero(self.table.field);

        self.table.matrix[0].push(ipa);
        for i in 0..self.table.length-1 {
            let weighted_sum = self.table.matrix[i as usize][Indices::InstructionPointer as usize].clone() * challenges[ChallengeIndices::A as usize]
                + self.table.matrix[i as usize][Indices::CurrentInstruction as usize].clone() * challenges[ChallengeIndices::B as usize]
                + self.table.matrix[i as usize][Indices::NextInstruction as usize].clone() * challenges[ChallengeIndices::C as usize] - challenges[ChallengeIndices::Alpha as usize];
            self.table.matrix[(i+1) as usize].push(ipa*weighted_sum); 
            ipa = ipa*weighted_sum;
    }
        self.table.matrix[0].push(mpa);
        for i in 0..self.table.length-1 {
            let weighted_sum = self.table.matrix[i as usize][Indices::Cycle as usize].clone() * challenges[ChallengeIndices::D as usize]
                + self.table.matrix[i as usize][Indices::MemoryPointer as usize].clone() * challenges[ChallengeIndices::E as usize]
                + self.table.matrix[i as usize][Indices::MemoryValue as usize].clone() * challenges[ChallengeIndices::F as usize] - challenges[ChallengeIndices::Beta as usize];
            self.table.matrix[(i+1) as usize].push(mpa*weighted_sum); 
            mpa = mpa*weighted_sum;
    }

        self.table.matrix[0].push(iea);
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, self.table.field) };
        for i in 0..self.table.length-1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f(',') == ci {
                iea = iea*challenges[ChallengeIndices::Gamma as usize] + self.table.matrix[i as usize][Indices::MemoryValue as usize];
                self.table.matrix[(i+1) as usize].push(iea); 
            }
            else{
                self.table.matrix[(i+1) as usize].push(iea);
            }
    }

        self.table.matrix[0].push(oea);
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, self.table.field) };
        for i in 0..self.table.length-1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('.') == ci {
                oea = oea*challenges[ChallengeIndices::Delta as usize] + self.table.matrix[i as usize][Indices::MemoryValue as usize];
                self.table.matrix[(i+1) as usize].push(oea); 
            }
            else{
                self.table.matrix[(i+1) as usize].push(oea);
            }
    }

    //define a selector polynomial for a specific instruction.
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
                acc *= indeterminate - f(c);
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
                    acc *= indeterminate - f(*c);
                }
            }
            deselectors.push((target_char, acc));
        }
    
        deselectors
    }
}
}


    //boundary constraints for the base coloumns
    // the values of instructionpermutaion ipa and mpa I am taking as 1
    // pub fn boundary_constraint(self)->Polynomial{

  

    

    

//  }











    

//}
    

// @todo test processor table padding