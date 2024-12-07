use crate::fields::{Field, FieldElement};
use super::{Table, roundup_npow2, derive_omicron};
pub struct InstructionTable {
    pub table: Table,
}

pub enum Indices {
    // Named indices for base columns
    Address,
    CurrentInstruction,
    NextInstruction,
    // Named indices for extension columns
    PermutationArg,
    EvaluationArg,
}

impl InstructionTable {
    pub fn new(field: Field, length:u128,  generator: FieldElement, order: u128, matrix: Vec<Vec<FieldElement>>) -> Self {
        let base_width = 3;
        let full_width = base_width + 2;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let mut gmatrix = vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        for i in 0..matrix.len() {
            gmatrix[i][Indices::Address as usize] = matrix[i][Indices::Address as usize];
            gmatrix[i][Indices::CurrentInstruction as usize] = matrix[i][Indices::CurrentInstruction as usize];
            gmatrix[i][Indices::NextInstruction as usize] = matrix[i][Indices::NextInstruction as usize];
        }
        let table = Table::new(field, base_width, full_width, length,  height, omicron, generator, order, gmatrix);
        Self { table }
    }

    pub fn get_table(&self) -> &Table {
        &self.table
    }

    // Note: Before padding initiate the matrix in table.
    // Add padding rows to convert the matrix derived from trace to a matrix of length of a power of 2
    pub fn pad(&mut self) {
        let zero = FieldElement::new(0, self.table.field);
        let length = self.table.length as usize;
        for i in 0..(self.table.height - self.table.length)as usize {
            let new_row = vec![self.table.matrix[length-1][Indices::Address as usize],zero, zero, zero, zero];
            self.table.matrix[self.table.length as usize + i] = new_row;
        }
    }

    // @todo add extension @soumyathakur44
}

// @todo test instruction table padding.
#[cfg(test)]
mod test_instruction {
    use super::Indices;
    use crate::vm::VirtualMachine;
    use crate::fields::Field;
    use crate::fields::FieldElement;
    #[test]
    fn test_padding() {
        let field = Field(18446744069414584321);
        let generator = field.generator();
        let order = 1<<32;
        let vm = VirtualMachine::new(field);
        let code2 = ">>[++-]<".to_string();
        let program = vm.compile(code2);
        let (runtime, _, _ ) = vm.run(&program, "".to_string());
        let (processor_matrix, _memory_matrix, instruction_matrix, _input_matrix, _output_matrix) = vm.simulate(&program, "".to_string());
        assert_eq!(instruction_matrix.len(), processor_matrix.len() + program.len());
        let ilen = instruction_matrix.len();
        let mut instruction_table = super::InstructionTable::new(field, instruction_matrix.len() as u128, generator, order, instruction_matrix.clone());
        instruction_table.pad();
        assert!(instruction_table.table.matrix.len().is_power_of_two());
        assert_eq!(instruction_matrix[ilen - 1][Indices::Address as usize], instruction_table.table.matrix.last().unwrap()[Indices::Address as usize]);
        assert_eq!(instruction_table.table.matrix.last().unwrap()[Indices::CurrentInstruction as usize], FieldElement::zero(field));
        assert_eq!(instruction_table.table.matrix.last().unwrap()[Indices::NextInstruction as usize], FieldElement::zero(field));
    }
}