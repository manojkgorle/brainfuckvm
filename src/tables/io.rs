use crate::fields::{Field, FieldElement};
use crate::tables::Table;
use crate::tables::{roundup_npow2, derive_omicron};
pub struct IOTable {
    pub table: Table,
}

pub enum Indices {
    // Named indices for base columns
    Column,
    // Named indices for extension columns
    EvaluationArg
}

impl IOTable {
    pub fn new(field: Field, length:u128,  generator: FieldElement, order: u128, matrix: Vec<Vec<FieldElement>>) -> Self {
        let base_width = 1;
        let full_width = base_width + 1;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let mut gmatrix = vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        for i in 0..matrix.len() {
            gmatrix[i][Indices::Column as usize] = matrix[i][Indices::Column as usize];
        
        }
        let table = Table::new(field, base_width, full_width, length, height, omicron, generator, order, gmatrix);
        Self { table }
    }

    // for padding we need to initialize the column column with zero. which is already done in the table.new.
    pub fn pad(&mut self) {}

    // add extension columns. @soumyathakur44
}

#[cfg(test)]
mod test_io {
    use super::{IOTable, Indices};
    use crate::fields::{Field,FieldElement};
    use crate::tables::instruction;
    use crate::vm::VirtualMachine;
    #[test]
    fn test_padding() {
        let field = Field(18446744069414584321);
        let generator = field.generator();
        let order = 1<<32;
        let vm = VirtualMachine::new(field);
        let code2 = ">>[++-]<,.,.".to_string();
        let program = vm.compile(code2);
        let (runtime, _, _ ) = vm.run(&program, "ab".to_string());
        let (processor_matrix, _memory_matrix, instruction_matrix, input_matrix, output_matrix) = vm.simulate(&program, "ab".to_string());
        let mut input_table = IOTable::new(field, input_matrix.len() as u128, generator, order, input_matrix);
        let mut output_table = IOTable::new(field, output_matrix.len() as u128, generator, order, output_matrix);
        input_table.pad();
        output_table.pad();
        assert_eq!(input_table.table.matrix.len(), 2);
        assert_eq!(output_table.table.matrix.len(), 2);
    }
}