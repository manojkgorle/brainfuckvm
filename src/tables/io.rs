#![allow(unused_variables)]
use crate::fields::{Field, FieldElement};
use crate::tables::Table;
use crate::tables::{derive_omicron, roundup_npow2};
pub struct IOTable {
    pub table: Table,
}

pub enum Indices {
    // Named indices for base columns
    Column,
    // Named indices for extension columns
    EvaluationArg,
}

impl IOTable {
    pub fn new(
        field: Field,
        length: u128,
        generator: FieldElement,
        order: u128,
        matrix: Vec<Vec<FieldElement>>,
    ) -> Self {
        let base_width = 1;
        let full_width = base_width + 1;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let mut gmatrix =
            vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        for i in 0..matrix.len() {
            gmatrix[i][Indices::Column as usize] = matrix[i][Indices::Column as usize];
        }
        let table = Table::new(
            field, base_width, full_width, length, height, omicron, generator, order, gmatrix,
        );
        Self { table }
    }

    pub fn pad(&mut self) {}

    pub fn extend_column_ea(
        &mut self,
        rand_field_elem: u128,
        challenge: FieldElement,
    ) -> Vec<FieldElement> {
        let mut ea = FieldElement::new(rand_field_elem, self.table.field); // take rand_field_elem as zero if no random secret implementation
        let mut terminal: Vec<FieldElement> = Vec::new();
        if(self.table.matrix.len()>0){
            self.table.matrix[0][1] = ea;
        for i in 0..self.table.length - 1 {
            ea = self.table.matrix[i as usize][1] * challenge
                + self.table.matrix[(i) as usize][0];
            self.table.matrix[(i + 1) as usize][1] = ea;
            //Tea = IOTable.matrix[length-1][1]
        }
        terminal.push(ea*challenge+self.table.matrix[self.table.length as usize-1][0]);
        }
        terminal
    }
}

#[cfg(test)]
mod test_io {
    use super::{IOTable, Indices};
    use crate::fields::{Field, FieldElement};
    use crate::tables::instruction;
    use crate::vm::VirtualMachine;
    #[test]
    fn test_padding() {
        let field = Field(18446744069414584321);
        let generator = field.generator();
        let order = 1 << 32;
        let vm = VirtualMachine::new(field);
        let code2 = ">>[++-]<,.,.".to_string();
        let program = vm.compile(code2);
        let (runtime, _, _) = vm.run(&program, "ab".to_string());
        let (processor_matrix, _memory_matrix, instruction_matrix, input_matrix, output_matrix) =
            vm.simulate(&program, "ab".to_string());
        let mut input_table = IOTable::new(
            field,
            input_matrix.len() as u128,
            generator,
            order,
            input_matrix,
        );
        let mut output_table = IOTable::new(
            field,
            output_matrix.len() as u128,
            generator,
            order,
            output_matrix,
        );
        input_table.pad();
        output_table.pad();

        for row in input_table.table.matrix.clone() {
            println!("{:?}", row)
        }
        for row in output_table.table.matrix.clone() {
            println!("{:?}", row)
        }
        assert_eq!(input_table.table.matrix.len(), 2);
        assert_eq!(output_table.table.matrix.len(), 2);
    }
}
