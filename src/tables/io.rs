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
    pub fn new(field: Field, length:u128, num_randomizers: u128, generator: FieldElement, order: u128) -> Self {
        let base_width = 1;
        let full_width = base_width + 1;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let matrix = vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        let table = Table::new(field, base_width, full_width, length, height, omicron, generator, order, matrix);
        Self { table: table }
    }

    pub fn pad(&mut self) {
        let matrix = &mut self.table.matrix;
        let field = self.table.field;
        let zero = FieldElement::zero(field);
        while matrix.len() & (matrix.len() - 1) != 0 {
            matrix.push(vec![zero;1]);
        }
    }

    // add extension columns.
}