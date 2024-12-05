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
    pub fn new(field: Field, length:u128,  generator: FieldElement, order: u128) -> Self {
        let base_width = 3;
        let full_width = base_width + 2;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let matrix = vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        let table = Table::new(field, base_width, full_width, length,  height, omicron, generator, order, matrix);
        Self { table: table }
    }

    pub fn get_table(&self) -> &Table {
        &self.table
    }

    // Note: Before padding initiate the matrix in table.
    // Add padding rows to convert the matrix derived from trace to a matrix of length of a power of 2
    pub fn pad(&mut self) {
        let zero = FieldElement::new(0, self.table.field);
        for _ in 0..(self.table.height - self.table.length)  {
            let new_row = vec![self.table.matrix.last().unwrap()[Indices::Address as usize],zero, zero, zero, zero];
            self.table.matrix.push(new_row);
        }
    }
}

// @todo test instruction table padding.