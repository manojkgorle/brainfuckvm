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

    pub fn pad(&mut self) {
        let matrix = &mut self.table.matrix;
        let field = self.table.field;
        let zero = FieldElement::zero(field);
        while matrix.len() & (matrix.len() - 1) != 0 {
            matrix.push(vec![zero;1]);
        }
    }

    pub fn extend_column_ea(&mut self, randFieldElem: u128, challenge: FieldElement){
        let mut ea = FieldElement::new(randFieldElem,self.table.field); // take randFieldElem as zero if no random secret implementation
        self.table.matrix[0][1] = ea; 
        for i in 0..self.table.length-1 {
            ea = self.table.matrix[i as usize][1]*challenge + self.table.matrix[(i+1) as usize][0];
            self.table.matrix[(i+1) as usize][1] = ea; 
        //Tea = IOTable.matrix[length-1][1] 
        }
    }

}