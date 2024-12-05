use crate::fields::{Field, FieldElement};
use crate::univariate_polynomial::Polynomial;
use super::processor::Indices as ProcessorIndices;
use super::Table;
use super::{roundup_npow2, derive_omicron};
pub struct Memory {
    pub table : Table,
}

pub enum Indices {
    // Named indices for base columns
    Cycle,
    MemoryPointer,
    MemoryValue,
    Dummy,
    // Named indices for extension columns
    PermutationArg,
}

impl Memory {

    pub fn new(field: Field, length:u128, num_randomizers: u128, generator: FieldElement, order: u128) -> Self {
        let base_width = 4;
        let full_width = base_width + 1;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let matrix = vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        let table = Table::new(field, base_width, full_width, length, num_randomizers, height, omicron, generator, order, matrix);
        Self { table: table }
    }

    pub fn pad(&mut self){
        let matrix = &mut self.table.matrix;
        let field = self.table.field;
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        // @todo scope to optimize.
        while matrix.len() & (matrix.len() - 1) != 0 {
            let last_row = matrix.last().unwrap();
            let mut new_row = vec![zero; 4];
            new_row[Indices::Cycle as usize] = last_row[Indices::Cycle as usize] + one;
            new_row[Indices::MemoryPointer as usize] = last_row[Indices::MemoryPointer as usize];
            new_row[Indices::MemoryValue as usize] = last_row[Indices::MemoryValue as usize];
            new_row[Indices::Dummy as usize] = zero;
            matrix.push(new_row);
        }
    }

    // processor matrix that goes here is unpadded
    pub fn derive_matrix(processor_matrix: &[Vec::<FieldElement>]) -> Vec<Vec<FieldElement>> {
        let field = processor_matrix[0][0].1;
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);

        let mut matrix: Vec<Vec<FieldElement>> = processor_matrix
        .iter()
        .filter(|pt| pt[ProcessorIndices::CurrentInstruction as usize] != zero)
        .map(|pt| vec![
            pt[ProcessorIndices::Cycle as usize].clone(),
            pt[ProcessorIndices::MemoryPointer as usize].clone(),
            pt[ProcessorIndices::MemoryValue as usize].clone(),
            zero, // Equivalent to 'zero' in Python
        ])
        .collect();

        // Sort by memory_pointer
        matrix.sort_by(|a, b| {
        a[Indices::MemoryPointer as usize]
            .partial_cmp(&b[Indices::MemoryPointer as usize])
            .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Insert dummy rows for smooth clock jumps
        //@todo dummy rows should not be added since the contraint that mv= mv* when cycle jumps more than 1 will be disturbed
        let mut i = 0;
        while i < matrix.len() - 1 {
            if matrix[i][Indices::MemoryPointer as usize] == matrix[i + 1][Indices::MemoryPointer as usize] && matrix[i + 1][Indices::Cycle as usize] != matrix[i][Indices::Cycle as usize] + one{
                matrix.insert(i + 1,vec![
                    matrix[i][Indices::Cycle as usize].clone() + one,
                    matrix[i][Indices::MemoryPointer as usize].clone(),
                    matrix[i][Indices::MemoryValue as usize].clone(),
                    one
                ],
            );
        }
        i += 1;
    }
    matrix
    }

    //the matrix taken here is padded
    pub fn extension_column_ppa(challenges: Vec<FieldElement>, challenge_indices: Vec<usize>, matrix: Vec<Vec<FieldElement>>)-> Vec<Vec<FieldElement>>{
        let e_matrix = matrix.clone();
        //@todo taking init 1 for now, change to random secret initial value which we check by difference constraint of Tmpa = Tppa
        let mut init = FieldElement::one(matrix.clone()[0][0].1);
        for mut row in e_matrix.clone() {
            let weighted_sum = row[Indices::Cycle as usize].clone() * challenges[challenge_indices[0]]
                + row[Indices::MemoryPointer as usize].clone() * challenges[challenge_indices[1]]
                + row[Indices::MemoryValue as usize] * challenges[challenge_indices[2]] - challenges[challenge_indices[3]];
            row.push(init*weighted_sum); 
            init = init*weighted_sum;
        }
        e_matrix
    }
    

    pub fn constraints(&self)-> Vec<Polynomial>{
        let interpolated = self.table.interpolate_columns(vec![0,1,2]);

        //Boundary constraints
        //interpolated[]
    }
}

// @todo test memory table padding

    

