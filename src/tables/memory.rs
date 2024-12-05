use crate::fields::{Field, FieldElement};
use crate::univariate_polynomial::Polynomial;
use crate::univariate_polynomial::interpolate_lagrange_polynomials;
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

pub enum IndicesPoly{
    Boundary,
    Transition,
    Terminal,
    Difference,
}

pub enum ChallengeIndices{
    A, B, C, D, E, F, Alpha, Beta, Delta, Gamma, Eta
}

impl Memory {
    pub fn new(field: Field, length:u128, generator: FieldElement, order: u128) -> Self {
        let base_width = 4;
        let full_width = base_width + 1;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let matrix= vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        let table = Table::new(field, base_width, full_width, length,  height, omicron, generator, order, matrix);
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
        //@todo i have written all extension and constraint functions assuming no dummy rows, i.e. same height as processor table
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
    pub fn extend_column_ppa(&mut self, challenges: Vec<FieldElement>){
        //@todo taking init 1 for now, change to random secret initial value which we check by difference constraint of Tmpa = Tppa
        let mut ppa = FieldElement::one(self.table.field);
        self.table.matrix[0].push(ppa); 
        for i in 0..self.table.length-1 {
            let weighted_sum = self.table.matrix[i as usize][Indices::Cycle as usize].clone() * challenges[ChallengeIndices::D as usize]
                + self.table.matrix[i as usize][Indices::MemoryPointer as usize].clone() * challenges[ChallengeIndices::E as usize]
                + self.table.matrix[i as usize][Indices::MemoryValue as usize].clone() * challenges[ChallengeIndices::F as usize] - challenges[ChallengeIndices::Beta as usize];
            self.table.matrix[(i+1) as usize].push(ppa*weighted_sum); 
            ppa = ppa*weighted_sum;
    }
    }
    
    //this is after padding and extension
    pub fn generate_AIR(&self, challenges: Vec<FieldElement>)-> Vec<Polynomial>{
        let interpolated = self.table.clone().interpolate_columns(vec![Indices::Cycle as u128, Indices::MemoryPointer as u128, Indices::MemoryValue as u128, Indices::PermutationArg as u128]);
        let CLK = interpolated[Indices::Cycle as usize].clone();
        let MP = interpolated[Indices::MemoryPointer as usize].clone();
        let MV = interpolated[Indices::MemoryValue as usize].clone();
        let PPA = interpolated[Indices::PermutationArg as usize].clone();

        let mut next_interpolated:Vec<Polynomial>=Vec::new();
        if self.table.height ==0{
            let poly=Polynomial::new_from_coefficients(vec![FieldElement::zero(Field::new(self.table.field.0))]);
            next_interpolated.push(poly);
           return next_interpolated;
        }
        
        let mut omicron_domain:Vec<FieldElement>=Vec::new();
        omicron_domain.push(self.table.omicron.pow(self.table.height-1));
        for i in 0..self.table.height-1{
            omicron_domain.push(self.table.omicron.pow(i));
        }

        for c in 0..self.table.matrix[0].len(){
            let mut trace:Vec<FieldElement>=Vec::new();
            for row in self.table.matrix.iter(){
                trace.push(row[c as usize]);
            }       
        let mut values:Vec<FieldElement>=Vec::new();
           
        values=trace.clone();
        if values.len()!=omicron_domain.len(){
            panic!("length of domain and values are unequal");
        };
        println!("domain ={:?}", omicron_domain);
        println!("values ={:?}", values);

        let poly= interpolate_lagrange_polynomials(omicron_domain.clone(), values);
        println!("poly ={:?}", poly);
            next_interpolated.push(poly);
        }

        let CLK_next = next_interpolated[Indices::Cycle as usize].clone();
        let MP_next = next_interpolated[Indices::MemoryPointer as usize].clone();
        let MV_next = next_interpolated[Indices::MemoryValue as usize].clone();
        let PPA_next = next_interpolated[Indices::PermutationArg as usize].clone();
        let one = Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        let mut AIR = vec![];

        //Boundary constraint: clk = mp = mv = 0, ppa = 1 (random secret using diff constr if we use it, not needed rn) 
        let boundaryAIR = CLK.clone()+MP.clone()+MV.clone()+PPA.clone()-Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        AIR.push(boundaryAIR);

        //Transition constraints: * == next
        //1. (mp+1-mp*).(mp-mp*)
        //2. (mp-mp*).mv*
        //3. (mp-mp*).(mv-mv*)
        //4. (clk - 1 - clk*).(mv*-mv)
        //5. ppa.(d.clk + e.mp + f.mv - beta) - ppa*
        let transitionAIR= (MP.clone()+one.clone()-MP_next.clone())*(MP.clone()-MP_next.clone()) 
        + (MP.clone()-MP_next.clone())*MV_next.clone() + (MP.clone()-MP_next.clone())*(MV.clone()-MV_next.clone())
        + (CLK.clone()-one.clone()-CLK_next)*(MV_next.clone()-MV.clone())
        + PPA.clone()*(CLK.scalar_mul(challenges[ChallengeIndices::D as usize]) + MP.scalar_mul(challenges[ChallengeIndices::E as usize]) + MV.scalar_mul(challenges[ChallengeIndices::F as usize]) - Polynomial::new_from_coefficients(vec![challenges[ChallengeIndices::Beta as usize]])) - PPA_next;
        AIR.push(transitionAIR);

        //Terminal constraints:
        //@todo Tppa = Tmpa --> include a constraint for this?
        //ppa.(d.clk+e.mp+f.mv-beta)-Tppa 
        //@todo Tppa and Tmpa given by prover, for now just taking it as empty polynomials to write constraint without error
        let Tppa = Polynomial::new_from_coefficients(vec![]);
        let Tmpa = Polynomial::new_from_coefficients(vec![]);
        let terminalAIR = PPA*(CLK.scalar_mul(challenges[ChallengeIndices::D as usize]) + MP.scalar_mul(challenges[ChallengeIndices::E as usize]) + MV.scalar_mul(challenges[ChallengeIndices::F as usize]) - Polynomial::new_from_coefficients(vec![challenges[ChallengeIndices::Beta as usize]])) - Tppa.clone()
        + Tppa - Tmpa;
        AIR.push(terminalAIR);

        AIR
    }

    pub fn generate_zerofier(&self)-> Vec<Polynomial>{  
        let mut zerofiers = vec![];
        let omicron = self.table.omicron;
        let x = Polynomial::new_from_coefficients(vec![FieldElement::zero(self.table.field), FieldElement::one(self.table.field)]);

        let boundary_zerofier = x - Polynomial::new_from_coefficients(vec![omicron.clone().pow(0)]);
        zerofiers.push(boundary_zerofier);

        //let transition_zerofier = Polynomial::new

        zerofiers
    }

}

//@todo test extend column ppa
//@todo test generate AIR
//@todo test generate zerofier
//@todo test generate quotient
//@todo test memory table padding

    

