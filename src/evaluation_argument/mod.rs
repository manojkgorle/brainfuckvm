use crate::fields::*;
use crate::univariate_polynomial::*;

#[derive(Debug, Clone)]

pub struct EvaluationArgument{
    field: Field,
    challenge_index: usize,
    terminal_index: usize,
    symbols: Vec<u128>,
}

impl EvaluationArgument{
    fn new(field:Field, challenge_index: usize, terminal_index: usize, symbols: Vec<u128>)->Self{
        Self {field, challenge_index, terminal_index, symbols}
    }
    
    fn compute_terminal(&self, challenges: &Vec<FieldElement>)-> FieldElement{
        let field = self.field;
        let iota = challenges[self.challenge_index];
        let mut acc = FieldElement::zero(field);

        for i in 0..self.symbols.len(){
            acc = iota.clone()*acc + FieldElement::new(self.symbols[i], field);
        }
        acc
    }

    fn select_terminal(&self, terminals: &Vec<FieldElement>)->FieldElement{
        terminals[self.terminal_index].clone()
    }
}

pub struct ProgramEvaluationArgument{
    field: Field,
    challenge_indices: Vec<usize>,
    terminal_index: usize,
    program: Vec<u128>,
}

impl ProgramEvaluationArgument{
    fn new(field: Field, challenge_indices: Vec<usize>, terminal_index: usize, program: Vec<u128>)->Self{
        Self {
            field,
            challenge_indices,
            terminal_index,
            program,
        }
    }

    fn compute_terminal(&self, challenges: &Vec<FieldElement>)->FieldElement{
        let field = self.field;
        let trimmed_challenges: Vec<FieldElement> = self
            .challenge_indices
            .iter()
            .map(|&i| challenges[i].clone())
            .collect();
        let [a, b, c, eta]: [FieldElement; 4] = trimmed_challenges.try_into().unwrap();
        let xfield = FieldElement::zero(field);
        let mut running_sum = FieldElement::zero(field);
        let mut previous_address = FieldElement::one(field);

        let padded_program: Vec<FieldElement> = self
            .program
            .iter()
            .map(|&p| FieldElement::new(p,field))
            .chain(std::iter::once(FieldElement::zero(field)))
            .collect();

        for i in 0..padded_program.len() - 1 {
            let address = FieldElement::new(i as u128, field);
            let current_instruction = padded_program[i].clone();
            let next_instruction = padded_program[i + 1].clone();

            if previous_address != address {
                running_sum = running_sum * eta.clone()
                    + a.clone() * address.clone()
                    + b.clone() * current_instruction
                    + c.clone() * next_instruction;
            }

            previous_address = address;
        }

        let index = padded_program.len() - 1;
        let address = FieldElement::new(index as u128, field);
        let current_instruction = padded_program[index].clone();
        let next_instruction = FieldElement::zero(field);

        running_sum = running_sum * eta.clone()
            + a.clone() * address
            + b.clone() * current_instruction
            + c.clone() * next_instruction;

        running_sum
    }

    fn select_terminal(&self, terminals: &[FieldElement]) -> FieldElement {
        terminals[self.terminal_index].clone()
    }
}



