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
    
    fn compute_terminal(&self, challenges: &[FieldElement])-> FieldElement{
        let field = self.field;
        let iota = challenges[self.challenge_index];
        let mut acc = FieldElement::zero(field);

        for i in 0..self.symbols.len(){
            acc = iota*acc + FieldElement::new(self.symbols[i], field);
        }
        acc
    }

    fn select_terminal(&self, terminals: &[FieldElement])->FieldElement{
        terminals[self.terminal_index]
    }
}

pub struct ProgramEvaluationArgument{
    field: Field,
    challenge_indices: Vec<usize>,
    terminal_index: usize,
    program: Vec<u128>,
    //@todo soumya- if program is passed as a vec of field elements directly, then can skip padded program part
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

    fn compute_terminal(&self, challenges: &[FieldElement])->FieldElement{
        let field = self.field;
        let trimmed_challenges: Vec<FieldElement> = self
            .challenge_indices
            .iter()
            .map(|&i| challenges[i])
            .collect();
        let [a, b, c, eta]: [FieldElement; 4] = trimmed_challenges.try_into().unwrap();
        let mut running_sum = FieldElement::zero(field);
        let mut previous_address = -1_isize;

        let padded_program: Vec<FieldElement> = self
            .program
            .iter()
            .map(|&p| FieldElement::new(p,field))
            .collect();
        
        //@todo ci goes till last element of program in last step of running sum, ni is zero for the last step
        
        for i in 0..padded_program.len()-1 {
            let address = i;
            let current_instruction = padded_program[i];
            let next_instruction = padded_program[i + 1];

            if previous_address != address as isize {
                running_sum = running_sum * eta
                    + a * FieldElement::new(address as u128, field)
                    + b * current_instruction
                    + c * next_instruction;
                println!("{}:{}:{}:{}:{}", i, address, current_instruction.0, next_instruction.0, running_sum.0 );
            }
            previous_address = address as isize;
        }

        let index = padded_program.len() - 1;
        let address = FieldElement::new(index as u128, field);
        let current_instruction = padded_program[index];
        let next_instruction = FieldElement::zero(field);

        running_sum = running_sum * eta
            + a * address
            + b * current_instruction
            + c * next_instruction;
        println!("{}:{}:{}:{}:{}", index, address.0, current_instruction.0, next_instruction.0, running_sum.0 );
        running_sum
    }

    fn select_terminal(&self, terminals: &[FieldElement]) -> FieldElement {
        terminals[self.terminal_index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_evaluation_argument() {
        // Create a mock field
        let field = Field::new(17); // Example prime field with modulus 17

        // Define challenge index, terminal index, and symbols
        let challenge_index = 0;
        let terminal_index = 1;
        let symbols = vec![3, 5, 7]; // Example coefficients for the polynomial

        // Create challenges
        let challenges = vec![
            FieldElement::new(2, field), // iota
            FieldElement::new(4, field), // Some other challenge
        ];

        // Create terminals
        let terminals = vec![
            FieldElement::new(10, field),
            FieldElement::new(12, field),
        ];

        // Instantiate EvaluationArgument
        let eval_arg = EvaluationArgument::new(field, challenge_index, terminal_index, symbols);

        // Compute terminal using challenges
        let computed_terminal = eval_arg.compute_terminal(&challenges);

        // Select terminal from the provided list
        let selected_terminal = eval_arg.select_terminal(&terminals);

        // Check results
        assert_eq!(computed_terminal, selected_terminal);
    }

    #[test]
    fn test_program_evaluation_argument() {
        // Create a mock field
        let field = Field::new(17); // Example prime field with modulus 17
        
        // Define challenge indices, terminal index, and program
        let challenge_indices = vec![0, 1, 2, 3];
        let terminal_index = 0;
        let program = vec![1, 2, 3, 4]; // Example program instructions

        // Create challenges
        let challenges = vec![
            FieldElement::new(2, field), // a
            FieldElement::new(3, field), // b
            FieldElement::new(5, field), // c
            FieldElement::new(7, field), // eta
        ];

        // Create terminals
        let terminals = vec![FieldElement::new(4, field)]; // Terminal at index 0

        // Instantiate ProgramEvaluationArgument
        let prog_eval_arg = ProgramEvaluationArgument::new(
            field,
            challenge_indices,
            terminal_index,
            program,
        );

        // Compute terminal using challenges
        let computed_terminal = prog_eval_arg.compute_terminal(&challenges);
        println!("{}",computed_terminal.0);
        let selected_terminal = prog_eval_arg.select_terminal(&terminals);

        assert_eq!(computed_terminal, selected_terminal); 
    }
}



