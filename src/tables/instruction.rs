use crate::fields::{Field, FieldElement};
use super::{Table, roundup_npow2, derive_omicron};
use crate::univariate_polynomial::*;
pub struct InstructionTable {
    pub table: Table,
}

pub enum Indices {
    // Named indices for base columns
    Address, //ip
    CurrentInstruction, //ci
    NextInstruction, //ni
    // Named indices for extension columns
    PermutationArg,
    EvaluationArg,
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
    
    pub fn extend_column(&mut self, randFieldElem: u128, challenges: Vec<FieldElement>){
        let mut ppa = FieldElement::new(randFieldElem,self.table.field); 
        //take randFieldElement = 1 when not implementing random secret diff constraint 
        let mut pea = FieldElement::zero(self.table.field);

        self.table.matrix[0 as usize][Indices::PermutationArg as usize] = ppa; 
        self.table.matrix[0 as usize][Indices::EvaluationArg as usize] = pea;
        //@todo set initial value of first row of ppa and pea

        for i in 0..self.table.length-1 {
            let weighted_sum = self.table.matrix[(i+1) as usize][Indices::Address as usize] * challenges[ChallengeIndices::A as usize]
                + self.table.matrix[(i+1) as usize][Indices::CurrentInstruction as usize] * challenges[ChallengeIndices::B as usize]
                + self.table.matrix[(i+1) as usize][Indices::NextInstruction as usize] * challenges[ChallengeIndices::C as usize];
            if self.table.matrix[(i+1) as usize][Indices::Address as usize]==self.table.matrix[i as usize][Indices::Address as usize]{
                ppa *= (weighted_sum - challenges[ChallengeIndices::Alpha as usize]);
                self.table.matrix[(i+1) as usize][Indices::PermutationArg as usize] = ppa; 
                self.table.matrix[(i+1) as usize][Indices::EvaluationArg as usize] = pea;
            }
            else{
                self.table.matrix[(i+1) as usize][Indices::PermutationArg as usize] = ppa; 
                self.table.matrix[(i+1) as usize][Indices::EvaluationArg as usize] = pea*challenges[ChallengeIndices::Eta as usize] + weighted_sum;
            }
        }
    }

    pub fn generate_air(&self, challenges: Vec<FieldElement>)-> Vec<Polynomial>{
        let interpolated = self.table.clone().interpolate_columns(vec![Indices::Address as u128, Indices::CurrentInstruction as u128, Indices::NextInstruction as u128, Indices::PermutationArg as u128, Indices::EvaluationArg as u128]);
        let IP = interpolated[Indices::Address as usize].clone();
        let CI = interpolated[Indices::CurrentInstruction as usize].clone();
        let NI = interpolated[Indices::NextInstruction as usize].clone();
        let PPA = interpolated[Indices::PermutationArg as usize].clone();
        let PEA = interpolated[Indices::EvaluationArg as usize].clone();

        let next_interpolated = self.table.clone().next_interpolate_columns(vec![Indices::Address as u128, Indices::CurrentInstruction as u128, Indices::NextInstruction as u128, Indices::PermutationArg as u128, Indices::EvaluationArg as u128]);

        let IP_next = next_interpolated[Indices::Address as usize].clone();
        let CI_next = next_interpolated[Indices::CurrentInstruction as usize].clone();
        let NI_next = next_interpolated[Indices::NextInstruction as usize].clone();
        let PPA_next = next_interpolated[Indices::PermutationArg as usize].clone();
        let PEA_next = next_interpolated[Indices::EvaluationArg as usize].clone();

        let one = Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        let mut AIR = vec![];

        //Boundary constraint: ip=0
        //@todo ppa and pea initial value from extended fn, see once

        let boundaryAIR = IP.clone()
        +PPA.clone()-Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)])
        +PEA.clone()-Polynomial::new_from_coefficients(vec![FieldElement::zero(self.table.field)]);;
        //@todo check this once!! initial value is not zero and one, set it to req value
        AIR.push(boundaryAIR);

        //Transition constraints: * == next
        //1. (ip-ip*).(ip*-ip-1)
        //2. (ip-ip*).(ni-ci*)
        //3. (ip*-ip-1).(ci*-ci)
        //4. (ip*-ip-1).(ni*-ni)
        //5. (ip+1-ip*).(ppa* - ppa.(a.ip*+b.ci*+c.ni*-alpha))
        //6. (ip-ip*).(ppa*-ppa)
        //7. (ip*-ip-1).(pea*-pea)
        //8. (ip*-ip).(pea* - pea.eta - (a.ip*+b.ci*+c.ni*))

        let transitionAIR= (IP.clone()-IP_next.clone())*(IP_next.clone()-IP.clone()-one.clone()) 
        + (IP.clone()-IP_next.clone())*(NI.clone()-CI_next.clone()) 
        + (IP_next.clone()-IP.clone()-one.clone())*(CI_next.clone()-CI.clone())
        + (IP_next.clone()-IP.clone()-one.clone())*(NI_next.clone()-NI.clone())
        + (IP.clone()+one.clone()-IP_next.clone())*(PPA_next.clone() - PPA.clone()*(IP_next.clone().scalar_mul(challenges[ChallengeIndices::A as usize])+CI_next.clone().scalar_mul(challenges[ChallengeIndices::B as usize])+NI_next.clone().scalar_mul(challenges[ChallengeIndices::C as usize])- Polynomial::new_from_coefficients(vec![challenges[ChallengeIndices::A as usize]])))
        + (IP.clone()-IP_next.clone())*(PPA_next.clone()-PPA.clone())
        + (IP_next.clone()-IP.clone()-one)*(PEA_next.clone()-PEA.clone())
        + (IP.clone()-IP_next.clone())*(PEA_next.clone() - PEA.clone().scalar_mul(challenges[ChallengeIndices::Eta as usize])-(IP_next.clone().scalar_mul(challenges[ChallengeIndices::A as usize])+CI_next.clone().scalar_mul(challenges[ChallengeIndices::B as usize])+NI_next.clone().scalar_mul(challenges[ChallengeIndices::C as usize])));
        AIR.push(transitionAIR);

        //Terminal constraints:
        //@todo Tppa = Tipa --> include a constraint for this?
        //ppa - Tppa
        //pea - Tpea
        //@todo Tppa and Tipa given by prover, for now just taking it as empty polynomials to write constraint without error
        //@todo Tpea is computed locally by verifier, taking empty polynomial for now
        
        let Tppa = Polynomial::new_from_coefficients(vec![]);
        let Tipa = Polynomial::new_from_coefficients(vec![]);
        let Tpea = Polynomial::new_from_coefficients(vec![]);
        let terminalAIR = PPA.clone() - Tppa.clone()
        + PEA.clone() - Tpea.clone()
        + Tppa.clone() - Tipa.clone();
        //@todo separate Tppa - Tipa term as it will cancel out 
        AIR.push(terminalAIR);

        AIR
    }

    pub fn generate_zerofier(&self)-> Vec<Polynomial>{  
        let mut zerofiers = vec![];
        let omicron = self.table.omicron;
        let x = Polynomial::new_from_coefficients(vec![FieldElement::zero(self.table.field), FieldElement::one(self.table.field)]);

        let boundary_zerofier = x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(0)]);
        zerofiers.push(boundary_zerofier);

        let mut transition_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length-1{
            transition_zerofier*=x.clone()-Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
        }
        zerofiers.push(transition_zerofier);

        let terminal_zerofier = x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(self.table.length-1)]);

        zerofiers
    }

    pub fn generate_quotients(&self, challenges: Vec<FieldElement>)->Vec<Polynomial>{
        let mut quotients = vec![];
        let AIR = self.generate_air(challenges);
        let zerofiers = self.generate_zerofier();

        for i in 0..AIR.len(){
            quotients.push(AIR[i].clone().q_div(zerofiers[i].clone()).0);
        }
        quotients
    }
    
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