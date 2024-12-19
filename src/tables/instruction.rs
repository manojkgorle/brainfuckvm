use super::{derive_omicron, roundup_npow2, Table};
use crate::fields::{Field, FieldElement};
use crate::univariate_polynomial::*;
pub struct InstructionTable {
    pub table: Table,
}

pub enum Indices {
    // Named indices for base columns
    Address,            //ip
    CurrentInstruction, //ci
    NextInstruction,    //ni
    // Named indices for extension columns
    PermutationArg,
    EvaluationArg,
}

pub enum IndicesPoly {
    Boundary,
    Transition,
    Terminal,
    Difference,
}
pub enum ChallengeIndices {
    A,
    B,
    C,
    D,
    E,
    F,
    Alpha,
    Beta,
    Delta,
    Gamma,
    Eta,
}

impl InstructionTable {
    pub fn new(
        field: Field,
        length: u128,
        generator: FieldElement,
        order: u128,
        matrix: &[Vec<FieldElement>],
    ) -> Self {
        let base_width = 3;
        let full_width = base_width + 2;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let mut gmatrix =
            vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        for i in 0..matrix.len() {
            gmatrix[i][Indices::Address as usize] = matrix[i][Indices::Address as usize];
            gmatrix[i][Indices::CurrentInstruction as usize] =
                matrix[i][Indices::CurrentInstruction as usize];
            gmatrix[i][Indices::NextInstruction as usize] =
                matrix[i][Indices::NextInstruction as usize];
        }
        let table = Table::new(
            field, base_width, full_width, length, height, omicron, generator, order, gmatrix,
        );
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
        for i in 0..(self.table.height - self.table.length) as usize {
            let new_row = vec![
                self.table.matrix[length - 1][Indices::Address as usize],
                zero,
                zero,
                zero,
                zero,
            ];
            self.table.matrix[self.table.length as usize + i] = new_row;
        }
    }
    #[warn(unused_variables)]
    pub fn extend_column(
        &mut self,
        _rand_field_elem: u128,
        challenges: Vec<FieldElement>,
    ) -> Vec<FieldElement> {
        let mut terminal: Vec<FieldElement> = Vec::new();
        let field = self.table.field;
        let mut ppa = FieldElement::one(field);
        let mut pea = self.table.matrix[(0) as usize][Indices::Address as usize]
            * challenges[ChallengeIndices::A as usize]
            + self.table.matrix[(0) as usize][Indices::CurrentInstruction as usize]
                * challenges[ChallengeIndices::B as usize]
            + self.table.matrix[(0) as usize][Indices::NextInstruction as usize]
                * challenges[ChallengeIndices::C as usize];

        self.table.matrix[(0) as usize][Indices::PermutationArg as usize] = ppa;
        self.table.matrix[(0) as usize][Indices::EvaluationArg as usize] = pea;

        for i in 0..self.table.length - 1 {
            let weighted_sum = self.table.matrix[(i + 1) as usize][Indices::Address as usize]
                * challenges[ChallengeIndices::A as usize]
                + self.table.matrix[(i + 1) as usize][Indices::CurrentInstruction as usize]
                    * challenges[ChallengeIndices::B as usize]
                + self.table.matrix[(i + 1) as usize][Indices::NextInstruction as usize]
                    * challenges[ChallengeIndices::C as usize];
            if self.table.matrix[(i + 1) as usize][Indices::Address as usize]
                == self.table.matrix[i as usize][Indices::Address as usize]
            {
                ppa *= weighted_sum - challenges[ChallengeIndices::Alpha as usize];
                self.table.matrix[(i + 1) as usize][Indices::PermutationArg as usize] = ppa;
                self.table.matrix[(i + 1) as usize][Indices::EvaluationArg as usize] = pea;
            } else {
                self.table.matrix[(i + 1) as usize][Indices::PermutationArg as usize] = ppa;
                pea = pea * challenges[ChallengeIndices::Eta as usize] + weighted_sum;
                self.table.matrix[(i + 1) as usize][Indices::EvaluationArg as usize] = pea;
            }
        }
        terminal.push(ppa);
        terminal.push(pea);
        terminal
    }

    pub fn generate_air(
        &self,
        challenges: Vec<FieldElement>,
        tppa: FieldElement,
        tpea: FieldElement,
    ) -> Vec<Polynomial> {
        let interpolated = self.table.clone().interpolate_columns(vec![
            Indices::Address as u128,
            Indices::CurrentInstruction as u128,
            Indices::NextInstruction as u128,
            Indices::PermutationArg as u128,
            Indices::EvaluationArg as u128,
        ]);
        let ip = interpolated[Indices::Address as usize].clone();
        let ci = interpolated[Indices::CurrentInstruction as usize].clone();
        let ni = interpolated[Indices::NextInstruction as usize].clone();
        let ppa = interpolated[Indices::PermutationArg as usize].clone();
        let pea = interpolated[Indices::EvaluationArg as usize].clone();

        let next_interpolated = self.table.clone().next_interpolate_columns(interpolated);

        let ip_next = next_interpolated[Indices::Address as usize].clone();
        let ci_next = next_interpolated[Indices::CurrentInstruction as usize].clone();
        let ni_next = next_interpolated[Indices::NextInstruction as usize].clone();
        let ppa_next = next_interpolated[Indices::PermutationArg as usize].clone();
        let pea_next = next_interpolated[Indices::EvaluationArg as usize].clone();

        let one = Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        let mut air = vec![];

        //Boundary constraint:
        //ip=0

        let boundaryair = ip.clone() + ppa.clone()
            - Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)])
            + pea.clone()
            - Polynomial::constant(
                self.table.matrix[(0) as usize][Indices::Address as usize]
                    * challenges[ChallengeIndices::A as usize]
                    + self.table.matrix[(0) as usize][Indices::CurrentInstruction as usize]
                        * challenges[ChallengeIndices::B as usize]
                    + self.table.matrix[(0) as usize][Indices::NextInstruction as usize]
                        * challenges[ChallengeIndices::C as usize],
            );

        air.push(boundaryair);

        //Transition constraints: * == next
        //1. (ip-ip*).(ip*-ip-1)
        //2. (ip-ip*).(ni-ci*)
        //3. (ip*-ip-1).(ci*-ci)
        //4. (ip*-ip-1).(ni*-ni)
        //5. (ip+1-ip*).(ppa* - ppa.(a.ip*+b.ci*+c.ni*-alpha))
        //6. (ip-ip*).(ppa*-ppa)
        //7. (ip*-ip-1).(pea*-pea)
        //8. (ip*-ip).(pea* - pea.eta - (a.ip*+b.ci*+c.ni*))

         let transitionair = (ip.clone() - ip_next.clone())
            * (ip_next.clone() - ip.clone() - one.clone())
            + (ip.clone() - ip_next.clone()) * (ni.clone() - ci_next.clone())
            + (ip_next.clone() - ip.clone() - one.clone()) * (ci_next.clone() - ci.clone())
            + (ip_next.clone() - ip.clone() - one.clone()) * (ni_next.clone() - ni.clone())
            + (ip.clone() + one.clone() - ip_next.clone())
                * (ppa_next.clone()
                    - ppa.clone()
                         * (ip_next
                             .clone()
                             .scalar_mul(challenges[ChallengeIndices::A as usize])
                            + ci_next
                                .clone()
                                .scalar_mul(challenges[ChallengeIndices::B as usize])
                            + ni_next
                                .clone()
                                .scalar_mul(challenges[ChallengeIndices::C as usize])
                            - Polynomial::new_from_coefficients(vec![
                                challenges[ChallengeIndices::Alpha as usize],
                            ])))
            + (ip.clone() - ip_next.clone()) * (ppa_next.clone() - ppa.clone())
            + (ip_next.clone() - ip.clone() - one) * (pea_next.clone() - pea.clone())
            + (ip.clone() - ip_next.clone())
                * (pea_next.clone()
                    - pea
                        .clone()
                        .scalar_mul(challenges[ChallengeIndices::Eta as usize])
                    - (ip_next
                        .clone()
                        .scalar_mul(challenges[ChallengeIndices::A as usize])
                        + ci_next
                            .clone()
                            .scalar_mul(challenges[ChallengeIndices::B as usize])
                        + ni_next
                            .clone()
                            .scalar_mul(challenges[ChallengeIndices::C as usize])));
        air.push(transitionair);

        //Terminal constraints:

        //ppa - Tppa
        //pea - tpea

        let terminalair =
            ppa.clone() - Polynomial::constant(tppa) + pea.clone() - Polynomial::constant(tpea);
        air.push(terminalair);

        air
    }

    pub fn generate_zerofier(&self) -> Vec<Polynomial> {
        let mut zerofiers = vec![];
        let omicron = self.table.omicron;
        let x = Polynomial::new_from_coefficients(vec![
            FieldElement::zero(self.table.field),
            FieldElement::one(self.table.field),
        ]);

        let boundary_zerofier =
            x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(0)]);
        zerofiers.push(boundary_zerofier);

        let mut transition_zerofier =
            Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        for i in 0..self.table.length - 1 {
            transition_zerofier *=
                x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
        }
        zerofiers.push(transition_zerofier);

        let terminal_zerofier = x.clone()
            - Polynomial::new_from_coefficients(vec![omicron.clone().pow(self.table.length - 1)]);
        zerofiers.push(terminal_zerofier);
        zerofiers
    }

    pub fn generate_quotients(
        &self,
        challenges: Vec<FieldElement>,
        tppa: FieldElement,
        tpea: FieldElement,
    ) -> Vec<Polynomial> {
        let mut quotients = vec![];
        let air = self.generate_air(challenges, tppa, tpea);
        let zerofiers = self.generate_zerofier();

        for i in 0..air.len() {
            assert_eq!(air[i].clone().q_div(zerofiers[i].clone()).1, Polynomial::constant(FieldElement::zero(self.table.field)));
            quotients.push(air[i].clone().q_div(zerofiers[i].clone()).0);
        }
        quotients
    }
}

#[cfg(test)]
mod test_instruction {
    use super::Indices;
    use crate::fields::Field;
    use crate::fields::FieldElement;
    use crate::tables::instruction::InstructionTable;
    use crate::tables::io::IOTable;
    use crate::tables::memory::MemoryTable;
    use crate::tables::processor::ProcessorTable;
    use crate::tables::roundup_npow2;
    use crate::univariate_polynomial::Polynomial;
    use crate::vm::VirtualMachine;

    #[test]
    fn test_padding() {
        let field = Field(18446744069414584321);
        let generator = field.generator();
        let order = 1 << 32;
        let vm = VirtualMachine::new(field);
        let code2 = ">>[++-]<".to_string();
        let program = vm.compile(code2);
        // let (runtime, _, _) = vm.run(&program, "".to_string());
        let (processor_matrix, _memory_matrix, instruction_matrix, _input_matrix, _output_matrix) =
            vm.simulate(&program, "".to_string());
        assert_eq!(
            instruction_matrix.len(),
            processor_matrix.len() + program.len()
        );
        let ilen = instruction_matrix.len();
        let mut instruction_table = super::InstructionTable::new(
            field,
            instruction_matrix.len() as u128,
            generator,
            order,
            &instruction_matrix,
        );
        instruction_table.pad();
        assert!(instruction_table.table.matrix.len().is_power_of_two());
        assert_eq!(
            instruction_matrix[ilen - 1][Indices::Address as usize],
            instruction_table.table.matrix.last().unwrap()[Indices::Address as usize]
        );
        assert_eq!(
            instruction_table.table.matrix.last().unwrap()[Indices::CurrentInstruction as usize],
            FieldElement::zero(field)
        );
        assert_eq!(
            instruction_table.table.matrix.last().unwrap()[Indices::NextInstruction as usize],
            FieldElement::zero(field)
        );
    }

    #[test]
    fn test_air() {
        let field = Field::new(18446744069414584321);
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        let two = one + one;
        let vm = VirtualMachine::new(field);
        let generator = field.generator().pow((1 << 32) - 1);
        let _omicron = generator.clone();
        let order = 1 << 32;
        // let code = "++>+++++[<+>-]++++++++[<++++++>-]<.".to_string();
        let code2 = "++>+-[+--]++.".to_string();
        let program = vm.compile(code2);
        println!("{:?}", program.clone());
        let (rt, _, _) = vm.run(&program, "".to_string());
        println!("{:?}", rt);
        // assert_eq!(program.len(), 2);
        let (processor_matrix, memory_matrix, instruction_matrix, input_matrix, output_matrix) =
            vm.simulate(&program, "".to_string());
        let challenges = vec![two; 11];
        let mut processor_table = ProcessorTable::new(
            field,
            processor_matrix.len() as u128,
            roundup_npow2(instruction_matrix.len() as u128),
            generator,
            order,
            &processor_matrix,
        );
        let mut memory_table = MemoryTable::new(
            field,
            memory_matrix.len() as u128,
            roundup_npow2(instruction_matrix.len() as u128),
            generator,
            order,
            &memory_matrix,
        );
        let mut instruction_table = InstructionTable::new(
            field,
            instruction_matrix.len() as u128,
            generator,
            order,
            &instruction_matrix,
        );
        let mut input_table = IOTable::new(
            field,
            input_matrix.len() as u128,
            generator,
            order,
            &input_matrix,
        );
        let mut output_table = IOTable::new(
            field,
            output_matrix.len() as u128,
            generator,
            order,
            &output_matrix,
        );

        processor_table.pad();
        memory_table.pad();
        instruction_table.pad();
        input_table.pad();
        output_table.pad();

        // println!("instruction table before extending columns");
        // for row in instruction_table.table.matrix.clone() {
        //     println!("{:?}", row);
        // }
        let terminal = instruction_table.extend_column(0, challenges.clone());
        let terminal2 = processor_table.extend_columns(challenges.clone());
        println!("instruction table after extending columns");
        for row in instruction_table.table.matrix.clone() {
            println!("{:?}", row);
        }
        // println!("processor table after extending columns");
        // for row in processor_table.table.matrix.clone() {
        //     println!("{:?}", row);
        // }
        println!("tppa: {:?}", terminal[0]);
        println!("tpea: {:?}", terminal[1]);
        println!("tipa: {:?}", terminal2[0]);

        let mut omicron_domain: Vec<FieldElement> = Vec::new();
        for i in 0..instruction_table.table.height {
            omicron_domain.push(instruction_table.table.omicron.pow(i));
            if i == 4 {
                println!("omicron_domain: {:?}", omicron_domain);
            }
        }
        let air = instruction_table.generate_air(challenges, terminal[0], terminal[1]);
        // println!("air");
        // for row in air.clone() {
        //     println!("{:?}", row);
        // }
        let zerofiers = instruction_table.generate_zerofier();
        let b = air[0].clone().evaluate(omicron_domain[0]);
        assert_eq!(b, zero);

        for v in 0..instruction_table.table.length - 1 {
            let transition = air[1].clone(); //.evaluate(omicron_domain[v as usize]);
            assert_eq!(air[1].clone().q_div(zerofiers[1].clone()).1, Polynomial::zero(field), "failed at {}", v);

            //println!("{:?}: {}: {}", transition, v, zerofiers[1].evaluate(omicron_domain[v as usize]));
            //assert_eq!(transition, zero);
        }

        let terminal =
            air[2].evaluate(omicron_domain[(instruction_table.table.length - 1) as usize]);
        assert_eq!(terminal, zero);

        // println!("Input table after extending columns");
        // for row in input_table.table.matrix.clone() {
        //     println!("{:?}", row);
        // }
        // println!("Output table after extending columns");
        // for row in output_table.table.matrix.clone() {
        //     println!("{:?}", row);
        // }
    }
}
