use super::processor::Indices as ProcessorIndices;
use super::Table;
use super::{derive_omicron, roundup_npow2};
use crate::fields::{Field, FieldElement};
use crate::univariate_polynomial::interpolate_lagrange_polynomials;
use crate::univariate_polynomial::Polynomial;
pub struct MemoryTable {
    pub table: Table,
}

pub enum Indices {
    // Named indices for base columns
    Cycle,
    MemoryPointer,
    MemoryValue,
    // Named indices for extension columns
    PermutationArg,
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

impl MemoryTable {
    pub fn new(
        field: Field,
        length: u128,
        height: u128,
        generator: FieldElement,
        order: u128,
        matrix: &[Vec<FieldElement>],
    ) -> Self {
        let base_width = 3;
        let full_width = base_width + 1;
        let omicron = derive_omicron(generator, order, height);

        let mut gmatrix =
            vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        for i in 0..matrix.len() {
            gmatrix[i][Indices::Cycle as usize] = matrix[i][Indices::Cycle as usize];
            gmatrix[i][Indices::MemoryPointer as usize] =
                matrix[i][Indices::MemoryPointer as usize];
            gmatrix[i][Indices::MemoryValue as usize] = matrix[i][Indices::MemoryValue as usize];
        }
        let table = Table::new(
            field, base_width, full_width, length, height, omicron, generator, order, gmatrix,
        );
        Self { table }
    }

    pub fn pad(&mut self) {
        let matrix = &mut self.table.matrix;
        let field = self.table.field;
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        let last_row = matrix[self.table.length as usize - 1].clone();
        for i in 0..(self.table.height - self.table.length) as usize {
            let mut curr_row = vec![zero; 4]; // 3 is number of columns.
            curr_row[Indices::Cycle as usize] =
                last_row[Indices::Cycle as usize] + FieldElement((i + 1) as u128, field) * one;
            curr_row[Indices::MemoryPointer as usize] = last_row[Indices::MemoryPointer as usize];
            curr_row[Indices::MemoryValue as usize] = last_row[Indices::MemoryValue as usize];
            matrix[self.table.length as usize + i] = curr_row;
        }
    }

    // processor matrix that goes here is unpadded
    pub fn derive_matrix(processor_matrix: &[Vec<FieldElement>]) -> Vec<Vec<FieldElement>> {
        let mut matrix: Vec<Vec<FieldElement>> = processor_matrix
            .iter()
            // .filter(|pt| pt[ProcessorIndices::CurrentInstruction as usize] != zero), we need the last processor row, right?
            .map(|pt| {
                vec![
                    pt[ProcessorIndices::Cycle as usize],
                    pt[ProcessorIndices::MemoryPointer as usize],
                    pt[ProcessorIndices::MemoryValue as usize],
                ]
            })
            .collect();

        // Sort by memory_pointer
        matrix.sort_by(|a, b| {
            a[Indices::MemoryPointer as usize]
                .partial_cmp(&b[Indices::MemoryPointer as usize])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        matrix
    }
    #[warn(non_snake_case)]

    //the matrix taken here is padded
    pub fn extend_column_ppa(
        &mut self,
        rand_field_elem: u128,
        challenges: Vec<FieldElement>,
    ) -> Vec<FieldElement> {
        let mut ppa = FieldElement::new(rand_field_elem, self.table.field);
        self.table.matrix[0][Indices::PermutationArg as usize] = ppa;
        for i in 0..self.table.length - 1 {
            let weighted_sum = self.table.matrix[i as usize][Indices::Cycle as usize]
                * challenges[ChallengeIndices::D as usize]
                + self.table.matrix[i as usize][Indices::MemoryPointer as usize]
                    * challenges[ChallengeIndices::E as usize]
                + self.table.matrix[i as usize][Indices::MemoryValue as usize]
                    * challenges[ChallengeIndices::F as usize]
                - challenges[ChallengeIndices::Beta as usize];
            ppa *= weighted_sum;
            self.table.matrix[(i + 1) as usize][Indices::PermutationArg as usize] = ppa;
        }
        let tppa = ppa
            * (self.table.matrix[(self.table.length - 1) as usize][Indices::Cycle as usize]
                * challenges[ChallengeIndices::D as usize]
                + self.table.matrix[(self.table.length - 1) as usize]
                    [Indices::MemoryPointer as usize]
                    * challenges[ChallengeIndices::E as usize]
                + self.table.matrix[(self.table.length - 1) as usize]
                    [Indices::MemoryValue as usize]
                    * challenges[ChallengeIndices::F as usize]
                - challenges[ChallengeIndices::Beta as usize]);
        let mut terminal: Vec<FieldElement> = Vec::with_capacity(1);
        terminal.push(tppa);
        terminal
    }

    //this is after padding and extension
    pub fn generate_air(
        &self,
        challenges: Vec<FieldElement>,
        tppa: FieldElement,
    ) -> Vec<Polynomial> {
        let interpolated = self.table.clone().interpolate_columns(vec![
            Indices::Cycle as u128,
            Indices::MemoryPointer as u128,
            Indices::MemoryValue as u128,
            Indices::PermutationArg as u128,
        ]);
        let clk = interpolated[Indices::Cycle as usize].clone();
        let mp = interpolated[Indices::MemoryPointer as usize].clone();
        let mv = interpolated[Indices::MemoryValue as usize].clone();
        let ppa = interpolated[Indices::PermutationArg as usize].clone();

        let next_interpolated = self.table.clone().next_interpolate_columns(interpolated);

        let clk_next = next_interpolated[Indices::Cycle as usize].clone();
        let mp_next = next_interpolated[Indices::MemoryPointer as usize].clone();
        let mv_next = next_interpolated[Indices::MemoryValue as usize].clone();
        let ppa_next = next_interpolated[Indices::PermutationArg as usize].clone();
        let one = Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        let mut air = vec![];

        //Boundary constraint: clk = mp = mv = 0, ppa = 1 (random secret using diff constr if we use it, not needed rn)
        let boundaryair = clk.clone() + mp.clone() + mv.clone() + ppa.clone()
            - Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        air.push(boundaryair);

        //Transition constraints: * == next
        //1. (mp+1-mp*).(mp-mp*)
        //2. (mp-mp*).mv*
        //4. (clk + 1 - clk*).(mv*-mv).(mp+1-mp*)
        //5. ppa.(d.clk + e.mp + f.mv - beta) - ppa*
        let mp_mp_next = mp.clone() - mp_next.clone();
        let mp_one_mp_next = mp_mp_next.clone() + one.clone();
        let transitionair = (mp_one_mp_next.clone()) * (mp_mp_next.clone())
            + (mp_mp_next.clone()) * mv_next.clone()
            + (clk.clone() + one.clone() - clk_next)
                * (mv_next.clone() - mv.clone())
                * (mp_one_mp_next)
            + ppa.clone()
                * (clk.scalar_mul(challenges[ChallengeIndices::D as usize])
                    + mp.scalar_mul(challenges[ChallengeIndices::E as usize])
                    + mv.scalar_mul(challenges[ChallengeIndices::F as usize])
                    - Polynomial::new_from_coefficients(vec![
                        challenges[ChallengeIndices::Beta as usize],
                    ]))
            - ppa_next;
        air.push(transitionair);

        //Terminal constraints:
        //ppa.(d.clk+e.mp+f.mv-beta)-Tppa

        let terminalair = ppa
            * (clk.scalar_mul(challenges[ChallengeIndices::D as usize])
                + mp.scalar_mul(challenges[ChallengeIndices::E as usize])
                + mv.scalar_mul(challenges[ChallengeIndices::F as usize])
                - Polynomial::new_from_coefficients(vec![
                    challenges[ChallengeIndices::Beta as usize],
                ]))
            - Polynomial::constant(tppa);

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
    ) -> Vec<Polynomial> {
        let mut quotients = vec![];
        let air = self.generate_air(challenges, tppa);
        let zerofiers = self.generate_zerofier();

        for i in 0..air.len() {
            assert_eq!(
                air[i].clone().q_div(zerofiers[i].clone()).1,
                Polynomial::constant(FieldElement::zero(self.table.field))
            );
            quotients.push(air[i].clone().q_div(zerofiers[i].clone()).0);
        }
        quotients
    }
}

#[cfg(test)]
mod test_memory_table {
    #![allow(unused_variables)]
    use super::MemoryTable;
    use crate::fields::{Field, FieldElement};
    use crate::tables::instruction::InstructionTable;
    use crate::tables::io::IOTable;
    use crate::tables::memory::{ChallengeIndices, Indices};
    use crate::tables::processor::ProcessorTable;
    use crate::tables::roundup_npow2;
    use crate::tables::Table;
    use crate::vm::VirtualMachine;

    #[test]
    fn test_padding() {
        let field = Field(18446744069414584321);
        let generator = field.generator();
        let order = 1 << 32;
        let vm = VirtualMachine::new(field);
        let code = "++>,<[>+.<-]".to_string();
        let program = vm.compile(code);
        let (processor_matrix, memory_matrix, instruction_matrix, input_matrix, output_matrix) =
            vm.simulate(&program, "a".to_string());
        let mlen = memory_matrix.len();
        let mut memory_table = MemoryTable::new(
            field,
            memory_matrix.len() as u128,
            roundup_npow2(instruction_matrix.len() as u128),
            generator,
            order,
            &memory_matrix,
        );
        memory_table.pad();
        assert!(memory_table.table.matrix.len().is_power_of_two());
        assert_eq!(
            memory_table.table.matrix.last().unwrap()[Indices::MemoryPointer as usize],
            memory_matrix.last().unwrap()[Indices::MemoryPointer as usize]
        );
    }
    #[test]
    pub fn test_memory_table_permutation_argument() {
        let field = Field(18446744069414584321);
        let vm = VirtualMachine::new(field);
        let code = "++>,<[>+.<-]".to_string();
        let program = vm.compile(code);
        let input = "a".to_string();
        // vm.run(&program, input.c);
        let (processor_matrix, memory_matrix, instruction_matrix, input_matrix, output_matrix) =
            vm.simulate(&program, input);
        let generator = field.generator();
        let order = 1 << 32;
        let zero = FieldElement::zero(field);
        let mut mem = MemoryTable::new(
            field,
            memory_matrix.len() as u128,
            roundup_npow2(instruction_matrix.len() as u128),
            generator,
            order,
            &&memory_matrix,
        );
        let mut challenges = vec![zero; 11];
        challenges[ChallengeIndices::Beta as usize] = FieldElement::new(3, field);
        challenges[ChallengeIndices::D as usize] = FieldElement::new(1, field);
        challenges[ChallengeIndices::E as usize] = FieldElement::new(1, field);
        challenges[ChallengeIndices::F as usize] = FieldElement::new(1, field);
        mem.extend_column_ppa(7, challenges);
    }

    #[test]
    fn test_extend_column_ppa() {
        // Mock field and FieldElement implementation (use actual implementation if available)
        let field = Field::new(101);
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);

        // Create a mock table structure
        let mut matrix = vec![
            vec![
                FieldElement::new(1, field),
                FieldElement::new(2, field),
                FieldElement::new(3, field),
                zero.clone(),
            ],
            vec![
                FieldElement::new(4, field),
                FieldElement::new(5, field),
                FieldElement::new(6, field),
                zero.clone(),
            ],
            vec![
                FieldElement::new(7, field),
                FieldElement::new(8, field),
                FieldElement::new(9, field),
                zero.clone(),
            ],
            vec![
                FieldElement::new(1, field),
                FieldElement::new(2, field),
                FieldElement::new(3, field),
                zero.clone(),
            ],
            vec![
                FieldElement::new(10, field),
                FieldElement::new(12, field),
                FieldElement::new(13, field),
                zero.clone(),
            ],
            vec![
                FieldElement::new(1, field),
                FieldElement::new(2, field),
                FieldElement::new(3, field),
                zero.clone(),
            ],
        ];
        // Create the challenges
        let challenges = vec![
            FieldElement::new(3, field), // D
            FieldElement::new(4, field), // E
            FieldElement::new(5, field), // F
            FieldElement::new(2, field), // Beta
        ];
        let mut ppa = FieldElement::new(1, field);
        println!("{}:0", ppa.0);
        matrix[0].push(FieldElement::one(field));
        for i in 0..matrix.len() - 1 {
            let weighted_sum = matrix[i as usize][0] * challenges[0]
                + matrix[i as usize][1] * challenges[1]
                + matrix[i as usize][2] * challenges[2]
                - challenges[3];
            matrix[(i + 1) as usize].push(ppa * weighted_sum);
            ppa *= weighted_sum;
            println!(
                "{} + {} + {} - {}",
                (matrix[i as usize][0] * challenges[0]).0,
                (matrix[i as usize][1] * challenges[1]).0,
                (matrix[i as usize][2] * challenges[2]).0,
                challenges[3].0
            );
            println!("{}:{}", ppa.0, i + 1);
        }
    }
    #[test]
    fn test_generate_air() {
        let field = Field::new((1 << 64) - (1 << 32) + 1);
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        let two = one + one;
        let vm = VirtualMachine::new(field);
        let generator = FieldElement::new(1753635133440165772, field);
        // let omicron = generator.clone();
        let order = 1 << 32;
        let code2 = "++>+++++[<+>-]++++++++[<++++++>-]<.".to_string();
        //let code2 = "++>+-[+--]++.".to_string();
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
        processor_table.table.generate_omicron_domain();
        memory_table.table.generate_omicron_domain();
        instruction_table.table.generate_omicron_domain();
        input_table.table.generate_omicron_domain();
        output_table.table.generate_omicron_domain();
        let terminal = memory_table.extend_column_ppa(1, challenges.clone());
        let terminal2 = processor_table.extend_columns(challenges.clone());
        let mut i = 0;
        println!("memory_table after extending columns");
        for row in memory_table.table.matrix.clone() {
            println!("{} : {:?}", i, row);
            i += 1;
        }
        println!("tppa: {:?}", terminal[0]);
        println!("tmpa: {:?}", terminal2[1]);
        let mut omicron_domain: Vec<FieldElement> = Vec::new();
        for i in 0..memory_table.table.height {
            omicron_domain.push(memory_table.table.omicron.pow(i));
        }

        let air = memory_table.generate_air(challenges, terminal[0]);

        let b = air[0].evaluate(omicron_domain[0]);

        assert_eq!(b, zero);
        for v in 0..rt - 1 {
            let t_all = air[1].evaluate(omicron_domain[v as usize]);
            //println!("{}", t_all);
            assert_eq!(t_all, zero, "failed at {}", v);
        }

        let v = omicron_domain[4];

        let t_air1 = air[2].evaluate(omicron_domain[rt as usize - 1]);
        //println!("{}",t_air1);
        assert_eq!(t_air1, zero);
    }

    #[test]
    fn test_not_generate_air() {
        let field = Field::new((1 << 64) - (1 << 32) + 1);
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        let two = one + one;
        let vm = VirtualMachine::new(field);
        let generator = FieldElement::new(1753635133440165772, field);
        // let omicron = generator.clone();
        let order = 1 << 32;
        // let code = "++>+++++[<+>-]++++++++[<++++++>-]<.".to_string();
        let code2 = ">>[++-]<".to_string();
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
        processor_table.table.generate_omicron_domain();
        memory_table.table.generate_omicron_domain();
        instruction_table.table.generate_omicron_domain();
        input_table.table.generate_omicron_domain();
        output_table.table.generate_omicron_domain();
        let terminal = memory_table.extend_column_ppa(1, challenges.clone());
        println!("memory_table after extending columns");
        for row in memory_table.table.matrix.clone() {
            println!("{:?}", row);
        }
        println!("tppa: {:?}", terminal[0]);
        let mut omicron_domain: Vec<FieldElement> = Vec::new();
        for i in 0..memory_table.table.height {
            omicron_domain.push(memory_table.table.omicron.pow(i));
            if i == 4 {
                println!("omicron_domain: {:?}", omicron_domain);
            }
        }

        let air = memory_table.generate_air(challenges, terminal[0]);

        let b = air[0].evaluate(omicron_domain[1]);

        assert_ne!(b, zero);
        for v in 0..rt - 1 {
            let t_all = air[1].evaluate(omicron_domain[v as usize]);
            // println!("{:?}", t_all);
            assert_eq!(t_all, zero);
        }

        let v = omicron_domain[3];

        let t_air1 = air[2].evaluate(v);
        assert_ne!(t_air1, zero);
    }
}
