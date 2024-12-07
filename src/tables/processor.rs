use super::Table;
use super::{derive_omicron, roundup_npow2};
use crate::fields::{Field, FieldElement};
use crate::fri::*;
use crate::univariate_polynomial::*;

pub struct ProcessorTable {
    table: Table,
}

pub enum Indices {
    // Named indices for registers
    Cycle,
    InstructionPointer,
    CurrentInstruction,
    NextInstruction,
    MemoryPointer,
    MemoryValue,
    MemoryValueInverse,
    // Named indices for extension columns
    InstructionPermutaion,
    MemoryPermuation,
    InputEvaluation,
    OutputEvaluation,
}
#[allow(non_camel_case_types)]
pub enum IndicesPoly {
    //i0 = [,
    //i1 = ],
    // i2 = <,
    // i3 = >,
    // i4 = +,
    // i5 = -,
    // i6 = ",",
    // i7 = ".",
    Boundary,
    Transition_i0,
    Transition_i1,
    Transition_i2,
    Transition_i3,
    Transition_i4,
    Transition_i5,
    Transition_i6,
    Transition_i7,
    Transition_all,
    //clk*-clk-1
    //inv.(1-inv.mv)
    //ci.(ipa.(a.ip+b.ci+c.ni-alpha)-ipa*)+(ipa*-ipa).deselector
    //mpa.(d.clk+e.mp+f.mv-beta)-mpa*
    //selector(,)(ci) . (iea.gamma + mv - iea*) + (ci -  “,”) . (iea - iea*)
    //selector(.)(ci) . (oea.delta + mv - oea*) + (ci -  “.”) . (oea - oea*)
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

impl ProcessorTable {
    pub fn new(
        field: Field,
        length: u128,
        generator: FieldElement,
        order: u128,
        matrix: Vec<Vec<FieldElement>>,
    ) -> Self {
        let base_width = 7;
        let full_width = base_width + 4;
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        let mut gmatrix =
            vec![vec![FieldElement::zero(field); full_width as usize]; height as usize];
        for i in 0..matrix.len() {
            gmatrix[i][Indices::Cycle as usize] = matrix[i][Indices::Cycle as usize];
            gmatrix[i][Indices::InstructionPointer as usize] =
                matrix[i][Indices::InstructionPointer as usize];
            gmatrix[i][Indices::CurrentInstruction as usize] =
                matrix[i][Indices::CurrentInstruction as usize];
            gmatrix[i][Indices::NextInstruction as usize] =
                matrix[i][Indices::NextInstruction as usize];
            gmatrix[i][Indices::MemoryPointer as usize] =
                matrix[i][Indices::MemoryPointer as usize];
            gmatrix[i][Indices::MemoryValue as usize] = matrix[i][Indices::MemoryValue as usize];
            gmatrix[i][Indices::MemoryValueInverse as usize] =
                matrix[i][Indices::MemoryValueInverse as usize];
        }
        let table = Table::new(
            field, base_width, full_width, length, height, omicron, generator, order, gmatrix,
        );
        Self { table }
    }

    // Note: Before padding initiate the matrix in table.
    // Add padding rows to convert the matrix derived from trace to a matrix of length of a power of 2
    pub fn pad(&mut self) {
        let matrix = &mut self.table.matrix;
        let field = self.table.field;
        let zero = FieldElement::zero(field);
        let last_row = matrix[self.table.length as usize - 1].clone();
        for i in 0..(self.table.height - self.table.length) as usize {
            // let last_row = &mut matrix.last().unwrap();
            let curr_row = &mut matrix[self.table.length as usize + i];
            // let mut new_row = vec![FieldElement::zero(self.table.field); self.table.full_width as usize];
            // Increment cycle by one for every new padded row.
            curr_row[Indices::Cycle as usize] = last_row[Indices::Cycle as usize]
                + FieldElement::new((i + 1) as u128, field) * FieldElement::one(self.table.field);
            // keep the instruction pointer same as the last row in every padded row.
            curr_row[Indices::InstructionPointer as usize] =
                last_row[Indices::InstructionPointer as usize];
            // keep the current instruction and next instruction as zero.
            curr_row[Indices::CurrentInstruction as usize] = zero;
            curr_row[Indices::NextInstruction as usize] = zero;
            // keep the memory pointer and memory value as last.
            curr_row[Indices::MemoryPointer as usize] = last_row[Indices::MemoryPointer as usize];
            curr_row[Indices::MemoryValue as usize] = last_row[Indices::MemoryValue as usize];
            curr_row[Indices::MemoryValueInverse as usize] =
                last_row[Indices::MemoryValueInverse as usize];
        }
    }

    //the matrix taken here is padded
    pub fn extend_columns(&mut self, challenges: Vec<FieldElement>) {
        //@todo Note: Taking init 1 for now, change to random secret initial value which we check by difference constraint of tmpa = Tppa
        let mut ipa = FieldElement::one(self.table.field);
        let mut mpa = FieldElement::one(self.table.field);
        let mut iea = FieldElement::zero(self.table.field);
        let mut oea = FieldElement::zero(self.table.field);

        self.table.matrix[0][Indices::InstructionPermutaion as usize] = ipa;
        self.table.matrix[0][Indices::MemoryPermuation as usize] = mpa;
        self.table.matrix[0][Indices::InputEvaluation as usize] = iea;
        self.table.matrix[0][Indices::OutputEvaluation as usize] = oea;
        // instruction permutation argument
        for i in 0..self.table.length - 1 {
            //ip * a + ci * b + ni * c
            let weighted_sum = self.table.matrix[i as usize][Indices::InstructionPointer as usize]
                * challenges[ChallengeIndices::A as usize]
                + self.table.matrix[i as usize][Indices::CurrentInstruction as usize]
                    * challenges[ChallengeIndices::B as usize]
                + self.table.matrix[i as usize][Indices::NextInstruction as usize]
                    * challenges[ChallengeIndices::C as usize]
                - challenges[ChallengeIndices::Alpha as usize];
            ipa *= weighted_sum;
            self.table.matrix[(i + 1) as usize][Indices::InstructionPermutaion as usize] = ipa;
        }

        for i in 0..self.table.length - 1 {
            let weighted_sum = self.table.matrix[i as usize][Indices::Cycle as usize]
                * challenges[ChallengeIndices::D as usize]
                + self.table.matrix[i as usize][Indices::MemoryPointer as usize]
                    * challenges[ChallengeIndices::E as usize]
                + self.table.matrix[i as usize][Indices::MemoryValue as usize]
                    * challenges[ChallengeIndices::F as usize]
                - challenges[ChallengeIndices::Beta as usize];
            mpa *= weighted_sum;
            self.table.matrix[(i + 1) as usize][Indices::MemoryPermuation as usize] = mpa;
        }

        let f =
            |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, self.table.field) };
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f(',') == ci {
                iea = iea * challenges[ChallengeIndices::Gamma as usize]
                    + self.table.matrix[i as usize][Indices::MemoryValue as usize];
                self.table.matrix[(i + 1) as usize][Indices::InputEvaluation as usize] = iea;
            } else {
                self.table.matrix[(i + 1) as usize][Indices::InputEvaluation as usize] = iea;
            }
        }

        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('.') == ci {
                oea = oea * challenges[ChallengeIndices::Delta as usize]
                    + self.table.matrix[i as usize][Indices::MemoryValue as usize];
                self.table.matrix[(i + 1) as usize][Indices::OutputEvaluation as usize] = oea;
            } else {
                self.table.matrix[(i + 1) as usize][Indices::OutputEvaluation as usize] = oea;
            }
        }
    }

    pub fn generate_zerofier(&self) -> Vec<Polynomial> {
        let mut zerofiers = vec![];
        let omicron = self.table.omicron;
        let x = Polynomial::new_from_coefficients(vec![
            FieldElement::zero(self.table.field),
            FieldElement::one(self.table.field),
        ]);
        let f =
            |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, self.table.field) };

        //boundary
        let boundary_zerofier =
            x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(0)]);
        zerofiers.push(boundary_zerofier);

        //i0
        let mut transition_i0_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('[') == ci {
                transition_i0_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i0_zerofier);

        //i1
        let mut transition_i1_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f(']') == ci {
                transition_i1_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i1_zerofier);

        //i2
        let mut transition_i2_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('<') == ci {
                transition_i2_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i2_zerofier);

        //i3
        let mut transition_i3_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('>') == ci {
                transition_i3_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i3_zerofier);

        //i4
        let mut transition_i4_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('+') == ci {
                transition_i4_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i4_zerofier);

        //i5
        let mut transition_i5_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('-') == ci {
                transition_i5_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i5_zerofier);

        //i6
        let mut transition_i6_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f(',') == ci {
                transition_i6_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i6_zerofier);

        //i7
        let mut transition_i7_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            let ci = self.table.matrix[i as usize][Indices::CurrentInstruction as usize];
            if f('.') == ci {
                transition_i7_zerofier *=
                    x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
            }
        }
        zerofiers.push(transition_i7_zerofier);

        //all
        let mut transition_all_zerofier = Polynomial::new_from_coefficients(vec![]);
        for i in 0..self.table.length - 1 {
            transition_all_zerofier *=
                x.clone() - Polynomial::new_from_coefficients(vec![omicron.clone().pow(i)]);
        }
        zerofiers.push(transition_all_zerofier);

        //terminal
        let terminal_zerofier = x.clone()
            - Polynomial::new_from_coefficients(vec![omicron.clone().pow(self.table.length - 1)]);

        zerofiers.push(terminal_zerofier);
        zerofiers
    }

    pub fn generate_quotients(&self, challenges: Vec<FieldElement>) -> Vec<Polynomial> {
        let mut quotients = vec![];
        let air = self.generate_air(challenges);
        let zerofiers = self.generate_zerofier();

        for i in 0..air.len() {
            quotients.push(air[i].clone().q_div(zerofiers[i].clone()).0);
        }
        quotients
    }

    //define a selector polynomial for a specific instruction.
    //this will return a non-zero value for instruction and zero for all other instructions
    pub fn selector_polynomial(instruction: char, ci: Polynomial, field: Field) -> Polynomial {
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, field) };
        let mut acc = Polynomial::constant(FieldElement::new(1, field)); // poly=1
        for c in "[]<>,.+-".chars() {
            if c != instruction {
                acc *= ci.clone() - Polynomial::constant(FieldElement::new(f(c).0, field));
            }
        }
        acc
    }

    // boundary constraints for the base coloumns
    // the values of instructionpermutaion ipa and mpa I am taking as 1
    pub fn generate_air(&self, challenges: Vec<FieldElement>) -> Vec<Polynomial> {
        let f =
            |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, self.table.field) };
        let interpolated = self.table.clone().interpolate_columns(vec![
            Indices::Cycle as u128,
            Indices::InstructionPointer as u128,
            Indices::CurrentInstruction as u128,
            Indices::NextInstruction as u128,
            Indices::MemoryPointer as u128,
            Indices::MemoryValue as u128,
            Indices::MemoryValueInverse as u128,
            Indices::InstructionPermutaion as u128,
            Indices::MemoryPermuation as u128,
            Indices::InputEvaluation as u128,
            Indices::OutputEvaluation as u128,
        ]);
        let clk = interpolated[Indices::Cycle as usize].clone();
        let ip = interpolated[Indices::InstructionPointer as usize].clone();
        let ci = interpolated[Indices::CurrentInstruction as usize].clone();
        let ni = interpolated[Indices::NextInstruction as usize].clone();
        let mp = interpolated[Indices::MemoryPointer as usize].clone();
        let mv = interpolated[Indices::MemoryValue as usize].clone();
        let inv_mv = interpolated[Indices::MemoryValueInverse as usize].clone();
        let ipa = interpolated[Indices::InstructionPermutaion as usize].clone();
        let mpa = interpolated[Indices::MemoryPermuation as usize].clone();
        let iea = interpolated[Indices::InputEvaluation as usize].clone();
        let oea = interpolated[Indices::OutputEvaluation as usize].clone();
        let next_interpolated = self.table.clone().next_interpolate_columns(vec![
            Indices::Cycle as u128,
            Indices::InstructionPointer as u128,
            Indices::CurrentInstruction as u128,
            Indices::NextInstruction as u128,
            Indices::MemoryPointer as u128,
            Indices::MemoryValue as u128,
            Indices::MemoryValueInverse as u128,
            Indices::InstructionPermutaion as u128,
            Indices::MemoryPermuation as u128,
            Indices::InputEvaluation as u128,
            Indices::OutputEvaluation as u128,
        ]);
        let clk_next = next_interpolated[Indices::Cycle as usize].clone();
        let ip_next = next_interpolated[Indices::InstructionPointer as usize].clone();
        let ci_next = next_interpolated[Indices::CurrentInstruction as usize].clone();
        let ni_next = next_interpolated[Indices::NextInstruction as usize].clone();
        let mp_next = next_interpolated[Indices::MemoryPointer as usize].clone();
        let mv_next = next_interpolated[Indices::MemoryValue as usize].clone();
        let inv_mv_next = next_interpolated[Indices::MemoryValueInverse as usize].clone();
        let ipa_next = next_interpolated[Indices::InstructionPermutaion as usize].clone();
        let mpa_next = next_interpolated[Indices::MemoryPermuation as usize].clone();
        let iea_next = next_interpolated[Indices::InputEvaluation as usize].clone();
        let oea_next = next_interpolated[Indices::OutputEvaluation as usize].clone();
        let one = Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        let mut air = vec![];
        //Boundary contsraints :clk=mp=mv=inv=ip=0
        //iea=oea=1 (we are using 1 instead of any random number)
        let poly_two =
            Polynomial::new_from_coefficients(vec![FieldElement::new(2, self.table.field)]);
        //@todo should we separate out ipa and mpa -2 into ipa -1 and mpa -1 so that malicious prover cant set ppa=2 and mpa =0 or something like that?
        let boundary_air = clk.clone()
            + ip.clone()
            + mp.clone()
            + mv.clone()
            + inv_mv.clone()
            + ipa.clone()
            + mpa.clone()
            + iea.clone()
            + oea.clone()
            - poly_two.clone();
        air.push(boundary_air);
        // Transition_i0,
        // Transition_i1,
        // Transition_i2,
        // Transition_i3,
        // Transition_i4,
        // Transition_i5,
        // Transition_i6,
        // Transition_i7,
        // Transition_all,
        let poly_one = Polynomial::new_from_coefficients(vec![FieldElement::one(self.table.field)]);
        let mv_is_zero = poly_one.clone() - mv.clone() * inv_mv.clone();
        //(ip⋆−ip−2)⋅mv+(ip⋆−ni)⋅iszero
        // mp⋆−mp
        // mv⋆−mv
        let trasition_i0 = (mv.clone() * (ip_next.clone() - ip.clone() - poly_two.clone())
            + mv_is_zero.clone() * (ip_next.clone() - ni.clone()))
            + (mp_next.clone() - mp.clone())
            + (mv_next.clone() - mv.clone());
        air.push(trasition_i0);
        // (ip⋆−ip−2)⋅iszero+(ip⋆−ni)⋅mv
        // mp⋆−mp
        // mv⋆−mv

        let trasition_i1 = mv_is_zero.clone() * (ip_next.clone() - ip.clone() - poly_two.clone())
            + (ip_next.clone() - ni.clone()) * mv.clone()
            + (mv_next.clone() - mv.clone())
            + (mp_next.clone() - mp.clone());
        air.push(trasition_i1);

        // ip⋆−ip−1
        // mp⋆−mp+1
        let trasition_i2 = (ip_next.clone() - ip.clone() - poly_one.clone())
            + (mp_next.clone() - mp.clone() + poly_one.clone());
        air.push(trasition_i2);
        // ip⋆−ip−1
        // mp⋆−mp+1
        let trasition_i3 = (ip_next.clone() - ip.clone() - poly_one.clone())
            + (mp_next.clone() - mp.clone() + poly_one.clone());
        air.push(trasition_i3);

        // ip⋆−ip−1
        // mp⋆−mp
        // mv⋆−mv−1
        let trasition_i4 = (ip_next.clone() - ip.clone() - poly_one.clone())
            + (mp_next.clone() - mp.clone())
            + (mv_next.clone() - mv.clone() - poly_one.clone());
        air.push(trasition_i4);
        // ip⋆−ip−1
        // mp⋆−mp
        // mv⋆−mv+1
        let trasition_i5 = (ip_next.clone() - ip.clone() - poly_one.clone())
            + (mp_next.clone() - mp.clone())
            + (mv_next.clone() - mv.clone() - poly_one.clone());
        air.push(trasition_i5);
        //  ip⋆−ip−1
        // mp⋆−mp

        let trasition_i6 =
            (ip_next.clone() - ip.clone() - poly_one.clone()) + (mp_next.clone() - mp.clone());
        air.push(trasition_i6);
        //  ip⋆−ip−1
        // mp⋆−mp
        // mv⋆−mv
        let trasition_i7 = (ip_next.clone() - ip.clone() - poly_one.clone())
            + (mp_next.clone() - mp.clone())
            + (mv_next.clone() - mv.clone());
        air.push(trasition_i7);
        //clk⋆−clk−1
        // inv⋅(1−inv⋅mv)
        //ci.(ipa.(a.ip+b.ci+c.ni-alpha)-ipa*) // this constrainst is become redundant becaue we are not using +(ipa*-ipa).deselector
        //mpa.(d.clk+e.mp+f.mv-beta)-mpa*
        //selector(,)(ci) . (iea.gamma + mv - iea*) + (ci -  “,”) . (iea - iea*)
        //selector(.)(ci) . (oea.delta + mv - oea*) + (ci -  “.”) . (oea - oea*)

        // /clk⋆−clk−1+inv⋅(1−inv⋅mv)
        //ci.(ipa.(a.ip+b.ci+c.ni-alpha)-ipa*)
        //mpa.(d.clk+e.mp+f.mv-beta)-mpa*
        // we are changing mv to mv*
        //selector(,)(ci) . (iea.gamma + mv* - iea*) + (ci -  “,”) . (iea - iea*)
        //selector(.)(ci) . (oea.delta + mv - oea*) + (ci -  “.”) . (oea - oea*)
        let trasition_all = (clk_next.clone() - clk.clone() - poly_one.clone())
            + (inv_mv.clone() * (mv_is_zero.clone()))
            + ci.clone()
                * (ipa.clone()
                    * (ip.scalar_mul(challenges[ChallengeIndices::A as usize])
                        + ci.scalar_mul(challenges[ChallengeIndices::B as usize])
                        + ni.scalar_mul(challenges[ChallengeIndices::C as usize])
                        - Polynomial::constant(challenges[ChallengeIndices::Alpha as usize]))
                    - ipa_next)
            + (mpa.clone()
                * (clk.scalar_mul(challenges[ChallengeIndices::D as usize])
                    + mp.scalar_mul(challenges[ChallengeIndices::E as usize])
                    + mv.scalar_mul(challenges[ChallengeIndices::F as usize])
                    - Polynomial::constant(challenges[ChallengeIndices::Beta as usize]))
                - mpa_next)
            + ProcessorTable::selector_polynomial(',', ci.clone(), self.table.field)
                * ci.clone()
                * (iea
                    .clone()
                    .scalar_mul(challenges[ChallengeIndices::Gamma as usize])
                    + mv_next.clone()
                    - iea_next.clone())
            + (ci.clone() - Polynomial::constant(f(','))) * (iea.clone() - iea_next.clone())
            + ProcessorTable::selector_polynomial('.', ci.clone(), self.table.field)
                * ci.clone()
                * (oea
                    .clone()
                    .scalar_mul(challenges[ChallengeIndices::Delta as usize])
                    + mv.clone()
                    - oea_next.clone())
            + (ci.clone() - Polynomial::constant(f('.'))) * (oea.clone() - oea_next.clone());
        air.push(trasition_all);
        //  Terminal constraints
        // tipa, tmpa- last row not accumulated so:
        //1.ipa.(a.ip+ b.ci+c.ni-alpha)-tipa
        //2.mpa.(d.clk+e.mp+f.mv-beta)-tmpa
        // tiea, toea- last element identical to terminal
        //3.iea-tiea   4. oea-toea
        let tipa = Polynomial::new_from_coefficients(vec![]);
        let tmpa = Polynomial::new_from_coefficients(vec![]);
        let tiea = Polynomial::new_from_coefficients(vec![]);
        let toea = Polynomial::new_from_coefficients(vec![]);
        let terminal_air = ipa.clone()
            * (ip
                .clone()
                .scalar_mul(challenges[ChallengeIndices::A as usize])
                + ci.clone()
                    .scalar_mul(challenges[ChallengeIndices::B as usize])
                + ni.clone()
                    .scalar_mul(challenges[ChallengeIndices::C as usize])
                - Polynomial::constant(challenges[ChallengeIndices::Beta as usize]))
            - tipa.clone()
            + mpa.clone()
                * (clk.scalar_mul(challenges[ChallengeIndices::D as usize])
                    + mp.clone()
                        .scalar_mul(challenges[ChallengeIndices::E as usize])
                    + mv.clone()
                        .scalar_mul(challenges[ChallengeIndices::F as usize])
                    - Polynomial::constant(challenges[ChallengeIndices::Beta as usize]))
            - tmpa
            + iea
            - tiea
            + oea
            - toea;
        air.push(terminal_air);
        air
    }
}
#[cfg(test)]
mod tests_processor_operations {
    use super::*;
    #[test]
    fn test_selector_poly() {
        let field = Field::new((1 << 64) - (1 << 32) + 1);
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, field) };
        let ci_array = vec!['[', ']', '<', '>', '-', '.', ']', ','];
        let mut values_array: Vec<FieldElement> = Vec::new();
        for c in ci_array.clone() {
            println!("{:?}", f(c));
            values_array.push(f(c));
        }
        let generator = FieldElement::new(1753635133440165772, field);
        let order = 1 << 32;
        let target_order = 8;
        //println!("generator ={:?}", generator);
        let omicron = derive_omicron(generator, order, target_order);
        //println!("omicron ={:?}", omicron);
        let domain = FriDomain::new(FieldElement::new(1, field), omicron, target_order);
        let ci = domain.interpolate(values_array);
        let poly = ProcessorTable::selector_polynomial('.', ci, field);
        //println!("{:?}", f(']'));
        for i in 0..target_order {
            println!(
                "Selector([) value at ci:{}, ci: {}",
                poly.evaluate(omicron.pow(i)),
                ci_array[i as usize]
            );
        }
    }
}

// @todo test processor table padding
#[cfg(test)]
mod test_processor {
    use std::time::Instant;

    use super::ProcessorTable;
    use crate::fields::{Field, FieldElement};
    use crate::tables::instruction::InstructionTable;
    use crate::tables::io::IOTable;
    use crate::tables::memory::Memory;
    use crate::tables::processor::Indices;
    use crate::vm::VirtualMachine;

    #[test]
    fn test_padding() {
        let field = Field(18446744069414584321);
        let generator = field.generator();
        let order = 1 << 32;
        let vm = VirtualMachine::new(field);
        let code2 = ">>[++-]<".to_string();
        let program = vm.compile(code2);
        let (runtime, _, _) = vm.run(&program, "".to_string());
        let (processor_matrix, _memory_matrix, _instruction_matrix, _input_matrix, _output_matrix) =
            vm.simulate(&program, "".to_string());
        assert!(runtime == processor_matrix.len() as i32);
        let mut processor_table = ProcessorTable::new(
            field,
            processor_matrix.len() as u128,
            generator,
            order,
            processor_matrix,
        );
        processor_table.pad();
        assert!(processor_table.table.matrix.len().is_power_of_two());
        assert_eq!(
            processor_table.table.matrix[(processor_table.table.matrix.len() - 1) as usize]
                [Indices::Cycle as usize]
                .0,
            processor_table.table.matrix.len() as u128 - 1
        );
    }

    #[test]
    fn test_air() {
        let field = Field::new(18446744069414584321);
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        let vm = VirtualMachine::new(field);
        let generator = field.generator();
        let order = 1 << 32;
        // let code = "++>+++++[<+>-]++++++++[<++++++>-]<.".to_string();
        let code2 = ">>[++-]<".to_string();
        let program = vm.compile(code2);
        vm.run(&program, "".to_string());
        // assert_eq!(program.len(), 2);
        let (processor_matrix, memory_matrix, instruction_matrix, input_matrix, output_matrix) =
            vm.simulate(&program, "".to_string());
        let challenges = vec![one; 11];
        let mut processor_table = ProcessorTable::new(
            field,
            processor_matrix.len() as u128,
            generator,
            order,
            processor_matrix,
        );
        let mut memory_table = Memory::new(
            field,
            memory_matrix.len() as u128,
            generator,
            order,
            memory_matrix,
        );
        let mut instruction_table = InstructionTable::new(
            field,
            instruction_matrix.len() as u128,
            generator,
            order,
            instruction_matrix,
        );
        let mut input_table = IOTable::new(
            field,
            input_matrix.len() as u128,
            generator,
            order,
            input_matrix,
        );
        let mut output_table = IOTable::new(
            field,
            output_matrix.len() as u128,
            generator,
            order,
            output_matrix,
        );

        processor_table.pad();
        memory_table.pad();
        instruction_table.pad();
        input_table.pad();
        output_table.pad();

        println!("processor table before extending columns");
        for row in processor_table.table.matrix.clone() {
            println!("{:?}", row);
        }
        processor_table.extend_columns(challenges.clone());
        println!("processor table after extending columns");
        for row in processor_table.table.matrix.clone() {
            println!("{:?}", row);
        }
        let air = processor_table.generate_air(challenges);
        println!("air");
        for row in air {
            println!("{:?}", row);
        }
    }
}
