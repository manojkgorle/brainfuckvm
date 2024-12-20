use crate::fields::{Field, FieldElement};
use crate::fri::*;
use crate::univariate_polynomial::{interpolate_lagrange_polynomials, Polynomial};
use chrono::Local;
use rand::*;
pub mod instruction;
pub mod io;
pub mod memory;
use rayon::{prelude::*, vec};
pub mod processor;
#[derive(Debug, Clone)]
// we are not going to use the randomizers in the table
pub struct Table {
    pub field: Field,
    base_width: u128,                  //number of base columns in the table.
    full_width: u128,                  //total no. of coloumn using extension and all
    length: u128,                      //Number of rows in the table.
    pub height: u128, //Represents the rounded-up next power of two of the table length
    omicron: FieldElement, //represent the generator eval_domain depends on the generator and the order of the subgroup
    generator: FieldElement, // A generator for the multiplicative group of the field
    order: u128,           //order of the generator.
    omicron_domain: Vec<FieldElement>, //The domain of the table
    pub matrix: Vec<Vec<FieldElement>>,
}

impl Table {
    // Constructor method to create a new instance of `table`
    pub fn new(
        field: Field,
        base_width: u128,
        full_width: u128,
        length: u128,
        height: u128,
        omicron: FieldElement,
        generator: FieldElement,
        order: u128,
        matrix: Vec<Vec<FieldElement>>,
    ) -> Self {
        Table {
            field,
            base_width,
            full_width,
            length,

            height,
            omicron,
            generator,
            order,
            omicron_domain: Vec::new(),
            matrix,
        }
    }
    pub fn new_2(
        field: Field,
        base_width: u128,
        full_width: u128,
        length: u128,
        generator: FieldElement,
        order: u128,
    ) -> Self {
        let height = roundup_npow2(length);
        let omicron = derive_omicron(generator, order, height);
        Table {
            field,
            base_width,
            full_width,
            length,
            height,
            omicron,
            generator,
            order,
            omicron_domain: Vec::new(),
            matrix: Vec::new(), // Initialize as empty
        }
    }

    pub fn get_interpolating_domain_length(&self) -> u128 {
        self.height
    }

    pub fn interpolate_degree(self) -> u128 {
        self.get_interpolating_domain_length() - 1
    }

    pub fn has_order_po2(order: u128) -> bool {
        (order & (order - 1)) == 0
    }

    pub fn generate_omicron_domain(&mut self) {
        let mut omicron_domain: Vec<FieldElement> = Vec::with_capacity(self.height as usize);
        for i in 0..self.height {
            omicron_domain.push(self.omicron.pow(i));
        }
        self.omicron_domain = omicron_domain;
    }

    // @todo optimize this
    pub fn interpolate_columns(self, column_indices: Vec<u128>) -> Vec<Polynomial> {
        let t = Local::now();
        if self.height == 0 {
            let poly = Polynomial::new_from_coefficients(vec![FieldElement::zero(self.field)]);
            return vec![poly];
        }

        let polynomial = column_indices
            .par_iter()
            .map(|c| {
                let mut trace: Vec<FieldElement> = Vec::new();
                let omicron_domain: Vec<FieldElement> = self.omicron_domain.clone();
                for row in self.matrix.iter() {
                    trace.push(row[*c as usize]);
                }
                if trace.len() != omicron_domain.len() {
                    panic!("length of domain and values are unequal");
                };
                let t2 = Local::now();
                let poly = interpolate_lagrange_polynomials(omicron_domain, trace);
                log::info!(
                    "Interpolating lagrange polynomials took: {}ms",
                    (Local::now() - t2).num_milliseconds()
                );
                poly
            })
            .collect();
        log::info!(
            "Interpolating columns took: {}ms",
            (Local::now() - t).num_milliseconds()
        );
        polynomial
    }

    pub fn next_interpolate_columns(self, interpolated: Vec<Polynomial>) -> Vec<Polynomial> {
        let t = Local::now();
        let mut next_interpolated: Vec<Polynomial> = Vec::new();
        for i in interpolated {
            let next_interpolate = i.compose(self.omicron);
            next_interpolated.push(next_interpolate)
        }
        log::debug!(
            "Next Interpolating columns took: {}ms",
            (Local::now() - t).num_milliseconds()
        );
        next_interpolated
    }

    // in the codewords matrix you are getting is coloumns*rows (rows and coloums reference to the intial table)
    pub fn lde(self, domain: FriDomain) -> Vec<Vec<FieldElement>> {
        let x = self.base_width;
        let polynomial = self.interpolate_columns((0..x).collect());
        let mut self_codewords = Vec::new();
        for p in polynomial {
            let codeword = domain.evaluate(p);
            self_codewords.push(codeword);
        }
        self_codewords
    }
}

pub fn roundup_npow2(len: u128) -> u128 {
    if len == 0 {
        return 0;
    } else if len == 1 {
        return 1;
    }
    // Calculate the next power of two
    let bit_representation = format!("{:b}", len - 1);
    1 << (bit_representation.len() as u128)
}

// mutable or clone doubt
pub fn derive_omicron(
    generator: FieldElement,
    generator_order: u128,
    target_order: u128,
) -> FieldElement {
    let mut t_generator = generator;
    let mut t_order = generator_order;
    while t_order != target_order {
        t_generator = t_generator.pow(2);
        t_order /= 2;
    }
    t_generator
}

pub fn has_order_po2(order: u128) -> bool {
    (order & (order - 1)) == 0
}

#[cfg(test)]
mod test_operations {
    #![allow(unused_variables)]
    use super::*;

    #[test]
    fn test_roundup_npow2() {
        let len: u128 = 2;
        let len2: u128 = 6;
        let len3 = 9;
        let round = roundup_npow2(len);

        let round2 = roundup_npow2(len2);
        let round3 = roundup_npow2(len3);
        println!("round2:{}", round2);
        assert_eq!(round, 2);
        assert_eq!(round2, 8);

        assert_eq!(round3, 16);
    }

    #[test]
    fn has_order_po2() {
        let order = 4_u128;
        let order2 = 5_u128;
        assert!(!Table::has_order_po2(order2));
        assert!(Table::has_order_po2(order));
    }

    #[test]
    fn test_interpolate_columns() {
        let field = Field::new(17);
        let offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        let domain = FriDomain::new(offset, omega, length);
        let generator = FieldElement::new(3, field);
        let order = 16;
        let mut table = Table::new_2(field, 3, 5, 4, generator, order);
        let mut matrix: Vec<Vec<FieldElement>> = Vec::new();
        matrix.push(vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
            FieldElement::new(5, field),
        ]);
        matrix.push(vec![
            FieldElement::new(6, field),
            FieldElement::new(7, field),
            FieldElement::new(8, field),
            FieldElement::new(9, field),
            FieldElement::new(10, field),
        ]);
        matrix.push(vec![
            FieldElement::new(11, field),
            FieldElement::new(12, field),
            FieldElement::new(13, field),
            FieldElement::new(14, field),
            FieldElement::new(15, field),
        ]);
        matrix.push(vec![
            FieldElement::new(16, field),
            FieldElement::new(17, field),
            FieldElement::new(18, field),
            FieldElement::new(19, field),
            FieldElement::new(20, field),
        ]);
        table.matrix = matrix;
        table.generate_omicron_domain();
        let column_indices = vec![0, 1, 2];
        let polynomials = table.interpolate_columns(column_indices);
        let expected_polynomials = [Polynomial::new_from_coefficients(vec![
            FieldElement::new(0, field),
            FieldElement::new(13, field),
            FieldElement::new(6, field),
            FieldElement::new(16, field),
        ])];
        println!("polynomials ={:?}", polynomials[0]);
        assert_eq!(polynomials[0], expected_polynomials[0]);
        let mut next_polynomials: Vec<Polynomial> = vec![];
        for i in polynomials {
            let next_polynomial = i.compose(omega);
            next_polynomials.push(next_polynomial);
        }
        println!("{:?}", next_polynomials[0]);
    }
    #[test]
    fn test_next_interpolate() {
        let field = Field::new(17);
        let offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        let domain = FriDomain::new(offset, omega, length);
        let generator = FieldElement::new(3, field);
        let order = 16;
        let mut table = Table::new_2(field, 3, 5, 4, generator, order);
        let mut matrix: Vec<Vec<FieldElement>> = Vec::new();
        matrix.push(vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
            FieldElement::new(5, field),
        ]);
        matrix.push(vec![
            FieldElement::new(6, field),
            FieldElement::new(7, field),
            FieldElement::new(8, field),
            FieldElement::new(9, field),
            FieldElement::new(10, field),
        ]);
        matrix.push(vec![
            FieldElement::new(11, field),
            FieldElement::new(12, field),
            FieldElement::new(13, field),
            FieldElement::new(14, field),
            FieldElement::new(15, field),
        ]);
        matrix.push(vec![
            FieldElement::new(16, field),
            FieldElement::new(17, field),
            FieldElement::new(18, field),
            FieldElement::new(19, field),
            FieldElement::new(20, field),
        ]);
        table.matrix = matrix;
        let column_indices = vec![0, 1, 2];
        // let polynomials = table.next_interpolate_columns(column_indices);
        // println!("polynomial_next{:?}",polynomials[0]);
    }

    #[test]
    fn test_lde() {
        let field = Field::new(17);
        let offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        let domain = FriDomain::new(offset, omega, length);
        let generator = FieldElement::new(3, field);
        let order = 16;
        let mut table = Table::new_2(field, 3, 5, 4, generator, order);
        let mut matrix: Vec<Vec<FieldElement>> = Vec::new();
        matrix.push(vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
            FieldElement::new(5, field),
        ]);
        matrix.push(vec![
            FieldElement::new(6, field),
            FieldElement::new(7, field),
            FieldElement::new(8, field),
            FieldElement::new(9, field),
            FieldElement::new(10, field),
        ]);
        matrix.push(vec![
            FieldElement::new(11, field),
            FieldElement::new(12, field),
            FieldElement::new(13, field),
            FieldElement::new(14, field),
            FieldElement::new(15, field),
        ]);
        matrix.push(vec![
            FieldElement::new(16, field),
            FieldElement::new(17, field),
            FieldElement::new(18, field),
            FieldElement::new(19, field),
            FieldElement::new(20, field),
        ]);
        table.matrix = matrix;
        table.generate_omicron_domain();
        let codewords = table.lde(domain);

        let expected_codewords = vec![
            vec![
                FieldElement::new(8, field),
                FieldElement::new(10, field),
                FieldElement::new(6, field),
                FieldElement::new(10, field),
            ],
            vec![
                FieldElement::new(9, field),
                FieldElement::new(11, field),
                FieldElement::new(7, field),
                FieldElement::new(11, field),
            ],
            vec![
                FieldElement::new(10, field),
                FieldElement::new(12, field),
                FieldElement::new(8, field),
                FieldElement::new(12, field),
            ],
        ];
        assert_eq!(codewords, expected_codewords);
    }
    #[test]
    fn test_derive_omicron() {
        let field = Field::new(1 + (1 << 64) - (1 << 32));
        let generator = FieldElement::new(1753635133440165772, field);
        let order = 1 << 32;
        let target_order = 8;
        let omicron = derive_omicron(generator, order, target_order);
        // let expected_omicron = FieldElement::new(49, field);
        // assert_eq!(omicron, expected_omicron);
    }
}
