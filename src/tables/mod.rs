use crate::fields::{Field, FieldElement};
use crate::multivariate_polynomial::*;
use crate::univariate_polynomial::{interpolate_lagrange_polynomials, Polynomial};
use rand::*;
use crate::fri::*;
use crate::ntt::*;
pub mod instruction;
pub mod memory;
pub mod processor;
pub mod io;
#[derive(Debug, Clone)]
// we are not going to use the randomizers in the table
pub struct Table {
    pub field: Field,
    base_width: u128,//number of base columns in the table.
    full_width: u128,//total no. of coloumn using extension and all
    length: u128,//Number of rows in the table.
     height: u128,//Represents the rounded-up next power of two of the table length 
    omicron: FieldElement,//represent the generator eval_domain depends on the generator and the order of the subgroup
    generator: FieldElement,// A generator for the multiplicative group of the field
    order: u128,//order of the generator.
    matrix: Vec<Vec<FieldElement>>,
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
            matrix,
        }
    }
    fn new_2(
        field: Field, 
        base_width: u128, 
        full_width: u128, 
        length: u128,
        generator:FieldElement, 
        order: u128
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
        matrix: Vec::new(), // Initialize as empty
    }
    }
    // dont know how to implement this method
    pub fn unit_distance(&self, omega_order:u128)->u128{
        if self.height ==0{
            return 0;
        }
        omega_order/self.height
    }

    pub fn get_interpolating_domain_length( &self)->u128{
        self.height
    }

    pub fn interpolate_degree(self)->u128{
        self.get_interpolating_domain_length()-1
    }

    pub fn has_order_po2(order: u128) -> bool {
        (order & (order - 1)) == 0
    }
    
    pub fn interpolate_columns(self, column_indices:Vec<u128>)->Vec<Polynomial>{
        let mut polynomial:Vec<Polynomial>=Vec::new();
        if self.height ==0{
            let poly=Polynomial::new_from_coefficients(vec![FieldElement::zero(Field::new(self.field.0))]);
            polynomial.push(poly);
           return  polynomial;
        }

        let mut omicron_domain:Vec<FieldElement>=Vec::new();
        for i in 0..self.height{
            omicron_domain.push(self.omicron.pow(i));
        }
        
        for c in column_indices{
            let mut trace:Vec<FieldElement>=Vec::new();
            for row in self.matrix.iter(){
                trace.push(row[c as usize]);
            }
            let values:Vec<FieldElement>=trace.clone();
            if values.len()!=omicron_domain.len(){
                panic!("length of domain and values are unequal");
            };
            let poly= interpolate_lagrange_polynomials(omicron_domain.clone(), values);
            println!("poly ={:?}", poly);
            polynomial.push(poly);  
        }
        polynomial
    }

    pub fn next_interpolate_columns(self, column_indices: Vec<u128>)->Vec<Polynomial>{
    let mut next_interpolated:Vec<Polynomial>=Vec::new();
        if self.height ==0{
            let poly=Polynomial::new_from_coefficients(vec![FieldElement::zero(Field::new(self.field.0))]);
            next_interpolated.push(poly);
           return next_interpolated;
        }
        
        let mut omicron_domain:Vec<FieldElement>=Vec::new();
        omicron_domain.push(self.omicron.pow(self.height-1));
        for i in 0..self.height-1{
            omicron_domain.push(self.omicron.pow(i));
        }

        for c in 0..self.matrix[0].len(){
            let mut trace:Vec<FieldElement>=Vec::new();
            for row in self.matrix.iter(){
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
        next_interpolated
    }

// in the codewords matrix you are getting is coloumns*rows (rows and coloums reference to the intial table)
pub fn lde(self,domain:FriDomain)->Vec<Vec<FieldElement>>{
    let x = self.base_width;
    let polynomial = self.interpolate_columns((0..x).collect());
    let mut self_codewords=Vec::new();
    for p in polynomial{
        let codeword=domain.evaluate(p);
        self_codewords.push(codeword);}
    self_codewords
}}



pub fn roundup_npow2( len:u128)->u128{
    if len==0{
        return 0;
    }else if len == 1 {
        return 1;
    }
    // Calculate the next power of two
    let bit_representation = format!("{:b}", len - 1);
    1 << (bit_representation.len() as u128)
}

// mutable or clone doubt
pub fn derive_omicron(generator:FieldElement,generator_order:u128,target_order:u128)->FieldElement{
    let mut t_generator=generator;
    let mut t_order=generator_order;
     while t_order!=target_order {
        t_generator=t_generator.pow(2);
            t_order/=2;
            println!("t_order ={:?}", t_order);
    }
    t_generator}
    
pub fn has_order_po2( order: u128) -> bool {
    (order & (order - 1)) == 0
}

#[cfg(test)]
mod test_operations{
    use super::*;



    #[test]
    fn test_roundup_npow2(){
        let len:u128 =2;
        let len2:u128 =6;
        let len3=9;
        let round= roundup_npow2(len);
     
        let round2= roundup_npow2(len2);
        let round3= roundup_npow2(len3);
        println!("round2:{}",round2);
        assert_eq!(round,2);
        assert_eq!(round2,8);

        assert_eq!(round3,16);
    }
    #[test]
    fn has_order_po2(){
        let order =4_u128;
        let order2 =5_u128;
        assert !(!Table::has_order_po2(order2));
        assert !(Table::has_order_po2(order));
    }
    #[test]
    fn test_interpolate_columns(){
        let field = Field::new(17);
        let   offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        let domain = FriDomain::new(offset, omega, length);
        let generator = FieldElement::new(3, field);
        let order = 16;
        let mut table = Table::new_2(field, 3, 5, 4, generator, order);
        let mut matrix: Vec<Vec<FieldElement>> = Vec::new();
        matrix.push(vec![FieldElement::new(1, field), FieldElement::new(2, field), FieldElement::new(3, field), FieldElement::new(4, field), FieldElement::new(5, field)]);
        matrix.push(vec![FieldElement::new(6, field), FieldElement::new(7, field), FieldElement::new(8, field), FieldElement::new(9, field), FieldElement::new(10, field)]);
        matrix.push(vec![FieldElement::new(11, field), FieldElement::new(12, field), FieldElement::new(13, field), FieldElement::new(14, field), FieldElement::new(15, field)]);
        matrix.push(vec![FieldElement::new(16, field), FieldElement::new(17, field), FieldElement::new(18, field), FieldElement::new(19, field), FieldElement::new(20, field)]);
        table.matrix = matrix;
        let column_indices = vec![0, 1, 2];
        let polynomials = table.interpolate_columns(column_indices);
        let expected_polynomials = [Polynomial::new_from_coefficients(vec![FieldElement::new(0, field), FieldElement::new(13, field), FieldElement::new(6, field), FieldElement::new(16, field)])];
            println!("polynomials ={:?}", polynomials[0]);
        assert_eq!(polynomials[0], expected_polynomials[0]);

    }
    #[test]
    fn test_lde(){
        let field = Field::new(17);
        let   offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        let domain = FriDomain::new(offset, omega, length);
        let generator = FieldElement::new(3, field);
        let order = 16;
        let mut table = Table::new_2(field, 3, 5, 4, generator, order);
        let mut matrix: Vec<Vec<FieldElement>> = Vec::new();
        matrix.push(vec![FieldElement::new(1, field), FieldElement::new(2, field), FieldElement::new(3, field), FieldElement::new(4, field), FieldElement::new(5, field)]);
        matrix.push(vec![FieldElement::new(6, field), FieldElement::new(7, field), FieldElement::new(8, field), FieldElement::new(9, field), FieldElement::new(10, field)]);
        matrix.push(vec![FieldElement::new(11, field), FieldElement::new(12, field), FieldElement::new(13, field), FieldElement::new(14, field), FieldElement::new(15, field)]);
        matrix.push(vec![FieldElement::new(16, field), FieldElement::new(17, field), FieldElement::new(18, field), FieldElement::new(19, field), FieldElement::new(20, field)]);
        table.matrix = matrix;
        let codewords = table.lde(domain);
        let expected_codewords = vec![
            vec![FieldElement::new(2, field), FieldElement::new(12, field), FieldElement::new(5, field), FieldElement::new(15, field)],
            vec![FieldElement::new(4, field), FieldElement::new(14, field), FieldElement::new(7, field), FieldElement::new(0, field)],
            vec![FieldElement::new(6, field), FieldElement::new(16, field), FieldElement::new(9, field), FieldElement::new(2, field)],
          
        ];
        assert_eq!(codewords, expected_codewords);

        
    }
    #[test]
    fn test_derive_omicron(){
        let field = Field::new( 1 + (1 << 64) - (1 << 32) );
        let generator = FieldElement::new(1753635133440165772, field);
        let order = 1 << 32 ;
        let target_order = 1 << 16;
        let omicron = derive_omicron(generator, order, target_order);
        println!("omicron ={:?}", omicron);
        let expected_omicron = FieldElement::new(49, field);
        assert_eq!(omicron, expected_omicron);
    }

}

