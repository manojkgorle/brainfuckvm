use crate::fields::{Field, FieldElement};
use crate::multivariate_polynomial::*;
use crate::univariate_polynomial::{interpolate_lagrange_polynomials, Polynomial};
use rand::*;
use crate::fri::*;
use crate::ntt::*;
pub mod instruction;
pub mod memory;
pub mod processor;
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
    fn new_2(field: Field, base_width: u128, full_width: u128, length: u128,
         generator:FieldElement, order: u128) -> Self {
     let height = roundup_npow2(length);
     let omicron = derive_omicron(generator, order, height);
     
     Table {
         field,
         base_width,
         full_width,
         length,
        //  num_randomizers,
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
        let mut domain:Vec<FieldElement>=vec![FieldElement::new(0,self.field);self.height as usize];
        for i in 0..self.height{
            domain[i as usize]=omicron_domain[i as usize];
        }

        for c in column_indices{
            let mut trace:Vec<FieldElement>=Vec::new();
     
            for row in self.matrix.iter(){
                trace.push(row[c as usize]);
            }
         let mut values:Vec<FieldElement>=Vec::new();
            let mut sum =FieldElement::zero(Field::new(self.field.0));
            let zero = FieldElement::zero(Field::new(self.field.0));
     for i in 0..self.height {
                if i < trace.len() as u128 {
                    // Push values from trace
                    values.push(trace[i as usize].clone());
                } else {
                    // Push the last value of trace or zero if trace is empty
                    let last_value = trace.last().unwrap_or(&zero).clone();
                    values.push(last_value);
                }
            }
            if values.len()!=omicron_domain.len(){
                panic!("length of domain and values are unequal");
            };
        let poly= interpolate_lagrange_polynomials(domain.clone(), values);
            polynomial.push(poly);  
        }
        polynomial
}

}

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

    while t_order!=target_order{
        t_generator=t_generator.pow(2);
            t_order/=2;
    }
    t_generator
}
    
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

}

