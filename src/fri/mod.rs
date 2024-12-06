use crate::fields::*;
use crate::merkle::*;
//we can use ntt fast fns for optimization, but rn we just implement using direct evaluate and multiply functions of polynomials
use crate::univariate_polynomial::*;
use blake2::{Blake2b512, Digest};
//some fns, structs and parameters are changed accordingly compared to fri.py and ntt.py, because we are not using extension fields
pub struct Fri{
    offset: FieldElement,
    omega: FieldElement,
    initial_domain_length: u128,
    domain: FriDomain,
    num_colinearity_tests: usize,
}

impl Fri{
    pub fn new(offset: FieldElement, omega: FieldElement, initial_domain_length: u128, num_colinearity_tests: usize)-> Self{
        let result = Fri{ offset: (offset), omega: (omega), initial_domain_length: (initial_domain_length), domain: (FriDomain::new(offset, omega, initial_domain_length)), num_colinearity_tests: (num_colinearity_tests) };
        assert!(result.num_rounds()>=1, "cannot do FRI with less than one round");
        result
    }
    pub fn num_rounds(&self)-> usize{
        let mut codeword_len = self.initial_domain_length;
        let mut num =0;
        while codeword_len>1_u128{
            codeword_len /= 2;
            num+=1;
        } 
        num
    }

}

pub struct FriDomain{
      pub  offset: FieldElement,
       pub  omega: FieldElement,
       pub length: u128,
}

impl FriDomain{
    pub fn new(offset: FieldElement,omega: FieldElement, length: u128)->Self{
        Self { offset, omega, length }
    }

    pub fn call(&self, index: usize)->FieldElement{
        self.omega.pow(index as u128)*self.offset
    }
    
    pub fn list(&self)->Vec<FieldElement>{
        let mut list: Vec<FieldElement>= vec![];
        for i in 0..self.length{
        list.push(self.omega.pow(i)*self.offset);
        }
        list
    }
    pub fn evaluate(&self, polynomial: Polynomial)-> Vec<FieldElement>{

        let polynomial = polynomial.scalar_mul(self.offset);
        let mut result: Vec<FieldElement> = vec![];
        // let mut  domain_gen:Vec<FieldElement>=Vec::new();
        for i in 0..self.length{
            // let x=self.omega.pow(i);
            result.push(polynomial.evaluate(self.omega.pow(i)));
            // domain_gen.push(x);
        }
        // println!("domain_gen ={:?}", domain_gen);
      
        result
    }
    //needed if we use extension field
    // pub fn xevaluate(&self, polynomial: Polynomial)-> Vec<FieldElement>{
    //     let polynomial = polynomial.scalar_mul(self.offset);
    //     let mut result: Vec<FieldElement> = vec![];
    //     for i in 0..polynomial.coefficients.len(){
    //         result.push(polynomial.evaluate(self.omega.pow(i as u128)));
    //     }
    //     result
    // }

    pub fn interpolate(&self, values: Vec<FieldElement>)->Polynomial{
        let mut list:Vec<FieldElement> =vec![];
        for i in 0..values.len(){
            list.push(self.omega.pow(i as u128));
        }
        
        interpolate_lagrange_polynomials(list, values).scalar_mul(self.offset.inverse())
    }
    //not written xinterpolate, as it is used for extension field
}

#[cfg(test)]
mod test_fri{
    use super::*;
    #[test]
    fn test_evaluate(){
        let field = Field::new(17);
        let   offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        println!("omega ={:?}", omega);
        let domain = FriDomain::new(offset, omega, length);
        let polynomial = Polynomial::new_from_coefficients(vec![FieldElement::new(1, field), FieldElement::new(2, field)]);
        let values = domain.evaluate(polynomial.clone());
        let finded=vec![FieldElement::new(6, field), FieldElement::new(3, field), FieldElement::new(15, field), FieldElement::new(1, field)];

        println!("values ={:?}", values);
        assert_eq!(values, finded);
        

    }
    #[test]
    fn test_interpolate(){
        let field = Field::new(17);
        let   offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        println!("omega ={:?}", omega);
        let domain = FriDomain::new(offset, omega, length);
        let polynomial = Polynomial::new_from_coefficients(vec![FieldElement::new(1, field), FieldElement::new(2, field)]);
        let values = domain.evaluate(polynomial.clone());
        let finded=vec![FieldElement::new(6, field), FieldElement::new(3, field)];
        let interpolated = domain.interpolate(finded);
        println!("interpolated ={:?}", interpolated);
        assert_eq!(interpolated.coefficients, polynomial.coefficients);
    }
}