use crate::fields::*;
use crate::merkle::*;
//we can use ntt fast fns for optimization, but rn we just implement using direct evaluate and multiply functions of polynomials
use crate::univariate_polynomial::*;

//some fns, structs and parameters are changed accordingly compared to fri.py and ntt.py, because we are not using extension fields
pub struct Fri{
    offset: FieldElement,
    omega: FieldElement,
    initial_domain_length: usize,
    domain: FriDomain,
    num_colinearity_tests: usize,
}

impl Fri{
    pub fn new(&self, offset: FieldElement, omega: FieldElement, initial_domain_length: usize, num_colinearity_tests: usize)-> Self{
        Self { offset: (offset), omega: (omega), initial_domain_length: (initial_domain_length), domain: (FriDomain::new(offset, omega, initial_domain_length)), num_colinearity_tests: (num_colinearity_tests) }
        assert!(self.num_rounds()>=1, "cannot do FRI with less than one round");
    }
    pub fn num_rounds(&self)-> usize{
        let mut codeword_len = self.initial_domain_length;
        let mut num =0;
        while codeword_len>1 as usize{
            codeword_len= codeword_len/2;
            num+=1;
        } 
        num
    }
    pub fn sample_index(byte_array: [u8], size:usize)->u8{
        
    }

}
pub struct FriDomain{
        offset: FieldElement,
        omega: FieldElement,
        length: usize,
}

impl FriDomain{
    pub fn new(offset: FieldElement,omega: FieldElement, length: usize)->Self{
        Self { offset, omega, length }
    }
    pub fn call(&self, index: usize)->FieldElement{
        let _x = self.omega.pow(index as u128)*self.offset;
        _x
    }
    pub fn list(&self)->Vec<FieldElement>{
        let mut list: Vec<FieldElement>= vec![];
        for i in 1..self.length{
        list.push(self.omega.pow(i as u128)*self.offset);
        }
        list
    }
    pub fn evaluate(&self, polynomial: Polynomial)-> Vec<FieldElement>{
        let polynomial = polynomial.scalar_mul(self.offset);
        let mut result: Vec<FieldElement> = vec![];
        for i in 1..polynomial.coefficients.len(){
            result.push(polynomial.evaluate(self.omega.pow(i as u128)));
        }
        result
    }
    //needed if we use extension field
    pub fn xevaluate(&self, polynomial: Polynomial)-> Vec<FieldElement>{
        let polynomial = polynomial.scalar_mul(self.offset);
        let mut result: Vec<FieldElement> = vec![];
        for i in 1..polynomial.coefficients.len(){
            result.push(polynomial.evaluate(self.omega.pow(i as u128)));
        }
        result
    }
    pub fn interpolate(&self, values: Vec<FieldElement>)->Polynomial{
        let mut list:Vec<FieldElement> =vec![];
        for i in 1..values.len(){
            list.push(self.omega.pow(i as u128));
        }
        let polynomial = interpolate_lagrange_polynomials(list, values).scalar_mul(self.offset.inverse());
        polynomial
    }
    //not written xinterpolate, as it is used for extension field
}