use crate::fields::*;
use crate::merkle::*;
//we can use ntt fast fns for optimization, but rn we just implement using direct evaluate and multiply functions of polynomials
use crate::univariate_polynomial::*;

pub struct Fri{
    
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
        let mut result: Vec<FieldElement> = vec![];
        for i in 1..polynomial.coefficients.len(){
            result.push(polynomial.evaluate(self.omega.pow(i as u128)));
        }
        result
    }
    pub fn xevaluate(&self, polynomial: Polynomial)-> Vec<FieldElement>{
        let polynomial = polynomial.scalar_mul(self.offset);
        let mut result: Vec<FieldElement> = vec![];
        for i in 1..polynomial.coefficients.len(){
            result.push(polynomial.evaluate(self.omega.pow(i as u128)));
        }
        result
    }
    

}