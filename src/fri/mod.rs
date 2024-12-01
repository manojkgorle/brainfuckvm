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
    fn new(&self, offset: FieldElement,omega: FieldElement, length: usize)->Self{
        Self { offset, omega, length }
    }
    fn call(&self, index: usize)->FieldElement{
        let _x = self.omega.pow(index as u128)*self.offset;
        _x
    }
}