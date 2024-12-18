use crate::fields::*;
use crate::univariate_polynomial::*;

//ntt is transforming coefficient form to point evaluation form, so that multiplication of polynomials can happen in O(nlogn), rather than O(n^2)
//We have already written functions like lagrange interpolation, multiply polynomial and evaluate polynomial, so we will write these functions mainly for optimization
//Will implement the ntt functions later on to improve time complexity of evaluation and multiplication

pub fn batch_inverse(array: &Vec<FieldElement>) -> Vec<FieldElement> {
    assert!(
        array.iter().all(|a| a.0 != 0_u128),
        "batch inverse does not work when input contains a zero"
    );

    let mut products = array.clone();

    for i in 1..products.len() {
        products[i] = products[i - 1] * array[i];
    }

    let mut acc = products[products.len() - 1].inverse();

    for i in (1..array.len()).rev() {
        products[i] = acc * products[i - 1];
        acc *= array[i];
    }

    products[0] = acc;

    products
}

#[cfg(test)]
mod test_ntt {
    use super::*;

    #[test]
    fn test_batch_inverse() {
        let field = Field::new(1000000007);
        let array = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
        ];

        let inverses = batch_inverse(&array);

        for (a, inv) in array.iter().zip(inverses.iter()) {
            let product = *a * *inv;
            assert_eq!(product, FieldElement::new(1, field));
        }
    }
}
