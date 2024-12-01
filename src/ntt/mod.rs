use crate::univariate_polynomial::*;
use crate:: fields::*;

pub fn batch_inverse(array: &Vec<FieldElement>)->Vec<FieldElement>{
    assert!(
        array.iter().all(|a| a.0!=0 as u128),
        "batch inverse does not work when input contains a zero"
    );

    let mut products = array.clone();

    for i in 1..products.len() {
        products[i] = products[i - 1].clone() * array[i].clone();
    }

    let mut acc = products[products.len()-1].inverse();

    for i in (1..array.len()).rev() {
        products[i] = acc * products[i - 1];
        acc *= array[i].clone();
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
            let product = a.clone()* inv.clone();
            assert_eq!(product, FieldElement::new(1, field));
        }
}
}
