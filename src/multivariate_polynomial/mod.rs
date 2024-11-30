
#[allow(dead_code)]
use std::collections::HashMap;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use crate::univariate_polynomial::*;
use crate::fields::{Field, FieldElement};
#[derive(Clone, Debug)]
struct MPolynomial {
    dictionary: HashMap<Vec<u128>, FieldElement>,  // Exponent vector -> coefficient
}
impl MPolynomial{
  // Constructor for creating a new MPolynomial
  pub fn new(dictionary: HashMap<Vec<u128>, FieldElement>) -> Self {
    MPolynomial { dictionary }
}

// Function to create a zero polynomial
pub fn zero() -> Self {
    MPolynomial {
        dictionary: HashMap::new(), // Empty dictionary represents a zero polynomial
    }
}

}
impl Add for MPolynomial{
    type Output = Self;
     fn add(self, other: MPolynomial) -> MPolynomial {
      let mut dictionary = HashMap::new();

      // Find the maximum number of variables
      let num_variables = std::cmp::max(
          self.dictionary.keys().map(|k| k.len()).max().unwrap_or(0),
          other.dictionary.keys().map(|k| k.len()).max().unwrap_or(0),
      );

      // Add the first polynomial's terms to the new dictionary
      for (k, v) in &self.dictionary {
          let mut pad = k.clone();
          pad.extend(std::iter::repeat(0).take(num_variables - k.len()));
          dictionary.insert(pad, v.clone());
      }

      // Add the second polynomial's terms to the new dictionary
      for (k, v) in &other.dictionary {
          let mut pad = k.clone();
          pad.extend(std::iter::repeat(0).take(num_variables - k.len()));
          if let Some(existing_coeff) = dictionary.get_mut(&pad) {
              // If the exponent vector already exists, add the coefficients
              *existing_coeff = existing_coeff.add(*v);
          } else {
              // If the exponent vector does not exist, insert the new one
              dictionary.insert(pad, v.clone());
          }
      }

      MPolynomial::new(dictionary)
  }}
  impl AddAssign for MPolynomial {
    fn add_assign(&mut self, other: Self) {
        // Find the maximum number of variables
        let num_variables = std::cmp::max(
            self.dictionary.keys().map(|k| k.len()).max().unwrap_or(0),
            other.dictionary.keys().map(|k| k.len()).max().unwrap_or(0),
        );

        // Add the second polynomial's terms to the first polynomial's dictionary
        for (k, v) in other.dictionary {
            let mut pad = k.clone();
            pad.extend(std::iter::repeat(0).take(num_variables - k.len()));
            if let Some(existing_coeff) = self.dictionary.get_mut(&pad) {
                // If the exponent vector already exists, add the coefficients
                *existing_coeff = existing_coeff.add(v);
            } else {
                // If the exponent vector does not exist, insert the new one
                self.dictionary.insert(pad, v);
            }
        }
    }}
    
  impl PartialEq for MPolynomial {
    fn eq(&self, other: &Self) -> bool {
        // Step 1: Check if the sizes of the dictionaries are equal
        if self.dictionary.len() != other.dictionary.len() {
            return false;
        }

        // Step 2: Check each key-value pair
        for (key, value) in &self.dictionary {
            // Check if the key exists in `other.dictionary`
            if let Some(other_value) = other.dictionary.get(key) {
                // Compare the values using `is_equal` or another comparison method
                if value.1 != other_value.1 {
                    return false;
                }
            } else {
                // Key not found in `other.dictionary`
                return false;
            }
        }

        // Step 3: If all checks pass, return true
        true
    }
}


#[cfg(test)]
mod test_MPolynomial_operation{
  use super::*;
  #[test]
  fn test_addition(){
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![1, 2], FieldElement::new(3, Field::new(5)));
    dictionary.insert(vec![2,1], FieldElement::new(3, Field::new(5)));
    let p1 = MPolynomial::new(dictionary);
    let mut dictionary2 = HashMap::new();
    dictionary2.insert(vec![2,1], FieldElement::new(2, Field::new(5)));
    dictionary2.insert(vec![1, 2], FieldElement::new(3, Field::new(5)));
  
    let p2 = MPolynomial::new(dictionary2);
    let mut dictionary3 = HashMap::new();
    dictionary3.insert(vec![1, 2], FieldElement::new(6, Field::new(5)));
    dictionary3.insert(vec![2, 1], FieldElement::new(5, Field::new(5)));
    let p3 = MPolynomial::new(dictionary3);
    assert_eq!(p1 + p2, p3);
  }
  #[test]
  fn test_addition_pad(){
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![1, 2,0,0], FieldElement::new(3, Field::new(7)));
    dictionary.insert(vec![2,1], FieldElement::new(3, Field::new(7)));
    let p1 = MPolynomial::new(dictionary);
    let mut dictionary2 = HashMap::new();
    dictionary2.insert(vec![2,1], FieldElement::new(2, Field::new(7)));
    dictionary2.insert(vec![1, 2], FieldElement::new(3, Field::new(7)));
  
    let p2 = MPolynomial::new(dictionary2);
    let mut dictionary3 = HashMap::new();
    dictionary3.insert(vec![1, 2,0,0], FieldElement::new(6, Field::new(7)));
    dictionary3.insert(vec![2, 1,0,0], FieldElement::new(5, Field::new(7)));
    let p3 = MPolynomial::new(dictionary3);
    assert_eq!(p1 + p2, p3);
  }
  #[test]
  fn test_eq(){
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![1, 2,3], FieldElement::new(3, Field::new(5)));
    let p1 = MPolynomial::new(dictionary);
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![1, 2,3], FieldElement::new(3, Field::new(5)));
    let p2 = MPolynomial::new(dictionary);
    assert_eq!(p1, p2);
  }
}

