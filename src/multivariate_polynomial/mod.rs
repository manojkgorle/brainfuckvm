
#[allow(dead_code)]
use std::collections::HashMap;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use crate::univariate_polynomial::*;
use crate::fields::{Field, FieldElement};
#[derive(Clone, Debug)]
struct MPolynomial {
    field: Field,
    dictionary: HashMap<Vec<u128>, FieldElement>,  // Exponent vector -> coefficient
}

impl MPolynomial{
  // Constructor for creating a new MPolynomial
    pub fn new(field: Field, dictionary: HashMap<Vec<u128>, FieldElement>) -> Self {
    MPolynomial { field, dictionary }
}

// Function to create a zero polynomial
pub fn zero(field: Field) -> Self {
    MPolynomial {
        field: field,
        dictionary: HashMap::new(), // Empty dictionary represents a zero polynomial
    }
}

//Function to calculate degree
pub fn degree(&self) -> FieldElement {
    if self.dictionary.is_empty() {
        return FieldElement::zero(self.field); // Return zero for an empty polynomial
    }

    let max_degree = self
        .dictionary
        .keys() // Get the keys (vectors of exponents)
        .map(|k| k.iter().sum::<u128>()) // Calculate the sum of exponents for each term
        .max() // Find the maximum sum
        .unwrap_or(0); // Default to 0 if no keys

    FieldElement::new(max_degree as u128, self.field) // Convert the maximum degree to a FieldElement
}

//Function to evaluate at a point
pub fn evaluate(&self, point: &Vec<FieldElement>)->FieldElement{
        let field = self.field;
        let mut acc = FieldElement::zero(field);

        for (exp, coeff) in &self.dictionary {
            let mut prod = *coeff;

            assert_eq!(
                point.len(),
                exp.len(),
                "Number of elements in point ({}) does not match number of variables ({}) for polynomial.",
                point.len(),
                exp.len()
            );

            for (i, &exp) in exp.iter().enumerate() {
                prod = prod * point[i].pow(exp as u128);
            }
            acc = acc + prod;
        }
        acc
    }

    pub fn evaluate_symbolic(
        &self,
        point: &Vec<MPolynomial>,
        memo: &mut HashMap<(u128, u128), MPolynomial>,
    ) -> MPolynomial {
        let field = self.field;
        let mut acc = MPolynomial::zero(self.field); // Accumulator for the result

        for (exp, coeff) in &self.dictionary {
            let mut prod = MPolynomial::constant(*coefficient, field.clone()); // Start with the coefficient as a constant

            for (i, &exp) in exponents.iter().enumerate() {
                let mut inner_acc = MPolynomial::one(&field); // Initialize to one

                let mut j = 0;
                while (1 << j) <= exp {
                    if !memo.contains_key(&(i, 1 << j)) {
                        if j == 0 {
                            memo.insert((i, 1 << j), point[i].clone());
                        } else {
                            let prev = memo.get(&(i, 1 << (j - 1))).unwrap().clone();
                            memo.insert((i, 1 << j), prev.clone() * prev);
                        }
                    }
                    let point_power = memo.get(&(i, 1 << j)).unwrap();

                    if (exp & (1 << j)) != 0 {
                        inner_acc = inner_acc * point_power.clone();
                    }
                    j += 1;
                }

                prod = prod * inner_acc; // Multiply by the result for this variable
            }

            acc = acc + prod; // Add the product for this term
        }

        acc
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
          let field = v.1;
          pad.extend(std::iter::repeat(0).take(num_variables - k.len()));
          if let Some(existing_coeff) = dictionary.get_mut(&pad) {
              // If the exponent vector already exists, add the coefficients
              *existing_coeff = existing_coeff.add(*v);
          } else {
              // If the exponent vector does not exist, insert the new one
              dictionary.insert(pad, v.clone());
          }
      }

      MPolynomial::new(self.field, dictionary)
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
    let field = Field::new(5);
    dictionary.insert(vec![1, 2], FieldElement::new(3, field));
    dictionary.insert(vec![2,1], FieldElement::new(3, field));
    let p1 = MPolynomial::new(field, dictionary);
    let mut dictionary2 = HashMap::new();
    dictionary2.insert(vec![2,1], FieldElement::new(2, field));
    dictionary2.insert(vec![1, 2], FieldElement::new(3, field));
  
    let p2 = MPolynomial::new(field, dictionary2);
    let mut dictionary3 = HashMap::new();
    dictionary3.insert(vec![1, 2], FieldElement::new(6, field));
    dictionary3.insert(vec![2, 1], FieldElement::new(5, field));
    let p3 = MPolynomial::new(field, dictionary3);
    assert_eq!(p1 + p2, p3);
  }
  #[test]
  fn test_addition_pad(){
    let mut dictionary = HashMap::new();
    let field = Field::new(7);
    dictionary.insert(vec![1, 2,0,0], FieldElement::new(3, field));
    dictionary.insert(vec![2,1], FieldElement::new(3, field));
    let p1 = MPolynomial::new(field, dictionary);
    let mut dictionary2 = HashMap::new();
    dictionary2.insert(vec![2,1], FieldElement::new(2, field));
    dictionary2.insert(vec![1, 2], FieldElement::new(3, field));
  
    let p2 = MPolynomial::new(field, dictionary2);
    let mut dictionary3 = HashMap::new();
    dictionary3.insert(vec![1, 2,0,0], FieldElement::new(6, field));
    dictionary3.insert(vec![2, 1,0,0], FieldElement::new(5, field));
    let p3 = MPolynomial::new(field, dictionary3);
    assert_eq!(p1 + p2, p3);
  }
  #[test]
  fn test_eq(){
    let mut dictionary = HashMap::new();
    let field = Field::new(5);
    dictionary.insert(vec![1, 2,3], FieldElement::new(3, field));
    let p1 = MPolynomial::new(field, dictionary);
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![1, 2,3], FieldElement::new(3, field));
    let p2 = MPolynomial::new(field, dictionary);
    assert_eq!(p1, p2);
  }

  #[test]
  fn test_mpoly_evaluate() {
    let mut dictionary = HashMap::new();
    let field = Field::new(5);
    dictionary.insert(vec![1, 2,3], FieldElement::new(3, field)); // 3xy^2z^3
    let p = MPolynomial::new(Field::new(5), dictionary);
    let point = vec![FieldElement::new(2, field), FieldElement::new(1, field), FieldElement::new(2, field)];
    assert_eq!(p.evaluate(&point), FieldElement::new(3, field));
  }

  #[test]
  fn test_mpoly_degree() {
    let mut dictionary = HashMap::new();
    let field = Field::new(17);
    dictionary.insert(vec![1,2,3], FieldElement::new(3, field)); // 3xy^2z^3
    dictionary.insert(vec![2,0,0], FieldElement::new(3, field)); // 3x^2
    let p = MPolynomial::new(Field::new(17), dictionary);

    assert_eq!(p.degree(), FieldElement::new(6, field));
  }
}

