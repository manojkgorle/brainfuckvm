
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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

    pub fn variables(num_variables: usize, field: Field) -> Vec<Self> {
      let mut variables = Vec::new();

      for i in 0..num_variables {
          let mut exponent = vec![0; num_variables];
          exponent[i] = 1; // Set the ith variable's exponent to 1

          let mut dictionary = HashMap::new();
          dictionary.insert(exponent, FieldElement::one(field));

          variables.push(MPolynomial::new(field, dictionary));
      }

      variables
  }
  // def lift(polynomial, variable_index):
  //       if polynomial.is_zero():
  //           return MPolynomial({})
  //       field = polynomial.coefficients[0].field
  //       variables = MPolynomial.variables(variable_index+1, field)
  //       x = variables[-1]
  //       acc = MPolynomial({})
  //       for i in range(len(polynomial.coefficients)):
  //           acc = acc + \
  //               MPolynomial.constant(polynomial.coefficients[i]) * (x ^ i)
  //       return acc

  pub fn lift(){
    
  }


    #[allow(dead_code)]
    pub fn neg(&self)->Self{
        let mut dictionary = HashMap::new();
        for (k, v) in &self.dictionary {
            dictionary.insert(k.clone(), v.neg());
        }
        MPolynomial::new(self.field, dictionary)
    }
    pub fn constant(x:FieldElement)->Self{
        let mut dictionary = HashMap::new();
        dictionary.insert(vec![0], x);
        MPolynomial::new(x.1, dictionary)
    }
    pub fn one (field:Field)->Self{
        MPolynomial::constant(FieldElement::one(field))
    }
    pub fn pow(&self, exp:u128)->Self{
      // write a testcase for a poly whose coefficients are zer0
      if self.is_zero() {
        return MPolynomial::new(self.field,HashMap::new());
    }
    let mut acc = MPolynomial::one(self.field); // Start with the multiplicative identity
        let mut base = self.clone();      // Clone the current polynomial as the base
        let mut exp = exp;          // Copy the exponent for manipulation

        while exp > 0 {
            if exp % 2 == 1 {
                acc = acc.mul(base.clone()); // Multiply by base if the current bit is 1
            }
            base = base.clone().mul(base); // Square the base
            exp /= 2;               // Shift to the next bit
        }

        acc
    }
    pub fn is_zero(&self)->bool{
        self.dictionary.is_empty()


   
    }}


  
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
          // let field = v.1;
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
impl Mul for MPolynomial{
  type Output = Self;
  fn mul(self, other: MPolynomial) -> MPolynomial {
    let mut dictionary = HashMap::new();
    
    // Calculate max number of variables
    let num_variables = self.dictionary.keys()
        .chain(other.dictionary.keys())
        .map(|k| k.len())
        .max()
        .unwrap_or(0);
    
    // Iterate through each term in both polynomials
    for (k0, v0) in &self.dictionary {
        for (k1, v1) in &other.dictionary {
            // Initialize exponent vector
            let mut exponent = vec![0; num_variables];
            
            // Add exponents from first polynomial
            for (k, &exp) in k0.iter().enumerate() {
                exponent[k] += exp;
            }
            
            // Add exponents from second polynomial
            for (k, &exp) in k1.iter().enumerate() {
                exponent[k] += exp;
            }
            
            // Convert to tuple for key
            let exponent_key = exponent.clone().into_iter().collect::<Vec<_>>();
            
            // Update or insert the coefficient
            dictionary
                .entry(exponent_key)
                .and_modify(|coeff| *coeff += *v0 * *v1)
                .or_insert(*v0 * *v1);
        }
    }
    
    // MPolynomial {dictionary }
    MPolynomial::new(self.field, dictionary)
}}
impl MulAssign for MPolynomial{
  fn mul_assign(&mut self, other: Self) {
    let mut dictionary = HashMap::new();
    
    // Calculate max number of variables
    let num_variables = self.dictionary.keys()
        .chain(other.dictionary.keys())
        .map(|k| k.len())
        .max()
        .unwrap_or(0);
    
    // Iterate through each term in both polynomials
    for (k0, v0) in &self.dictionary {
        for (k1, v1) in &other.dictionary {
            // Initialize exponent vector
            let mut exponent = vec![0; num_variables];
            
            // Add exponents from first polynomial
            for (k, &exp) in k0.iter().enumerate() {
                exponent[k] += exp;
            }
            
            // Add exponents from second polynomial
            for (k, &exp) in k1.iter().enumerate() {
                exponent[k] += exp;
            }
            
            // Convert to tuple for key
            let exponent_key = exponent.clone().into_iter().collect::<Vec<_>>();
            
            // Update or insert the coefficient
            dictionary
                .entry(exponent_key)
                .and_modify(|coeff| *coeff += *v0 * *v1)
                .or_insert(*v0 * *v1);
        }
    }
    
    // Update the dictionary
    self.dictionary = dictionary;
}}

impl Sub for MPolynomial{
  type Output=Self;
  fn sub(self, other :MPolynomial)->Self{
    let other2=other.neg();
    self+other2
  }
}
impl SubAssign for MPolynomial{
  fn sub_assign(&mut self, other:Self){
    let other2=other.neg();
    *self+=other2;
  }
}


#[cfg(test)]
mod test_mpolynomial_operation {
  use super::*;

  #[test]
  fn test_addition() {
    let mut dictionary = HashMap::new();
    let field = Field::new(5);
    dictionary.insert(vec![1, 2], FieldElement::new(3, field));
    dictionary.insert(vec![2, 1], FieldElement::new(3, field));
    let p1 = MPolynomial::new(field, dictionary);
    let mut dictionary2 = HashMap::new();
    dictionary2.insert(vec![2, 1], FieldElement::new(2, field));
    dictionary2.insert(vec![1, 2], FieldElement::new(3, field));
  
    let p2 = MPolynomial::new(field, dictionary2);
    let mut dictionary3 = HashMap::new();
    dictionary3.insert(vec![1, 2], FieldElement::new(6, field));
    dictionary3.insert(vec![2, 1], FieldElement::new(5, field));
    let p3 = MPolynomial::new(field, dictionary3);
    assert_eq!(p1 + p2, p3);
  }

  #[test]
  fn test_addition_pad() {
    let mut dictionary = HashMap::new();
    let field = Field::new(7);
    dictionary.insert(vec![1, 2, 0, 0], FieldElement::new(3, field));
    dictionary.insert(vec![2, 1], FieldElement::new(3, field));
    let p1 = MPolynomial::new(field, dictionary);
    let mut dictionary2 = HashMap::new();
    dictionary2.insert(vec![2, 1], FieldElement::new(2, field));
    dictionary2.insert(vec![1, 2], FieldElement::new(3, field));
  
    let p2 = MPolynomial::new(field, dictionary2);
    let mut dictionary3 = HashMap::new();
    dictionary3.insert(vec![1, 2, 0, 0], FieldElement::new(6, field));
    dictionary3.insert(vec![2, 1, 0, 0], FieldElement::new(5, field));
    let p3 = MPolynomial::new(field, dictionary3);
    assert_eq!(p1 + p2, p3);
  }

  #[test]
  fn test_eq() {
    let mut dictionary = HashMap::new();
    let field = Field::new(5);
    dictionary.insert(vec![1, 2, 3], FieldElement::new(3, field));
    let p1 = MPolynomial::new(field, dictionary);
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![1, 2, 3], FieldElement::new(3, field));
    let p2 = MPolynomial::new(field, dictionary);
    assert_eq!(p1, p2);
  }

  #[test]
  fn test_mpoly_evaluate() {
    let mut dictionary = HashMap::new();
    let field = Field::new(5);
    dictionary.insert(vec![1, 2, 3], FieldElement::new(3, field)); // 3xy^2z^3
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

  #[test]
  fn test_mul() {
    let mut dictionary = HashMap::new();
    let field = Field::new(5);
    dictionary.insert(vec![1, 2, 3], FieldElement::new(3, field));
    let p1 = MPolynomial::new(field, dictionary);
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![1, 2, 3], FieldElement::new(3, field));
    let p2 = MPolynomial::new(field, dictionary);
    let mut dictionary = HashMap::new();
    dictionary.insert(vec![2, 4, 6], FieldElement::new(4, field));
    let p3 = MPolynomial::new(field, dictionary);
    assert_eq!(p1.mul(p2), p3);
  }
}


