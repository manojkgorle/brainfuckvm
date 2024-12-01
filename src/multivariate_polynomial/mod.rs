use crate::fields::{Field, FieldElement};
use crate::univariate_polynomial::*;
use std::collections::HashMap;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
#[derive(Clone, Debug)]
struct MPolynomial {
    field: Field,
    dictionary: HashMap<Vec<u128>, FieldElement>, // Exponent vector -> coefficient
}

impl MPolynomial {
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
    pub fn degree(&self) -> i128 {
        if self.dictionary.is_empty() {
            return -1; // Return zero for an empty polynomial
        }

        let max_degree = self
            .dictionary
            .iter() // Iterate over the key-value pairs in the dictionary
            .filter(|(_, coefficient)| coefficient.0 != 0 as u128) // Filter out terms with zero coefficients
            .map(|(exponents, _)| exponents.iter().sum::<u128>()) // Sum the exponents for non-zero terms
            .max() // Find the maximum sum of exponents
            .unwrap_or(0); // Default to 0 if no non-zero terms

        max_degree as i128 // Convert the maximum degree to a FieldElement
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

    //Function to evaluate at a point
    #[allow(dead_code)]
    pub fn evaluate(&self, point: &Vec<FieldElement>) -> FieldElement {
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
        let mut acc = MPolynomial::zero(field); // Accumulator for the result

        for (exp, coeff) in &self.dictionary {
            let c = FieldElement::new(1, field);
            let mut prod = MPolynomial::constant(c); // Start with the coefficient as a constant

            for (i, &exp) in exp.iter().enumerate() {
                let mut inner_acc = MPolynomial::one(field); // Initialize to one

                let mut j = 0;
                while (1 << j) <= exp {
                    if !memo.contains_key(&(i as u128, 1 << j)) {
                        if j == 0 {
                            memo.insert((i as u128, 1 << j), point[i].clone());
                        } else {
                            let prev = memo.get(&(i as u128, 1 << (j - 1))).unwrap().clone();
                            memo.insert((i as u128, 1 << j), prev.clone() * prev);
                        }
                    }
                    let point_power = memo.get(&(i as u128, 1 << j)).unwrap();

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


    pub fn symbolic_degree_bound(&self, max_degrees: Vec<u128>) -> i128 {
        // Check if the polynomial is empty
        let field = self.field;
        if self.dictionary.is_empty() {
            return -1;
        }
        // Ensure the max_degrees vector is correctly sized
        let num_variables = self.dictionary.keys().next().unwrap_or(&vec![]).len();
        if max_degrees.len() < num_variables {
            panic!("max_degrees length does not match the number of variables in the polynomial.");
        }

        // Make sure all degrees in max_degrees are the same
        if !max_degrees.iter().all(|&degree| degree == max_degrees[0]) {
            panic!("max_degrees must contain the same integer repeated.");
        }

        let mut total_degree_bound = 0;

        for (exp, coeff) in &self.dictionary {
            if coeff.0 == 0 as u128 {
                continue;
            }

            let mut term_degree_bound = 0;
            for (exp, max_degree) in exp.iter().zip(max_degrees.iter()) {
                term_degree_bound += (*exp as u128) * (*max_degree);
            }

            total_degree_bound = total_degree_bound.max(term_degree_bound as i128);
        }

        total_degree_bound
    }

    pub fn partial_evaluate(&self, partial_assignment: HashMap<usize, u128>) -> MPolynomial {
        let field = self.field;
        let num_variables = self.dictionary.keys().len();
        let mut complete_assignment = MPolynomial::variables(num_variables, field);
        
        for (index, value) in partial_assignment {
            complete_assignment[index] = MPolynomial::constant(FieldElement::new(value, field));
        }

        let mut polynomial = MPolynomial::zero(field);
        for (key, value) in &self.dictionary {
            let mut term = MPolynomial::constant(*value);
            for (i, &exp) in key.iter().enumerate() {
                term = term.mul(complete_assignment[i].pow(exp));
            }
            polynomial = polynomial.add(term);
        }

        polynomial
    }


    

  pub fn lift(polynomial:&Polynomial, variable_index:usize)->MPolynomial{
    let field = Field::new(polynomial.coefficients[0].modulus());
    if polynomial.is_all_zeros(){
      return MPolynomial::zero(field);
    }
    
    let variables = MPolynomial::variables(variable_index+1, field);
    let x = variables[variable_index].clone();
    let mut acc = MPolynomial::zero(field);
    for i in 0..polynomial.coefficients.len(){
      acc = acc + MPolynomial::constant(polynomial.coefficients[i].clone()) * x.pow(i as u128);
    }
    acc

    
  }


    #[allow(dead_code)]
    pub fn neg(&self) -> Self {
        let mut dictionary = HashMap::new();
        for (k, v) in &self.dictionary {
            dictionary.insert(k.clone(), v.neg());
        }
        MPolynomial::new(self.field, dictionary)
    }
    pub fn constant(x: FieldElement) -> Self {
        let mut dictionary = HashMap::new();
        dictionary.insert(vec![0], x);
        MPolynomial::new(x.1, dictionary)
    }
    pub fn one(field: Field) -> Self {
        MPolynomial::constant(FieldElement::one(field))
    }
    pub fn pow(&self, exp: u128) -> Self {
        // write a testcase for a poly whose coefficients are zer0
        if self.is_zero() {
            return MPolynomial::new(self.field, HashMap::new());
        }
        let mut acc = MPolynomial::one(self.field); // Start with the multiplicative identity
        let mut base = self.clone(); // Clone the current polynomial as the base
        let mut exp = exp; // Copy the exponent for manipulation

        while exp > 0 {
            if exp % 2 == 1 {
                acc = acc.mul(base.clone()); // Multiply by base if the current bit is 1
            }
            base = base.clone().mul(base); // Square the base
            exp /= 2; // Shift to the next bit
        }

        acc
    }
    pub fn is_zero(&self) -> bool {
        self.dictionary.is_empty()
    }
    pub fn str(&self)->String{
        let mut terms = Vec::new();
        for (k, v) in &self.dictionary {
            let mut term = String::new();
            for (i, &exp) in k.iter().enumerate() {
                if exp == 0 {
                    continue;
                }
                if term.len() > 0 {
                    term.push_str("*");
                }
                term.push_str(&format!("x_{}^{}", i, exp));
            }
            term.push_str(&format!("*{:?}", v.0));
            terms.push(term);
        }
        terms.join(" + ")
    }
}

impl Add for MPolynomial {
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
    }
}

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
    }
}

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
impl Mul for MPolynomial {
    type Output = Self;
    fn mul(self, other: MPolynomial) -> MPolynomial {
        let mut dictionary = HashMap::new();

        // Calculate max number of variables
        let num_variables = self
            .dictionary
            .keys()
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
    }
}
impl MulAssign for MPolynomial {
    fn mul_assign(&mut self, other: Self) {
        let mut dictionary = HashMap::new();

        // Calculate max number of variables
        let num_variables = self
            .dictionary
            .keys()
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
    }
}

impl Sub for MPolynomial {
    type Output = Self;
    fn sub(self, other: MPolynomial) -> Self {
        let other2 = other.neg();
        self + other2
    }
}
impl SubAssign for MPolynomial {
    fn sub_assign(&mut self, other: Self) {
        let other2 = other.neg();
        *self += other2;
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
        let point = vec![
            FieldElement::new(2, field),
            FieldElement::new(1, field),
            FieldElement::new(2, field),
        ];
        assert_eq!(p.evaluate(&point), FieldElement::new(3, field));
    }

    #[test]
    fn test_mpoly_degree() {
        let mut dictionary = HashMap::new();
        let field = Field::new(17);
        dictionary.insert(vec![1, 2, 3], FieldElement::new(3, field)); // 3xy^2z^3
        dictionary.insert(vec![2, 0, 0], FieldElement::new(3, field)); // 3x^2
        let p = MPolynomial::new(Field::new(17), dictionary);

        assert_eq!(p.degree(), 6 as i128);
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

    #[test]
    fn test_mpoly_variables() {
        let field = Field::new(5);
        let variables = MPolynomial::variables(4, field);
        let mut dictionary = HashMap::new();
        dictionary.insert(vec![0, 0, 1, 0], FieldElement::new(3, field));
        let p1 = MPolynomial::new(field, dictionary);
        assert_eq!(variables[2], p1);
    }

    #[test]
    fn test_poly_eval_symbolic() {
        let field = Field::new(17);
        // Define variables x, y, z, t
        let variables = MPolynomial::variables(4, field);
        let z = variables[2].clone();
        let t = variables[3].clone();

        // Define polynomial 3x^2y + 2x
        let mut dictionary = HashMap::new();
        dictionary.insert(vec![2, 1], FieldElement::new(3, field)); // 3x^2y
        dictionary.insert(vec![1, 0], FieldElement::new(2, field)); // 2x
        let poly = MPolynomial::new(field, dictionary);

        let mut dictionary = HashMap::new();
        dictionary.insert(vec![0, 0, 2, 1], FieldElement::new(3, field)); // 3z^2t
        dictionary.insert(vec![0, 0, 1, 0], FieldElement::new(2, field)); // 2z
        let poly2 = MPolynomial::new(field, dictionary);

        // Memoization cache
        let mut memo = HashMap::new();

        // Evaluate symbolically with x, y
        let result = poly.evaluate_symbolic(&vec![z, t], &mut memo);

        // Expected result: 3x^2y + 2x
        assert_eq!(result, poly2);
        // Check if symbolic result matches
    }

    #[test]
    fn test_symbolic_degree_bound() {
        let field = Field::new(17);
        let mut dictionary = HashMap::new();
        dictionary.insert(vec![1, 2, 0], FieldElement::new(3, field)); // term 3x^1 * y^2
        dictionary.insert(vec![0, 0, 3], FieldElement::new(3, field)); // term 3z^3
        let p = MPolynomial::new(field, dictionary);

    // Test with max_degrees as [3, 3, 3] for each variable
    let max_degrees = vec![4, 4, 4];
    let degree_bound = p.symbolic_degree_bound(max_degrees);
    assert_eq!(degree_bound, 12); 
}

#[test]
    fn test_partial_evaluate() {
        let field = Field::new(1000000007);

        // Create a polynomial: 3x^2y + 5y^2z + 7
        let mut dictionary = HashMap::new();
        dictionary.insert(vec![2, 1, 0], FieldElement::new(3, field)); // 3x^2y
        dictionary.insert(vec![0, 2, 1], FieldElement::new(5, field)); // 5y^2z
        dictionary.insert(vec![0, 0, 0], FieldElement::new(7, field)); // 7 (constant)

        let polynomial = MPolynomial::new(field, dictionary);

        // Define the partial assignment: x = 2, z = 4
        let mut partial_assignment = HashMap::new();
        partial_assignment.insert(1, 2); // y = 2
        partial_assignment.insert(2, 4); // z = 44

        // Evaluate the polynomial
        let evaluated = polynomial.partial_evaluate(partial_assignment);


        // Define expected results as symbolic terms
        let mut expected_dict = HashMap::new();
        expected_dict.insert(vec![1, 0,0], FieldElement::new(6, field)); // 12x
        expected_dict.insert(vec![0, 0,0], FieldElement::new(80, field)); // 80
        expected_dict.insert(vec![0, 0,0], FieldElement::new(7, field));  // 7 (constant)
        let expected_polynomial = MPolynomial::new(field, expected_dict);

        println!("{}", evaluated.str());
        println!("{}", expected_polynomial.str());
        //assert_eq!(evaluated.dictionary, expected_polynomial.dictionary);
    }

}
