use chrono::Local;

use crate::fields::{Field, FieldElement};
use crate::fields::{NEG_ORDER, P};
use rayon::prelude::*;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};
// How should we be interpolating polynomials?
// -> lagrange interpolation
// -> iNTT inverse Number Theoretic Transformation
// How should we be evaluating polynomials?
// -> Normal evaluation
// -> NTT Number Theoretic Transformation

// @todo time to redefine algebra?

/// a_0 + a_1 * x + a_2 * x^2 + a_3 * x^3 + ... + a_n * x^n
#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coefficients: Vec<FieldElement>,
    pub evaluation_form: Vec<(FieldElement, FieldElement)>,
}

impl Polynomial {
    pub fn new_from_coefficients(coefficients: Vec<FieldElement>) -> Self {
        Polynomial {
            coefficients,
            evaluation_form: Vec::new(),
        }
    }

    pub fn new_from_evaluation(evaluation_form: Vec<(FieldElement, FieldElement)>) -> Self {
        Polynomial {
            coefficients: Vec::new(),
            evaluation_form,
        }
    }

    pub fn evaluate(&self, x: FieldElement) -> FieldElement {
        let mut result = FieldElement::zero(Field::new(x.modulus()));

        for i in 0..self.coefficients.len() {
            result += self.coefficients[i] * x.pow(i as u128);
        }
        result
    }
    // new function added
    pub fn evaluate_domain(&self, x: Vec<FieldElement>) -> Vec<FieldElement> {
        let len = x.len();
        let mut result = Vec::with_capacity(len);
        for i in 0..len {
            result.push(self.evaluate(x[i]));
        }
        result
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    pub fn scalar_mul(&self, scalar: FieldElement) -> Self {
        let len = self.coefficients.len();
        let mut result = Vec::with_capacity(len);
        for i in 0..len {
            result.push(self.coefficients[i] * scalar);
        }
        Polynomial::new_from_coefficients(result)
    }

    pub fn is_all_zeros(&self) -> bool {
        for i in 0..self.coefficients.len() {
            if self.coefficients[i].0 != 0 {
                return false;
            }
        }
        true
    }

    pub fn scalar_div(&self, scalar: FieldElement) -> Self {
        let len = self.coefficients.len();
        let mut result = Vec::with_capacity(len);

        for i in 0..len {
            result.push(self.coefficients[i] / scalar);
        }
        // self.coefficients = result.clone();
        Polynomial::new_from_coefficients(result)
    }
    pub fn zero(field: Field) -> Self {
        Polynomial::new_from_coefficients(vec![FieldElement::zero(Field::new(field.0))])
    }

    pub fn q_div(self, poly2: Self) -> (Self, Self) {
        let field = Field::new(self.coefficients[0].modulus());
        let n = self.coefficients.len();
        let m = poly2.coefficients.len();
        let zero = FieldElement::zero(field);
        if m == 0 {
            return (self, Polynomial::zero(field));
        }
        if n == 0 {
            return (self, Polynomial::zero(field));
        }
        if n < m {
            return (Polynomial::new_from_coefficients(vec![zero]), self);
        }
        let mut q = Vec::with_capacity(n - m + 1);
        let mut poly1_coeff = self.clone().coefficients;
        let mut poly2_coeff = poly2.clone().coefficients;
        poly1_coeff.reverse();
        poly2_coeff.reverse();
        for i in 0..n - m + 1 {
            let mut other_coeff = poly2_coeff.clone();
            other_coeff.append(&mut vec![zero; n - m - i]);
            let q_temp = poly1_coeff[0] / other_coeff[0];
            let other_poly = Polynomial::new_from_coefficients(other_coeff).scalar_mul(q_temp);
            poly1_coeff = (Polynomial::new_from_coefficients(poly1_coeff) - other_poly.clone())
                .coefficients[1..]
                .to_vec();
            q.push(q_temp);
        }
        q.reverse();
        poly1_coeff.reverse();
        let poly1 = Polynomial::new_from_coefficients(poly1_coeff);
        if poly1.is_all_zeros() {
            (
                Polynomial::new_from_coefficients(q),
                Polynomial::new_from_coefficients(vec![zero]),
            )
        } else {
            (Polynomial::new_from_coefficients(q), poly1)
        }
    }

    pub fn compose(mut self, mono: FieldElement) -> Self {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] *= mono.pow(i as u128);
        }
        self
    }
    // new function added
    pub fn pow(self, exp: u128) -> Self {
        let mut res = Polynomial::new_from_coefficients(vec![FieldElement::one(Field::new(
            self.coefficients[0].modulus(),
        ))]);

        // Identity polynomial (constant 1)
        let mut base = self.clone(); // Clone the base polynomial
        let mut exp = exp;

        while exp > 0 {
            if exp % 2 == 1 {
                res *= base.clone(); // Multiply res by base polynomial if current bit is 1
            }

            base = base.clone() * base; // Square the base polynomial
            exp /= 2; // Right shift the exponent (equivalent to dividing by 2)
        }

        res
    }
    // new function added
    pub fn scale(&self, factor: u128) -> Self {
        let len = self.coefficients.len();
        let mut result = Vec::with_capacity(len);

        for i in 0..len {
            result.push(FieldElement::new(
                self.coefficients[i].0 * factor.pow(i as u32),
                self.coefficients[i].1,
            ));
        }
        Polynomial::new_from_coefficients(result)
    }
    // new function added
    pub fn collinearity(evaluation_form: Vec<(FieldElement, FieldElement)>) -> bool {
        let input: Vec<FieldElement> = evaluation_form.iter().map(|(x, _)| *x).collect();
        let output: Vec<FieldElement> = evaluation_form.iter().map(|(_, y)| *y).collect();
        let poly = interpolate_lagrange_polynomials(input.clone(), output.clone());
        poly.degree() == 1
    }
    pub fn constant(constant: FieldElement) -> Self {
        Polynomial::new_from_coefficients(vec![constant])
    }
}

impl Add for Polynomial {
    type Output = Polynomial;

    fn add(self, other: Polynomial) -> Polynomial {
        let mut i = 0;
        let field = Field::new(other.coefficients[0].modulus());
        let coeff = self.coefficients.iter().map(|x| x.0).collect::<Vec<u128>>();
        let coeff2 = other
            .coefficients
            .iter()
            .map(|x| x.0)
            .collect::<Vec<u128>>();
        let len = coeff.len().max(coeff2.len());
        let mut result = Vec::with_capacity(len);
        while i < coeff.len() && i < coeff2.len() {
            result.push(coeff[i] + coeff2[i]);
            i += 1;
        }

        while i < coeff.len() {
            result.push(coeff[i]);
            i += 1;
        }

        while i < coeff2.len() {
            result.push(coeff2[i]);
            i += 1;
        }
        let res = result
            .iter()
            .map(|x| FieldElement::new(*x, field))
            .collect::<Vec<FieldElement>>();
        Polynomial::new_from_coefficients(res)
    }
}

impl AddAssign for Polynomial {
    fn add_assign(&mut self, other: Polynomial) {
        let mut i = 0;
        let field = Field::new(self.coefficients[0].modulus());
        let coeff = self.coefficients.iter().map(|x| x.0).collect::<Vec<u128>>();
        let coeff2 = other
            .coefficients
            .iter()
            .map(|x| x.0)
            .collect::<Vec<u128>>();
        let len = coeff.len().max(coeff2.len());
        let mut result = Vec::with_capacity(len);
        while i < coeff.len() && i < coeff2.len() {
            result.push(coeff[i] + coeff2[i]);
            i += 1;
        }

        while i < coeff.len() {
            result.push(coeff[i]);
            i += 1;
        }

        while i < coeff2.len() {
            result.push(coeff2[i]);
            i += 1;
        }

        self.coefficients = result
            .iter()
            .map(|x| FieldElement::new(*x, field))
            .collect::<Vec<FieldElement>>();
    }
}

impl Sub for Polynomial {
    type Output = Polynomial;

    fn sub(self, other: Polynomial) -> Polynomial {
        let len = self.coefficients.len().max(other.coefficients.len());
        let mut result = Vec::with_capacity(len);
        let mut i = 0;

        while i < self.coefficients.len() && i < other.coefficients.len() {
            result.push(self.coefficients[i] - other.coefficients[i]);
            i += 1;
        }

        while i < self.coefficients.len() {
            result.push(self.coefficients[i]);
            i += 1;
        }

        while i < other.coefficients.len() {
            result.push(-other.coefficients[i]);
            i += 1;
        }

        Polynomial::new_from_coefficients(result)
    }
}

impl SubAssign for Polynomial {
    fn sub_assign(&mut self, other: Polynomial) {
        let len = self.coefficients.len().max(other.coefficients.len());
        let mut result = Vec::with_capacity(len);
        let mut i = 0;

        while i < self.coefficients.len() && i < other.coefficients.len() {
            result.push(self.coefficients[i] - other.coefficients[i]);
            i += 1;
        }

        while i < self.coefficients.len() {
            result.push(self.coefficients[i]);
            i += 1;
        }

        while i < other.coefficients.len() {
            result.push(-other.coefficients[i]);
            i += 1;
        }

        self.coefficients = result;
    }
}

impl Mul for Polynomial {
    type Output = Polynomial;

    fn mul(self, other: Polynomial) -> Polynomial {
        let field = Field::new(self.coefficients[0].modulus());
        let mut result = vec![0; self.coefficients.len() + other.coefficients.len() - 1];
        let coeff = self.coefficients.iter().map(|x| x.0).collect::<Vec<u128>>();
        let coeff2 = other
            .coefficients
            .iter()
            .map(|x| x.0)
            .collect::<Vec<u128>>();
        for i in 0..coeff.len() {
            for j in 0..coeff2.len() {
                result[i + j] += (coeff[i] * coeff2[j]) % field.0;
            }
        }
        let res = result
            .iter()
            .map(|x| FieldElement::new(*x, field))
            .collect::<Vec<FieldElement>>();
        Polynomial::new_from_coefficients(res)
    }
}

impl MulAssign for Polynomial {
    fn mul_assign(&mut self, other: Polynomial) {
        if self.coefficients.len() > 0 {
            let field = Field::new(self.coefficients[0].modulus());
            let mut result = vec![0; self.coefficients.len() + other.coefficients.len() - 1];
            let coeff = self.coefficients.iter().map(|x| x.0).collect::<Vec<u128>>();
            let coeff2 = other
                .coefficients
                .iter()
                .map(|x| x.0)
                .collect::<Vec<u128>>();
            for i in 0..coeff.len() {
                for j in 0..coeff2.len() {
                    result[i + j] += (coeff[i] * coeff2[j]) % field.0;
                }
            }
            self.coefficients = result
                .iter()
                .map(|x| FieldElement::new(*x, field))
                .collect::<Vec<FieldElement>>();
        }
    }
}

impl Div for Polynomial {
    type Output = Polynomial;

    fn div(self, other: Polynomial) -> Polynomial {
        let (q, r) = self.q_div(other);
        if r.is_all_zeros() {
            q
        } else {
            panic!("Division not possible");
        }
    }
}

impl DivAssign for Polynomial {
    fn div_assign(&mut self, other: Polynomial) {
        let (q, r) = self.clone().q_div(other);
        if r.is_all_zeros() {
            self.coefficients = q.coefficients;
        } else {
            panic!("Division not possible");
        }
    }
}
// zerofier domain
pub fn gen_polynomial_from_roots(roots: Vec<FieldElement>) -> Polynomial {
    let field = Field::new(roots[0].modulus());
    let mut result = vec![FieldElement::new(1, field); 1];
    for i in 0..roots.len() {
        let mut new_result = vec![FieldElement::new(0, field); 2];
        new_result[0] = -roots[i];
        new_result[1] = FieldElement::new(1, field);
        let new_polynomial = Polynomial::new_from_coefficients(new_result);
        result = (Polynomial::new_from_coefficients(result) * new_polynomial).coefficients;
    }
    Polynomial::new_from_coefficients(result[..roots.len() + 1].to_vec())
}
//interpolate_domain

pub fn gen_lagrange_polynomials(x: Vec<FieldElement>) -> Vec<Polynomial> {
    let n = x.len();
    let mut lagrange_polynomials = Vec::with_capacity(n);
    let one = FieldElement::one(x[0].1);
    for i in 0..n {
        let mut denominator = Vec::with_capacity(n);

        let roots = &x[..i];

        let numerator = gen_polynomial_from_roots([roots, &x[i + 1..]].concat());
        for j in 0..n {
            if i == j {
                continue;
            }
            denominator.push(x[i] - x[j]);
        }
        let den_sum = denominator.iter().fold(one, |acc, x| acc * *x);
        let lagrange_polynomial = numerator.scalar_div(den_sum);
        lagrange_polynomials.push(lagrange_polynomial);
    }

    lagrange_polynomials
}

pub fn gen_lagrange_polynomials_parallel(x: Vec<FieldElement>) -> Vec<Polynomial> {
    let n = x.len();
    let lagrange_polynomials = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut denominator = Vec::with_capacity(n);

            let roots = &x[..i];

            let numerator: Polynomial = gen_polynomial_from_roots([roots, &x[i + 1..]].concat());

            for j in 0..n {
                if i == j {
                    continue;
                }
                denominator.push(x[i] - x[j]);
            }

            let den_sum = denominator.iter().fold(
                FieldElement::new(1, Field::new(x[i].modulus())),
                |acc, x| acc * *x,
            );

            numerator.scalar_div(den_sum)
        })
        .collect();

    lagrange_polynomials
}

pub fn interpolate_lagrange_polynomials(x: Vec<FieldElement>, y: Vec<FieldElement>) -> Polynomial {
    let n = x.len();

    let lagrange_polynomials = gen_lagrange_polynomials_parallel(x.clone());
    let field = Field::new(x[0].modulus());
    let mut result = Polynomial::new_from_coefficients(vec![FieldElement::new(0, field); n]);

    for i in 0..n {
        result += lagrange_polynomials[i].scalar_mul(y[i]);
    }

    result
}

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        if self.coefficients.len() != other.coefficients.len() {
            return false;
        }
        for i in 0..self.coefficients.len() {
            if self.coefficients[i] != other.coefficients[i] {
                return false;
            }
        }
        true
    }
}

// fn karatsuba_multiply(a: &[u128], b: &[u128]) -> Vec<u128> {
//     let n = a.len();
//     let m = b.len();

//     // Base case: if the polynomials are small, use the naive multiplication
//     if n < 32 || m < 32 {
//         return naive_multiply(a, b);
//     }

//     // Ensure both polynomials have the same length by padding with zeros
//     let max_len = n.max(m);
//     let mut a_padded = vec![0; max_len];
//     let mut b_padded = vec![0; max_len];
//     a_padded[..n].copy_from_slice(a);
//     b_padded[..m].copy_from_slice(b);

//     // Split the polynomials into two halves
//     let mid = max_len / 2;
//     let (a_low, a_high) = a_padded.split_at(mid);
//     let (b_low, b_high) = b_padded.split_at(mid);

//     // Recursive calls to Karatsuba
//     let z0 = karatsuba_multiply(a_low, b_low);
//     let z1 = karatsuba_multiply(&vec_add(a_low, a_high), &vec_add(b_low, b_high));
//     let z2 = karatsuba_multiply(a_high, b_high);

//     // Combine the results
//     let mut result = vec![0; 2 * max_len - 1];
//     for (i, &val) in z0.iter().enumerate() {
//         result[i] += val;
//     }
//     for (i, &val) in z1.iter().enumerate() {
//         result[i + mid] += val - z0[i] - z2[i];
//     }
//     for (i, &val) in z2.iter().enumerate() {
//         result[i + 2 * mid] += val;
//     }

//     result
// }

// // Helper function to add two vectors element-wise
// fn vec_add(a: &[u128], b: &[u128]) -> Vec<u128> {
//     a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
// }

// // Naive polynomial multiplication (for small inputs)
// fn naive_multiply(a: &[u128], b: &[u128]) -> Vec<u128> {
//     let n = a.len();
//     let m = b.len();
//     let mut result = vec![0; n + m - 1];

//     for i in 0..n {
//         for j in 0..m {
//             result[i + j] += a[i] * b[j];
//         }
//     }

//     result
// }

#[cfg(test)]
mod test_polynomials {
    use super::*;

    #[test]
    fn test_evaluate() {
        let field = Field::new(7);
        let coefficients = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let polynomial = Polynomial::new_from_coefficients(coefficients);
        let x = FieldElement::new(4, field);
        let result = polynomial.evaluate(x);
        assert_eq!(result.0, 2);
    }

    #[test]
    fn test_add() {
        let field = Field::new(7);
        let coefficients1 = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let coefficients2 = vec![FieldElement::new(3, field), FieldElement::new(4, field)];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let result = polynomial1 + polynomial2;
        assert_eq!(result.coefficients[0].0, 4);
        assert_eq!(result.coefficients[1].0, 6);
    }

    #[test]
    fn test_add_poly2_larger() {
        let field = Field::new(7);
        let coefficients1 = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let coefficients2 = vec![
            FieldElement::new(3, field),
            FieldElement::new(4, field),
            FieldElement::new(5, field),
        ];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let result = polynomial1 + polynomial2;
        assert_eq!(result.coefficients[0].0, 4);
        assert_eq!(result.coefficients[1].0, 6);
        assert_eq!(result.coefficients[2].0, 5);
    }

    #[test]
    fn test_add_poly1_larger() {
        let field = Field::new(7);
        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(4, field),
            FieldElement::new(3, field),
        ];
        let coefficients2 = vec![FieldElement::new(3, field), FieldElement::new(4, field)];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let result = polynomial1 + polynomial2;
        assert_eq!(result.coefficients[0].0, 4);
        assert_eq!(result.coefficients[1].0, 1);
        assert_eq!(result.coefficients[2].0, 3);
    }

    #[test]
    fn test_sub() {
        let field = Field::new(7);
        let coefficients1 = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let coefficients2 = vec![FieldElement::new(3, field), FieldElement::new(4, field)];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let result = polynomial1 - polynomial2;
        assert_eq!(result.coefficients[0].0, 5);
        assert_eq!(result.coefficients[1].0, 5);
    }

    #[test]
    fn test_sub_poly2_larger() {
        let field = Field::new(7);
        let coefficients1 = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let coefficients2 = vec![
            FieldElement::new(3, field),
            FieldElement::new(4, field),
            FieldElement::new(5, field),
        ];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let result = polynomial1 - polynomial2;
        assert_eq!(result.coefficients[0].0, 5);
        assert_eq!(result.coefficients[1].0, 5);
        assert_eq!(result.coefficients[2].0, 2);
    }

    #[test]
    fn test_sub_poly1_larger() {
        let field = Field::new(7);
        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(4, field),
            FieldElement::new(3, field),
        ];
        let coefficients2 = vec![FieldElement::new(3, field), FieldElement::new(4, field)];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let result = polynomial1 - polynomial2;
        assert_eq!(result.coefficients[0].0, 5);
        assert_eq!(result.coefficients[1].0, 0);
        assert_eq!(result.coefficients[2].0, 3);
    }

    #[test]
    fn test_mul() {
        let field = Field::new(7);
        let coefficients1 = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let coefficients2 = vec![FieldElement::new(3, field), FieldElement::new(4, field)];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let result = polynomial1 * polynomial2;
        assert_eq!(result.coefficients[0].0, 3);
        assert_eq!(result.coefficients[1].0, 3);
        assert_eq!(result.coefficients[2].0, 1);
    }

    #[test]
    fn test_scalar_mul() {
        let field = Field::new(7);
        let coefficients = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let polynomial = Polynomial::new_from_coefficients(coefficients);
        let scalar = FieldElement::new(3, field);
        let result = polynomial.scalar_mul(scalar);
        assert_eq!(result.coefficients[0].0, 3);
        assert_eq!(result.coefficients[1].0, 6);
    }

    #[test]
    fn test_scalar_div() {
        let field = Field::new(7);
        let coefficients = vec![FieldElement::new(1, field), FieldElement::new(2, field)];
        let polynomial = Polynomial::new_from_coefficients(coefficients);
        let scalar = FieldElement::new(4, field);
        let result = polynomial.scalar_div(scalar);
        assert_eq!(result.coefficients[0].0, 2);
        assert_eq!(result.coefficients[1].0, 4);
    }

    #[test]
    fn test_polynomial_from_roots() {
        let field = Field::new(7);
        let roots = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];
        let polynomial = gen_polynomial_from_roots(roots);
        assert_eq!(polynomial.coefficients[0].0, 1);
        assert_eq!(polynomial.coefficients[1].0, 4);
        assert_eq!(polynomial.coefficients[2].0, 1);
        assert_eq!(polynomial.coefficients[3].0, 1);
    }

    #[test]
    fn test_gen_lagrange_poly() {
        let field = Field::new(7);
        let x = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(5, field),
        ];
        let lagrange_polynomials = gen_lagrange_polynomials(x);
        assert_eq!(lagrange_polynomials.len(), 3);
        assert_eq!(lagrange_polynomials[0].coefficients[0].0, 5);
        assert_eq!(lagrange_polynomials[0].coefficients[1].0, 2);
        assert_eq!(lagrange_polynomials[0].coefficients[2].0, 5);
        assert_eq!(lagrange_polynomials[1].coefficients[0].0, 2);
        assert_eq!(lagrange_polynomials[1].coefficients[1].0, 0);
        assert_eq!(lagrange_polynomials[1].coefficients[2].0, 3);
        assert_eq!(lagrange_polynomials[2].coefficients[0].0, 1);
        assert_eq!(lagrange_polynomials[2].coefficients[1].0, 5);
        assert_eq!(lagrange_polynomials[2].coefficients[2].0, 6);
    }

    #[test]
    fn test_gen_lagrange_poly2() {
        let field = Field::new(7);
        let x = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(5, field),
            FieldElement::new(6, field),
        ];
        let lagrange_polynomials = gen_lagrange_polynomials(x);
        assert_eq!(lagrange_polynomials.len(), 4);
        assert_eq!(lagrange_polynomials[0].coefficients[0].0, 4);
        assert_eq!(lagrange_polynomials[0].coefficients[1].0, 0);
        assert_eq!(lagrange_polynomials[0].coefficients[2].0, 0);
        assert_eq!(lagrange_polynomials[0].coefficients[3].0, 4);
        assert_eq!(lagrange_polynomials[1].coefficients[0].0, 4);
        assert_eq!(lagrange_polynomials[1].coefficients[1].0, 4);
        assert_eq!(lagrange_polynomials[1].coefficients[2].0, 6);
        assert_eq!(lagrange_polynomials[1].coefficients[3].0, 6);
        assert_eq!(lagrange_polynomials[2].coefficients[0].0, 6);
        assert_eq!(lagrange_polynomials[2].coefficients[1].0, 1);
        assert_eq!(lagrange_polynomials[2].coefficients[2].0, 3);
        assert_eq!(lagrange_polynomials[2].coefficients[3].0, 1);
        assert_eq!(lagrange_polynomials[3].coefficients[0].0, 1);
        assert_eq!(lagrange_polynomials[3].coefficients[1].0, 2);
        assert_eq!(lagrange_polynomials[3].coefficients[2].0, 5);
        assert_eq!(lagrange_polynomials[3].coefficients[3].0, 3);
    }

    #[test]
    fn test_interpolate_lagrange_polynomials() {
        let field = Field::new(7);
        let x = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(5, field),
        ];
        let y = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];
        let result = interpolate_lagrange_polynomials(x, y);
        println!("{:?}", result.coefficients);
        assert_eq!(result.coefficients[0].0, 5);
        assert_eq!(result.coefficients[1].0, 3);
        assert_eq!(result.coefficients[2].0, 1);
    }

    #[test]
    fn test_qdiv() {
        let field = Field::new(7);
        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
        ];
        let coefficients2 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];
        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let (q, r) = polynomial1.q_div(polynomial2);
        assert_eq!(q.coefficients.len(), 2);
        assert_eq!(r.coefficients.len(), 2);
        assert_eq!(q.coefficients[0].0, 4);
        assert_eq!(q.coefficients[1].0, 6);
        assert_eq!(r.coefficients[0].0, 4);
        assert_eq!(r.coefficients[1].0, 2);
    }

    #[test]
    fn test_qdiv2() {
        let field = Field::new(7);

        let coefficients1 = vec![
            FieldElement::new(2, field),
            FieldElement::new(1, field),
            FieldElement::new(3, field),
            FieldElement::new(0, field),
            FieldElement::new(6, field),
        ];

        let coefficients2 = vec![
            FieldElement::new(3, field),
            FieldElement::new(5, field),
            FieldElement::new(0, field),
            FieldElement::new(4, field),
        ];

        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);
        let (q, r) = polynomial1.q_div(polynomial2);

        assert_eq!(r.coefficients.len(), 3);
        assert_eq!(q.coefficients.len(), 2);
        assert_eq!(r.coefficients[0].0, 2);
        assert_eq!(r.coefficients[1].0, 0);
        assert_eq!(r.coefficients[2].0, 6);
        assert_eq!(q.coefficients[0].0, 0);
        assert_eq!(q.coefficients[1].0, 5);
    }

    #[test]
    fn test_qdiv_same_n_m() {
        let field = Field::new(7);
        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];
        let coefficients2 = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
        ];

        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);

        let (q, r) = polynomial1.q_div(polynomial2);
        assert_eq!(q.coefficients.len(), 1);
        assert_eq!(r.coefficients.len(), 2);
        assert_eq!(q.coefficients[0].0, 6);
        assert_eq!(r.coefficients[0].0, 3);
        assert_eq!(r.coefficients[1].0, 5);
    }

    #[test]
    fn test_qdiv_m_larger_than_n() {
        let field = Field::new(7);
        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];
        let coefficients2 = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
            FieldElement::new(5, field),
        ];

        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);

        let (q, r) = polynomial1.q_div(polynomial2);
        assert_eq!(q.coefficients.len(), 1);
        assert_eq!(r.coefficients.len(), 3);
        assert_eq!(q.coefficients[0].0, 0);
        assert_eq!(r.coefficients[0].0, 1);
        assert_eq!(r.coefficients[1].0, 2);
        assert_eq!(r.coefficients[2].0, 3);
    }

    #[test]
    fn test_qdiv_no_remainder() {
        let field = Field::new(7);

        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let coefficients2 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);

        let (q, r) = polynomial1.q_div(polynomial2);

        assert_eq!(q.coefficients.len(), 1);
        assert_eq!(r.coefficients.len(), 1);
        assert_eq!(q.coefficients[0].0, 1);
        assert_eq!(r.coefficients[0].0, 0);
    }

    #[test]
    fn test_div() {
        let field = Field::new(7);

        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let coefficients2 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);

        let result = polynomial1 / polynomial2;

        assert_eq!(result.coefficients.len(), 1);
        assert_eq!(result.coefficients[0].0, 1);
    }

    #[test]
    #[should_panic]
    fn test_div_not_possible() {
        let field = Field::new(7);

        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let coefficients2 = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
        ];

        let polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);

        let _ = polynomial1 / polynomial2;
    }

    #[test]
    fn test_div_assign() {
        let field = Field::new(7);

        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let coefficients2 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let mut polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);

        polynomial1 /= polynomial2;

        assert_eq!(polynomial1.coefficients.len(), 1);
        assert_eq!(polynomial1.coefficients[0].0, 1);
    }

    #[test]
    #[should_panic]
    fn test_div_assign_not_possible() {
        let field = Field::new(7);

        let coefficients1 = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let coefficients2 = vec![
            FieldElement::new(2, field),
            FieldElement::new(3, field),
            FieldElement::new(4, field),
        ];

        let mut polynomial1 = Polynomial::new_from_coefficients(coefficients1);
        let polynomial2 = Polynomial::new_from_coefficients(coefficients2);

        polynomial1 /= polynomial2;
    }

    #[test]
    fn test_mono_compose() {
        let field = Field::new(7);

        let coefficients = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];

        let polynomial = Polynomial::new_from_coefficients(coefficients);

        let mono = FieldElement::new(4, field);

        let polynomial = polynomial.compose(mono);

        assert_eq!(polynomial.coefficients[0].0, 1);
        assert_eq!(polynomial.coefficients[1].0, 1);
        assert_eq!(polynomial.coefficients[2].0, 6);
    }
    #[test]
    fn test_pow() {
        let field = Field::new(7);
        let coefficients = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];
        let polynomial = Polynomial::new_from_coefficients(coefficients);
        let result = polynomial.pow(2);
        assert_eq!(result.coefficients[0].0, 1);
        assert_eq!(result.coefficients[1].0, 4);
        assert_eq!(result.coefficients[2].0, 3);
        assert_eq!(result.coefficients[3].0, 5);
        assert_eq!(result.coefficients[4].0, 2);
    }
    #[test]
    fn test_scale() {
        let field = Field::new(7);
        let coefficients = vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ];
        let polynomial = Polynomial::new_from_coefficients(coefficients);
        let result = polynomial.scale(2);
        assert_eq!(result.coefficients[0].0, 1);
        assert_eq!(result.coefficients[1].0, 4);
        assert_eq!(result.coefficients[2].0, 5);
    }
    #[test]
    fn test_collinearity() {
        let field = Field::new(7);
        let evaluation_form = vec![
            (FieldElement::new(1, field), FieldElement::new(1, field)),
            (FieldElement::new(2, field), FieldElement::new(2, field)),
            (FieldElement::new(3, field), FieldElement::new(3, field)),
        ];
        let result = Polynomial::collinearity(evaluation_form);
        assert!(!result);
    }
    #[test]
    fn test_collinearity2() {
        let field = Field::new(7);
        let evaluation_form = vec![
            (FieldElement::new(1, field), FieldElement::new(1, field)),
            (FieldElement::new(2, field), FieldElement::new(2, field)),
        ];
        let result = Polynomial::collinearity(evaluation_form);
        assert!(result);
    }

    // #[test]
    // fn bench_mul() {
    //     let field = Field::new(P as u128);
    //     let a = vec![1, 2, 3, 4, 5,6, 7, 8, 0, 1, 2, 3, 4, 5,6, 7, 8, 0, 1, 2, 3, 4, 5,6, 7, 8, 0,1, 2, 3, 4, 5,6, 7, 8, 9];
    //     let b  = vec![1, 2, 3, 4, 5,6, 7, 8, 0, 1, 2, 3, 4, 5,6, 7, 8, 0, 1, 2, 3, 4, 5,6, 7, 8, 0,1, 2, 3, 4, 5,6, 7, 8, 9];
    //     let t = Local::now();
    //     let poly = Polynomial::new_from_coefficients(a.iter().map(|x| FieldElement::new(*x, field)).collect());
    //     let poly2 = Polynomial::new_from_coefficients(b.iter().map(|x| FieldElement::new(*x, field)).collect());
    //     for _ in 0..100 {
    //         karatsuba_multiply(&a, &b);
    //     }
    //     println!("Time: {:?}", Local::now().signed_duration_since(t));
    //     let t = Local::now();
    //     for _ in 0..100 {
    //        let _ = poly.clone() * poly2.clone();
    //     }
    //     println!("Time: {:?}", Local::now().signed_duration_since(t));

    // }
}
