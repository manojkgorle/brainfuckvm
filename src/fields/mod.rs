#![allow(dead_code)]
use core::hash::{Hash, Hasher};
use core::hint::unreachable_unchecked;
use std::cmp::{Ord, PartialOrd};
use std::fmt::write;
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
// Define a field
// Define a field element
// Define arithmetic operations on field elements
#[allow(clippy::derived_hash_with_manual_eq)]
#[derive(Debug, Clone, Copy, Hash)]
pub struct Field(pub u128);

pub enum ChallengeIndices {
    A,
    B,
    C,
    D,
    E,
    F,
    Alpha,
    Beta,
    Delta,
    Gamma,
    Eta,
}
/// The Goldilocks prime
pub const P: u64 = 0xFFFF_FFFF_0000_0001;
/// Two's complement of `ORDER`, i.e. `2^64 - ORDER = 2^32 - 1`.
pub const NEG_ORDER: u64 = P.wrapping_neg();
impl Field {
    #[inline(always)]
    pub fn new(x: u128) -> Field {
        Field(x)
    }
    #[inline(always)]
    pub fn primitive_nth_root(self, n: u128) -> FieldElement {
        if self.0 == 1 + (1 << 64) - (1 << 32) {
            assert!(
                n <= 1 << 32 && (n & (n - 1)) == 0,
                "Field does not have nth root of unity where n > 2^32 or not power of two."
            );
            let mut root = FieldElement(1753635133440165772, self);

            let mut order = 1 << 32;

            while order != n {
                root = root.pow(2); // Square the root
                order /= 2; // Halve the order
            }

            root
        } else {
            panic!("Unknown field, can't return root of unity.");
        }
    }
    #[inline(always)]
    pub fn generator(self) -> FieldElement {
        assert!(
            self.0 == 1 + (1 << 64) - (1 << 32),
            "Do not know generator for other fields beyond 2^64 - 2^32 + 1"
        );
        FieldElement(7, self)
    }
}

impl PartialEq for Field {
    #[inline(always)]
    fn eq(&self, other: &Field) -> bool {
        self.0 == other.0
    }
}

#[derive(Clone, Copy)]
pub struct FieldElement(pub u128, pub Field);

impl FieldElement {
    #[inline(always)]
    pub fn new(x: u128, field: Field) -> FieldElement {
        FieldElement(x % field.0, field)
    }

    #[inline(always)]
    pub fn zero(field: Field) -> FieldElement {
        FieldElement(0, field)
    }

    #[inline(always)]
    pub fn one(field: Field) -> FieldElement {
        FieldElement(1, field)
    }

    #[inline(always)]
    pub fn modulus(&self) -> u128 {
        self.1 .0
    }

    #[inline(always)]
    pub fn inverse(&self) -> FieldElement {
        let mut inv = 1;
        let mut base = self.0;
        let mut exp = self.1 .0 - 2;
        while exp > 0 {
            if exp % 2 == 1 {
                inv = (inv * base) % self.1 .0;
            }
            base = (base * base) % self.1 .0;
            exp /= 2;
        }
        FieldElement(inv, self.1)
    }

    #[inline(always)]
    pub fn pow(&self, exp: u128) -> FieldElement {
        let mut res = 1;
        let mut base = self.0;
        let mut exp = exp;
        while exp > 0 {
            if exp % 2 == 1 {
                res = (res * base) % self.1 .0;
            }
            base = (base * base) % self.1 .0;
            exp /= 2;
        }
        FieldElement(res % self.1 .0, self.1)
    }

    #[inline(always)]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut e = self.0.to_be_bytes().to_vec();
        let mut f = self.1 .0.to_be_bytes().to_vec();
        e.append(&mut f);
        e
    }

    #[inline(always)]
    pub fn from_bytes(bytes: &[u8]) -> FieldElement {
        let mut x = [0u8; 16];
        let mut y = [0u8; 16];
        x.copy_from_slice(&bytes[..16]);
        y.copy_from_slice(&bytes[16..]);
        FieldElement(
            u128::from_be_bytes(x) % u128::from_be_bytes(y),
            Field(u128::from_be_bytes(y)),
        )
    }
}

impl Add for FieldElement {
    type Output = FieldElement;
    #[inline(always)]
    fn add(self, other: FieldElement) -> FieldElement {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        FieldElement((self.0 + other.0) % self.1 .0, self.1)
    }
}

impl AddAssign for FieldElement {
    #[inline(always)]
    fn add_assign(&mut self, other: FieldElement) {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        self.0 = (self.0 + other.0) % self.1 .0;
    }
}

impl Sub for FieldElement {
    type Output = FieldElement;
    #[inline(always)]
    fn sub(self, other: FieldElement) -> FieldElement {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        if self.0 < other.0 {
            FieldElement((self.0 + self.1 .0 - other.0) % self.1 .0, self.1)
        } else {
            FieldElement((self.0 - other.0) % self.1 .0, self.1)
        }
    }
}

impl SubAssign for FieldElement {
    #[inline(always)]
    fn sub_assign(&mut self, other: FieldElement) {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        if self.0 < other.0 {
            self.0 = (self.0 + self.1 .0 - other.0) % self.1 .0;
        } else {
            self.0 = (self.0 - other.0) % self.1 .0;
        }
    }
}

impl Mul for FieldElement {
    type Output = FieldElement;
    #[inline(always)]
    fn mul(self, other: FieldElement) -> FieldElement {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        FieldElement((self.0 * other.0) % self.1 .0, self.1)
    }
}

impl MulAssign for FieldElement {
    #[inline(always)]
    fn mul_assign(&mut self, other: FieldElement) {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        self.0 = (self.0 * other.0) % self.1 .0;
    }
}

impl Div for FieldElement {
    type Output = FieldElement;
    #[inline(always)]
    fn div(self, other: FieldElement) -> FieldElement {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        let mut inv = 1;
        let mut base = other.0;
        let mut exp = self.1 .0 - 2;
        while exp > 0 {
            if exp % 2 == 1 {
                inv = (inv * base) % self.1 .0;
            }
            base = (base * base) % self.1 .0;
            exp /= 2;
        }
        FieldElement((self.0 * inv) % self.1 .0, self.1)
    }
}

impl DivAssign for FieldElement {
    #[inline(always)]
    fn div_assign(&mut self, other: FieldElement) {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        let mut inv = 1;
        let mut base = other.0;
        let mut exp = self.1 .0 - 2;
        while exp > 0 {
            if exp % 2 == 1 {
                inv = (inv * base) % self.1 .0;
            }
            base = (base * base) % self.1 .0;
            exp /= 2;
        }
        self.0 = (self.0 * inv) % self.1 .0;
    }
}

impl Neg for FieldElement {
    type Output = FieldElement;
    #[inline(always)]
    fn neg(self) -> FieldElement {
        FieldElement(self.1 .0 - self.0, self.1)
    }
}

impl PartialEq for FieldElement {
    #[inline(always)]
    fn eq(&self, other: &FieldElement) -> bool {
        if self.1 != other.1 {
            return false;
        }
        self.0 == other.0
    }
}

impl Eq for FieldElement {}

impl Hash for FieldElement {
    #[inline(always)]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
        self.1.hash(state);
    }
}

impl PartialOrd for FieldElement {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.1 != other.1 {
            return None;
        }
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for FieldElement {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.1 != other.1 {
            panic!("Fields must be same");
        }
        self.0.cmp(&other.0)
    }
}

impl Display for FieldElement {
    #[inline(always)]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Debug for FieldElement {
    #[inline(always)]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod test_field_operations {
    use chrono::Local;

    #[allow(arithmetic_overflow)]
    use super::*;
    use std::primitive;

    #[test]
    fn test_field_add() {
        let field = Field::new(7);
        let a = FieldElement::new(1, field);
        let b = FieldElement::new(2, field);
        let c = a + b;
        assert_eq!(c.0, 3);
    }

    #[test]
    fn test_field_sub() {
        let field = Field::new(7);
        let a = FieldElement::new(1, field);
        let b = FieldElement::new(2, field);
        let c = a - b;
        assert_eq!(c.0, 6);
    }

    #[test]
    fn test_field_mul() {
        let field = Field::new(7);
        let a = FieldElement::new(1, field);
        let b = FieldElement::new(2, field);
        let c = a * b;
        assert_eq!(c.0, 2);
    }

    #[test]
    fn test_field_div() {
        let field = Field::new(7);
        let a = FieldElement::new(1, field);
        let b = FieldElement::new(2, field);
        let c = a / b;
        assert_eq!(c.0, 4);
    }

    #[test]
    fn test_field_inverse() {
        let field = Field::new(7);
        let a = FieldElement::new(2, field);
        let b = a.inverse();
        assert_eq!(b.0, 4);
    }

    #[test]
    fn test_field_pow() {
        let field = Field::new(7);
        let a = FieldElement::new(2, field);
        let b = a.pow(3);
        assert_eq!(b.0, 1);
    }

    #[test]
    #[should_panic]
    fn test_diff_field() {
        let field1 = Field::new(7);
        let field2 = Field::new(8);
        let a = FieldElement::new(1, field1);
        let b = FieldElement::new(2, field2);
        let _ = a + b;
    }

    #[test]
    fn test_negative_number() {
        let field = Field::new(7);
        let a = FieldElement::new(2, field);
        assert_eq!((-a).0, 5);
    }
    #[test]
    fn test_primitive_root_of_unity() {
        let field = Field::new(1 + (1 << 64) - (1 << 32));
        let x = FieldElement::new(1753635133440165772, field);
        let y = x.pow(2);
        let b = field.primitive_nth_root(1 << 32);
        assert_eq!(b.0, x.0);
        let c = field.primitive_nth_root(1 << 31);
        assert_eq!(c.0, y.0);
    }

    #[test]
    fn test_encoding() {
        let field = Field::new(100007);
        let a = FieldElement::new(256, field);
        let b = a.to_bytes();
        let c = FieldElement::from_bytes(&b);
        assert_eq!(a.0, c.0);
    }

    #[test]
    fn bench_field_mul() {
        let field = Field::new(7);
        let a = FieldElement::new(1, field);
        let b = FieldElement::new(2, field);
        let c = a * b;
        assert_eq!(c.0, 2);

        let field = Field::new(P.into());
        let a = FieldElement::new(1687837, field);
        let b = FieldElement::new(1687837, field);
        let t = Local::now();
        for _ in 0..100000000 {
            let _ = a * b;
        }
        println!("Time taken: {:?}", Local::now().signed_duration_since(t));
    }
}
