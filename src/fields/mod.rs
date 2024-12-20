#![allow(dead_code)]
use core::hash::{Hash, Hasher};
use core::hint::unreachable_unchecked;
use std::cmp::{Ord, PartialOrd};
use std::fmt::write;
use std::fmt::{Debug, Display};
use std::iter::Product;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// The Goldilocks prime
pub const P: u64 = 0xFFFF_FFFF_0000_0001;
/// Two's complement of `ORDER`, i.e. `2^64 - ORDER = 2^32 - 1`.
pub const NEG_ORDER: u64 = P.wrapping_neg();

#[derive(Clone, Copy)]
pub struct FieldElement(pub u64);

impl FieldElement {
    pub const ORDER_U64: u64 = P;
    pub const NEG_ORDER: u64 = Self::ORDER_U64.wrapping_neg();
    pub const TWO_ADICITY: usize = 32;
    pub const ZERO: Self = Self::new(0);
    pub const ONE: Self = Self::new(1);
    pub const TWO: Self = Self::new(2);
    pub const NEG_ONE: Self = Self::new(Self::ORDER_U64 - 1);
    pub const GENERATOR: Self = Self::new(7);

    #[inline(always)]
    pub const fn new(x: u64) -> Self {
        if x >= Self::ORDER_U64 {
            Self(x - Self::ORDER_U64)
        } else {
            Self(x)
        }
    }

    pub fn is_zero(&self) -> bool {
        self.0 == 0 || self.0 == Self::ORDER_U64
    }

    pub fn primitive_nth_root(bits: usize) -> Self {
        assert!(bits <= Self::TWO_ADICITY);
        let base = Self::new(1_753_635_133_440_165_772); // generates the whole 2^TWO_ADICITY group
        base.exp_power_of_2(Self::TWO_ADICITY - bits)
    }
    #[inline]
    pub fn as_canonical_u64(&self) -> u64 {
        let mut c = self.0;
        // We only need one condition subtraction, since 2 * ORDER would not fit in a u64.
        if c >= Self::ORDER_U64 {
            c -= Self::ORDER_U64;
        }
        c
    }

    #[inline(always)]
    pub fn from_canonical_u64(n: u64) -> Self {
        Self::new(n)
    }
    /// The elementary function `double(a) = 2*a`.
    ///
    /// This function should be preferred over calling `a + a` or `TWO * a` as a faster implementation may be available for some algebras.
    /// If the field has characteristic 2 then this returns 0.
    #[must_use]
    pub fn double(&self) -> Self {
        self.clone() + self.clone()
    }

    /// The elementary function `square(a) = a^2`.
    ///
    /// This function should be preferred over calling `a * a`, as a faster implementation may be available for some algebras.
    #[must_use]
    pub fn square(&self) -> Self {
        self.clone() * self.clone()
    }

    /// The elementary function `cube(a) = a^3`.
    ///
    /// This function should be preferred over calling `a * a * a`, as a faster implementation may be available for some algebras.
    #[must_use]
    pub fn cube(&self) -> Self {
        self.square() * self.clone()
    }

    /// Compute self^{2^power_log} by repeated squaring.
    #[must_use]
    pub fn exp_power_of_2(&self, power_log: usize) -> Self {
        let mut res = self.clone();
        for _ in 0..power_log {
            res = res.square();
        }
        res
    }

    /// code directly borowwed from plonky3. thanks for these sexy curves, ahh field ops, oops.
    #[inline(always)]
    pub fn inverse(&self) -> FieldElement {
        // if self.is_zero() {
        //     return None;
        // }

        // From Fermat's little theorem, in a prime field `F_p`, the inverse of `a` is `a^(p-2)`.
        //
        // compute a^(p - 2) using 72 multiplications
        // The exponent p - 2 is represented in binary as:
        // 0b1111111111111111111111111111111011111111111111111111111111111111
        // Adapted from: https://github.com/facebook/winterfell/blob/d238a1/math/src/field/f64/mod.rs#L136-L164

        // compute base^11
        let t2 = self.square() * *self;

        // compute base^111
        let t3 = t2.square() * *self;

        // compute base^111111 (6 ones)
        // repeatedly square t3 3 times and multiply by t3
        let t6 = exp_acc::<3>(t3, t3);
        let t60 = t6.square();
        let t7 = t60 * *self;

        // compute base^111111111111 (12 ones)
        // repeatedly square t6 6 times and multiply by t6
        let t12 = exp_acc::<5>(t60, t6);

        // compute base^111111111111111111111111 (24 ones)
        // repeatedly square t12 12 times and multiply by t12
        let t24 = exp_acc::<12>(t12, t12);

        // compute base^1111111111111111111111111111111 (31 ones)
        // repeatedly square t24 6 times and multiply by t6 first. then square t30 and
        // multiply by base
        let t31 = exp_acc::<7>(t24, t7);

        // compute base^111111111111111111111111111111101111111111111111111111111111111
        // repeatedly square t31 32 times and multiply by t31
        let t63 = exp_acc::<32>(t31, t31);

        // compute base^1111111111111111111111111111111011111111111111111111111111111111
        t63.square() * *self
    }

    // @todo optimize this. or correct this.
    #[inline(always)]
    pub fn pow(&self, exp: u64) -> FieldElement {
        let mut res = 1;
        let mut base = self.0;
        let mut exp = exp;
        while exp > 0 {
            if exp % 2 == 1 {
                res = res * base;
            }
            base = base * base;
            exp /= 2;
        }
        FieldElement(res)
    }

    #[inline(always)]
    pub fn to_bytes(&self) -> Vec<u8> {
        self.0.to_be_bytes().to_vec()
    }

    #[inline(always)]
    pub fn from_bytes(bytes: &[u8]) -> FieldElement {
        let mut x = [0u8; 8];
        x.copy_from_slice(&bytes[..8]);
        FieldElement(u64::from_be_bytes(x))
    }
}

impl Add for FieldElement {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        let (sum, over) = self.0.overflowing_add(rhs.0);
        let (mut sum, over) = sum.overflowing_add(u64::from(over) * Self::NEG_ORDER);
        if over {
            // NB: self.value > Self::ORDER && rhs.value > Self::ORDER is necessary but not
            // sufficient for double-overflow.
            // This assume does two things:
            //  1. If compiler knows that either self.value or rhs.value <= ORDER, then it can skip
            //     this check.
            //  2. Hints to the compiler how rare this double-overflow is (thus handled better with
            //     a branch).
            assume(self.0 > Self::ORDER_U64 && rhs.0 > Self::ORDER_U64);
            branch_hint();
            sum += Self::NEG_ORDER; // Cannot overflow.
        }
        Self::new(sum)
    }
}

impl AddAssign for FieldElement {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for FieldElement {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let (diff, under) = self.0.overflowing_sub(rhs.0);
        let (mut diff, under) = diff.overflowing_sub(u64::from(under) * Self::NEG_ORDER);
        if under {
            // NB: self.value < NEG_ORDER - 1 && rhs.value > ORDER is necessary but not
            // sufficient for double-underflow.
            // This assume does two things:
            //  1. If compiler knows that either self.value >= NEG_ORDER - 1 or rhs.value <= ORDER,
            //     then it can skip this check.
            //  2. Hints to the compiler how rare this double-underflow is (thus handled better
            //     with a branch).
            assume(self.0 < Self::NEG_ORDER - 1 && rhs.0 > Self::ORDER_U64);
            branch_hint();
            diff -= Self::NEG_ORDER; // Cannot underflow.
        }
        Self::new(diff)
    }
}

impl SubAssign for FieldElement {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for FieldElement {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        reduce128(u128::from(self.0) * u128::from(rhs.0))
    }
}

impl MulAssign for FieldElement {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Product for FieldElement {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|x, y| x * y).unwrap_or(Self::ONE)
    }
}

impl Div for FieldElement {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl Neg for FieldElement {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(Self::ORDER_U64 - self.as_canonical_u64())
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.as_canonical_u64() == other.as_canonical_u64()
    }
}

impl Eq for FieldElement {}

impl Hash for FieldElement {
    #[inline(always)]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl Ord for FieldElement {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.as_canonical_u64().cmp(&other.as_canonical_u64())
    }
}

impl PartialOrd for FieldElement {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
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

// code directly borrowed from plonky3 goldilocks crate.

/// Squares the base N number of times and multiplies the result by the tail value.
#[inline(always)]
fn exp_acc<const N: usize>(base: FieldElement, tail: FieldElement) -> FieldElement {
    base.exp_power_of_2(N) * tail
}

#[inline(always)]
pub fn assume(p: bool) {
    debug_assert!(p);
    if !p {
        unsafe {
            unreachable_unchecked();
        }
    }
}

/// Try to force Rust to emit a branch. Example:
///
/// ```no_run
/// let x = 100;
/// if x > 20 {
///     println!("x is big!");
///     p3_util::branch_hint();
/// } else {
///     println!("x is small!");
/// }
/// ```
///
/// This function has no semantics. It is a hint only.
#[inline(always)]
pub fn branch_hint() {
    // NOTE: These are the currently supported assembly architectures. See the
    // [nightly reference](https://doc.rust-lang.org/nightly/reference/inline-assembly.html) for
    // the most up-to-date list.
    #[cfg(any(
        target_arch = "aarch64",
        target_arch = "arm",
        target_arch = "riscv32",
        target_arch = "riscv64",
        target_arch = "x86",
        target_arch = "x86_64",
    ))]
    unsafe {
        core::arch::asm!("", options(nomem, nostack, preserves_flags));
    }
}

/// Reduces to a 64-bit value. The result might not be in canonical form; it could be in between the
/// field order and `2^64`.
#[inline]
pub(crate) fn reduce128(x: u128) -> FieldElement {
    let (x_lo, x_hi) = split(x); // This is a no-op
    let x_hi_hi = x_hi >> 32;
    let x_hi_lo = x_hi & FieldElement::NEG_ORDER;

    let (mut t0, borrow) = x_lo.overflowing_sub(x_hi_hi);
    if borrow {
        branch_hint(); // A borrow is exceedingly rare. It is faster to branch.
        t0 -= FieldElement::NEG_ORDER; // Cannot underflow.
    }
    let t1 = x_hi_lo * FieldElement::NEG_ORDER;
    let t2 = unsafe { add_no_canonicalize_trashing_input(t0, t1) };
    FieldElement::new(t2)
}

#[inline]
#[allow(clippy::cast_possible_truncation)]
const fn split(x: u128) -> (u64, u64) {
    (x as u64, (x >> 64) as u64)
}

/// Fast addition modulo ORDER for x86-64.
/// This function is marked unsafe for the following reasons:
///   - It is only correct if x + y < 2**64 + ORDER = 0x1ffffffff00000001.
///   - It is only faster in some circumstances. In particular, on x86 it overwrites both inputs in
///     the registers, so its use is not recommended when either input will be used again.
#[inline(always)]
#[cfg(target_arch = "x86_64")]
unsafe fn add_no_canonicalize_trashing_input(x: u64, y: u64) -> u64 {
    let res_wrapped: u64;
    let adjustment: u64;
    core::arch::asm!(
        "add {0}, {1}",
        // Trick. The carry flag is set iff the addition overflowed.
        // sbb x, y does x := x - y - CF. In our case, x and y are both {1:e}, so it simply does
        // {1:e} := 0xffffffff on overflow and {1:e} := 0 otherwise. {1:e} is the low 32 bits of
        // {1}; the high 32-bits are zeroed on write. In the end, we end up with 0xffffffff in {1}
        // on overflow; this happens be NEG_ORDER.
        // Note that the CPU does not realize that the result of sbb x, x does not actually depend
        // on x. We must write the result to a register that we know to be ready. We have a
        // dependency on {1} anyway, so let's use it.
        "sbb {1:e}, {1:e}",
        inlateout(reg) x => res_wrapped,
        inlateout(reg) y => adjustment,
        options(pure, nomem, nostack),
    );
    assume(x != 0 || (res_wrapped == y && adjustment == 0));
    assume(y != 0 || (res_wrapped == x && adjustment == 0));
    // Add NEG_ORDER == subtract ORDER.
    // Cannot overflow unless the assumption if x + y < 2**64 + ORDER is incorrect.
    res_wrapped + adjustment
}

#[inline(always)]
#[cfg(not(target_arch = "x86_64"))]
unsafe fn add_no_canonicalize_trashing_input(x: u64, y: u64) -> u64 {
    let (res_wrapped, carry) = x.overflowing_add(y);
    // Below cannot overflow unless the assumption if x + y < 2**64 + ORDER is incorrect.
    res_wrapped + Goldilocks::NEG_ORDER * u64::from(carry)
}

#[cfg(test)]
mod test_field_operations {
    use chrono::Local;

    #[allow(arithmetic_overflow)]
    use super::*;
    use std::primitive;

    #[test]
    fn test() {
        let f = FieldElement::new(100);
        assert_eq!(f.as_canonical_u64(), 100);

        // Over the Goldilocks field, the following set of equations hold
        // p               = 0
        // 2^64 - 2^32 + 1 = 0
        // 2^64            = 2^32 - 1
        let f = FieldElement::new(u64::MAX);
        assert_eq!(f.as_canonical_u64(), u32::MAX as u64 - 1);

        let f = FieldElement::from_canonical_u64(u64::MAX);
        assert_eq!(f.as_canonical_u64(), u32::MAX as u64 - 1);

        let f = FieldElement::from_canonical_u64(0);
        assert!(f.is_zero());

        let f = FieldElement::from_canonical_u64(FieldElement::ORDER_U64);
        assert!(f.is_zero());

        assert_eq!(FieldElement::GENERATOR.as_canonical_u64(), 7_u64);

        let f_1 = FieldElement::new(1);
        let f_1_copy = FieldElement::new(1);

        let expected_result = FieldElement::ZERO;
        assert_eq!(f_1 - f_1_copy, expected_result);

        let expected_result = FieldElement::new(2);
        assert_eq!(f_1 + f_1_copy, expected_result);

        let f_2 = FieldElement::new(2);
        let expected_result = FieldElement::new(3);
        assert_eq!(f_1 + f_1_copy * f_2, expected_result);

        let expected_result = FieldElement::new(5);
        assert_eq!(f_1 + f_2 * f_2, expected_result);

        let f_p_minus_1 = FieldElement::from_canonical_u64(FieldElement::ORDER_U64 - 1);
        let expected_result = FieldElement::ZERO;
        assert_eq!(f_1 + f_p_minus_1, expected_result);

        let f_p_minus_2 = FieldElement::from_canonical_u64(FieldElement::ORDER_U64 - 2);
        let expected_result = FieldElement::from_canonical_u64(FieldElement::ORDER_U64 - 3);
        assert_eq!(f_p_minus_1 + f_p_minus_2, expected_result);

        let expected_result = FieldElement::new(1);
        assert_eq!(f_p_minus_1 - f_p_minus_2, expected_result);

        let expected_result = f_p_minus_1;
        assert_eq!(f_p_minus_2 - f_p_minus_1, expected_result);

        let expected_result = f_p_minus_2;
        assert_eq!(f_p_minus_1 - f_1, expected_result);

        let expected_result = FieldElement::new(3);
        assert_eq!(f_2 * f_2 - f_1, expected_result);

        // Generator check
        let expected_multiplicative_group_generator = FieldElement::new(7);
        assert_eq!(
            FieldElement::GENERATOR,
            expected_multiplicative_group_generator
        );

        // Check on `reduce_u128`
        let x = u128::MAX;
        let y = reduce128(x);
        // The following equalitiy sequence holds, modulo p = 2^64 - 2^32 + 1
        // 2^128 - 1 = (2^64 - 1) * (2^64 + 1)
        //           = (2^32 - 1 - 1) * (2^32 - 1 + 1)
        //           = (2^32 - 2) * (2^32)
        //           = 2^64 - 2 * 2^32
        //           = 2^64 - 2^33
        //           = 2^32 - 1 - 2^33
        //           = - 2^32 - 1
        let expected_result = -FieldElement::new(2_u64.pow(32)) - FieldElement::new(1);
        assert_eq!(y, expected_result);

        // assert_eq!(f.exp_u64(10540996611094048183).exp_const_u64::<7>(), f);
        // assert_eq!(y.exp_u64(10540996611094048183).exp_const_u64::<7>(), y);
        // assert_eq!(f_2.exp_u64(10540996611094048183).exp_const_u64::<7>(), f_2);
    }
    #[test]
    fn test_field_add() {
        let f1 = FieldElement::new(100);
        let f2 = FieldElement::new(200 + FieldElement::ORDER_U64);
        let f3 = f1 + f2;
        assert_eq!(f3.as_canonical_u64(), 300);
    }

    #[test]
    fn test_field_sub() {}

    #[test]
    fn test_field_mul() {}

    #[test]
    fn test_field_div() {}

    #[test]
    fn test_field_inverse() {}

    #[test]
    fn test_field_pow() {}

    #[test]
    fn test_negative_number() {}
    #[test]
    fn test_primitive_root_of_unity() {}

    #[test]
    fn test_encoding() {}
}
