use crate::fields::*;
use crate::ntt::*;
use crate::tables::*;
use crate::univariate_polynomial::*;

pub struct PermutationArgument {
    all_tables: Table,
    lhs: (usize, usize),
    rhs: (usize, usize),
}

//this method has been implemented in all the tables respectively since it has conditions specific for each table.
