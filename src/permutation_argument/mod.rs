use crate::fields::*;
use crate::ntt::*;
use crate::tables::*;
use crate::univariate_polynomial::*;

pub struct PermutationArgument {
    all_tables: Table,
    lhs: (usize, usize),
    rhs: (usize, usize),
}

// impl PermutationArgument{
//     pub fn new(all_tables: table, lhs: (usize, usize), rhs: (usize, usize)) -> Self  {
//         Self { all_tables, lhs, rhs }
//     }
//     check if fri domain needs to be a different struct or vvector of field elements simply
//     pub fn quotient(&self, fri_domain: Vec<FieldElement>)->Vec<FieldElement>{
//         let field = self.all_tables.field;

//     }
//}
