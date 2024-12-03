use crate::fields::{Field, FieldElement};
use crate::multivariate_polynomial::*;
use crate::univariate_polynomial::{interpolate_lagrange_polynomials, Polynomial};
use rand::*;
use crate::fri::*;
use crate::ntt::*;
// we will implement abstract methods in rust using the traits.

#[derive(Debug, Clone)]
pub struct Table {
    pub field: Field,
    base_width: u128,//number of base columns in the table.
    full_width: u128,//total no. of coloumn using extension and all
    length: u128,//Number of rows in the table.
    num_randomizers: u128,//number of randomizer columns added for cryptographic or redundancy purposes.
    height: u128,//Represents the rounded-up next power of two of the table length 
    omicron: FieldElement,//represent the generator eval_domain depends on the generator and the order of the subgroup
    generator: FieldElement,// A generator for the multiplicative group of the field
    order: u128,//order of the generator.
    matrix: Vec<Vec<FieldElement>>,
}
impl Table {
    // Constructor method to create a new instance of `table`
    pub fn new(
        field: Field,
        base_width: u128,
        full_width: u128,
        length: u128,
        num_randomizers: u128,
        height: u128,
        omicron: FieldElement,
        generator: FieldElement,
        order: u128,
        matrix: Vec<Vec<FieldElement>>,
    ) -> Self {
        Table {
            field,
            base_width,
            full_width,
            length,
            num_randomizers,
            height,
            omicron,
            generator,
            order,
            matrix,
        }
    }
    fn new_2(field: Field, base_width: u128, full_width: u128, length: u128,
        num_randomizers: u128, generator:FieldElement, order: u128) -> Self {
     let height = roundup_npow2(length);
     let omicron = derive_omicron(generator, order, height);
     
     Table {
         field,
         base_width,
         full_width,
         length,
         num_randomizers,
         height,
         omicron,
         generator,
         order,
         matrix: Vec::new(), // Initialize as empty
     }
 }
    // dont know how to implement this method
    pub fn unit_distance(&self, omega_order:u128)->u128{
        if self.height ==0{
            return 0;
        }
        omega_order/self.height


    }
    // wrong implementation in py
    pub fn get_interpolating_domain_length( &self)->u128{
        self.height
    }
    pub fn interpolate_degree(self)->u128{
        self.get_interpolating_domain_length()-1

    }
      pub fn has_order_po2(order: u128) -> bool {
        (order & (order - 1)) == 0
    }
    
    pub fn interpolate_columns(self,omega:FieldElement,omega_order:u128, column_indices:Vec<u128>)->Vec<Polynomial>{
        if !has_order_po2(omega_order){
            panic!("omega does not have claimed order");

        }
        let mut polynomial:Vec<Polynomial>=Vec::new();
        if self.height ==0{
            let poly=Polynomial::new_from_coefficients(vec![FieldElement::zero(Field::new(self.field.0))]);
            polynomial.push(poly);
           return  polynomial;
        }

        let mut omicron_domain:Vec<FieldElement>=Vec::new();
        let mut randomizer_domain:Vec<FieldElement>=Vec::new();
        for i in 0..self.height{
            omicron_domain.push(self.omicron.pow(i));
        }
        for i in 0..self.height{
        randomizer_domain.push(omega.pow(2*i+1));}
        let mut domain:Vec<FieldElement>=vec![FieldElement::new(0,self.field);self.height as usize];
        for i in 0..self.height{
            domain[i as usize]=omicron_domain[i as usize]+randomizer_domain[i as usize];
        }

        for c in column_indices{
            let mut trace:Vec<FieldElement>=Vec::new();
     
            for row in self.matrix.iter(){
                trace.push(row[c as usize]);
            }
            let mut randomizers:Vec<FieldElement>=Vec::new();
            // doubt
            for _i in 0..self.num_randomizers{
                randomizers.push(FieldElement::new(random::<u128>(),Field::new(self.field.0)));
            }
            let mut values:Vec<FieldElement>=Vec::new();
            let mut sum =FieldElement::zero(Field::new(self.field.0));
            let zero = FieldElement::zero(Field::new(self.field.0));
      
            for i in 0..self.height {
         
                let x = trace.get(i as usize).unwrap_or(&zero);
                sum =*x + randomizers[i as usize];
                values.push(sum);
            }
            if values.len()!=omicron_domain.len(){
                panic!("length of domain and values are unequal");
            };
        let poly= interpolate_lagrange_polynomials(domain.clone(), values);
            polynomial.push(poly);  
        }
        polynomial
    }
//  fn lde(self,domain:FriDomain)->Vec<FieldElement>{
//     let polynomials = self.interpolate_columns(domain.omega, self.height, (0..self.full_width).collect());
//     for p in polynomials{

//     }


// }

pub fn boundary_constraints_ext(self,challeneges:Vec<FieldElement>){
}
// pub fn boundary_quotients(self,fri_domain:FriDomain,codewords:Vec<Vec<FieldElement>>,challenges:Vec<FieldElement>)->Vec<Vec<FieldElement>>{
//     if (codewords.len()==0){
//         println!("panic! because codewords' argument must have nonzero length")
// }
// let mut quotient_codewords=Vec::new();
// let boundary_constraints= self.boundary_constraints_ext(challenges);
// let mut zeroifier = Vec::new();
// let fri_values = fri_domain.list();
// for i in 0..fri_domain.length{
//     zeroifier.push(fri_values[i ]-FieldElement(1,self.field ))


// }
// let zerofier_inverse=batch_inverse(&zeroifier);
// for l in 0..boundary_constraints.len(){
//     // considering mpo as the univariate polynomial
//   let  mpo:Polynomial=boundary_constraints[l as usize];
//     for i in 0..fri_domain.length{
//         let x=mpo.evaluate()


//     }
// }






}







    pub fn roundup_npow2( len:u128)->u128{
        if len==0{
            return    0;
         
        }
        else if len == 1 {
            return 1;
        }
        // Calculate the next power of two
        let bit_representation = format!("{:b}", len - 1);
      
       
      
       1 <<
    (bit_representation.len() as u128)

    }
    // mutable or clone doubt
    pub fn derive_omicron(generator:FieldElement,generator_order:u128,target_order:u128)->FieldElement{
        let mut t_generator=generator;
        let mut t_order=generator_order;

        while t_order!=target_order{
          t_generator=t_generator.pow(2);
            t_order/=2;
        }
        t_generator



    }
    pub fn has_order_po2( order: u128) -> bool {
        (order & (order - 1)) == 0
    }



    

#[cfg(test)]
mod test_operations{
    use super::*;


    #[test]
    fn test_roundup_npow2(){
        let len:u128 =2;
        let len2:u128 =6;
        let len3=9;
        let round= roundup_npow2(len);
     
        let round2= roundup_npow2(len2);
        let round3= roundup_npow2(len3);
        println!("round2:{}",round2);
        assert_eq!(round,2);
        assert_eq!(round2,8);

        assert_eq!(round3,16);
    }
    #[test]
    fn has_order_po2(){
        let order =4_u128;
        let order2 =5_u128;
        assert !(!Table::has_order_po2(order2));
        assert !(Table::has_order_po2(order));
    }

}



