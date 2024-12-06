use crate::fields::*;
use crate::merkle::*;
use crate::channel::*;
use crate::univariate_polynomial::*;
//we can use ntt fast fns for optimization, but rn we just implement using direct evaluate and multiply functions of polynomials
//some fns, structs and parameters are changed accordingly compared to fri.py and ntt.py, because we are not using extension fields

//pub fn combination_polynomial(boundary:Vec<Polynomials, quotient, ) ->  {
//     a1.pow(2) + a0.pow(2)
// }

// pub fn generate_trace(a0: FieldElement, a1: FieldElement, n: usize) -> Vec<FieldElement> {
//     let mut trace = vec![a0, a1];
//     for i in 2..n {
//         let next = fibonnaci_sq(trace[i - 2], trace[i - 1]);
//         trace.push(next);
//     }
//     trace
// }

// /// Generates the evaluation domain for the trace polynomial,
// /// given trace and generator g.
// /// trace is used to find the number of elements in trace.
// pub fn generate_eval_domain_for_trace(
//     trace: &[FieldElement],
//     g: FieldElement,
// ) -> Vec<FieldElement> {
//     trace
//         .iter()
//         .enumerate()
//         .map(|(i, _)| (g.pow(i as u64)))
//         .collect()
// }

// /// Generates a new evaluation domain used after applying the fri operator.
// /// Eval domain len is a power of 2.
// pub fn next_eval_domain(eval_domain: &[FieldElement]) -> Vec<FieldElement> {
//     let new_eval_domain = eval_domain[..eval_domain.len() / 2].to_vec();
//     new_eval_domain
//         .into_iter()
//         .enumerate()
//         .map(|(_, x)| x.pow(2))
//         .collect()
// }

// /// Applies fri operator.
// /// old_polynomial is the polynomial the fri operator should be applied.
// /// beta is a field element(random) given by verifier.
// pub fn next_fri_polynomial(old_polynomial: &Polynomial, beta: FieldElement) -> Polynomial {
//     // old_polynomial = g(x^2) + x * h(x^2)

//     // check the len of the polynomial
//     let len = old_polynomial.coefficients.len();
//     // h(y)
//     let mut odd_poly: Vec<FieldElement> = Vec::new();
//     // g(y)
//     let mut even_poly: Vec<FieldElement> = Vec::new();
//     for i in 0..len {
//         if i % 2 == 0 {
//             even_poly.push(old_polynomial.coefficients[i]);
//         } else {
//             odd_poly.push(old_polynomial.coefficients[i]);
//         }
//     }

//     // g(y) + beta * h(y)
//     Polynomial::new_from_coefficients(even_poly)
//         + Polynomial::new_from_coefficients(odd_poly).scalar_mul(beta)
// }

// /// Generates the next fri layer, evaluation domain and evaluations.
// pub fn next_fri_layer(
//     old_polynomial: &Polynomial,
//     beta: FieldElement,
//     domain: &[FieldElement],
// ) -> (Polynomial, Vec<FieldElement>, Vec<FieldElement>) {
//     let new_eval_domain = next_eval_domain(domain);
//     let new_polynomial = next_fri_polynomial(old_polynomial, beta);
//     let new_evaluations = new_eval_domain
//         .iter()
//         .map(|x| new_polynomial.evaluate(*x))
//         .collect();
//     (new_polynomial, new_eval_domain, new_evaluations)
// }

// /// takes compostion polynomial, which is the first FRI polynomial
// /// eval_domain or coset, the first fri domain.
// /// evaluations, the evaluations of the composition polynomial at eval_domain.
// /// merkle_tree constructed from the evaluations.
// /// channel to send the data to the verifier and get random numbers.
// ///
// /// returns the fri polynomials, fri domains, fri evaluations and fri merkle trees.
// pub fn fri_commit(
//     composition_poynomial: Polynomial,
//     eval_domain: &[FieldElement],
//     evaluations: &[FieldElement],
//     merkle_root: MerkleTree,
//     channel: &mut Channel,
// ) -> (
//     Vec<Polynomial>,
//     Vec<Vec<FieldElement>>,
//     Vec<Vec<FieldElement>>,
//     Vec<MerkleTree>,
// ) {
//     let mut fri_polys = vec![composition_poynomial.clone()];
//     let mut fri_domains = vec![eval_domain.to_vec()];
//     let mut fri_layers = vec![evaluations.to_vec()];
//     let mut fri_merkle = vec![merkle_root];

//     let field = composition_poynomial.coefficients[0].1;
//     while (fri_polys[fri_polys.len() - 1]).degree() > 0 {
//         let beta = channel.receive_random_field_element(field);
//         let (next_poly, next_eval_domain, next_layer) = next_fri_layer(
//             &fri_polys[fri_polys.len() - 1],
//             beta,
//             &fri_domains[fri_domains.len() - 1],
//         );
//         let next_merkle_tree = MerkleTree::new(&next_layer);
//         fri_polys.push(next_poly);
//         fri_domains.push(next_eval_domain);
//         fri_layers.push(next_layer);
//         // send the next merkle tree root to the verifier
//         channel.send(next_merkle_tree.inner.root().unwrap().to_vec());
//         fri_merkle.push(next_merkle_tree);
//     }

//     // send the last layers free term to the verifier
//     channel.send(fri_layers[fri_layers.len() - 1][0].to_bytes());
//     (fri_polys, fri_domains, fri_layers, fri_merkle)
// }

// /// function takes index and channel along with fri_layers and fri_merkles, sends necessary data over the channel that is used for verifying correctness of fri layers.
// /// It iterates over fri layers and fri merkles and in each iteration it sends:
// /// i.   The element of fri layer at the given index(using fri layer)
// /// ii.  authentication path from the corresponding fri merkle tree.
// /// iii. elements fri sibling. for x it sends -x. if element is cp_i(x), then its sibling is cp_i(-x).
// /// iv.  authentication path of the elements sibling.
// pub fn decommit_fri_layers(
//     idx: usize,
//     fri_layers: &[Vec<FieldElement>],
//     fri_merkle: &[MerkleTree],
//     channel: &mut Channel,
// ) {
//     log::debug!("Decommitting on fri layers for query {}", idx);
//     // we dont send authentication path for element in last layer, as all elements are equal, regardless of query, as they are evaluations of a constant polynomial
//     for (layer, merkle) in fri_layers[..fri_layers.len() - 1]
//         .iter()
//         .zip(fri_merkle[..fri_merkle.len() - 1].iter())
//     {
//         log::debug!("sending elements and merkle proofs for layer");
//         let length = layer.len();
//         let elem_idx = idx % length;
//         channel.send(layer[elem_idx].to_bytes());
//         let proof = merkle.get_authentication_path(elem_idx);
//         channel.send(proof);
//         let sibling_idx = (idx + length / 2) % length;
//         channel.send(layer[sibling_idx].to_bytes());
//         let sibling_proof = merkle.get_authentication_path(sibling_idx);
//         channel.send(sibling_proof);
//     }

//     // send the last layer element.
//     log::debug!("sending element of last layer");
//     channel.send(fri_layers.last().unwrap()[0].to_bytes());
// }

// /// sends
// pub fn decommit_on_query(
//     idx: usize,
//     blow_up_factor: usize,
//     f_eval: &[FieldElement],
//     f_merkle: &MerkleTree,
//     fri_layers: &[Vec<FieldElement>],
//     fri_merkle: &[MerkleTree],
//     channel: &mut Channel,
// ) {
//     log::debug!("Decommitting on query {}", idx);
//     // trace domain is blew up 8 times to get the evaluation domain(i.e. the coset). so x, g*x, g^2 * x are seperated by 8 indices in between.
//     // get f(x), f(g*x), f(g^2 * x) and send them over the channel, along with the merkle proofs.
//     assert!(idx + 2 * blow_up_factor < f_eval.len());
//     channel.send(f_eval[idx].to_bytes()); // f(x)
//     channel.send(f_merkle.get_authentication_path(idx)); // merkle proof for f(x)
//     channel.send(f_eval[idx + blow_up_factor].to_bytes()); // f(g*x)
//     channel.send(f_merkle.get_authentication_path(idx + blow_up_factor)); // merkle proof for f(g*x)
//     channel.send(f_eval[idx + 2 * blow_up_factor].to_bytes()); // f(g^2 * x)
//     channel.send(f_merkle.get_authentication_path(idx + 2 * blow_up_factor)); // merkle proof for f(g^2 * x)
//     decommit_fri_layers(idx, fri_layers, fri_merkle, channel)
// }

// /// decommits for each query x and provides consistency proof for the fri layers.
// /// cp(x) = g(x^2) + x*h(x^2)
// /// cp(-x) = g(x^2) - x*h(x^2)
// /// g(x^2) = (cp(x) + cp(-x))/2
// /// h(x^2) = (cp(x) - cp(-x))/2x
// /// so cp_i(x^2) = g(x^2) + beta*h(x^2)
// /// we only need cp_i(x) and cp_i(-x) to compute cp_i+1(x^2)
// ///
// /// for every query x: x belongs to evaluation domain.
// /// we send f_eval(x), f_eval(g*x), f_eval(g^2 * x) and their merkle proofs.
// /// using f_eval(x), f_eval(g*x), f_eval(g^2*x), we compute cp(x)
// /// and send cp(x) and cp(-x) to the verifier along with their merkle proofs.
// /// continue till the last layer has been reached.
// /// if prover is successful in convencing that, the fri_layer is consistent with the evaluations, then the verifier will accept the proof.
// pub fn decommit_fri(
//     num_of_queries: usize,
//     blow_up_factor: usize,
//     maximum_random_int: u64,
//     f_eval: &[FieldElement],
//     f_merkle: &MerkleTree,
//     fri_layers: &[Vec<FieldElement>],
//     fri_merkle: &[MerkleTree],
//     channel: &mut Channel,
// ) {
//     for _ in 0..num_of_queries {
//         let idx = channel.receive_random_int(0, maximum_random_int, true);
//         decommit_on_query(
//             idx as usize,
//             blow_up_factor,
//             f_eval,
//             f_merkle,
//             fri_layers,
//             fri_merkle,
//             channel,
//         );
//     }
// }

// #[cfg(test)]
// mod test_fri_layer {
//     use crate::{field::Field, utils::*};
//     #[test]
//     fn test_fri() {
//         let field = Field::new(7);
//         let poly = Polynomial::new_from_coefficients(vec![
//             FieldElement(2, field),
//             FieldElement(3, field),
//             FieldElement(0, field),
//             FieldElement(1, field),
//         ]);
//         let domain = vec![FieldElement::new(2, field), FieldElement::new(5, field)];
//         let beta = FieldElement(3, field);
//         let (next_poly, next_eval_domain, next_evaluations) = next_fri_layer(&poly, beta, &domain);
//         assert_eq!(next_poly.coefficients.len(), 2);
//         assert_eq!(next_poly.coefficients[0].0, 4);
//         assert_eq!(next_poly.coefficients[1].0, 3);
//         assert_eq!(next_eval_domain.len(), 1);
//         assert_eq!(next_eval_domain[0].0, 4);
//         assert_eq!(next_evaluations.len(), 1);
//         assert_eq!(next_evaluations[0].0, 2);
//     }
// }

pub struct Fri{
    offset: FieldElement,
    omega: FieldElement,
    initial_domain_length: u128,
    domain: FriDomain,
    num_colinearity_tests: usize,
}

impl Fri{
    pub fn new(offset: FieldElement, omega: FieldElement, initial_domain_length: u128, num_colinearity_tests: usize)-> Self{
        let result = Fri{ offset: (offset), omega: (omega), initial_domain_length: (initial_domain_length), domain: (FriDomain::new(offset, omega, initial_domain_length)), num_colinearity_tests: (num_colinearity_tests) };
        assert!(result.num_rounds()>=1, "cannot do FRI with less than one round");
        result
    }
    pub fn num_rounds(&self)-> usize{
        let mut codeword_len = self.initial_domain_length;
        let mut num =0;
        while codeword_len>1_u128{
            codeword_len /= 2;
            num+=1;
        } 
        num
    }

}

pub struct FriDomain{
      pub  offset: FieldElement,
       pub  omega: FieldElement,
       pub length: u128,
}

impl FriDomain{
    pub fn new(offset: FieldElement,omega: FieldElement, length: u128)->Self{
        Self { offset, omega, length }
    }

    pub fn call(&self, index: usize)->FieldElement{
        self.omega.pow(index as u128)*self.offset
    }
     
    pub fn list(&self)->Vec<FieldElement>{
        let mut list: Vec<FieldElement>= vec![];
        for i in 0..self.length{
        list.push(self.omega.pow(i)*self.offset);
        }
        list
    }

    pub fn evaluate(&self, polynomial: Polynomial)-> Vec<FieldElement>{
        let polynomial = polynomial.scalar_mul(self.offset);
        let mut result: Vec<FieldElement> = vec![];
        // let mut  domain_gen:Vec<FieldElement>=Vec::new();
        for i in 0..self.length{
            // let x=self.omega.pow(i);
            result.push(polynomial.evaluate(self.omega.pow(i)));
            // domain_gen.push(x);
        }
        // println!("domain_gen ={:?}", domain_gen);
      
        result
    }
    //needed if we use extension field
    // pub fn xevaluate(&self, polynomial: Polynomial)-> Vec<FieldElement>{
    //     let polynomial = polynomial.scalar_mul(self.offset);
    //     let mut result: Vec<FieldElement> = vec![];
    //     for i in 0..polynomial.coefficients.len(){
    //         result.push(polynomial.evaluate(self.omega.pow(i as u128)));
    //     }
    //     result
    // }

    pub fn interpolate(&self, values: Vec<FieldElement>)->Polynomial{
        let mut list:Vec<FieldElement> =vec![];
        for i in 0..values.len(){
            list.push(self.omega.pow(i as u128));
        }
        
        interpolate_lagrange_polynomials(list, values).scalar_mul(self.offset.inverse())
    }
    //not written xinterpolate, as it is used for extension field
}

#[cfg(test)]
mod test_fri{
    use super::*;
    #[test]
    fn test_evaluate(){
        let field = Field::new(17);
        let   offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        println!("omega ={:?}", omega);
        let domain = FriDomain::new(offset, omega, length);
        let polynomial = Polynomial::new_from_coefficients(vec![FieldElement::new(1, field), FieldElement::new(2, field)]);
        let values = domain.evaluate(polynomial.clone());
        let finded=vec![FieldElement::new(6, field), FieldElement::new(3, field), FieldElement::new(15, field), FieldElement::new(1, field)];

        println!("values ={:?}", values);
        assert_eq!(values, finded);
        

    }
    #[test]
    fn test_interpolate(){
        let field = Field::new(17);
        let   offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        println!("omega ={:?}", omega);
        let domain = FriDomain::new(offset, omega, length);
        let polynomial = Polynomial::new_from_coefficients(vec![FieldElement::new(1, field), FieldElement::new(2, field)]);
        let values = domain.evaluate(polynomial.clone());
        let finded=vec![FieldElement::new(6, field), FieldElement::new(3, field)];
        let interpolated = domain.interpolate(finded);
        println!("interpolated ={:?}", interpolated);
        assert_eq!(interpolated.coefficients, polynomial.coefficients);
    }
}