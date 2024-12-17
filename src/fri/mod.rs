use alloy::primitives::ruint::BaseConvertError;

use crate::channel::*;
use crate::fields::*;
use crate::merkle::*;
use crate::univariate_polynomial::*;
//we can use ntt fast fns for optimization, but rn we just implement using direct evaluate and multiply functions of polynomials
//some fns, structs and parameters are changed accordingly compared to fri.py and ntt.py, because we are not using extension fields

//@todo pass difference quotient polynomial as a parameter too if using random secret initials
//challenges = [alpha, beta], given by fiat shamir
//@todo this has to be inside prover/main?, check once
//boundary_q includes boundary constraints for all tables together, similarly for others

pub fn combination_polynomial(
    processor_q: Vec<Polynomial>,
    memory_q: Vec<Polynomial>,
    instruction_q: Vec<Polynomial>,
    challenges: Vec<FieldElement>,
    height: usize,
    field: Field,
) -> Polynomial {
    let mut combination = Polynomial::new_from_coefficients(vec![]);
    let alpha = challenges[0];
    let beta = challenges[1];
    let x = Polynomial::new_from_coefficients(vec![
        FieldElement::zero(field),
        FieldElement::one(field),
    ]);
    let degree = height - 1;
    //@todo what should be degree here since processor and instruction table can have different heights
    //@todo we can also pass a single vector of all quotient

    for i in 0..processor_q.clone().len() {
        if processor_q[i].degree()<degree {
            let d = degree - processor_q[i].clone().degree();
        combination += Polynomial::new_from_coefficients(vec![alpha]) * processor_q[i].clone()
            + Polynomial::new_from_coefficients(vec![beta])
                * x.clone().pow(d as u128)
                * processor_q[i].clone();
        }
        else {println!("processor quotient {} degree greater than degree max", i);}
    }

    for i in 0..memory_q.clone().len() {
        if memory_q[i].degree()<degree {
        let d = degree - memory_q[i].clone().degree();
        combination += Polynomial::new_from_coefficients(vec![alpha]) * memory_q[i].clone()
            + Polynomial::new_from_coefficients(vec![beta])
                * x.clone().pow(d as u128)
                * memory_q[i].clone();
        }
        else {println!("memory quotient {} degree greater than degree max", i);}
    }

    for i in 0..instruction_q.clone().len() {
        if instruction_q[i].degree()<degree {
        let d = degree - instruction_q[i].clone().degree();
        combination += Polynomial::new_from_coefficients(vec![alpha]) * instruction_q[i].clone()
            + Polynomial::new_from_coefficients(vec![beta])
                * x.clone().pow(d as u128)
                * instruction_q[i].clone();
        }
        else {println!("instruction quotient {} degree greater than degree max", i);}
    }

    combination
}

/// Generates the evaluation domain given offset, height and expansion factor of lde
pub fn generate_eval_domain(
    height: usize,
    expansion_f: usize,
    offset: FieldElement,
    field: Field,
) -> FriDomain {
    let n = height * expansion_f;
    let omicron = field.primitive_nth_root(n as u128);

    FriDomain::new(offset, omicron, n as u128)
}

/// Generates a new evaluation domain used after applying the fri operator.
/// Eval domain len is a power of 2.
pub fn next_eval_domain(eval_domain: FriDomain) -> FriDomain {
    FriDomain::new(
        eval_domain.offset,
        eval_domain.omega.pow(2),
        eval_domain.length / 2,
    )
}

/// Applies fri operator.
/// old_polynomial is the polynomial the fri operator should be applied.
/// beta is a field element(random) given by verifier.
pub fn next_fri_polynomial(old_polynomial: &Polynomial, beta: FieldElement) -> Polynomial {
    //old_polynomial = g(x^2) + x * h(x^2);

    // check the len of the polynomial
    let len = old_polynomial.coefficients.len();
    // h(y)
    let mut odd_poly: Vec<FieldElement> = Vec::new();
    // g(y)
    let mut even_poly: Vec<FieldElement> = Vec::new();
    for i in 0..len {
        if i % 2 == 0 {
            even_poly.push(old_polynomial.coefficients[i]);
        } else {
            odd_poly.push(old_polynomial.coefficients[i]);
        }
    }
    // g(y) + beta * h(y)
    Polynomial::new_from_coefficients(even_poly)
        + Polynomial::new_from_coefficients(odd_poly).scalar_mul(beta)
}

/// Generates the next fri layer, evaluation domain and evaluations.
pub fn next_fri_layer(
    old_polynomial: Polynomial,
    beta: FieldElement,
    domain: FriDomain,
) -> (Polynomial, FriDomain, Vec<FieldElement>) {
    let new_eval_domain = next_eval_domain(domain);
    let new_polynomial = next_fri_polynomial(&old_polynomial, beta);
    let mut new_evaluations = vec![];
    for i in 0..new_eval_domain.length {
        new_evaluations.push(new_polynomial.evaluate(new_eval_domain.omega.pow(i)));
    }
    (new_polynomial, new_eval_domain, new_evaluations)
}

/// takes compostion polynomial, which is the first FRI polynomial
/// eval_domain or coset, the first fri domain.
/// evaluations, the evaluations of the composition polynomial at eval_domain.
/// merkle_tree constructed from the evaluations.
/// channel to send the data to the verifier and get random numbers.
///
/// returns the fri polynomials, fri domains, fri evaluations and fri merkle trees.
pub fn fri_commit(
    composition_poynomial: Polynomial,
    eval_domain: FriDomain,
    evaluations: Vec<FieldElement>,
    merkle_root: MerkleTree,
    channel: &mut Channel,
) -> (
    Vec<Polynomial>,
    Vec<FriDomain>,
    Vec<Vec<FieldElement>>,
    Vec<MerkleTree>,
) {
    let mut fri_polys = vec![composition_poynomial.clone()];
    let mut fri_domains = vec![eval_domain];
    let mut fri_layers = vec![evaluations];
    let mut fri_merkle = vec![merkle_root];

    let field = composition_poynomial.coefficients[0].1;
    while (fri_polys[fri_polys.len() - 1]).clone().degree() > 0 {
        let beta = channel.receive_random_field_element(field);
        let (next_poly, next_eval_domain, next_layer) = next_fri_layer(
            fri_polys[fri_polys.len() - 1].clone(),
            beta,
            fri_domains[fri_domains.len() - 1].clone(),
        );
        let next_merkle_tree = MerkleTree::new(&next_layer);
        fri_polys.push(next_poly);
        fri_domains.push(next_eval_domain);
        fri_layers.push(next_layer.clone());
        // send the next merkle tree root to the verifier
        channel.send(next_merkle_tree.inner.root().unwrap().to_vec());
        println!("fri root sent to compressed proof, layer len: {}", next_layer.clone().len());
        fri_merkle.push(next_merkle_tree);
    }

    // send the last layers free term to the verifier
    channel.send(fri_layers[fri_layers.len() - 1][0].to_bytes());
    println!("fri root sent to compressed proof, layer len: {}", fri_layers[fri_layers.len()-1].clone().len());
    (fri_polys, fri_domains, fri_layers, fri_merkle)
}

/// function takes index and channel along with fri_layers and fri_merkles, sends necessary data over the channel that is used for verifying correctness of fri layers.
/// It iterates over fri layers and fri merkles and in each iteration it sends:
/// i.   The element of fri layer at the given index(using fri layer)
/// ii.  authentication path from the corresponding fri merkle tree.
/// iii. elements fri sibling. for x it sends -x. if element is cp_i(x), then its sibling is cp_i(-x).
/// iv.  authentication path of the elements sibling.
pub fn decommit_fri_layers(
    idx: usize,
    fri_layers: &[Vec<FieldElement>],
    fri_merkle: &[MerkleTree],
    channel: &mut Channel,
) {
    log::debug!("Decommitting on fri layers for query {}", idx);
    // we dont send authentication path for element in last layer, as all elements are equal, regardless of query, as they are evaluations of a constant polynomial
    //println!("frilayer len: {}", fri_layers.len());
    for (layer, merkle) in fri_layers[..fri_layers.len() - 1]
        .iter()
        .zip(fri_merkle[..fri_merkle.len() - 1].iter())
    {
        //println!("fri layer lengths: {}", layer.len());
        log::debug!("sending elements and merkle proofs for layer");
        let length = layer.len();
        let elem_idx = idx % length;
        channel.send(layer[elem_idx].to_bytes());
        let proof = merkle.get_authentication_path(elem_idx);
        channel.send(proof);
        let sibling_idx = (idx + length / 2) % length;
        channel.send(layer[sibling_idx].to_bytes());
        let sibling_proof = merkle.get_authentication_path(sibling_idx);
        channel.send(sibling_proof);
        println!("4 vec sent to compressed proof");
    }

    // send the last layer element.
    log::debug!("sending element of last layer");
    channel.send(fri_layers.last().unwrap()[0].to_bytes());
    println!("1 vec sent to compressed proof");
}

/// sends
pub fn decommit_on_query(
    idx: usize,
    blow_up_factor: usize,        //expansion_f
    f_eval: Vec<&[FieldElement]>, //this contains basecodewords zipped, and extension codewords zipped
    f_merkle: Vec<&MerkleTree>, //this contains MerkleTree of base codewords zipped, and extension codewords zipped
    fri_layers: &[Vec<FieldElement>],
    fri_merkle: &[MerkleTree],
    channel: &mut Channel,
) {
    log::debug!("Decommitting on query {}", idx);
    //we need basecodewords zipped and extension codewords zipped evaluations at x and gx, which will be separated by blowupfactor
    // at x -> clk, ip, ci, ni, mp, mv, inv: (in zipped basecodeword); ipa, mpa, iea, oea: (in zipped extension codeword)
    // at gx -> clk*, ip*, ci*, ni*, mp*, mv*, inv*: (""); ipa*, mpa*, iea*, oea*: ("")
    // get basecodeword[idx], basecodewords[idx+blowupfactor] and extensioncodeword[idx], extensioncodeword[idx+blowupfactor] and send them over the channel, along with the merkle proofs.
    assert!(idx + blow_up_factor < f_eval[0].len());
    let base_x=(f_eval[0][idx].to_bytes().clone());
    channel.send(f_eval[0][idx].to_bytes()); //basecodeword[idx] or f(x)
    let base_x_auth=f_merkle[0].get_authentication_path(idx).clone();
    
    channel.send(f_merkle[0].get_authentication_path(idx)); // merkle proof for basecodeword[idx] or f(x)
    // for i in 0..base_x_auth.len(){
    //     print!("{} ", base_x_auth[i]);
    // } println!("\n base_x_auth of prover");
    // println!("{} of prover",idx);
    // for i in 0..base_x.len(){
    //     print!("{}  ", base_x[i]);
    // } println!("\nbase_x of prover");
    channel.send(f_eval[0][idx + blow_up_factor].to_bytes()); //basecodeword[idx+blowupfactor] or f(g*x)
    channel.send(f_merkle[0].get_authentication_path(idx + blow_up_factor)); // merkle proof for basecodeword[idx+blowupfactor] or f(g*x)

    channel.send(f_eval[1][idx].to_bytes()); //extensioncodeword[idx] or f(x)
    channel.send(f_merkle[1].get_authentication_path(idx)); // merkle proof for extensioncodeword[idx] or f(x)
    channel.send(f_eval[1][idx + blow_up_factor].to_bytes()); //extensioncodeword[idx+blowupfactor] or f(g*x)
    channel.send(f_merkle[1].get_authentication_path(idx + blow_up_factor)); // merkle proof for extensioncodeword[idx+blowupfactor] or f(g*x)
    println!("8 vec sent to compressed proof");
    decommit_fri_layers(idx, fri_layers, fri_merkle, channel)
}

pub fn decommit_fri(
    num_of_queries: usize,
    blow_up_factor: usize,
    maximum_random_int: u64,
    f_eval: Vec<&[FieldElement]>,
    f_merkle: Vec<&MerkleTree>,
    fri_layers: &[Vec<FieldElement>],
    fri_merkle: &[MerkleTree],
    channel: &mut Channel,
) {
    for _ in 0..num_of_queries {
        let idx = channel.receive_random_int(0, maximum_random_int, true);
        println!("Prover index received from channel: {}", idx);
        decommit_on_query(
            idx as usize,
            blow_up_factor,
            f_eval.clone(),
            f_merkle.clone(),
            fri_layers,
            fri_merkle,
            channel,
        );
    }
}

pub struct Fri {
    offset: FieldElement,
    omega: FieldElement,
    initial_domain_length: u128,
    domain: FriDomain,
    num_colinearity_tests: usize,
    expansion_f: usize,
}

impl Fri {
    pub fn new(
        offset: FieldElement,
        omega: FieldElement,
        initial_domain_length: u128,
        num_colinearity_tests: usize,
        expansion_f: usize,
    ) -> Self {
        let result = Fri {
            offset: (offset),
            omega: (omega),
            initial_domain_length: (initial_domain_length),
            domain: (FriDomain::new(offset, omega, initial_domain_length)),
            num_colinearity_tests: (num_colinearity_tests),
            expansion_f: (expansion_f),
        };
        assert!(
            result.num_rounds() >= 1,
            "cannot do FRI with less than one round"
        );
        result
    }
    pub fn num_rounds(&self) -> usize {
        let mut codeword_len = self.initial_domain_length;
        let mut num = 0;
        while codeword_len > self.expansion_f as u128 {
            codeword_len /= 2;
            num += 1;
        }
        num
    }
}

#[derive(Debug, Clone)]
pub struct FriDomain {
    pub offset: FieldElement,
    pub omega: FieldElement,
    pub length: u128,
}

impl FriDomain {
    pub fn new(offset: FieldElement, omega: FieldElement, length: u128) -> Self {
        Self {
            offset,
            omega,
            length,
        }
    }

    pub fn call(&self, index: usize) -> FieldElement {
        self.omega.pow(index as u128) * self.offset
    }

    pub fn list(&self) -> Vec<FieldElement> {
        let mut list: Vec<FieldElement> = vec![];
        for i in 0..self.length {
            list.push(self.omega.pow(i) * self.offset);
        }
        list
    }

    pub fn evaluate(&self, polynomial: Polynomial) -> Vec<FieldElement> {
        let polynomial = polynomial.scale(self.offset.0);
        let mut result: Vec<FieldElement> = vec![];
        // let mut  domain_gen:Vec<FieldElement>=Vec::new();
        for i in 0..self.length {
            // let x=self.omega.pow(i);
            result.push(polynomial.evaluate(self.omega.pow(i)));
            // domain_gen.push(x);
        }
        // println!("domain_gen ={:?}", domain_gen);
        result
    }
    //needed if we use extension field
    pub fn xevaluate(&self, polynomial: Polynomial) -> Vec<FieldElement> {
        let polynomial = polynomial.scalar_mul(self.offset);
        let mut result: Vec<FieldElement> = vec![];
        for i in 0..polynomial.coefficients.len() {
            result.push(polynomial.evaluate(self.omega.pow(i as u128)));
        }
        result
    }
    // interpolate with the given offset
    pub fn interpolate(&self, values: Vec<FieldElement>) -> Polynomial {
        let mut list: Vec<FieldElement> = vec![];
        for i in 0..values.len() {
            list.push(self.omega.pow(i as u128));
        }

        interpolate_lagrange_polynomials(list, values).scalar_mul(self.offset.inverse())
    }
    //interpolate without offset
    pub fn real_interpolate(&self, values: Vec<FieldElement>) -> Polynomial {
        let mut list: Vec<FieldElement> = vec![];
        for i in 0..values.len() {
            list.push(self.omega.pow(i as u128));
        }

        interpolate_lagrange_polynomials(list, values)
    }

    //not written xinterpolate, as it is used for extension field
}

#[cfg(test)]
//@todo write test for fri layer, in stark 101, eval domain was a vector, here it is FriDomain struct, thus some code and values will change accordingly

mod test_fri_layer {
    use super::*;
    //use crate::{field::Field, utils::*};
    #[test]
    fn test_fri() {
        let field = Field::new(17);
        let poly = Polynomial::new_from_coefficients(vec![
            FieldElement(2, field),
            FieldElement(3, field),
            FieldElement(0, field),
            FieldElement(1, field),
        ]);
        let domain = FriDomain {
            offset: (FieldElement::one(field)),
            omega: (FieldElement(4, field)),
            length: (4),
        };

        let beta = FieldElement(3, field);

        let (next_poly, next_eval_domain, next_evaluations) = next_fri_layer(poly, beta, domain);

        assert_eq!(next_poly.coefficients.len(), 2);
        assert_eq!(next_poly.coefficients[0].0, 11);
        assert_eq!(next_poly.coefficients[1].0, 3);
        assert_eq!(next_eval_domain.length, 2);
        assert_eq!(next_eval_domain.omega, FieldElement::new(16, field));
        assert_eq!(next_evaluations.len(), 2);
        assert_eq!(next_evaluations[0].0, 14);
        assert_eq!(next_evaluations[1].0, 8);
    }
}
mod test_fri_domain {
    #![allow(unused_variables)]
    use super::*;
    #[test]
    fn test_evaluate() {
        let field = Field::new(17);
        let offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        println!("omega ={:?}", omega);
        let domain = FriDomain::new(offset, omega, length);
        let polynomial = Polynomial::new_from_coefficients(vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
        ]);
        let values = domain.evaluate(polynomial.clone());
        let finded = vec![
            FieldElement::new(5, field),
            FieldElement::new(2, field),
            FieldElement::new(14, field),
            FieldElement::new(0, field),
        ];
        assert_eq!(values, finded);
    }
    #[test]
    fn test_interpolate() {
        let field = Field::new(17);
        let offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        println!("omega ={:?}", omega);
        let domain = FriDomain::new(offset, omega, length);
        let polynomial = Polynomial::new_from_coefficients(vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
        ]);
        let values = domain.evaluate(polynomial.clone());
        let finded = vec![FieldElement::new(6, field), FieldElement::new(3, field)];
        let interpolated = domain.interpolate(finded);
        println!("interpolated ={:?}", interpolated);
        assert_eq!(interpolated.coefficients, polynomial.coefficients);
    }
    #[test]
    fn test_evaluate2() {
        let field = Field::new(17);
        let offset = FieldElement::new(2, field);
        let length = 4_u128;
        let omega = FieldElement::new(13, field);
        let domain = FriDomain::new(offset, omega, length);
        let polynomial = Polynomial::new_from_coefficients(vec![
            FieldElement::new(1, field),
            FieldElement::new(2, field),
            FieldElement::new(3, field),
        ]);
        let values = domain.evaluate(polynomial.clone());
        let finded = vec![
            FieldElement::new(0, field),
            FieldElement::new(7, field),
            FieldElement::new(9, field),
            FieldElement::new(5, field),
        ];
        assert_eq!(values, finded)
    }
}
