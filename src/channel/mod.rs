use alloy::{hex::ToHexExt, primitives::U256};

use crate::fields::{Field, FieldElement};
pub struct Channel {
    pub state: String,
    // proof contains all the messages that are exchanged between prover and verifier.
    pub proof: Vec<Vec<u8>>,
    // compressed proof contains only the messages that are sent from prover to verifier.
    pub compressed_proof: Vec<Vec<u8>>,
}

impl Channel {
    pub fn new() -> Channel {
        Channel {
            state: String::new(),
            proof: Vec::new(),
            compressed_proof: Vec::new(),
        }
    }

    pub fn send(&mut self, s: Vec<u8>) {
        log::debug!("sending to channel");
        let data_for_digest = format!("{}:{:?}", self.state, s.clone());
        self.state = sha256::digest(data_for_digest);
        self.proof.push(s.clone());
        self.compressed_proof.push(s);
    }

    pub fn receive_random_field_element(&mut self, field: Field) -> FieldElement {
        let received_int = self.receive_random_int(0, (field.0 - 1) as u64, true);
        log::debug!("received_int: {}", received_int);
        FieldElement::new(received_int as u128, field)
    }

    pub fn receive_random_int(&mut self, min: u64, max: u64, show_in_proof: bool) -> u64 {
        // convert state to hexadecimal number
        let num = (U256::from(min) + U256::from_str_radix(&self.state, 16).unwrap())
            % U256::from(max - min + 1);
        let t_num = min + num.into_limbs()[0];
        let state = self.state.clone() + &t_num.to_be_bytes().to_vec().encode_hex();
        self.state = sha256::digest(state);
        if show_in_proof {
            self.proof.push(num.into_limbs()[0].to_be_bytes().to_vec());
        }
        min + num.into_limbs()[0]
    }

    pub fn proof_size(&self) -> usize {
        let mut size = 0;
        for proof in &self.proof {
            size += proof.len();
        }
        size
    }

    pub fn compressed_proof_size(&self) -> usize {
        let mut size = 0;
        for proof in &self.compressed_proof {
            size += proof.len();
        }
        size
    }
}