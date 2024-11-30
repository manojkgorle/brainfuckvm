use crate::fields::{Field, FieldElement};

#[derive(Debug, Clone)]
pub struct table{
    pub field: Field,
    base_width: u128,
    full_width: u128,
    length: u128,
    num_randomizers: Vec<u128>,
    height: u128,
    omicron: FieldElement,
    generator: FieldElement,
    //order: u128,
    matrix: Vec<Vec<FieldElement>>
}

impl table{
    pub fn new(self, field: Field, base_width:u128, full_width:u128, length:u128, height:u128){
        //num_randomizers
        self.omicron = field.primitive_nth_root(height);
        self.generator = field.generator();
        self.matrix = Vec<Vec<FieldElement>>::new();
    }
}

