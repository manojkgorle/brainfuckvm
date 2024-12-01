use crate::fields::{Field, FieldElement};

#[derive(Debug, Clone)]
pub struct table {
    pub field: Field,
    base_width: u128,
    full_width: u128,
    length: u128,
    num_randomizers: Vec<u128>,
    height: u128,
    omicron: FieldElement,
    generator: FieldElement,
    //order: u128,
    matrix: Vec<Vec<FieldElement>>,
}

impl table {
    pub fn new(
        field: Field,
        base_width: u128,
        full_width: u128,
        length: u128,
        height: u128,
    ) -> Self {
        //num_randomizers
        let omicron = field.primitive_nth_root(height);
        let generator = field.generator();
        let matrix = Vec::new();
        let num_randomizers = Vec::new();

        table {
            field: field,
            base_width: base_width,
            full_width: full_width,
            length: length,
            num_randomizers: num_randomizers,
            height: height,
            omicron: omicron,
            generator: generator,
            matrix: matrix,
        }
    }
}
