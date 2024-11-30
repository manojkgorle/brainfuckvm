use core::panic;
use std::collections::HashMap;

use crate::fields::{Field, FieldElement};

// Struct holding all the registers.
pub struct Register {
    pub field: Field,
    pub cycle: FieldElement,
    pub instruction_pointer: FieldElement,
    pub current_instruction: FieldElement,
    pub next_instruction: FieldElement,
    pub memory_pointer: FieldElement,
    pub memory_value: FieldElement,
    pub memory_value_inverse: FieldElement,
}

impl Register {
    pub fn new(field: Field) -> Register {
        Register {
            field,
            cycle: FieldElement::zero(field),
            instruction_pointer: FieldElement::zero(field),
            current_instruction: FieldElement::zero(field),
            next_instruction: FieldElement::zero(field),
            memory_pointer: FieldElement::zero(field),
            memory_value: FieldElement::zero(field),
            memory_value_inverse: FieldElement::zero(field),
        }
    }
}

// Struct holding the virtual machine.
pub struct VirtualMachine {
    pub field: Field,
}

impl VirtualMachine {
    pub fn new(field: Field) -> VirtualMachine {
        VirtualMachine { field }
    }

    // compiles and executes the program.
    pub fn execute(&self, code: String) -> (i32, String, String) {
        let program = self.compile(code);
        self.run(&program, "".to_string())
    }

    // compiles the program and generates instruction set.
    pub fn compile(&self, code: String) -> Vec<FieldElement> {
        let field = self.field;
        let zero = FieldElement::zero(field);
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, field) };

        let mut program: Vec<FieldElement> = Vec::new();
        let mut stack: Vec<usize> = Vec::new();

        for c in code.chars() {
            program.push(f(c));

            if c == '[' {
                program.push(zero);
                stack.push(program.len() - 1);
            } else if c == ']' {
                program.push(FieldElement::new(
                    (stack[stack.len() - 1] + 1) as u128,
                    field,
                ));
                program[stack[stack.len() - 1]] = FieldElement::new(program.len() as u128, field);
                stack.pop();
            }
        }

        program
    }

    // runs the program and returns the output
    pub fn run(&self, program: &[FieldElement], input_data: String) -> (i32, String, String) {
        let field = self.field;
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, field) };

        // intialize the registers and state.
        let mut instruction_pointer = 0;
        let mut memory_pointer = FieldElement::zero(field);
        let mut memory: HashMap<FieldElement, FieldElement> = HashMap::new();
        let mut output_data: String = String::new();
        let mut input_counter = 0;
        let mut running_time = 1;
        while instruction_pointer < program.len() {
            if program[instruction_pointer] == f('+') {
                instruction_pointer += 1;
                let v = *memory.get(&memory_pointer).unwrap_or(&zero) + one;
                memory.insert(memory_pointer, v);
            } else if program[instruction_pointer] == f('-') {
                instruction_pointer += 1;
                let v = *memory.get(&memory_pointer).unwrap_or(&zero) - one;
                memory.insert(memory_pointer, v);
            } else if program[instruction_pointer] == f('>') {
                instruction_pointer += 1;
                memory_pointer = memory_pointer + one;
            } else if program[instruction_pointer] == f('<') {
                instruction_pointer += 1;
                memory_pointer = memory_pointer - one;
            } else if program[instruction_pointer] == f('[') {
                if *memory.get(&memory_pointer).unwrap_or(&zero) == zero {
                    instruction_pointer = program[instruction_pointer + 1].0 as usize;
                } else {
                    instruction_pointer += 2;
                }
            } else if program[instruction_pointer] == f(']') {
                if *memory.get(&memory_pointer).unwrap_or(&zero) != zero {
                    instruction_pointer = program[instruction_pointer + 1].0 as usize;
                } else {
                    instruction_pointer += 2;
                }
            } else if program[instruction_pointer] == f('.') {
                instruction_pointer += 1;
                let result_char = char::from_u32((memory[&memory_pointer].0 % 256) as u32)
                    .expect("Value out of valid Unicode range");
                output_data.push(result_char);
            } else if program[instruction_pointer] == f(',') {
                instruction_pointer += 1;
                let mut c = char::default();
                if input_counter < input_data.len() {
                    c = input_data.chars().nth(input_counter).unwrap();
                    input_counter += 1;
                } else {
                    // @todo implement getch handler?
                }
                memory.insert(memory_pointer, f(c));
            } else {
                panic!(
                    "unrecognized instruction at {}, {:?}",
                    instruction_pointer,
                    char::from_u32((program[instruction_pointer].0 % 256) as u32)
                );
            }

            running_time += 1;
        }

        (running_time, input_data, output_data)
    }

    pub fn simulate(&self, program: &[FieldElement], input_data: String) {
        let field = self.field;
        let zero = FieldElement::zero(field);
        let one = FieldElement::one(field);
        let two = FieldElement::new(2, field);
        let f = |x: char| -> FieldElement { FieldElement::new((x as u32) as u128, field) };
        let mut register = Register::new(field);
        register.current_instruction = program[0];
        if program.len() == 1 {
            register.next_instruction = zero;
        } else {
            register.next_instruction = program[1];
        }

        let mut memory: HashMap<FieldElement, FieldElement> = HashMap::new();
        let input_counter = 0;
        let output_data = String::new();

        // @todo create tables.
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compile() {
        let vm = VirtualMachine::new(Field(18446744069414584321));
        let code = "++++++++[>++++[>++>+++>+++>+<<<<-]>+>+>->>+[<]<-]>>.>---.+++++++..+++.>>.<-.<.+++.------.--------.>>+.>++.".to_string();
        let program = vm.compile(code);
        assert_eq!(program.len(), 112);

        let code2 = ">>[++-]<".to_string();
        let program2 = vm.compile(code2);
        assert_eq!(program2.len(), 10);
    }

    #[test]
    fn test_execute() {
        let vm = VirtualMachine::new(Field(18446744069414584321));
        let code = "++++++++[>++++[>++>+++>+++>+<<<<-]>+>+>->>+[<]<-]>>.>---.+++++++..+++.>>.<-.<.+++.------.--------.>>+.>++.".to_string();
        let code2 = ">>[++-]<".to_string();
        let (running_time, input_data, output_data) = vm.execute(code);
        let expected_output = "Hello World!\n";
        assert_eq!(output_data, expected_output);
    }
}
