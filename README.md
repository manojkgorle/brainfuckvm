#  OOTY-VM
![alt text](IMG_3080.JPG)

The OOTY-VM contains the Rust version of Brainstark.

We have implemented all the cryptography libraries needed for this VM from scratch. Before going through this VM, you must have a strong understanding of the STARK protocol. To help with that, we recommend going through the [stark101](https://starkware.co/stark-101/) by Starkware and reviewing the [stark101](https://starkware.co/stark-101/) implementation in Rust.

We have made some modifications to the [Brainstark](https://aszepieniec.github.io/stark-brainfuck/index), such as using univariate polynomials and applying the FRI protocol on the combination of quotient polynomials. Additionally, we have implemented optimizations to reduce the proving time, such as parallel computing for Lagrange polynomials, among other improvements.


## Structure
- [main.rs](./src/main.rs)- test for prover and verifier 
- [vm](./src/vm/mod.rs)- for compiling the program and creating all the matrix needed for proving the vm
- [stark](./src/stark/mod.rs)- prover and verifier function
- [tables/mod.rs](./src/tables/mod.rs)- all the useful function and structure needed for creating the tables.
- [tables/instruction.rs](./src/tables/instruction.rs)- A module in which written all function for creating the instruction table, all the costraints and quotient poly of the instruction table.
- [tables/memory.rs](./src/tables/memory.rs)- A module in which written all function for creating the memory table, extension columns, all the constraints and quotient poly of the memory table.
- [tables/processor.rs](./src/tables/processor.rs)- A module in which written all function for creating the processor table, extension columns, all the constraints and quotient poly of the processor table.
- [tables/io.rs](./src/tables/io.rs)- A module in which written all function for creating the input and output table, extension columns all the constraints and quotient poly of the memory table.
- [fri](./src/fri/mod.rs)- A module for the fri protcol
- [field](./src/fields/mod.rs) - A module for field operations.
- [polynomial](./src/univariate_polynomial/mod.rs) - A module for univariate polynomials in a field.
- [merkle tree](./src/merkle/mod.rs) - A wrapper around, rs_merkle crate.
- [channel](./src/channel/mod.rs) - Simulates prover and verifier interaction.

## Resources:

- [Brainstark part 1](https://aszepieniec.github.io/stark-brainfuck/engine)
- [Brainstark part 2](https://aszepieniec.github.io/stark-brainfuck/brainfuck)
- [Brainstark part 3](https://aszepieniec.github.io/stark-brainfuck/arithmetization)
- [Brainfuck ISA](https://en.wikipedia.org/wiki/Brainfuck#:~:text=The%20language%20takes%20its%20name,the%20boundaries%20of%20computer%20programming.)
- [Writing your own VM](https://www.jmeiners.com/lc3-vm/#:registers)
- [How to code fri from scratch](https://blog.lambdaclass.com/how-to-code-fri-from-scratch/)
- [A summary on the fri low degree testing](https://eprint.iacr.org/2022/1216.pdf)




