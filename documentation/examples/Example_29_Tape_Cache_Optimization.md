Example 29 - Tape cache optimization {#Example_29_Tape_cache_optimization}
=======

**Goal:** Demonstrate the use of std::complex with CoDiPack.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD, \ref Example_02_Custom_adjoint_vector_evaluation

**Function:**
\snippet examples/Example_28_Complex_numbers.cpp Function implementations

**Full code:**
\snippet examples/Example_29_Tape_Cache_Optimization.cpp Example 29 - Tape cache optimization

**Additional information:**
The cache optimizer performs a lifetime analysis of the identifiers on the tape. The identifiers are redistributed
for more cache performance during the reverse evaluation.

Since identifiers are redistributed, identifiers from values like `x` should not be used. Instead, the identifier of
`x` should be stored and then the stored value should be given to the optimizer.

The cache optimization is only meaningfully if the tape is evaluated at least 10 times, e.g., as in a reverse
accumulation process.
