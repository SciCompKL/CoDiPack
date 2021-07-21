Example 9 - OpenMP reverse evaluation {#Example_09_OpenMP_reverse_evaluation}
=======

**Goal:** Use OpenMP to evaluate different adjoint vectors at the same time.

**Prerequisite:** \ref Example_02_Custom_adjoint_vector_evaluation

**Function:** \ref func_simpleNto2
\snippet examples/Example_09_OpenMP_reverse_evaluation.cpp Function

**Full code:**
\snippet examples/Example_09_OpenMP_reverse_evaluation.cpp Example 9 - OpenMP reverse evaluation

`tape.evaluate()` can not be used with OpenMP since both threads would use the internal adjoint vector of the tape
concurrently.
