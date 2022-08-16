Example 21 - Special handling of linear system solvers {#Example_21_Special_handling_of_linear_system_solvers}
=======

**Goal:** Add a special handling of linear system solvers to your code.

**Prequesties:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet examples/Example_21_Special_handling_of_linear_system_solvers.cpp Function

- The example demonstrates how the already implemented Eigen wrapper for the linear system handling can be used.
- In case you are using Eigen, then you only need to extend from [codi::EigenLinearSystem](@ref codi::EigenLinearSystem)
  and define which solver you are using:
\snippet examples/Example_21_Special_handling_of_linear_system_solvers.cpp Specialization of Eigen solver
- For an abitrary lineary system solver, you have to implement the interface [codi::LinearSystemInterface](@ref codi::LinearSystemInterface).
  Please have a look at [codi::EigenLinearSystem](@ref codi::EigenLinearSystem) and the class documentation of codi::LinearSystemInterface.
  - Most of the methods are iterators for matrices and vectors as well as creation and deletion methods for the matrices and
    vectors.


**Full code:**
\snippet examples/Example_21_Special_handling_of_linear_system_solvers.cpp Example 21 - Special handling of linear system solvers

