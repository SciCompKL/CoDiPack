Example 14 - RefrenceActiveType {#Example_14_ReferenceActiveType}
=======

**Goal:** Optimize statements where the same argument occurs multiple times.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:** \ref func_poly1D
\snippet examples/Example_14_ReferenceActiveType.cpp Function

**Full code:**
\snippet examples/Example_14_ReferenceActiveType.cpp Example 14 - ReferenceActiveType

The codi:ReferenceActiveType can be used to optimize statements where the same argument is used multiple times. This
works only with a Jacobian taping approach and 12 bytes are saved per extra occurrence of the variable.

An alternative is to use the flag `-DCODI_RemoveDuplicateJacobianArguments=1`. This will enable an extra pass in
Jacobian tapes that removes such duplicated arguments but introduces a slight performance reduction.
