Example 10 - External function helper {#Example_10_External_function_helper}
=======

**Goal:** Add external functions to the tape via a helper structure.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:** \ref func_linearSystemSolve
\snippet examples/Example_10_External_function_helper.cpp Function

**Full code:**
\snippet examples/Example_10_External_function_helper.cpp Example 10 - External function helper

The forward and reverse equations for the linear system solve are defined in \ref func_linearSystemSolve. The function
`solve2_rev` implements the reverse mode equation according to this definition. The reverse mode function has to adhere
to the header definition [ReverseFunc](@ref codi::ExternalFunctionHelper::ReverseFunc). All arguments are derived from
the usual source transformation notation where `_b` corresponds to the bar value of the primal. The sizes of `x` and `y`
are given by `m` and `n` respectively. The parameter `d` holds user defined data. This data is provided by the helper
and can be retrieved with a call to
[getExternalFunctionUserData](@ref codi::ExternalFunctionHelper::getExternalFunctionUserData). For details on how to use
this structure please see \ref Example_11_External_function_user_data.

The steps for using the
codi::ExternalFUnctionHelper are then quite simple. The user has to:
 - Create the helper
 - Add the inputs
 - Add the outputs
 - Call the primal function
 - Add the reverse function to the tape
 
For the fourth step the user has two choices. Either a primal evaluation function with the header
[PrimalFunc](@ref codi::ExternalFunctionHelper::PrimalFunc) can be implemented or the original function can be
evaluated with the CoDiPack type. In the second option nothing is recorded during this evaluation.

After the function is registered on the tape, the helper is reset so that it can be used to push another external
function.
