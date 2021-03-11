Example 2 - Custom adjoint vector evaluation {#Example_02_Custom_adjoint_vector_evaluation}
=======

**Goal:** Learn to use a custom adjoint vector in a CoDiPack reverse mode tape.

**Prequesties:** \ref Tutorial_4_Vector_mode_AD, [Identifier management](@ref IdentifierManagement).

**Function:** \ref func_simpleNto2
\snippet examples/Example_02_Custom_adjoint_vector_evaluation.cpp Function

**Full code:**
\snippet examples/Example_02_Custom_adjoint_vector_evaluation.cpp Example 2 - Custom adjoint vector helper

Notes on using the codi::CustomAdjointVectorHelper:
 - Internal adjoint vector is automatically created.
 - The default tape is the global tape. (Can be changed with [setTape](@ref codi::CustomAdjointVectorInterface::setTape)
 - Identifiers need to be used from the variables.
 - Multiple instances of the CustomAdjointVectorHelper are thread safe. (Each thread can have its own instance.)
