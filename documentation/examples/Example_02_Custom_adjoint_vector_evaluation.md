Example 2 - Custom adjoint vector evaluation {#Example_02_Custom_adjoint_vector_evaluation}
=======

**Goal:** Learn to use a custom type in an adjoint vector for a CoDiPack reverse mode tape.

**Prerequisite:** \ref Tutorial_04_Vector_mode_AD, [Identifier management](@ref IdentifierManagement).

**Function:** \ref func_simpleNto2
\snippet examples/Example_02_Custom_adjoint_vector_evaluation.cpp Function

**Full code:**
\snippet examples/Example_02_Custom_adjoint_vector_evaluation.cpp Example 2 - Custom adjoint vector helper

The custom adjoint vector helper structure ([CustomAdjointVectorHelper](@ref codi::CustomAdjointVectorHelper)) allows to
use a custom type for the adjoint vector, by providing a simple tape like interface. The structure creates its own
adjoint vector of the correct size, which is independent of the adjoint vector created by the tape.

The biggest change is that all operations for the reverse mode have to be done on the helper structure. That is seeding,
evaluation and retrieving.

Notes on using the codi::CustomAdjointVectorHelper:
 - Internal adjoint vector is automatically created.
 - The default tape is the global tape. (Can be changed with [setTape](@ref codi::CustomAdjointVectorInterface::setTape))
 - Identifiers need to be used from the variables.
 - Multiple instances of the CustomAdjointVectorHelper are thread safe. (Each thread can have its own instance.)
