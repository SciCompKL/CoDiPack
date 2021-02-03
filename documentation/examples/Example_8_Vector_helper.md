Example 8 - Custom adjoint vector helper {#Example_8_Vector_helper}
=======

**Goal:** Vector mode with a custom adjoint vector.

**Prequesties:** \ref Tutorial_4_Vector_mode_AD

**Function:** \ref func_simpleNto2
\snippet examples/Example_8_Vector_helper.cpp Function

**Full code:**
\snippet examples/Example_8_Vector_helper.cpp Example 8 - Vector helper

The custom adjoint vector helper structure ([CustomAdjointVectorHelper](@ref codi::CustomAdjointVectorHelper) allows to
use a custom type for the adjoint vector, by providing a simple tape like interface. Custom adjoint vectors can also be
used via the codi::CustomAdjointVectorEvaluationTapeInterface, but the helper hides some details from the user.

The biggest change is that all operations for the reverse mode have to be done on the helper structure. That is seeding,
evaluation and retrieving.

Since the helper is hardcoded to the adjoint vector type, there is also a codi::CustomAdjointVectorInterface definition.
This interface allows to program for an arbitrary vector dimension. In order to access the vector with this interface the
codi::VectorAccessInterface needs to be used. It can be obtained via the method
[getVectorInterface](@ref codi::CustomAdjointVectorInterface::getVectorInterface()).
