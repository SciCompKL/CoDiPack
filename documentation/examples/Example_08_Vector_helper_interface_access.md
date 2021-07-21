Example 8 - Generalized adjoint vector interface {#Example_08_Vector_helper_interface_access}
=======

**Goal:** Use the generialized interface for the adjoint vector access.

**Prerequisite:** \ref Example_02_Custom_adjoint_vector_evaluation

**Function:** \ref func_simpleNto2
\snippet examples/Example_08_Vector_helper_interface_access.cpp Function

**Full code:**
\snippet examples/Example_08_Vector_helper_interface_access.cpp Example 8 - Vector helper interface access

Since the ([CustomAdjointVectorHelper](@ref codi::CustomAdjointVectorHelper)) is hardcoded to the adjoint vector type,
there is also a codi::CustomAdjointVectorInterface definition. This interface allows to program for an arbitrary
vector dimension. In order to access the vector with this interface the codi::VectorAccessInterface needs to be used.
It can be obtained via the method [getVectorInterface](@ref codi::CustomAdjointVectorInterface::getVectorInterface()).
