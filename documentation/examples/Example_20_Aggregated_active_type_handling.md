Example 20 - Aggregated active type handling {#Example_20_Aggregated_active_type_handling}
=======

**Goal:** Learn how to generalize the data extraction for external functions.

**Prequesties:** \ref Example_11_External_function_user_data

**Function:** \ref func_simple1to1
\snippet examples/Example_20_Aggregated_active_type_handling.cpp Function

**Full code:**
\snippet examples/Example_20_Aggregated_active_type_handling.cpp Example 20 - Aggregated active type handling


The example shows how a function that can be called with `double` and `std::complex<double>` can be differentiated with
external functions. The implementation for the differentiation is generalized for the template parameter of the
function. In the recording process, the helper structure `codi::RealTraits` is used for the generalization. For the
reverse handling in the external function, the `codi::AggregatedTypeVectorAccessWrapperFactory` is used to create a
wrapped version of the `codi::VectorAccessInterface`. In addition, `codi::ComputationTraits` are used for a
generalization of the transpose. The advantage of using these traits and the wrapper is that aggregated types can be
used in a similar fashion to standard CoDiPack types. In this case, the same code covers `codi::RealReverse` and
`std::complex<codi::RealReverse>`.
