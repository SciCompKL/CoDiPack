Example 20 - Aggregated active type handling {#Example_20_Aggregated_active_type_handling}
=======

**Goal:** Learn how to generalize the data extraction for external functions.

**Prequesties:** \ref Example_11_External_function_user_data

**Function:** \ref func_simple1to1
\snippet examples/Example_20_Aggregated_active_type_handling.cpp Function

**Full code:**
\snippet examples/Example_20_Aggregated_active_type_handling.cpp Example 20 - Aggregated active type handling


In the primal application the helper structure from codi::RealTraits are used. For the reverse handling in the external
function the codi::AggregatedTypeVectorAccessWrapperFactory is used to create a wrapped version of the
codi::VectorAccessInterface. In addition codi::ComputationTraits are used for a generlization of the transpose.
