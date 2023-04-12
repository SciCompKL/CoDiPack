Example 11 - External function user data {#Example_11_External_function_user_data}
=======

**Goal:** Learn how to use the external function user data structure.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:** \ref func_simple1to1
\snippet examples/Example_11_External_function_user_data.cpp Function

**Full code:**
\snippet examples/Example_11_External_function_user_data.cpp Example 11 - External function user data

The data in codi::ExternalFunctionUserData is accessed in a toroidal pattern. After all added data items have been accessed the next call will start
again at the first data item. This allows for multiple reverse interpretations.
