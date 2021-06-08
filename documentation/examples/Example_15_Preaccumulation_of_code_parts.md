Example 15 - Preaccumulation of code parts {#Example_15_Preaccumulation_of_code_parts}
=======

**Goal:** Reduce the memory consumption of code regions by storing the Jacobian for that region.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:** 
\snippet examples/Example_15_Preaccumulation_of_code_parts.cpp Function

**Full code:**
\snippet examples/Example_15_Preaccumulation_of_code_parts.cpp Example 15 - Preaccumulation of code parts

The example demonstrates on an ODE how preaccumulation works. Without preaccumulation 165170 Byte are used. If
preaccumulation is enabled, then only 261 Byte are required.

Preaccumulation reduces the memory if the code region has only a few inputs and output, but requires a lot of
computational effort.
