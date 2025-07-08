Example 28 - Complex numbers {#Example_28_Complex_numbers}
=======

**Goal:** Demonstrate the use of std::complex with CoDiPack.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet examples/Example_28_Complex_numbers.cpp Function implementations

**Full code:**
\snippet examples/Example_28_Complex_numbers.cpp Example 28 - Complex numbers

**Additional information:**
std::complex is automatically specialized for CoDiPack types. With some compilers this can lead to problems.
Use the option `-DCODI_SpecializeStdComplex=0` to disable this behaviour.
Without the specialization, complex types can be defined by using `codi::ActiveComplex<CoDiType>`, e.g. `codi::ActiveComplex<codi::RealReverse>`.
