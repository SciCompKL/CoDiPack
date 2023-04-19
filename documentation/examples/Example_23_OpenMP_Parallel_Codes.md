Example 23 - OpenMP Parallel Codes {#Example_23_OpenMP_Parallel_Codes}
=======

**Goal:** Learn how to differentiate OpenMP parallel codes with CoDiPack and OpDiLib.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Full code:**
\snippet examples/Example_23_OpenMP_Parallel_Codes.cpp Example 23 - OpenMP Parallel Codes

OpenMP parallel codes are differentiated with the help of the OpDiLib AD tool add-on (https://www.scicomp.uni-kl.de/software/opdi).
Please refer to the OpDiLib documentation for further examples and modes of operation.

The usual AD workflow is the same as for a serial code, except for the initialization and finalization of OpDiLib. In the example
at hand, OpDiLib's macro backend is used, where pragmas such as `#pragma omp parallel` are replaced by macros such as
`OPDI_PARALLEL()` together with end macros such as `OPDI_END_PARALLEL`. If you have a compiler with OMPT support, you can use
OpDiLib's OMPT backend that supports the usage of unmodified OpenMP pragmas (see https://github.com/SciCompKL/OpDiLib
for examples).

OpDiLib can only be used with CoDiPack types that are thread-safe for applications in OpenMP parallel codes. Right now, the dedicated
type that supports this is codi::RealReverseIndexOpenMP.
