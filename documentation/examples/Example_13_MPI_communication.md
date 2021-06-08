Example 13 - MPI communication {#Example_13_MPI_communication}
=======

**Goal:** Learn how MPI communication is differentiated with CoDiPack and MeDiPack.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Full code:**
\snippet examples/Example_13_MPI_communication.cpp Example 13 - MPI communication

The differentiation of MPI communication is done via the MeDiPack library (https://www.scicomp.uni-kl.de/software/medi).
As the example shows all 'MPI_*" routines and functions need to be replaced with the MeDiPack wrapper functions 'AMPI_*'.
Otherwise the MPI functions can be used as usual. The only difference is, that the CoDiPack MPI Datatype needs to be
used in all communications where CoDiPack types are involved.

If different CoDiPack types are used (e.g. codi::RealReverse and codi::RealReverseIndex), then for each one a CoDiPack
MPI Datatype needs to be created.
