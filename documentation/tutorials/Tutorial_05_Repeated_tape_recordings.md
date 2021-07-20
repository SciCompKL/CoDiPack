Tutorial 5 - Repeated tape recordings {#Tutorial_05_Repeated_tape_recordings}
=======

**Goal:** Be aware of possible side effects in multiple reverse evaluations

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet tutorials/Tutorial_05_Repeated_tape_recordings.cpp Function

**Full code:**
\snippet tutorials/Tutorial_05_Repeated_tape_recordings.cpp Tutorial 5 - Repeated tape recordings

The example in [Tutorial 2](@ref Tutorial_02_Reverse_mode_AD) already covers everything that is required for multiple
tape recordings. Most importantly, it explains that only a call to [reset](@ref codi::ReverseTapeInterface::reset) will
completely delete the last recorded function evaluation ([setActive](@ref codi::ReverseTapeInterface::setActive) does
not reset a CoDiPack tape). In the context of functions with no side effects this will work without any problems.
Though problems can occur if the functions have side effects and left over values from old computations are used during
the tape recording.

The problem is caused by the linear index management of the CoDiPack codi::RealReverse type. Identifiers from this
index managers are not global and reused during each new tape recording (after a reset of the tape). CoDiPack can not
check whether an identifier is generated during the current recording or if it is an _old_ one. A variable contains an
_old_ identifier if the variable has not been overwritten in the current recording and was active in the previous one.

This is demonstrated in the code above. The global variable in `func` is only updated if the parameter is true.
In **Step 1** the global variable is updated which creates the side effect.
**Step 2** still uses, but does not update the global variable and it is should be treated as a constant in the CoDiPack
taping process. However, since the variable has not been reset, it still holds the old identifier and the result is
wrong (\f$4^3(4+4^2)\neq 756\f$). In **Step 3** the identifier of the global variable is reset with a call to
[deactivateValue](@ref codi::IdentifierInformationTapeInterface::deactivateValue). Afterwards the derivative result is
correct (576).

The same evaluation is done with the codi::RealReverseIndex type. With this type the result is correct without handling
the side effect. The codi::RealReverseIndex type uses a reuse index manager and these identifiers are global. Therefore,
on a tape reset the identifiers are not invalidated. Nevertheless, it is still a good idea to deactivate unused values,
such that the run time activity analysis of CoDiPack only records values, that are dependent on the registered input
variables.

Have a look at \ref ActiveTypeOverview for the properties of the default CoDiPack types.
