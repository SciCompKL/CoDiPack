Tutorial 5 - Multiple tape recordings {#Tutorial_5_Multiple_tape_recordings}
=======

**Goal:** Be aware of possible side effects in multiple reverse evaluations

**Prequesties:** \ref Tutorial_2_Reverse_mode_AD

**Function:**
\snippet tutorials/Tutorial_5_Multiple_tape_recordings.cpp Function

**Full code:**
\snippet tutorials/Tutorial_5_Multiple_tape_recordings.cpp Tutorial 5 - Multiple tape recordings

The example in [Tutorial 2](@ref Tutorial_2_Reverse_mode_AD) shows already everything that is required for multiple 
tape recordings. The main thing to know is that a [setActive](@ref codi::ReverseTapeInterface::setActive) call does not
reset a CoDiPack tape. Only the call to [reset](@ref codi::ReverseTapeInterface::reset) will completely delete the last
recorded function evaluation. In the context of functions with no side effects this will work without any problems. The
problems start to occur if the functions have side effects and left over values from old computations are used during
the tape recording.

The problem is caused by the linear index management of the CoDiPack codi::RealReverse type. Identifiers from this
index managers are not global and reused during each new tape recording (after a reset of the tape). If an old
identifier is still used in a variable (e.g. the variable has not been overwritten in the current recording and was
active in the previous one), then CoDiPack can not check if this identifier is an _old_ one.

This is demonstrated in the code above. The global variable in `func` is only updated if the parameter is true.
In **Step 1** the global variable is updated which creates the side effect.
**Step 2** does not update the global variable and therefore it is constant with respect to the CoDiPack taping process.
Since the variable has not been reset, it still holds the old identifier and the result is wrong. In **Step 3** the
identifier of the global variable is rest with a call to
[deactivateValue](@ref codi::IdentifierInformationTapeInterface::deactivateValue). Afterwards the derivative result is
correct.

The same evaluation is done with the codi::RealReverseIndex type. With this type the result is correct without handling
the side effect. The codi::RealReverseIndex type uses a reuse index manager and these identifiers are global. Therefore,
on a tape reset the identifiers are not invalidated. Nevertheless, it is still a good idea to deactivate unused values,
such that the run time activity analysis of CoDiPack only records values, that are dependent on the registered input
variables.

Have a look at \ref ActiveTypeOverview for the properties of the default CoDiPack types.
