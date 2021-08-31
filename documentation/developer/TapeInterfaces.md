Tape Interface Design {#TapeInterfaces}
=======

The full tape interface of CoDiPack is broken down into several smaller interfaces. This has the advantages that
 - related functions are clustered together,
 - documentation can be provided in more concise format,
 - other classes can specify more clearly which functionality is required,
 - new tape implementations can start with the basic interfaces

A disadvantage is that functionality is spread over more files, but the advantages should outweigh this.

The four basic CoDiPack tapes codi::JacobianLinearTape, codi::JacobianReuseTape, codi::PrimalValueLinearTape and
codi::PrimalValueReuseTape implement the codi::FullTapeInterface. This interface consists of all interfaces defined
in codi/tapes/interfaces and defines all functionality that is required for a tape implementation in order to work with
all helper structures of CoDiPack.

An implementation of the full tape interface requires significant effort for a new tape development. It is most of the time enough
to just implement a subset of the interfaces for a working tape implementation. For reverse tapes the minimum required
functionality is defined in the codi::ReverseTapeInterface. A development tutorial on how to implement a simple reverse
tape with the functionality available in CoDiPack is provided in \ref Developer_Simple_Tape.

Here is the list of all interfaces currently defined in CoDiPack.
 - Internal:
   - codi::InternalStatementRecordingTapeInterface: Used in codi::ActiveType to initialize the storing of expressions.
 - Basic (should be implemented by a tape):
   - codi::GradientAccessTapeInterface: Allow access to the tangent and/or adjoint values.
   - codi::IdentifierInformationTapeInterface: Check for activity of values.
   - codi::ReverseTapeInterface: Basic reverse tape implementation.
 - Advanced (can be implemented by a tape):
   - codi::CustomAdjointVectorEvaluationTapeInterface: Use user-defined adjoint vectors in the tape evaluation.
   - codi::DataManagementTapeInterface: Tape file IO, direct access to data allocation.
   - codi::ExternalFunctionTapeInterface: Allow user-defined functions in the tape evaluation.
   - codi::ForwardEvaluationTapeInterface: Allow a forward AD mode evaluation of the tape.
   - codi::ManualStatementPushTapeInterface: Push user-computed Jacobian data to the tape.
   - codi::PositionalEvaluationTapeInterface: Partly evaluation of the tape.
   - codi::PreaccumulationEvaluationTapeInterface: Tape evaluations without a state change.
   - codi::PrimalEvaluationTapeInterface: Allow a primal (re)evaluation of the tape.

For a detailed explanation of the functionality, see the respective interface documentation.

