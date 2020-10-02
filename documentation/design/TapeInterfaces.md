Tape Interface Design {#TapeInterfaces}
=======

The full tape interface of CoDiPack is broken down into several smaller interfaces. This has the advantage that:
 - Function with the same topic are clustered together
 - Documentation can be provided in more concise format
 - Other classes can specify more clearly which functionality is required
 - New tape implementations can start with the basic interfaces. (Especially with regard to development guideline TODO: ref)

The disadvantages are that all functionality is spread over more files, but the advantages should outweigh this.

The four basic CoDiPack tapes codi::JacobianLinearTape, codi::JacobianReuseTape, codi::PrimalValueLinearTape and
codi::PrimalValueReuseTape implement the codi::FullTapeInterface. This interface consists of all interfaces defined
in codi/tapes/interfaces and defines all functionality that is required for a tape implementation in order to work with
all helper structures of CoDiPack.

Since the full tape interface is quite cumbersome to implement for a new tape development. It is most of the time enough
to just implement a subset of the interfaces for a working tape implementation. For reverse tapes the minimum required
functionality is defined in the codi::ReverseTapeInterface. A development tutorial on how to implement a simple reverse
tape with the functionality available in CoDiPack is provided int TODO: ref.

The list of all interfaces currently defined in CoDiPack is:
 - Internal:
   - codi::InternalStatementRecordingInterface: Used in codi::ActiveType to initialize the storing of expressions.
 - Basic: (Should be implemented by a tape)
   - codi::GradientAccessTapeInterface: Allow access to the tangent and/or adjoint values.
   - codi::IdentifierInformationTapeInterface: Check for activity of values.
   - codi::ReverseTapeInterface: Basic reverse tape implementation basics.
 - Advanced: (Can be implemented by a tape)
   - codi::CustomAdjointVectorEvaluationTapeInterface: Use user defined adjoint vectors in the tape evaluation.
   - codi::DataManagementTapeInterface: Tape file IO, direct access to data allocation.
   - codi::ExternalFunctionTapeInterface: Allow user defined function in the tape evaluation.
   - codi::ForwardEvaluationTapeInterface: Allow a forward AD mode evaluation of the tape.
   - codi::ManualStatementPushTapeInterface: Push user computed Jacobian data to the tape.
   - codi::PositionalEvaluationTapeInterface: Partly evaluation of the tape.
   - codi::PreaccumulationEvaluationTapeInterface: Tape evaluations without a state change.
   - codi::PrimalEvaluationTapeInterface: Allow a primal (re)evaluation of the tape.

For a detailed explanation of the defined functionality see the interface documentation.

