Tutorial 2 - Reverse mode AD {#Tutorial_02_Reverse_mode_AD}
=======

**Goal:** Introduction to reverse mode AD with CoDiPack.

**Prerequisite:** AD reverse mode. See \ref sec_reverseAD

**Function:** \ref func_simple1to1
\snippet tutorials/Tutorial_02_Reverse_mode_AD.cpp Function

**Full code:**
\snippet tutorials/Tutorial_02_Reverse_mode_AD.cpp Tutorial 2 - Reverse mode AD

The derivative computation with the reverse mode of CoDiPack needs in total nine steps:
  - Start the recording
  - A definition of the input variables.
  - Record the function evaluation on the tape.
  - A definition of the output variables.
  - Stop the recording.
  - Setting the direction of the derivative.
  - Evaluate the tape recorded by CoDiPack.
  - Get the results of the derivatives.
  - Clear the tape and adjoints.
  
Most of these steps need to call functions on the global tape structure of CoDiPack. This structure can be accessed with
the function [getTape](@ref codi::LhsExpressionInterface::getTape).

#### Step 1: Start the recording

For this step, it is sufficient to call the method [setActive](@ref codi::ReverseTapeInterface::setActive). This will
enable the recording in CoDiPack for all statements evaluated after the `setActive` call.

#### Step 2: Defining the input variables

In the reverse AD equation, the variable \f$x\f$ describes the vector of input variables. CoDiPack needs to know about
all these values. On each one, the method [registerInput](@ref codi::ReverseTapeInterface::registerInput) needs to be
called.

#### Step 3: Recording the function evaluation

In this step, CoDiPack is only indirectly involved. The function \f$f\f$ needs to be evaluated in the program and CoDiPack
needs to record the statements that are called during the evaluation. It is therefore necessary to write the function
\f$f\f$ such that it uses the CoDiPack type. How this is done depends on the program that is differentiated.
The best option is to write the function as a template function such that the calculation type is flexible. The second
option is most of the time used when software with a large code base is differentiated. Here, a global typedef like
`using Real = codi::RealReverse` is used and all doubles in the program are changed to this typedef. The calculation
type can then be changed during compile time and different executables can be generated.

#### Step 4: Defining the output variables

In the reverse AD equation, the variable \f$y\f$ describes the vector of output variables. CoDiPack needs to know about
all these values. On each one, the method [registerOutput](@ref codi::ReverseTapeInterface::registerOutput) needs to be
called.

#### Step 5: Stop the recording

For this step it is sufficient to call the method [setPassive](@ref codi::ReverseTapeInterface::setPassive). This will
disable the recording in CoDiPack for all statements evaluated after the `setPassive` call.

#### Step 6: Setting the direction of the derivative

In the reverse AD equation, the variable \f$\bar y\f$ is the vector of adjoint output variables. Since these
variables serve in the reverse mode of AD as inputs, they need to be defined before the reverse evaluation is started.
The ''bar'' values of the variables can be accessed via the functions [gradient](@ref codi::LhsExpressionInterface::gradient) and
[setGradient](@ref codi::LhsExpressionInterface::setGradient). The default value is zero for the gradients.

For an access to the gradient data for variables which run out of scope, see \ref IdentifierManagement.

#### Step 7: Evaluating the tape

The simplest way to evaluate the full tape is to call the method [evaluate](@ref codi::ReverseTapeInterface::evaluate).
This will evaluate the full tape from the last recorded statement to the first one.

For a partial evaluation of the tape see \ref Example_03_Positional_tape_evaluations. There are also other evaluation possibilities available, which are
described in other [Tutorials](@ref TutorialsAndExamples).

#### Step 8: Get the resulting derivatives

In the reverse AD equation, the variable \f$\bar x\f$ is the vector of adjoint input variables. After an
evaluation of the tape, these values are populated with the derivative data. The ''bar'' values of the variables can be
accessed via the functions [gradient](@ref codi::LhsExpressionInterface::gradient() const) and
[getGradient](@ref codi::LhsExpressionInterface::getGradient).

CoDiPack will clear these values only when the tape is reset or the adjoints are cleared. If these calls are forgotten
and multiples evaluations are performed, then the derivative data will be accumulated.

#### Step 9: Clear the tape and adjoints

The last step is very important when multiple tape evaluations or recordings are performed during one program execution.
A tape can be reset via the function [reset](@ref codi::ReverseTapeInterface::reset).

Depending on the use case, a full tape reset might not be necessary. See \ref Tutorial_03_Full_jacobian_computation and
\ref Tutorial_05_Repeated_tape_recordings for more details.





