Tutorial 1 - Forward mode AD {#Tutorial_01_Forward_mode_AD}
=======

**Goal:** Introduction to forward mode AD with CoDiPack.

**Prerequisite:** AD forward mode. See \ref sec_forwardAD

**Function:** \ref func_simple1to1
\snippet tutorials/Tutorial_01_Forward_mode_AD.cpp Function

**Full code:**
\snippet tutorials/Tutorial_01_Forward_mode_AD.cpp Tutorial 1 - Forward mode AD

The derivative computation with the forward mode of CoDiPack is quite simple and needs only three steps:
  - Set the direction of the derivative on the input variables.
  - Evaluate the function.
  - Get the direction of the derivative from the output variables.

The simplicity comes from the implementation. The tangent data is stored in the active types and no tape is recorded.
All tangent evaluations are directly done in the statement evaluation. This is what ADOL-C calls the ''tapeless'' mode.

If you need to record a tape and evaluate it in a forward manner, please have a look at example \ref Example_06_Forward_tape_evaluation.

#### Step 1: Set the tangent seeding.

In the forward AD equation, the variable \f$x\f$ describes the vector of input variables. On these values the tangent
direction \f$\dot x\f$ needs to be set. For a single variable this can be done with the functions
[gradient](@ref codi::LhsExpressionInterface::gradient) and [setGradient](@ref codi::LhsExpressionInterface::setGradient).

#### Step 2: Tangent function evaluation

In this step, CoDiPack is only indirectly involved. The function \f$f\f$ needs to be evaluated in the program and CoDiPack
needs to evaluate the forward AD mode equation for each statement that are called during the evaluation. It is therefore
necessary to write the function\f$f\f$ such that it uses the CoDiPack type. How this is done depends on the program that
is differentiated. The best option is to write the function as a template function such that the calculation type is
flexible. The second option is most of the time used when software with a large code base is differentiated. Here, a
global typedef like `using Real = codi::RealForward` is used and all doubles in the program are changed to this typedef.
The calculation type can then be changed during compile time and different executables can be generated.

#### Step 3: Get the directional derivative

In the forward AD equation, the variable \f$y\f$ describes the vector of output variables. During the function
evaluation, CoDiPack computed the directional derivative for these variables which is the vector \f$\dot y\f$. For a
single variable, the tangent information can be extracted with the functions
[gradient](@ref codi::LhsExpressionInterface::gradient() const) and
[getGradient](@ref codi::LhsExpressionInterface::getGradient).

#### Notes on multiple tangent computations ####

The forward mode is very simple to use and multiple tangents evaluations do not require any additional effort. The only
think to keep in mind is that tangent values are only reset by CoDiPack if the value is overwritten. Tangent seedings,
that are set on the input values, are not reset and need to be reset by the user.

Left over tangent seedings can also happen if the computational path changes and values from an old evaluation are used.
See \ref Example_01_Old_tangent_leftovers_forward_mode for an example.
