Tutorial 3 - Full Jacobian computation (Multiple reverse evaluations) {#Tutorial_03_Full_jacobian_computation}
=======

**Goal:** 
 - Compute several direction with the AD forward mode to create the Jacobian.
 - Evaluate the reverse tape multiple times to compute the Jacobian.

**Prerequisite:** \ref Tutorial_01_Forward_mode_AD, \ref Tutorial_02_Reverse_mode_AD

**Function:** \ref func_simpleNto2
\snippet tutorials/Tutorial_03_Full_jacobian_computation.cpp Function

**Full code:**
\snippet tutorials/Tutorial_03_Full_jacobian_computation.cpp Tutorial 3 - Full Jacobian computation

With the forward and reverse AD modes only directional derivatives or reverse directional derivatives can be computed.
For a computation of the full Jacobian matrix \f$\frac{\d f}{\d x}\f$ the tangent direction \f$\dot x\f$ or adjoint
direction \f$\bar y\f$ need to be set to the unit vectors and the corresponding AD mode has to be evaluated for each
unit vector.

In general, the Jacobian computation for both modes, forward and reverse, uses all things that have been explained in
the basic tutorials.

#### Forward mode ####

For a Jacobian evaluation in the forward mode an iteration over the input dimension is required (**Step 1**). Then a
regular forwad mode evaluation follows (**Step 2** and **Step 3**). The new element is **Step 4**. Here, the tangent
seeding from **Step 1** is cleared so that only unit vectors are set as tangent direction.

#### Reverse mode ####

The Jacobian evaluation in the reverse mode differs a little bit form the procedure in
\ref Tutorial_02_Reverse_mode_AD. The tape is recorded as usual, but only once (**Step 1**). Afterwards, a regular seeding
of the output values is performed and the tape is evaluated (**Step 3** and **Step 4**). The additional function call is
required in **Step 5**. Here, the adjoint vector of CoDiPack is cleared, such that no adjoint values from the previous 
tape evaluation disturb the result. The function [clearAdjoints](@ref codi::ReverseTapeInterface::clearAdjoints) only 
clears the adjoint in contrast to a tape [reset](@ref codi::ReverseTapeInterface::reset) which also clears the tape data.
The optional parameter on the tape reset function allows to skip the clear of the adjoints.

#### Performance improvement ####

CoDiPack clears the adjoints of all intermediate variables automatically, only the ones from input variables are left in
place. (Adjoints of output variables are also cleared.) Therefore a complete reset of the adjoint vector is not necessary.
After the adjoint of an output value is accessed, the adjoint value can be cleared manually with `x[j].gradient() = 0.0`
in **Step 4**. Then **Step 5** is no longer necessary. If not all input variables are cleared this way,
subsequent results will be wrong.







