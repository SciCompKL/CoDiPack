Example 3 - Positional tape evaluations {#Example_03_Positional_tape_evaluations}
=======

**Goal:** Evaluate only parts of CoDiPack tapes.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet examples/Example_03_Positional_tape_evaluations.cpp Function

**Full code:**
\snippet examples/Example_03_Positional_tape_evaluations.cpp Example 3 - Positional tape evaluation

In this example, some part of the function `func()` can be replaced by another function `funcInner`. The tape recording
and evaluation for func() can either be done in the usual way or by using funcInner through storing its position during
the recording process and manually providing its derivatives (u1_d, u2_d).

In general, positional evaluation can have multiple use cases:
 - Record independent parts of the application and access them on a need basis
 - Directly evaluate recorded functions to compute the Jacobians and store them (See also
    \ref Example_15_Preaccumulation_of_code_parts)
 - etc.

For a partly evaluated/reset tape the following should always be kept in mind:
 - Is the adjoint vector cleared properly
 - Are all active variables still valid after the reset.
