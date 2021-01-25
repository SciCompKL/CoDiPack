Example 3 - Positional tape evaluations {#Example_3_Positional_tape_evaluations}
=======

**Goal:** Evaluate only parts of CoDiPack tapes.

**Prequesties:** \ref Tutorial_2_Reverse_mode_AD

**Function:**
\snippet tutorials/Example_3_Positional_tape_evaluations.cpp Function

**Full code:**
\snippet tutorials/Example_3_Positional_tape_evaluations.cpp Example 3 - Positional tape evaluation

Positional evaluation can have multiple use cases:
 - Record independent parts of the application and access them on a need basis
 - Directly evaluate recorded functions to compute the Jacobians and store them (See also TODO Preaccumulation)
 - etc.

For a partly evaluated/reset tape the following should always be kept in mind:
 - Is the adjoint vector cleared properly
 - Are all active variables still valid after the reset.
