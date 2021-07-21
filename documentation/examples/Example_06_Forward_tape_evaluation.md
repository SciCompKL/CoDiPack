Example 6 - Forward mode tape evaluation {#Example_06_Forward_tape_evaluation}
=======

**Goal:** Evaluate a tape with forward AD mode.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet examples/Example_06_Forward_tape_evaluation.cpp Function

**Full code:**
\snippet examples/Example_06_Forward_tape_evaluation.cpp Example 6 - Forward tape evaluation

Forward mode tape evaluation works the same as a reverse mode tape evaluation. Only the inputs need to be seeded and the
derivatives are obtained from the outputs.
