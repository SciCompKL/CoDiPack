Example 25 - Tape Writers {#Example_25_Tape_Writers}
=======

**Goal:** Write a tape to storage in a binary, text, graphical or math format. 

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet examples/Example_25_Tape_Writers.cpp Function

**Generate Tape**
\snippet examples/Example_25_Tape_Writers.cpp Generate Tape

**Full code:**
\snippet examples/Example_25_Tape_Writers.cpp Example 25 - Tape Writers

The tape writers are used to store a tape. The writers can create a binary or text file which can later be used in a new context to restore a tape. The writers can also create a graphical (for both Jacobian and primal value tapes) or a math representation of the tape (only for primal value tapes). These formats are useful for visualization and debugging.