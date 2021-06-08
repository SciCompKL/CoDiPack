Example 7 - Primal tape evaluation {#Example_07_Primal_tape_evaluation}
=======

**Goal:** Reevaluate a tape at a different position.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet examples/Example_07_Primal_tape_evaluation.cpp Function

**Full code:**
\snippet examples/Example_07_Primal_tape_evaluation.cpp Example 7 - Primal tape evaluation

A primal reevaluation follows the same rules as an evaluation of the tape. There are a few things to consider:
 - A primal value taping approach needs to be used, e.g. codi::RealReversePrimal
 - The primal value needs to be changed on the tape with [setPrimal](@ref codi::PrimalEvaluationTapeInterface::setPrimal()).
 - New primal values of the outputs have to be received with [getPrimal](@ref codi::PrimalEvaluationTapeInterface::getPrimal()).
 - CoDiPack does not record branching statements. These changes are not considered in a primal reevaluation.
