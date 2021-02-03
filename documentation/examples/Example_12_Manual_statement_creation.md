Example 12 - Manual statement creation {#Example_12_Manual_statement_creation}
=======

**Goal:** Add a statement with manually computed primal and Jacobian values to the tape.

**Prequesties:** \ref Tutorial_2_Reverse_mode_AD

**Function:** \ref func_poly2D
\snippet examples/Example_12_Manual_statement_creation.cpp Function

**Full code:**
\snippet examples/Example_12_Manual_statement_creation.cpp Example 12 - Manual statement creation


The statement helper can be used to push multiple statements on the tape. It resets itself after each
[endPushStatement](@ref codi::StatementPushHelper::endPushStatement) call. Other calls that push a complete statement
will also reset the structure.

Otherwise the use is quite simple. First, the values need to be computed outside of the taping process. Either by
disabling the tape or by computing with non CoDiPack values. Afterwards, the helper can be used to push the statement:
 - Start with a call to [startPushStatement](@ref codi::StatementPushHelper::startPushStatement)
 - Add all arguments and their Jacobians with [pushArgument](@ref codi::StatementPushHelper::pushArgument)
 - Finish the statement with [endPushStatement](@ref codi:StatementPushHelper::endPushStatement)
 
The other [pushStatement](@ref codi::StatementPushHelper::pushStatement) functions will always push a complete statement.
