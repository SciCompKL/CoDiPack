Changelog {#Changelog}
===========================

### v 2.?.? - ???-??-??
 - Features:
  * New helper for adding Enzyme-generated derivative functions to the tape. See \ref Example_24_Enzyme_external_function_helper.
  * Recover primal values from primal values tapes in ExternalFunctionHelper.

 - Bugfix:
  * Uninitialized values in external function helper.
  * External function outputs in Jacobian tapes no longer use unused indices.

### v 2.1.0 - 2023-05-10
 - Features:
   * Helpers for linear system solvers. See Example 21 for details.
   * Event system. See Example 22 for details.
   * Support for shared-memory parallel settings, specifically OpenMP.
   * OpDiLib bindings.
   * CMake support.

 - Bugfix:
   * Initialization fixes.
   * Initialization of passive values.
   * Preaccumulation disables the tape during the algorithm.
   * Proper assigning of 1st order derivatives to 2nd order derivatives.
   * Proper reset of reuse index mangers.
   * Swap of primal value tapes resizes primal value vectors.
   * Tapes handle small chunk sizes correctly.

 - Other:
   * Expression trait for the number of operations.
   * Forward tape evaluation tests for primal value tapes.
   * Gcc warning of no return statement is now an error.
   * Make std::fmin and std::fmax available in namespace.
   * clangd compatibility with CODI_IDE.

### v 2.0.2 - 2022-06-14
 - Bugfix:
   * Prevent copy constrution of tapes.
   * Proper activity tracking for linear index tapes in MeDiPack.

### v 2.0.1 - 2021-11-15
 - Bugfix:
   * Fixed a stack corruption in the reverse interpretation of tapes with a linear index manager.
   * Renamed folders aux -> misc for windows compatibility.

### v 2.0.0 - 2021-09-01
 - Complete rewrite of CoDiPack. See \ref Migration_1_9_2_0 for some migration help.
   * Drop of modular class architecture.
   * Iterators for expression trees.
   * First support for aggregated types in external functions (e.g. std::complex).
     See \ref Example_20_Aggregated_active_type_handling.
   * Autocompletion of template arguments in IDEs. See \ref TemplateDeclaration.
   * Overhaul of tutorials and examples. See \ref TutorialsAndExamples.

### v 1.9.3 - 2020-05-18
 - Bufix:
   * PreaccumulationHelper with changing sizes could give a segmentation fault
   * PreaccumulationHelper with changing zero patterns gave wrong results
   * Dirty adjoint vector after computeJacobian call in Algorithms with a forward evaluation

### v 1.9.2 - 2020-04-28
 - Core functionality:
    * Support for remainder and round function.
 - MeDiPack bindings:
    * Updated MeDiPack bindings to MeDiPack 1.2 (not backwards compatible).
 - MSVC compatibility:
    * Renaming interface -> inter.

### v 1.9.1 - 2020-01-13
 - Bugfix:
   * Missing declaration of MaxDerivativeOrder in UnaryOp type traits.

### v 1.9.0 - 2019-10-30
 - Helper structures:
   * EvaluationHelper:
     + Simplifies the computation of Jacobian and Hessian matrices for functions, function objects and lambda functions.
       Only the function needs to be provided and the helper will compute all derivatives.
     + See Tutorial B1, B1.1 and B1.2.
   * TapeHelper:
     + Provides a more convenient handling of the tape recording and derivative computation process. All CoDiPack
       specifics are hidden and the user can compute the full Jacobian or Hessian matrix with a simple function call.
     + See Tutorial B2
 - Expression rework:
   * Binary and unary expressions are now defined via logic objects that provide the derivative and primal evaluation
     functions.
   * Users can now change the derivative logic for custom types by the specialization of these logic objects.
 - Test suite:
   * The test suite checks now also primal and second order derivatives.
 - Bugfix:
   * Overflow check for linear index handlers

### v 1.8.0 - 2019-01-07
 - Interface:
    - Added function to disable active variables
 - Feature: On the fly combination of entries for the same argument in Jacobian tapes.
   - See CODI_EnableCombineJacobianArguments for details
 - New tutorials:
   - Tutorial for recording of several different tapes in one application at different times
 - Internal:
   - Removed intermediate lambda functions in tape evaluation functions
   - Tape modules are now implemented as structures and no longer included as super macros

### v 1.7.0 - 2018-10-30
 - Support for forward evaluation of tapes finalized
 - Support for primal evaluation of tapes
 - MeDiPack interfaces are now provided with CoDiPack
 - Removed AdjointMPI interfaces

### v 1.6.0 - 2018-03-15
 - First support for forward evaluation of tapes
   * Used in the preaccumulation helper for code sections
     that have more output values than input values.
 - Default Change: Adjoint values are reset by default
   * 10% to 20% performance improvement possible
   * Adjoint values in the RealReverse types are now automatically reset
     during the reverse interpretation. A call to clearAdjoints is no longer
     necessary, only the input values need to be set manually to zero.
 - Bugfix: Input values after a computation for RealReverseIndex tapes
   * If an index for an input value was used before the input value is
     registered, then the result for the input was zero and other
     input values received wrong updates.
 - Bugfix: Memory leak in external function helper.
 - Generalized function for stack evaluation
    * Intermediate function are no longer required for the stack
      evaluation.

### v 1.5.0 - 2017-12-18

 - Support for custom adjoint vectors
   * User provided adjoint vectors can now be used in the reverse evaluation
   * Support added for external functions
 - External function interface addition
   * Extra parameter which provides a general access to the adjoint vector
 - Several helper structures:
   * TapeVectorHelper
     + Easy API for the new custom adjoint vector features
     + See Tutorial A4
   * ExternalFunctionHelper:
     + Helper structure for handling libraries which can not be differentiated or
       to optimize the reverse evaluation for larger code parts
     + See Tutorial A1
   * PreaccumulationHelper:
     + Helper structure for the memory optimization of small code parts
     + If the code part has only few input and output values but is computational
       expensive, the Jacobian is stored instead of the evaluation trace
     + See Tutorial A2
   * StatementPushHelper:
     + Helper structure for the memory optimization of single statements
     + If the Jacobian of the statement can be computed efficiently manually
       this information can be directly provided to CoDiPack
     + See Tutorial A3

### v 1.4.1 - 2017-09-26

 - Bugfix for index managers
   * The statement 'a = a;' would deactivate a as a variable

### v 1.4.0 - 2017-02-16

 - Generalized types: RealForwardGen, RealReverseGen, etc.
 - Removed Float types (Use generalized types)
 - Higher order derivative helper
   * Select single derivatives by there order and number
   * Set all derivatives of a specific order
   * See tutorials 7, 7.1 and 7.2
 - Write/read tapes to and from files
 - Management of multiple tapes
 - External functions interface changed
   * User functions have now a pointer to the calling tape.

### v 1.3.0 - 2016-09-14

 - Primal value taping
   * New chunk types: RealReversePrimal, RealReversePrimalIndex
   * New unchecked types: RealReversePrimalUnchecked, RealReversePrimalIndexUnchecked
   * New vector types: RealReversePrimalVec, RealReversePrimalIndexVec
 - Specialization for numeric_limits
 - Added erf and erfc

### v 1.2.3 - 2016-06-13

 - Bugfix for higher order derivatives:
   * The check for zero Jacobians and adjoints was wrong
 - Added fmin and fmax as expressions

### v 1.2.2 - 2016-05-26

 - Bugfix for prefix increment operator

### v 1.2.1 - 2016-04-20

 - Introduction of the macro CODI_INLINE. This can be used to force the inlining of CoDiPack.

### v 1.2 - 2016-04-07

 - ReferenceActiveReal type for the optimizations of arguments
   * e.g. w = x*x*x + sin(x);
   * stores 1 Jacobian instead of 4
 - Assign optimization index manger (new default one)
 - Full vector mode support
 - Module based implementation of the tapes
 - Own assert command

### v 1.1 - 2015-12-23

 - Index reuse chunk and unchecked tape
 - Template files for binary and unary expressions
 - Support for manual push operations on the tape

### v 1.0 - 2015-06-19
 - Chunk tape
 - Unchecked tape
 - External functions
