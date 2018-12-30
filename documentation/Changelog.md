Changelog      {#Changelog}
===========================

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
