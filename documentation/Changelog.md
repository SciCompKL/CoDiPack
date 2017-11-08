Changelog      {#Changelog}
============

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
