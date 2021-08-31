Example tape implementation with CoDiPack {#Developer_Simple_Tape}
=======

**Goal:** Get to know how a simple operator taping approach can be implemented with CoDiPack

**Prerequisite:** AD reverse mode, see \ref sec_reverseAD; [Identifier management](@ref IdentifierManagement)

**Full code:**
\snippet developer/simpleTape.cpp Simple Tape

The implementation of CoDiPack tape starts with two things,
 - deciding on a taping strategy and the data layout as well as
 - the codi:ReverseTapeInterface.

We start with the discussion of the data layout. In order to keep things simple we want to implement an
operator taping approach that stores the primal input values of each operation. The
identifiers are not managed, we will simply generate a new one whenever required.

An operator taping approach will break down all statements in the program to their operations for the AD tool. A statement
like `w = sqrt(a * a + b * b)` will be broken down into the single assignment code
```
  t1 = a * a;
  t2 = b * b;
  t3 = t1 + t2;
  w = sqrt(t3);
```

With this consideration we only need to handle binary and unary operations. For each of these operations we need to
store the necessary data in order to be able to evaluate the reverse AD equation. For a binary operation
\f$ y = f_{bin}(x, z)\f$ the required data consists of
 - the primal values of \f$x\f$ and \f$z\f$,
 - the identifiers of \f$x\f$, \f$z\f$ and \f$y\f$ as well as
 - an enumerator to identify the operation \f$f_{bin}\f$.

The analogous data is required for a unary operation \f$ y = f_{un}(x)\f$ where the data items with respect to \f$z\f$
are omitted, of course. Each type of data in the above list will be stored in a separate data stream (codi::DataInterface), that
is, we have one data stream for primal values, one for identifiers and one for the operators.

With the defined data layout it is then possible to evaluate the reverse AD equation for each operation, therefore we
can have a look at the codi::ReverseTapeInterface and start the implementation.

#### ReverseTapeInterface

The codi::ReverseTapeInterface defines the basic functionality for the reverse evaluation of a tape and the management
of the taping process. It inherits from the two interfaces codi::GradientAccessTapeInterface and
codi::InternalStatementRecordingTapeInterface. The first one defines the basic access to the adjoint variables. The second
one defines the functionality which is called by the CoDiPack expressions to signal to the tape that a statement needs
to be recorded. Both interfaces are required such that the basic active type (codi::ActiveType) is able to work with the
tape implementation.

In the following the implementation will be discussed step by step.

#### Data definition

Type definitions:
\snippet developer/simpleTape.cpp Data stream - Type definition

Member definitions:
\snippet developer/simpleTape.cpp Data stream - Member definition

Member creation:
\snippet developer/simpleTape.cpp Data stream - Member creation

The three steps from above will create the data management defined at the beginning. In the first step, we define the
three data streams. In this example we use the chunked data stream. It allocates chunks of data every time the already
allocated memory is used up. The chunk definition tells each stream what kind of data should be stored. In our example
each stream consists only of one data item per entry. The `OperatorCode` type is an enum class and will be
covered in more detail later. Note that the type definitions are nested. In order to define
the data stream for the identifier data, the data stream for the operator data is provided as the nested stream.
The design decision for this nesting will become clear when we discuss the reverse interpretation of the data.

The second step defines the members and the third step initializes them. The important thing here is the initialization
of the nested data stream with a call to codi::DataInterface::setNested.

#### Identifier management and adjoint vector

As defined in the taping strategy, identifiers are not managed. Every time an identifier is required, a new one is
generated. Therefore we only keep track of the largest one and define the adjoint vector as a standard vector:
\snippet developer/simpleTape.cpp Identifiers - Member definition

The basic size of the adjoint vector is 1 so that we can always safely access the first entry:
\snippet developer/simpleTape.cpp Identifiers - Member initialization

This allows us to implement the codi::GradientAccessTapeInterface:
\snippet developer/simpleTape.cpp Adjoint - Access

The helper function `checkAndResizeAdjoints` checks the size of the vector and only resizes it if necessary. The
generation of a new identifier is also quite simple but standardized in `generateIdentifier`:
\snippet developer/simpleTape.cpp Identifiers - Helper

Since we can now generate identifiers, it is possible to implement the functions for registering input and output
variables. An input registration will just generate a new identifier and set it as the identifier for the value. We do
not need to create a statement on the tape here since the identifiers and statements are not interlocked (removing the
left hand side identifier from the data stream in a linear index management scheme would interlock them). In the
`registerOutput` implementation, nothing needs to be done since all identifiers are assigned to a unique value
(copy statement optimizations would require to assign a unique identifier here which would store a copy operation).
\snippet developer/simpleTape.cpp Identifiers - Registration

The identifiers are stored in the AD type provided by CoDiPack. The initialization of the identifier in
the AD value is done by the function `initIdentifier` required by the codi::InternalStatementRecordingTapeInterface. We
implement an online activity analysis in this tape. Therefore, all identifiers in the AD values can be initialized with
zero. The zero identifier is used in our implementation to track _passive_ values. These are values that do not depend
on the input values. How this is done is explained in the next section.
\snippet developer/simpleTape.cpp Identifiers - Initialization

#### Storing of expressions/operators

The entry point for the storing of expressions is defined in the codi::InternalStatementRecordingTapeInterface. The method
[store](@ref codi::InternalStatementRecordingTapeInterface::store) is called every time a new value is assigned to a
CoDiPack type. The implementation in the tape performs the basic activity checking of the tape. If the tape is currently
in recording mode, then it forwards the data to the method `storeOperator`:
\snippet developer/simpleTape.cpp Storing - Entry

`storeOperator` is the entry point for storing the expressions. Since CoDiPack usually employs a statement taping
strategy with expressions (\ref Expressions), the storing of pure operators is a little bit more complicated. We need to walk
recursively through the expression tree and store all operators in this tree. In order to do that we define a helper
class `StoreOperator_Impl`. This class is then specialized for all the different nodes in a CoDiPack expression. In our
case these will be
 - codi::UnaryExpression,
 - codi::BinaryExpression,
 - codi::ConstantExpression,
 - codi::LhsExpressionInterface.

The detailed process will be explained for the unary operator, the others are similar:
\snippet developer/simpleTape.cpp Storing - Unary operator

The arguments to the function are the expression `exp`, the tape, the result value reference `resultValue`, the result
identifier reference `resultIdentifier` and the boolean `copy`. `copy` is used to indicate a copy operation only if the
first level of the expression is a codi::LhsExpression. This is explained later. `tape` is just used to access the data
streams and to call `storeOperator` recursively. `exp` ist the current node in the expression tree we want to store. The
two reference arguments `resultValue` and `resultIdentifier` define the value and the identifier of the current node.
The method has to populate these two values and store everything for the reverse interpretation.

The first step is to call `storeOperator` on the nested expression. This will populate the arguments `argValue` and
`argIdentifier` with the result of the nested expression and the identifier for the result. Since this is at least the
second level in the expression tree, we do not want to generate copy operations.

Afterwards, it is checked if the result of the nested expression is active. This is the case if at least one argument of
the operator is active, that is, it has a nonzero identifier. If one argument is active, we store the operation
otherwise the result identifier of this operation is also set to passive (zero).

The storing of the operation starts with the reservation of the data items. All reservations for one operation have to
be done in one block and from the most nested stream to the root stream. This ensures that all data is stored in the
same evaluation block (see also codi::DataInterface). Afterwards the new identifier for the result of the operator is generated and
the data is stored. As described in the introduction we store the primal values of the arguments and the identifiers for
the argument as well as the result. The operator entry is retrieved via a lookup function. This function is specialized
for each operator node in a CoDiPack expression tree. The definition of the enumeration and the lookup function is:
\snippet developer/simpleTape.cpp Storing - Operator codes


The implementation for the other operators is nearly the same:
\snippet developer/simpleTape.cpp Storing - Other operators

The specialization for the constant expressions returns just the value with a passive identifier. The left hand side
expression specialization checks for the copy flag. This will only be true if the method is directly called from the
root node of the expression. In this case we have the assignment `c = a` which needs to be stored. Otherwise it is a
node in the expression tree, e.g. `a` or `b` in `c = a * b`. In this case we just provide the identifier and the primal
value.

### Reverse evaluation

The implementation of the reverse evaluation procedure is done in two functions. The entry point is the
[evaluate](@ref codi::ReverseTapeInterface::evaluate) function from the codi::ReverseTapeInterface definition:
\snippet developer/simpleTape.cpp Evaluation - Entry

Here, we use the [evaluateReverse](@ref codi::DataInterface::evaluateReverse) function for the stack evaluation. This
function will call the provided function object `SimpleTape::evaluateStack` for each valid region in the stored data.
A valid region is defined as a contiguous set of memory for all data streams that is not interrupted by any chunk boundaries.
It is therefore safe to access all memory pointers in the specified ranges. For each data stream, the range and the
pointers to the data are added to the argument list, that is, `size_t& curPos`, `size_t endPos`, `Data1* data1`, `Data2*
data2`, etc.. First the data of the root vector is appended, then the one of the nested vector and so on. Additional
arguments are added at the beginning of the argument list. The implementation of `SimpleTape::evaluateStack` shows the final
result of the argument aggregation:
\snippet developer/simpleTape.cpp Evaluation - Stack

The current position offsets are proved by reference, it is expected that they are decremented in the function and reach
the end offset. CoDiPack performs an assertion on this assumption in order to detect tape corruptions during the
evaluation.

The body of the function consists of three major parts:
 - the while loop goes through all operator codes,
 - the first switch detects if a unary or binary operation was stored and gets the data from the other two data streams,
 - the second switch implements for each operator code the associated reverse AD formulation.

### Remaining interface functions

The other functions from the codi::ReverseTapeInterface are the activity functions, the reset functions and the
statistics function. Most of them do not require a detailed discussion, only the reset of the data stream needs to be
explained. It is sufficient to call reset on the root data stream, which will recursively call `reset` on all nested
data streams.
\snippet developer/simpleTape.cpp Other - Activity
\snippet developer/simpleTape.cpp Other - Misc


### Example and analysis

The tape is now ready to use with the codi::ActiveType. For demonstration purposes the example is once evaluated with
the `SimpleTape` and once with the default `codi::RealReverse` type. The example itself consists only of a single
statement that is differentiated.
\snippet developer/simpleTape.cpp Example
Since the example is quite simple, it is possible to understand the taping process by
steping through it with the debugger.

The output and data consumption of both tapes is:
```
Simple tape:
c = 0.354971
d c/d a = 0.96017
d c/d b = -0.1455
-------------------------------------
Example tape
-------------------------------------
  Total memory used      :     200.00 B
  Total memory allocated :      16.06 KB
-------------------------------------
Adjoint vector
-------------------------------------
  Number of adjoints     :          8
  Memory allocated       :      64.00 B
-------------------------------------
Index manager
-------------------------------------
  Max. live indices      :          8
-------------------------------------
Primal data
-------------------------------------
  Total number           :          8
  Number of chunks       :          1
  Memory used            :      64.00 B
  Memory allocated       :       8.00 KB
-------------------------------------
Identifier data
-------------------------------------
  Total number           :         13
  Number of chunks       :          1
  Memory used            :      52.00 B
  Memory allocated       :       4.00 KB
-------------------------------------
Operator data
-------------------------------------
  Total number           :          5
  Number of chunks       :          1
  Memory used            :      20.00 B
  Memory allocated       :       4.00 KB
-------------------------------------



codi::RealReverse:
c = 0.354971
d c/d a = 0.96017
d c/d b = -0.1455
-------------------------------------
CoDi Tape Statistics ( JacobianLinearTape )
-------------------------------------
  Total memory used      :       96.00 B
  Total memory allocated :       28.50 MB
-------------------------------------
Adjoint vector
-------------------------------------
  Number of adjoints     :           4
  Memory allocated       :       32.00 B
-------------------------------------
Index manager
-------------------------------------
  Max. live indices      :           4
-------------------------------------
Statement entries
-------------------------------------
  Total number           :           4
  Number of chunks       :           1
  Memory used            :        4.00 B
  Memory allocated       :        2.00 MB
-------------------------------------
Jacobian entries
-------------------------------------
  Total number           :           5
  Number of chunks       :           1
  Memory used            :       60.00 B
  Memory allocated       :       24.00 MB
-------------------------------------
External function entries
-------------------------------------
  Total number           :           0
  Number of chunks       :           1
  Memory used            :        0.00 B
  Memory allocated       :        2.50 MB
-------------------------------------
```
As expected, they yield the same result, but the required memory is different. The allocated memory does not interest
us since this memory is defined by the default chunk size, that is different for the two tapes.The codi::RealReverse
tape needs 96 bytes of memory and the simple tape requires 200 bytes of memory. In order to keep it simple, no advanced
optimizations have been implemented in the SimpleTape. The data management is also quite simple and does not result in
the minimal memory footprint.
