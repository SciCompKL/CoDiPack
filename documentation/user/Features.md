Features {#Features}
=======

#### Direct data access (Identifier management) {#IdentifierManagement}

CoDiPack allows the user to have direct access to the tape and therefore also to the adjoint values. This design
decision makes it at some places a little bit more complicated to work with CoDiPack but it allows the user to use
their own data structures for handling the adjoint and tangent values.

Each variable that uses a reverse CoDiPack type (e.g. codi::RealReverse) has an associated identifier which is in the
default configuration an integer. This identifier can be obtained by calling the function
[getIdentifier](@ref codi::LhsExpressionInterface::getIdentifier()). It can then be used in the tape functions
to get the adjoint value (e.g. [getGradient](@ref codi::GradientAccessTapeInterface::getGradient())).

In CoDiPack the identifier for a variable can change. Each time the variable is assigned a new value, the variable
will get a new identifier. (Under some circumstances the identifier can stay the same but this is not the rule.) The
identifier of the variable will therefore only represent the last assignment and not all values (and therefore also derivatives)
that have been assigned to this variable. This also means, that if a variable is overwritten during the computation and
the derivative with respect to the original value is required, then the identifier for this value needs to be stored.

An example is:
\snippet examples/gradientAccessTapeInterface.cpp Gradient Access

In this case x is used both as an input and output variable.

#### Statement level recording (Expression templates) {#ExpressionTemplates}

CoDiPack uses expression template to do a statement level recording approach for all of its tapes. The statement
```
  w = sqrt((a + b) * (c - d));
```
would be split into four statements by an operator overloading approach. That is:
```
 t1 = a + b;
 t2 = c - d;
 t3 = t1 * t2;
 w = sqrt(t3);
```
The reverse mode of AD requires a certain amount of memory per statement and per argument to that statement. For an
Jacobian taping approach that uses an operator level taping approach this is 4.125 byte (int + boolean) per statement and
12 byte (int + double) per argument. For the above example this would be a total of 100.5 byte (4 statements,
7 arguments). For an Jacobian taping approach that uses statement level taping the memory requirements are 7 byte
(int + uint8) per statement and 12 byte (int + double) per argument. For the above example this would be a total of
55 byte (1 statement, 4 arguments). The switch from an operator level taping approach to a statement level taping
approach gives for the above example a memory reduction by 45%.

Since not all statements contain multiple operations, the total memory gain is usually lower but still significant.
On the other side, if the program is written such that only one operation is done per statement, then the statement
level taping will increase the memory slightly. Nevertheless, the statement level taping approach is usually benifical.

Expression template are the key to enable the AD tool implementation to access the whole statement instead of just
seeing single operators. The basic idea is to change the return type of an operator. Let
\f$\text{Real} \circ \text{Real} \rightarrow \text{Real}\f$ be the default definition of a binary operator, where
\f$\text{Real}\f$ is the basic computation type (e.g. double). The layout of the operator is now changed to
\f$\text{Expr}_A \circ \text{Expr}_B \rightarrow \text{Expr}_{A \circ B}\f$ where \f$\text{Expr}\f$ is a general
expression class. The definition of the operator is now such that it expects two expressions as arguments. The subscripts
\f$A\f$ and \f$B\f$ represent the operations evaluated by these expressions. The return value of the operator is now
not the value, but a new expression that incorporates the operation of the operator. The implementation of the
expressions is done such that for each operator an associated class is implemented. This class knows how to
evaluate the operation and captures the arguments. It is then used as the return value of the operator implementation.
For the example above the generated expression would look like
```
SQRT<MUL<ADD<Real, Real>, SUB<Real, Real>>>
```
where `Real` stands for an actual argument to the statement. The assign operation on the `Real` implementation will
now see the whole expression and can use it to store the appropriate data for the reverse AD mode.


#### Online activity analysis {#ActivityAnalysis}

The taping process in CoDiPack is organized such that it can automatically detect the dependency relation with respect to
the input variable. Computations, that do not depend on the input variables, are not recorded. This is done by giving a
special role to the zero index. This index is used to identify all _passive_ variables. This means that all values
with an non zero index are _active_ variables.

The initial state of the program is that all variables have a zero index. A call to
[registerInput](@ref codi::ReverseTapeInterface::registerInput) will tell CoDiPack that this variable is an input an it
will gain a non zero index. The variable is now active. When storing a statement, CoDiPack inspects all arguments and
checks if there is a least one that has a non zero index. If this is the case, then the result of the statement is also
active, gains a non zero index, and the statement is recorded. Otherwise, the computation in the statement does not
depend on the input variables and the result will be passive, that is it receives a zero index and the statement is not
recorded.

With this tracking technique only the active parts of the computations are recorded.
