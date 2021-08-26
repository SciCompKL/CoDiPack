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
will get a new identifier (under special circumstances, the identifier can stay the same but this is not the default). The
identifier of the variable will therefore only represent the last assignment and not all values (and therefore also derivatives)
that have been assigned to this variable. This also means that if a variable is overwritten during the computation and
the derivative with respect to the original value is required, then the identifier for this value needs to be stored.

An example is:
\snippet examples/gradientAccessTapeInterface.cpp Gradient Access

In this case x is used both as an input and output variable.

#### Statement level recording (Expression templates) {#StatementLevelRecording}

See \subpage Expressions.

#### Online activity analysis {#ActivityAnalysis}

The taping process in CoDiPack is organized such that it can automatically detect the dependency relation with respect to
the input variable. Computations that do not depend on the input variables are not recorded. This is done by giving a
special role to the zero index. This index is used to identify all _passive_ variables. This means that all values
with an nonzero index are _active_ variables.

The initial state of the program is that all variables have a zero index. A call to
[registerInput](@ref codi::ReverseTapeInterface::registerInput) will tell CoDiPack that this variable is an input and it
will receive a nonzero index. The variable is now active. When storing a statement, CoDiPack inspects all arguments and
checks if there is a least one that has a nonzero index. If this is the case, then the result of the statement is also
active, receives a nonzero index, and the statement is recorded. Otherwise, the computation in the statement does not
depend on the input variables and the result will be passive, that is, it receives the zero index and the statement is not
recorded.

With this tracking technique, only the active parts of the computations are recorded.
