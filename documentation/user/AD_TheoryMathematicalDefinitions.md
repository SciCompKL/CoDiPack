AD theory and mathematical definitions {#AD_TheoryMathematicalDefinitions}
=======

On this page we want to give a brief introduction into the algorithmic differentiation (AD) theory. For a detailed introduction please see the
book [_Evaluating Derivatives_](https://doi.org/10.1137/1.9780898717761) from Griewank and Walther or
[_The Art of Differentiating Computer Programs: An Introduction to Algorithmic Differentiation_](https://doi.org/10.1137/1.9781611972078) from Naumann.

With AD the derivative of computer programs can be computed. It can either be applied to the full program or to a
part of the program. Therefore we define a general function
 \f[ y = f(x) \f]
that represents the computer program which should be differentiated.
In the definition, \f$ x \in \R^n \f$ is the input vector with the size \f$ n \f$ and \f$y \in \R^m \f$ is the output vector
with the size \f$ m \f$.

Forward AD Equation {#sec_forwardAD}
-------
The forward AD mode computes now the directional derivative of \f$f\f$ in the direction \f$\dot x\f$. The equation reads
\f[ \dot y = \frac{df}{dx} \dot x \f]
where \f$\dot y\f$ represents the directional derivative. This equation describes what AD will compute but it does not
describe how it is done (the matrix \f$df/dx\f$ is not setup directly).

The _how_ can be explained by looking at the program evaluation on the CPU. During the compilation process and
evaluation, the program is broken down in so called elemental operations \f$\phi_i : \R^{n_i} \rightarrow \R \f$,
\f$ w = \phi_i(u) \f$ where \f$ u \f$ is a vector that contains the \f$ n_i \f$ arguments of the operation. Usually these are the binary operations like \f$+\f$, \f$*\f$, etc. and unary functions like
\f$sin\f$, \f$exp\f$, etc.. For all of these elemental operations, the derivative is known and the program can be seen
as big concatenated evaluation of these elemental operations. By applying the chain rule and the directional derivative
to the chain of elemental operations, the forward mode AD theory is established. The result is that alongside each
elemental operation \f$\phi_i\f$ the forward AD equation \f[ \dot w = \frac{d\phi_i}{du} \dot u \f] needs to be
evaluated.

The computation is done such that for each primal variable (e.g. \f$w\f$) a tangent variable (e.g. \f$\dot w\f$) is also
allocated. The user can then set the seeding of the tangent variables before the program is run. During the program
evaluation, for each elemental operation the corresponding forward AD equation is evaluated. After the program evaluation
is finished, the user can get the directional derivative by accessing the tangent values of the outputs.

Reverse AD Equation {#sec_reverseAD}
-------

The reverse AD mode computes the adjoint directional derivative of \f$f\f$ in the adjoint direction \f$\bar y\f$. The
equation reads
\f[ \bar x = \frac{df}{dx}^T \bar y \f]
where \f$\bar x\f$ represents the adjoint directional derivative. This equation describes what AD will compute but, again, it
does not describe how it is done (the matrix \f$df/dx^T\f$ is not setup directly).

The _how_ can be explained by the identity \f$ \scalar{\dot y}{\bar y} = \scalar{\frac{df}{dx} \dot x}{\bar y} = \scalar{\dot x}{\frac{df}{dx}^T \bar y} = \scalar{\dot x}{\bar x}\f$.
It describes that the reverse AD mode is just the discrete adjoint of the forward AD mode. How the discrete adjoint
looks is derived in the books above. The result is that for each elemental operation a slightly modified version  of the
reverse AD equation
\f[ \bar u \aeq \frac{d\phi_i}{du}^T \bar w; \quad \bar w = 0 \f]
needs to be evaluated. First the adjoint variables of the inputs are updated, afterwards the adjoint value of the output
is reset to zero. This equation can no longer be evaluated alongside the primal computation. It has to be evaluated for
every elemental operation in reverse order, that is from \f$k\f$ to \f$1\f$ where \f$k \in \N\f$ is the last elemental
operation of the primal program.

Since the elemental operations need to be evaluated in reverse order, for each operation the necessary information for
the reversal needs to be stored. Therefore, the computation is done such that the user has to declare (register) all
input variables. Based on this declaration the AD tool records the necessary information. After the computation is
finished the user can set the seeding on \f$\bar y\f$ and start the reverse interpretation of the recorded data. This
will populate the adjoint sensitivity information in \f$\bar x\f$.

Mathematical naming conventions {#sec_namingConventions}
-------

The Jacobian of \f$ f \f$ is defined by:
 \f[ J = \frac{\d f}{\d x} \in \R^{m \times n} \f]
The number of rows (\f$ m \f$) represents the number of output variables and the number of columns (\f$ n \f$)
represents the number of input variable. The derivative for the i-th output with respect to the j-th input is
represented by \f$ J_{i,j} \f$.

The Hessian of \f$ f \f$ is defined by:
 \f[ H = \frac{\d^2 f}{\d^2 x} \in \R^{m \times n \times n} \f]
The first dimension (\f$ m \f$) represents the number of output variables, the second and third dimension (\f$ n \f$) represents the
first and second derivative with respect to the input variables.
The second derivative for the i-th output with respect to the j-th and k-th input is
represented by \f$ H_{i,j,k} \f$.
