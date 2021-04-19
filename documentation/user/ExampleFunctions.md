Examples functions {#ExampleFunctions}
=======

These are functions that are used in tutorials and examples.

#### Simple real valued function #### {#func_simple1to1}
Definition:
\f[
  y = f(x) = x * x * x = x^3
\f]
Jacobian:
\f[
  \frac{\d f}{\d x}(x) = 3 * x * x = 3 * x^2
\f]

#### Simple vector valued function #### {#func_simpleNto2}
\f[
  y = f(x) = \left(\sum_{i = 1}^n x_i, \prod_{i = 1}^n x_i \right)
\f]
Jacobian:
\f[
  \frac{\d f}{\d x_j}(x) = \left( 1.0, \prod_{i = 1 \ldots n \wedge i \not = j } x_i \right)
\f]

#### Simple real valued function for higher order derivatives #### {#func_simple1to1_higher}
\f[
  y = f(x) = 3*x^7 \eqdot
\f]
Jacobian:
\f[
  \frac{\d f}{\d x}(x) = 21 * x^6, \quad
\f]
Higher order derivatives:
\f[
  \frac{\d^2 f}{\d^2 x}(x) = 126 * x^5, \quad
  \frac{\d^3 f}{\d^3 x}(x) = 630 * x^4, \quad
  \frac{\d^4 f}{\d^4 x}(x) = 2520 * x^3, \quad
  \frac{\d^5 f}{\d^5 x}(x) = 7560 * x^2, \quad
  \frac{\d^6 f}{\d^6 x}(x) = 15120 * x \eqdot
\f]

#### Linear system solve #### {#func_linearSystemSolve}
\f[
  x = A^{-1}b
\f]
Forward mode:
\f[
  \dot x = A^{-1}(\dot b - \dot Ax)
\f]
Reverse mode:
\f[
  \begin{aligned}
      s = & A^{-T}\bar x\\
      \bar A \aeq & -s \cdot x^T \\
      \bar b \aeq & s \\
      \bar x = & 0 \eqdot
  \end{aligned}
\f]

#### 1D polynomial #### {#func_poly1D}
\f[
  w = f(x) = 3x^4 + 5x^3 - 3x^2 + 2x -4
\f]

Derivatives:
\f[
  \frac{\d f}{\d x}(x) = 12x^3 + 15x^2 - 6x + 2
\f]

#### 2D polynomial #### {#func_poly2D}
\f[
  w = f(x, y) = \sum_{i = 1 \ldots n, j = 1 \ldots n} A_{i,j} x^{i - 1} y^{j - 1}
\f]

Derivatives:
\f[
  \frac{\d f}{\d x}(x, y) = \sum_{i = 2 \ldots n, j = 1 \ldots n} (i - 1) A_{i,j} x^{i - 2}  y^{j - 1}
\f]
\f[
  \frac{\d f}{\d y}(x, y) = \sum_{i = 1 \ldots n, j = 2 \ldots n} (j - 1) A_{i,j} x^{i - 1}  y^{j - 2}
\f]
