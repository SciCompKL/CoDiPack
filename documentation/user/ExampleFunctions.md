Examples functions {#ExampleFunctions}
=======

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
