Tutorial 6 - Higher order derivatives {#Tutorial_06_Higher_order_types_helper_acces}
============

**Goal:** Learn how to define higher order types and how to compute higher order derivatives.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:** \ref func_simple1to1_higher
\snippet tutorials/Tutorial_06_Higher_order_types_helper_access.cpp Function

**Full code:**
\snippet tutorials/Tutorial_06_Higher_order_types_helper_access.cpp Tutorial 6 - Higher order derivatives

#### Defining higher order types ####

Higher order derivative types in CoDiPack can be defined by nesting the provided CoDiPack types. The simplest on is the
nesting of forward types:
~~~~{.cpp}
using t1s = codi::RealForwardGen<double>;
using t2s = codi::RealForwardGen<t1s>;
using t3s = codi::RealForwardGen<t2s>;
using t4s = codi::RealForwardGen<t3s>;
using t5s = codi::RealForwardGen<t4s>;
using t6s = codi::RealForwardGen<t5s>;
~~~~
The types are created by using the general types where the user can specify the type of the value and gradient. The
nesting is done in the
example up to the 6-th order. This nesting is called forward-over-forward since the forward AD mode is applied on the
forward AD mode.

The reverse types can also be used to create higher oder types. The general recomendation is to create first a nesing
with forward mode types and then apply once the reverse type:
~~~~{.cpp}
using r6s = codi::RealReverseGen<t5s>;
~~~~
Here, a 6-th order type is constructed by applying five times the forward mode and then once the reverse mode. This
nesting is called forward-over-reverse.

With these two techniques arbitray higher order types can be constructed. There are also the two other posibilites to
use a reverse-over-forward nesting or a reverse-over-reverse nesting. The first one is mathematical
identical to the forward-over-reverse approach but it would create a very large reverse tape since all higher order
computations are recorded on the tape. The second approach is theoretical possible but the seccond application of the
reversal transforms the reverse code into a _forward_ derivative code. Therefore it is more appropritate to use a
forward mode type in the first place.

#### Setting and getting derivatives ####

In order to know, which directions need to be set to get e.g. second order derivatives, the forward AD mode equation
\ref sec_forwardAD
needs to be modified. (We use the notations from Uwe Naumans book "The Art of Differentiating Computer Programs".) Each
application of the forward mode introduces a superscript with the number of the application. For the first application
this is
\f[
  \dot w^{(1)} = d\phi/du * \dot u^{(1)} \eqdot
\f]
The second application gets now the superscript \f$\cdot^{(2)}\f$. Two superscirpts of the same mode are merged
together. The total number of entries shows the derivative order of the value. The second application of the forward
mode to the above equation yields (including all primal and derivative equations)
\f{align*}{
  w =& \phi(u) \\
  \dot w^{(1)} =& \frac{\d \phi}{\d u}(u) \dot u^{(1)} \\
  \dot w^{(2)} =& \frac{\d \phi}{\d u}(u) \dot u^{(2)} \\
  \dot w^{(1,2)} =& \frac{\d^2 \phi}{\d^2 u}(u) \dot u^{(1)} \dot u^{(2)} + \frac{\d \phi}{\d u}(u) \dot u^{(1,2)}\eqdot
\f}

From these equations we learn, that all first order tangent directions \f$\dot u^{(1)}\f$ and \f$\dot u^{(2)}\f$ need to
be set in order to get the second order derivative \f$\frac{\d^2 \phi}{\d^2 u}(u)\f$. The second order tangent direction
\f$\dot u^{(1,2)}\f$ needs to be zero, otherwise the result will contain additional information.

This can be extended to arbitrary derivative orders. All first order derivatives need to be seeded in order to
compute the highest order derivative. For forward-over-reverse higher order types this is also true, but here n-1
direction need to be set during the recording and 1 direction during the reversal.

This tutorial uses the [DerivativeAccess](@ref codi::DerivativeAccess) helper for the management of the derivative
directions. It provides convinience functions that allow to set all derivatives on a specific order. For a compile time
interface see the example \ref Example_04_Higher_order_types_helper_access_compile_time. Example
\ref Example_05_Higher_order_types_direct_access shows how the gradients can be accessd without the helper.
