Tutorial B2: Tape helper {#TutorialB2}
============

This tutorial describes the basic use of the codi::TapeHelper which gives a simpler access to the CoDiPack tapes.
The user does not need to care about which variables are registered as inputs or outputs and does not need to remember
them. Also, all the tape management - like resets or adjoint clears - is handled by the codi::TapeHelper.

We want to differentiate a function that computes the angle between two vectors \f$a\f$ and \f$b\f$.
Mathematical this is done by computing the dot product of the normalized vectors and taking the arcus cosine function of
the result:
\f[
  \alpha = f(a, b) = \arccos\left(\frac{\scalar{a}{b}}{\norm{a} \norm{b}}\right)
\f]

This function is implemented such that the angle and the two norms of the vectors are returned:
~~~~{.cpp}
const size_t n = 10;

void dotWithNorms(double const* a, double const* b, size_t n, double& alpha, double& aNorm, double& bNorm) {
  alpha = double(); // Dot product is accumulated in alpha
  aNorm = double();
  bNorm = double();

  for(size_t i = 0; i < n; i += 1) {
    alpha += a[i] * b[i];
    aNorm += a[i] * a[i];
    bNorm += b[i] * b[i];
  }

  aNorm = sqrt(aNorm);
  bNorm = sqrt(bNorm);
  alpha = acos(alpha / (aNorm * bNorm));
}
~~~~

In order to record the derivative for this function we need to create a codi::TapeHelper instance. The only difficulty
here is the choice of the CoDiPack type. By default all `Reverse` CoDiPack types in codi.hpp are supported by the helper.
Some features are not available for some types which will be explained later. codi.hpp provides two default type definitions
codi::JacobianComputationType and codi::HessianComputationType which can also be used in the helper class. The first one
supports Jacobian computations and the second one supports Hessian and Jacobian computations. So a basic setup of a
codi::TapeHelper is:
~~~~{.cpp}
using TH = codi::TapeHelper<codi::HessianComputationType>;
TH th;
~~~~
which will create the tape helper. For the recording of the tape the helper provides four functions:
[startRecording](@ref codi::TapeHelper::startRecording),
[registerInput](@ref codi::TapeHelper::registerInput)
[registerOutput](@ref codi::TapeHelper::registerOutput) and
[stopRecording](@ref codi::TapeHelper::stopRecording). `startRecording` and `stopRecording` define the borders for the
code that will be recorded on the tape. `registerInput` and `registerOutput` will tell CoDiPack which are independent and
dependent variables of the source code that is being taped. In our case the recording of `dotWithNorms` will look
like:
~~~~{.cpp}
const size_t n = 10;

std::vector<codi::HessianComputationType> a(n);
std::vector<codi::HessianComputationType> b(n);
for(size_t i = 0; i < n; i += 1) {
  a[i] = i;
  b[i] = pow(-1, i);
}

th.startRecording();
for(size_t i = 0; i < n; i += 1) {
  th.registerInput(a[i]);
}
for(size_t i = 0; i < n; i += 1) {
  th.registerInput(b[i]);
}

codi::HessianComputationType alpah, aNorm, bNorm;
dotWithNorms(a.data(), b.data(), n, alpha, aNorm, bNorm);

th.registerOutput(alpha);
th.registerOutput(aNorm);
th.registerOutput(bNorm);

th.stopRecording();
~~~~
The initialization of the vectors `a` and `b` can be inside the recorded section or outside. Since, CoDiPack does not
track these initializations, it makes not difference. The more important matter is the order in which the values of `a`
and `b` are registered as inputs. The codi::TapeHelper class identifies all variables in the gradient vectors, Jacobians
and Hessians by the order in which they are registered on the tape. Therefore, the first variable corresponds to the
first entry, the second to the second etc.. For the output values this is the same. In the code above, vector `a` is
identified by the first 10 entries in the input vector which are the indices 0 to 9. For `b` the indices 10 to 19 are
used. We could also have registered the vectors `a` and `b` directly after there values have been initialized. Then the
values of `a` and `b` would have been interleaved in the CoDiPack representation. Vector `a` taking every even index and `b`
every uneven index.

After the tape is recorded, the next step is the evaluation of the Hessians, Jacobians and gradients. For the evaluation
of these objects the codi::TapeHelper provides several functions:
[evalJacobian](@ref codi::TapeHelper::evalJacobian),
[evalHessian](@ref codi::TapeHelper::evalHessian),
[evalForward](@ref codi::TapeHelper::evalForward) and
[evalReverse](@ref codi::TapeHelper::evalReverse).
We will first have a look at `evalJacobian` and `evalHessian` and then look at `evalForward` and `evalReverse`.

The Jacobian and Hessian evaluation function will use the recorded information to compute the respective derivative
information of the recorded evaluation. If the recorded evaluation can be described by the mathematical function
  \f[ F: \R^n \rightarrow \R^m \f]
then `evalJacobian` will compute
  \f[ J = \frac{\d F}{\d x} \in \R^{m \times n} \f]
and `evalHessian` will compute
  \f[ H = \frac{\d^2 F}{\d^2 x} \in \R^{m \times n \times n} \eqdot \f]
The helper will choose the best way to evaluate the desired object. The code for evaluating the Hessian and Jacobian
will look like:
~~~~{.cpp}
TH::JacobianType& jac = th.createJacobian();
TH::HessianType& hes = th.createHessian();

th.evalJacobian(jac);
th.evalHessian(hes);

// use values

th.deleteJacobian(jac);
th.deleteHessian(hes);
~~~~

Since the first order derivatives for the Jacobian are also computed during the hessian evaluation, it is also possible
to evaluate both at the same time with `th.evalHessian(hes, jac);`. If we would like to compute the Hessian and Jacobian
at a different position, there are two options. The first option is to use the `eval...At` methods of the
codi::TapeHelper. If the CoDiPack type supports this reevaluation of the tape, it is usually the fastest and simplest
method. The second option is to use the Helper structure for a second tape recording with the new values. The first
approach is supported only by `RealReversePrimal...` types of CoDiPack. The second approach is supported by all CoDiPack
reverse types.

In our case the first option with the reevaluation would look like:
~~~~{.cpp}
TH::Real* x = th.createPrimalVectorInput();
TH::Real* y = th.createPrimalVectorOutput();

for(size_t i = 0; i < n; i += 1) {
  x[0 + i] = i * i;
  x[n + i] = pow(-1, i + 1);
}

th.evalHessianAt(x, hes, y, jac);

// use values

th.deletePrimalVector(x);
th.deletePrimalVector(y);
~~~~
For the reevaluation we have to create an input vector `x` which holds the values for the the new position where the
tape is evaluated. Otherwise the call to `evalHessianAt` is the same as `evalHessian`, the difference is only that we
provide a new input position with `x` and also retrieve the new output values with `y`.

The second option requires the same calls as for the recording and therefore we will not repeat it here.

The codi::TapeHelper will always choose the most appropriate way to evaluate the Jacobian or Hessian. It also allows, with
the [evalForward](@ref codi::TapeHelper::evalForward) and [evalReverse](@ref codi::TapeHelper::evalReverse) methods, to
directly perform a forward mode or reverse mode AD evaluation. `evalForward` will compute the AD forward mode which
computes
\f[ \dot y = \frac{\d F}{\d x}(x)\dot x \eqdot \f]
`evalReverse` will compute the AD reverse mode which computes
\f[ \bar x = \frac{\d F}{\d x}^T(x)\bar y \eqdot \f]
An example reverse mode evaluation would look like:
~~~~{.cpp}
TH::GradientValue* x_b = th.createGradientVectorInput();
TH::GradientValue* y_b = th.createGradientVectorOutput();

y_b[0] = {1.0, 0.0, 0.0, 0.0};
y_b[1] = {0.0, 1.0, 0.0, 0.0};
y_b[2] = {0.0, 0.0, 1.0, 0.0};

th.evalReverse(y_b, x_b);

// use values

th.deleteGradientVector(x_b);
th.deleteGradientVector(y_b);
~~~~
Since the default [Hessian  computation type](@ref codi::HessianComputationType) uses a vector mode of size four, it is possible to
compute the derivative with respect to all three output values simultaneously. In the example, each output variable is
seeded with a different unit vector. The result in `x_b` contains then the derivative in the same dimension
as it is chosen for the output value. For a more indep tutorial about the vector mode please see [tutorial 6](@ref Tutorial6).
There are also definitions of the codi::HessianComputationType and codi::JacobianComputationType for scalar values
available: codi::HessianComputationScalarType, codi::JacobianComputationScalarType.

As already seen in the examples above, the codi::TapeHelper provides several ease of access functions that can be used
to create all vectors and derivative objects which are used in the evaluation functions. For a proper release of the
acquired resources from the `create` method, the matching `delete` method needs to be called on the same object. Each
call of the `create` methods will create a new objects, such that e.g. multiple Jacobians can be used at the same time.
All creation methods are:
[createGradientVectorInput](@ref codi::TapeHelper::createGradientVectorInput),
[createGradientVectorOutput](@ref codi::TapeHelper::createGradientVectorOutput),
[createGradientPrimalInput](@ref codi::TapeHelper::createPrimalVectorInput),
[createGradientPrimalOutput](@ref codi::TapeHelper::createPrimalVectorOutput),
[createJacobian](@ref codi::TapeHelper::createJacobian) and
[createHessian](@ref codi::TapeHelper::createHessian).
All deletion methods are:
[deleteGradientVector](@ref codi::TapeHelper::deleteGradientVector),
[deleteGradientPrimal](@ref codi::TapeHelper::deletePrimalVector),
[deleteJacobian](@ref codi::TapeHelper::deleteJacobian) and
[deleteHessian](@ref codi::TapeHelper::deleteHessian).
There are two additional functions to retrive the number of inputs and outputs, they are
[getInputSize](@ref codi::TapeHelper::getInputSize) and [getOutputSize](@ref codi::TapeHelper::getOutputSize).

The full code for the tutorial is:
~~~~{.cpp}
#include <codi.hpp>

#include <iostream>

template<typename Real>
void dotWithNorms(Real const* a, Real const* b, size_t n, Real& alpha, Real& aNorm, Real& bNorm) {
  alpha = Real(); // Dot product is accumulated in alpha
  aNorm = Real();
  bNorm = Real();

  for(size_t i = 0; i < n; i += 1) {
    alpha += a[i] * b[i];
    aNorm += a[i] * a[i];
    bNorm += b[i] * b[i];
  }

  aNorm = sqrt(aNorm);
  bNorm = sqrt(bNorm);
  alpha = acos(alpha / (aNorm * bNorm));
}

template<typename Vec>
void printVector(std::string const& name, Vec const& v, size_t length, size_t offset) {
  std::cout << "Vector " << name << ": {";
  for(size_t i = 0; i < length; i += 1) {
    if(i != 0) {
      std::cout << ", ";
    }
    std::cout << v[offset + i];
  }
  std::cout << "}" << std::endl;
}

template<typename Vec>
void printVectorDim(std::string const& name, Vec const& v, size_t length, size_t offset, size_t dim) {
  std::cout << "Vector " << name << ": {";
  for(size_t i = 0; i < length; i += 1) {
    if(i != 0) {
      std::cout << ", ";
    }
    std::cout << v[offset + i][dim];
  }
  std::cout << "}" << std::endl;
}


template<typename Jac>
void printJacCol(std::string const& text, Jac const &jac, size_t col) {
  std::cout << text <<": {";
  for(size_t j = 0; j < jac.getN(); j += 1) {
    if(j != 0) {
      std::cout << ", ";
    }
    std::cout << jac(col, j);
  }
  std::cout << "}" << std::endl;
}

template<typename Hes>
void printHesForOutput(std::string const& text, Hes const &hes, size_t output) {
  std::cout << text <<": {\n";
  for(size_t j = 0; j < hes.getN(); j += 1) {
    std::cout << "  ";
    for(size_t k = 0; k < hes.getN(); k += 1) {
      if(k != 0) {
        std::cout << ", ";
      }
      std::cout << hes(output, j, k);
    }
    std::cout << "\n";
  }
  std::cout << "}" << std::endl;
}

int main(int nargs, char** args) {

  int mode = 1;

  if(2 <= nargs) {
    mode = std::stoi(args[1]);

    if(mode < 1 || 2 < mode) {
      std::cerr << "Error: Please enter a mode from 1 to 2, it was '" << mode << "'." << std::endl;
      std::cerr << "  Mode  1: separate evaluation of Hessian and Jacobian" << std::endl;
      std::cerr << "  Mode  2: combined evaluation of Hessian and Jacobian" << std::endl;

      exit(-1);
    }
  }

  using TH = codi::TapeHelper<codi::HessianComputationType>;
  TH th;

  const size_t n = 10;

  std::vector<codi::HessianComputationType> a(n);
  std::vector<codi::HessianComputationType> b(n);
  for(size_t i = 0; i < n; i += 1) {
    a[i] = i;
    b[i] = pow(-1, i);
  }

  th.startRecording();
  for(size_t i = 0; i < n; i += 1) {
    th.registerInput(a[i]);
  }
  for(size_t i = 0; i < n; i += 1) {
    th.registerInput(b[i]);
  }

  codi::HessianComputationType alpha, aNorm, bNorm;
  dotWithNorms(a.data(), b.data(), n, alpha, aNorm, bNorm);

  th.registerOutput(alpha);
  th.registerOutput(aNorm);
  th.registerOutput(bNorm);

  th.stopRecording();

  TH::JacobianType& jac = th.createJacobian();
  TH::HessianType& hes = th.createHessian();

  if(1 == mode) {
    th.evalJacobian(jac);
    th.evalHessian(hes);
  } else {
    th.evalHessian(hes, jac);
  }

  printVector("a", a, n, 0);
  printVector("b", b, n, 0);
  std::cout << std::endl;
  printJacCol("Jacobian with respect to alpha: ", jac, 0);
  printJacCol("Jacobian with respect to aNorm: ", jac, 1);
  printJacCol("Jacobian with respect to bNorm: ", jac, 2);
  std::cout << std::endl;
  printHesForOutput("Hessian with respect to alpha: ", hes, 0);
  printHesForOutput("Hessian with respect to aNorm: ", hes, 1);
  printHesForOutput("Hessian with respect to bNorm: ", hes, 2);


  // Evaluate at different position
  TH::Real* x = th.createPrimalVectorInput();
  TH::Real* y = th.createPrimalVectorOutput();

  for(size_t i = 0; i < n; i += 1) {
    x[0 + i] = i * i;
    x[n + i] = pow(-1, i + 1);
  }

  if(1 == mode) {
    th.evalJacobianAt(x, jac, y);
    // Jacobian evaluation already shifted the point for the evaluation. No second ...At call is required here
    th.evalHessian(hes);
  } else {
    th.evalHessianAt(x, hes, y, jac);
  }

  printVector("a", a, n, 0);
  printVector("b", b, n, 0);
  std::cout << std::endl;
  printJacCol("Jacobian with respect to alpha: ", jac, 0);
  printJacCol("Jacobian with respect to aNorm: ", jac, 1);
  printJacCol("Jacobian with respect to bNorm: ", jac, 2);
  std::cout << std::endl;
  printHesForOutput("Hessian with respect to alpha: ", hes, 0);
  printHesForOutput("Hessian with respect to aNorm: ", hes, 1);
  printHesForOutput("Hessian with respect to bNorm: ", hes, 2);

  // Evaluate gradient
  TH::GradientValue* x_b = th.createGradientVectorInput();
  TH::GradientValue* y_b = th.createGradientVectorOutput();

  y_b[0] = {1.0, 0.0, 0.0, 0.0};
  y_b[1] = {0.0, 1.0, 0.0, 0.0};
  y_b[2] = {0.0, 0.0, 1.0, 0.0};

  th.evalReverse(y_b, x_b);

  std::cout << "Reverse evaluation for alpha_b:" << std::endl;
  printVectorDim("a_b", x_b, n, 0, 0);
  printVectorDim("b_b", x_b, n, n, 0);
  std::cout << std::endl;
  std::cout << "Reverse evaluation for aNorm_b:" << std::endl;
  printVectorDim("a_b", x_b, n, 0, 1);
  printVectorDim("b_b", x_b, n, n, 1);
  std::cout << std::endl;
  std::cout << "Reverse evaluation for bNorm_b:" << std::endl;
  printVectorDim("a_b", x_b, n, 0, 2);
  printVectorDim("b_b", x_b, n, n, 2);

  // Clean up vectors
  th.deleteGradientVector(x_b);
  th.deleteGradientVector(y_b);

  th.deletePrimalVector(x);
  th.deletePrimalVector(y);

  th.deleteJacobian(jac);
  th.deleteHessian(hes);

  return 0;
}
~~~~
