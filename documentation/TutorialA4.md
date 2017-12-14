Tutorial A4: Advanced vector mode usage {#TutorialA4}
============

[Tutorial 6](@ref Tutorial6) introduces the vector mode of CoDiPack.
The approach there is to change the derivative type such that the adjoint
data consists of vectors with the size \f$d\f$. The disadvantage of this
approach is that each time a new vector size is required, the whole
application needs to be recompiled.

The reverse evaluation from tutorial 6 is the following:
~~~~{.cpp}
  codi::RealReverseVec<2> xR[5];
  codi::RealReverseVec<2> yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverseVec<2>::TapeType& tape = codi::RealReverseVec<2>::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  yR[0].gradient()[0] = 1.0;
  yR[1].gradient()[1] = 1.0;
  tape.evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobiR[i][0] = xR[i].getGradient()[0];
    jacobiR[i][1] = xR[i].getGradient()[1];
  }
~~~~
A close look at the code shows that the recording part is independent of
the chosen vector size. Only after the tape is set to passive, the
chosen vector starts to matter. If the evaluation of the tape can be
generalized only the translation unit with the evaluation needs to be
recompiled on a vector size change.

The codi::TapeVectorHelper is designed to provide the required abstraction.
The above code example with the vector helper looks like:
~~~~{.cpp}
  codi::RealReverse xR[5];
  codi::RealReverse yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  codi::TapeVectorHelper<codi::RealReverse, codi::Direction<double, 2> > vh;
  vh.gradient(yR[0].getGradientData())[0] = 1.0;
  vh.gradient(yR[1].getGradientData())[1] = 1.0;
  vh.evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobiR[i][0] = vh.getGradient(xR[i].getGradientData())[0];
    jacobiR[i][1] = vh.getGradient(xR[i].getGradientData())[1];
  }
~~~~
The first major difference is that instead of using the codi::RealReverseVec
type the codi::RealReverse type is used in the recording section of the code.
After the recording is finished the vector helper is created and used to
perform the reverse evaluation. All calls are the same with the only
difference that the [getGradientData](@ref codi::ActiveReal::getGradientData)
method of the codi::ActiveReal type is used in the vector helper functions
to provide the information which adjoint variable is accessed.

The [getGradientData](@ref codi::ActiveReal::getGradientData) method returns
the internal identifier of CoDiPack for the adjoint values. The vector helper
uses the same identification to mange its own adjoint vector.

In some programs it might occur that the number of vector directions is
arbitrary and needs to be determined at runtime. This cases can be handled
by using the codi::TapeVectorHelperInterface in conjunction with the
codi::TapeVectorHelperInterface::getAdjointInterface method.
The above example with this interface looks like:
~~~~{.cpp}
  codi::TapeVectorHelperInterface<codi::RealReverse>* vh = new codi::TapeVectorHelper<codi::RealReverse, codi::Direction<double, 2> >();
  codi::AdjointInterface<codi::RealReverse::Real>* ai = vh->getAdjointInterface();

  for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
    ai->updateAdjoint(yR[dim].getGradientData(), dim, 1.0);
  }
  vh->evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
      jacobiR[i][dim] = ai->getAdjoint(xR[i].getGradientData(), dim);
    }
  }

  delete vh;
~~~~
The adjoint interface is a more general abstraction for the access of the
adjoint vector and also used by the external function mechanism to provide
access to the adjoint vector.

With the codi::TapeVectorHelper it is now possible to program a generalized
way for the evaluation of the reverse AD mode. This works for all tapes which
are implemented in CoDiPack. For the primal value tapes the code needs to
be compiled with the flag `CODI_EnableVariableAdjointInterfaceInPrimalTapes`.
This flag changes the implementation of the primal value tapes such that
they support the generalized vector evaluation.

The full code for the tutorial is:
~~~~{.cpp}
#include <codi.hpp>
#include <iostream>

template<typename Real>
void func(const Real* x, size_t l, Real* y) {
  y[0] = 0.0;
  y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}

void vectorType() {
  std::cout << "codi::RealReverse( vector type ):" << std::endl;
  // Reverse vector mode
  codi::RealReverseVec<2> xR[5];
  codi::RealReverseVec<2> yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverseVec<2>::TapeType& tape = codi::RealReverseVec<2>::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  yR[0].gradient()[0] = 1.0;
  yR[1].gradient()[1] = 1.0;
  tape.evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobiR[i][0] = xR[i].getGradient()[0];
    jacobiR[i][1] = xR[i].getGradient()[1];
  }

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }
}

void vectorHelper() {
  std::cout << "codi::RealReverse( vector helper):" << std::endl;
  // Reverse vector mode

  codi::RealReverse xR[5];
  codi::RealReverse yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  codi::TapeVectorHelper<codi::RealReverse, codi::Direction<double, 2> > vh;
  vh.gradient(yR[0].getGradientData())[0] = 1.0;
  vh.gradient(yR[1].getGradientData())[1] = 1.0;
  vh.evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobiR[i][0] = vh.getGradient(xR[i].getGradientData())[0];
    jacobiR[i][1] = vh.getGradient(xR[i].getGradientData())[1];
  }

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }
}

void vectorHelperInterface() {
  std::cout << "codi::RealReverse( vector helper interface):" << std::endl;
  // Reverse vector mode

  codi::RealReverse xR[5];
  codi::RealReverse yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  codi::TapeVectorHelperInterface<codi::RealReverse>* vh = new codi::TapeVectorHelper<codi::RealReverse, codi::Direction<double, 2> >();
  codi::AdjointInterface<codi::RealReverse::Real>* ai = vh->getAdjointInterface();

  for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
    ai->updateAdjoint(yR[dim].getGradientData(), dim, 1.0);
  }
  vh->evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
      jacobiR[i][dim] = ai->getAdjoint(xR[i].getGradientData(), dim);
    }
  }

  delete vh;

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }
}

int main(int nargs, char** args) {

  vectorType();

  codi::RealReverse::getGlobalTape().reset();
  vectorHelper();

  codi::RealReverse::getGlobalTape().reset();
  vectorHelperInterface();

  return 0;
}
~~~~
