
#include <iostream>
#include <codi.hpp>

//! [External function]
using Real = codi::RealReverse;
using Tape = typename Real::Tape;
using VAI = codi::VectorAccessInterface<typename Real::Real, typename Real::Identifier>;

void printSomething(Tape* tape, void* data, VAI* vai) {

  std::cout << "Hello from the reverse run." << std::endl;

  int index = *((int*)data);
  std::cout << "Adjoint of x is " << vai->getAdjoint(index, 0) << std::endl;
}

void deleteSomething(Tape* tape, void* data) {

  std::cout << "Hello from the cleanup crew." << std::endl;

  int* index = (int*)data;
  delete index;
}

int main(int nargs, char** args) {

  Tape& tape = Real::getTape();

  // Recording
  Real x = 10.0;

  tape.setActive();
  tape.registerInput(x);

  tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(printSomething, new int(x.getIdentifier()), deleteSomething));
  Real y = 42.0 * x * x;
  tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(printSomething, new int(x.getIdentifier()), deleteSomething));

  tape.registerOutput(y);
  tape.setPassive();

  // Reverse evaluation
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "Gradient of dy/dx: " << x.getGradient() << std::endl;

  tape.reset();
}
//! [External function]
