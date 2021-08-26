
#include <iostream>
#include <codi.hpp>

int main(int nargs, char** args) {
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Tape& tape = Real::getTape();
  tape.setActive();

  Real x = 10.0;
  tape.registerInput(x);

  Real y[4];

//! [Statement push helper]
  y[0] = x * x;

  codi::StatementPushHelper<Real> sh;

  // Regular use
  sh.startPushStatement();
  sh.pushArgument(x, 2.0 * x.value());
  sh.endPushStatement(y[1], x.value() * x.value());

  // Iterator based push
  std::vector<Real> valuesIter;
  std::vector<double> jacobiansIter;
  valuesIter.push_back(x);
  jacobiansIter.push_back(2.0 * x.value());

  sh.pushStatement(y[2], x.value() * x.value(), valuesIter.begin(), valuesIter.end(), jacobiansIter.begin());

  // Array based push

  Real valuesArray[] = {x};
  double jacobiansArray[] = {2.0 * x.value()};

  sh.pushStatement(y[3], x.value() * x.value(), valuesArray, jacobiansArray, 1);
//! [Statement push helper]

  for(int i = 0; i < 4; i += 1) {
    tape.registerOutput(y[i]);
  }
  tape.setPassive();

  for(int i = 0; i < 4; i += 1) {
    tape.clearAdjoints();
    y[i].setGradient(1.0);
    tape.evaluate();
    std::cout << "Gradient of dy[" << i << "]/dx: " << x.getGradient() << std::endl;
  }
}
