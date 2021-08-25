
#include <iostream>
#include <codi.hpp>

int main(int nargs, char** args) {
//! [Identifier Activity]
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Tape& tape = Real::getTape();

  tape.setActive();

  Real x = 200.0;
  std::cout << "Passive x: " << x.getIdentifier() << " " << tape.isIdentifierActive(x.getIdentifier()) << std::endl;
  tape.registerInput(x);
  std::cout << "Active x: " << x.getIdentifier() << " " << tape.isIdentifierActive(x.getIdentifier()) << std::endl;
  tape.deactivateValue(x);
  std::cout << "Passive x: " << x.getIdentifier() << " " << tape.isIdentifierActive(x.getIdentifier()) << std::endl;
//! [Identifier Activity]
}
