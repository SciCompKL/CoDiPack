//! [Tutorial 2 - Reverse mode AD]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealReverse;
using Tape = typename Real::Tape;

//! [Function]
Real func(const Real& x) {
  return x * x * x;
}
//! [Function]

int main(int nargs, char** args) {
  Real x = 4.0;

  Tape& tape = Real::getTape();
  tape.setActive();        // Step 1: Start recording

  tape.registerInput(x);   // Step 2: Register inputs
  Real y = func(x);        // Step 3: Call function
  tape.registerOutput(y);  // Step 4: Register outputs

  tape.setPassive();       // Step 5: Stop recording
  y.setGradient(1.0);      // Step 6: Set seeding
  tape.evaluate();         // Step 7: Perform reverse evaluation

  // Step 8: Access gradients
  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

  tape.reset();            // Step 9: Clean tape and adjoints

  return 0;
}
//! [Tutorial 2 - Reverse mode AD]
