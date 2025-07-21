//! [Example 29 - Tape cache optimization]
#include <codi.hpp>
#include <iostream>

//! [Function]
template<typename Real>
void func(const Real* x, size_t l, Real* y) {
  y[0] = 0.0;
  y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}
//! [Function]

int main(int nargs, char** args) {

  using Real = codi::RealReverseIndex;
  using Identifier = typename Real::Identifier;
  using Tape = typename Real::Tape;

  Real x[5];
  Real y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  // Step 1: Record the tape.
  Tape& tape = Real::getTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }

  func(x, 5, y);

  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();

  // Step 2: Gather the input and output identifiers.
  Identifier xIds[5];
  Identifier yIds[2];
  for(int i = 0; i < 5; i += 1) {
    xIds[i] = x[i].getIdentifier();
  }
  for(int i = 0; i < 2; i += 1) {
    yIds[i] = y[i].getIdentifier();
  }

  // Step 3: Define the input and output iterators.
  auto iterX = [&xIds](auto&& func) {
    for(size_t i = 0; i < 5; ++i) {
      func(xIds[i]);
    }
  };
  auto iterY = [&yIds](auto&& func) {
    for(size_t i = 0; i < 2; ++i) {
      func(yIds[i]);
    }
  };

  // Step 4: Apply the optimization.
  codi::IdentifierCacheOptimizerHotCold<Tape> co{tape};
  co.eval(iterX, iterY);

  // Step 5: Do a tape evaluation with the translated ids.
  codi::Jacobian<double> jacobian(2,5);
  for(size_t curY = 0; curY < 2; curY += 1) {
    tape.gradient(yIds[curY]) = 1.0;
    tape.evaluate();

    for(size_t curX = 0; curX < 5; curX += 1) {
      jacobian(curY,curX) = tape.gradient(xIds[curX]);
      tape.gradient(xIds[curX]) = 0.0;
    }
  }

  std::cout << "Reverse Jacobian:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;

  tape.reset();

  return 0;
}
//! [Example 29 - Tape cache optimization]
