//! [Example 25 - Tape Writers]
#include <codi.hpp>
#include <iostream>

//! [Function]
template<typename Real>
Real func(const std::vector<Real>& x) {
  return x[0] * x[1] * x[0];
}
//! [Function]

//! [Generate Tape]
template<typename Real>
void generateTape(std::string const& fileName) {
  using Tape = typename Real::Tape;
  using Identifier = typename Real::Identifier;

  std::vector<Identifier> x_id;
  std::vector<Identifier> y_id;

  // Step 1: Do a normal recording
  std::vector<Real> x = {4.0, 3.0};

  Tape& tape = Real::getTape();
  tape.setActive();

  tape.registerInput(x[0]);
  tape.registerInput(x[1]);

  Real y = func(x);
  tape.registerOutput(y);

  // Step 2: Record the inputs and outputs in std::vectors
  x_id.push_back(x[0].getIdentifier());
  x_id.push_back(x[1].getIdentifier());
  y_id.push_back(y.getIdentifier());

  tape.setPassive();

  // Step 3: Write the tape to storage. Select between a text, binary or graphical file type. For primal value tapes, a
  // math representation can also be selected.

  // Text format
  tape.writeTape(codi::createWriter<Real>(fileName + "_text.txt", x_id, y_id, codi::FileType::Text));
  // Graphical format
  tape.writeTape(codi::createWriter<Real>(fileName + "_graph.dot", x_id, y_id, codi::FileType::Graph));

  // The tape can be still be evaluated as before
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << fileName << ":" << std::endl;
  std::cout << "df/dx[0](4.0) = " << x[0].getGradient() << std::endl;
  std::cout << "df/dx[1](3.0) = " << x[1].getGradient() << std::endl << std::endl;

  tape.reset();
}
//! [Generate Tape]

int main(int nargs, char** args) {
  // Example with two different tape types
  generateTape<codi::RealReverse>("jacobian_linear");
  generateTape<codi::RealReversePrimalIndex>("primal_reuse");

  return 0;
}
//! [Example 25 - Tape Writers]
