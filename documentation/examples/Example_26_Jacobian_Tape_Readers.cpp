//! [Example 26 - Jacobian Tape Readers]
#include <codi.hpp>
#include <iostream>

int main(int nargs, char** args) {
  // Step 1: Match the Real type with the stored type
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;
  using Identifier = typename Real::Identifier;

  const std::string directory = "documentation/generated_files/";

  // Step 2: Use the file name of the .txt or .dat file to restore the tape
  auto reader = codi::readTapeFile<Real>(directory + "jacobian_linear_text.txt");

  // Step 3: Get the restored tape and the IO
  Tape& tape = reader->getTape();
  std::vector<Identifier> const& x_id = reader->getInputs();
  std::vector<Identifier> const& y_id = reader->getOutputs();

  // Step 4: Seed the restored IO and evaluate the tape
  tape.gradient(y_id[0]) = 1.0;

  tape.evaluate();

  // Step 5: View the results
  std::cout << "df/dx[0] = " << tape.getGradient(x_id[0]) << std::endl;
  std::cout << "df/dx[1] = " << tape.getGradient(x_id[1]) << std::endl;

  return 0;
}
//! [Example 26 - Jacobian Tape Readers]
