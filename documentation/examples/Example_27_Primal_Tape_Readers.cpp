//! [Example 27 - Primal Tape Readers]
#include <codi.hpp>
#include <iostream>

// Step 1: Include the .hpp file that was generated during the writing process
#include "../generated_files/primal_reuse_text.hpp"

int main(int nargs, char** args) {
  // Step 2: Match the Real type with the stored type
  using Real = codi::RealReversePrimalIndex;
  using Tape = typename Real::Tape;
  using Identifier = typename Real::Identifier;

  const std::string directory = "documentation/generated_files/";

  // Step 3: Get the evaluation handles from the included .hpp
  std::vector<typename Tape::EvalHandle> evalHandles = primal_reuse_textCreateEvalHandles<Tape>();

  // Step 4: Use the file name of the .txt or .dat file to restore the tape. Include the evaluation handles.
  auto reader = codi::readTapeFile<Real>(directory + "primal_reuse_text.txt", evalHandles);

  // Step 5: Get the restored tape and the IO
  Tape& tape = reader->getTape();
  std::vector<Identifier> const& x_id = reader->getInputs();
  std::vector<Identifier> const& y_id = reader->getOutputs();

  // Step 6: Seed the restored IO and evaluate the tape
  tape.gradient(y_id[0]) = 1.0;

  tape.evaluate();

  // Step 7: View the results
  std::cout << "df/dx[0] = " << tape.getGradient(x_id[0]) << std::endl;
  std::cout << "df/dx[1] = " << tape.getGradient(x_id[1]) << std::endl;

  return 0;
}
//! [Example 27 - Primal Tape Readers]
