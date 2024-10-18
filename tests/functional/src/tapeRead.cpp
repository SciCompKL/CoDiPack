#include <codi.hpp>

#include "../results/tapeWrite/primal_linearBinary.hpp"
#include "../results/tapeWrite/primal_linearText.hpp"
#include "../results/tapeWrite/primal_multiuseBinary.hpp"
#include "../results/tapeWrite/primal_multiuseText.hpp"
#include "../results/tapeWrite/primal_reuseBinary.hpp"
#include "../results/tapeWrite/primal_reuseText.hpp"

template<typename Real, typename Id, typename Tape>
void func(Tape& tape, std::vector<Real>& x, std::vector<Real>& y, std::vector<Id>& x_id, std::vector<Id>& y_id) {
  Real sum = 0.0;
  Real mul = 1.0;

  for (size_t i = 0; i < x.size(); i += 1) {
    x[i] = (i + 1);
    tape.registerInput(x[i]);
    x_id[i] = x[i].getIdentifier();

    y[i] = sin(x[i]);

    // t1 and t2 will always have the same index but with different meanings.
    {  // Force delete of t1 after scope.
      Real t1 = sum + y[i];
      sum += t1;
    }
    {  // Force delete of t1 after scope.
      Real t2 = mul * y[i];
      mul *= t2;
    }
  }

  for (size_t i = 0; i < x.size(); i += 1) {
    if (0 == (i % 2)) {
      y[i] += sum;
    } else {
      y[i] *= mul;
    }
  }

  for (size_t i = 0; i < x.size(); i += 1) {
    y_id[i] = y[i].getIdentifier();
    tape.registerOutput(y[i]);
  }
}

template<typename Tape, typename Id, typename Grad>
void evalTape(Tape& tape, std::vector<Id> const& x_id, std::vector<Id> const& y_id, std::vector<Grad> const& seed,
              std::vector<Grad>& grad) {
  for (size_t i = 0; i < x_id.size(); i += 1) {
    tape.gradient(y_id[i]) = seed[i];
  }

  tape.evaluate();

  for (size_t i = 0; i < x_id.size(); i += 1) {
    grad[i] = tape.gradient(x_id[i]);
  }
}

template<typename Grad>
void roundOff(Grad& N, int n) {
  Grad scale = pow(10, n);
  N = std::round(N * scale) / scale;
}

void openFile(FILE*& fileHandle, std::string const& name, std::string const& mode) {
  fileHandle = fopen(name.c_str(), mode.c_str());
  if (nullptr == fileHandle) {
    CODI_EXCEPTION("Could not open file %s", name.c_str());
  }
}

template<typename Grad>
void compare(std::vector<Grad> const& base, std::vector<Grad> const& other, FILE*& resultHandle) {
  bool hasDiff = false;
  for (size_t i = 0; i < base.size(); i += 1) {
    Grad diff = std::abs((other[i] - base[i]));
    roundOff(diff, 10);
    if (diff != 0.0) {
      fprintf(resultHandle, "diff at %zu: %0.12e (%0.12e %0.12e)\n", i, diff, base[i], other[i]);
      hasDiff = true;
    }
  }

  if (!hasDiff) {
    fprintf(resultHandle, "No differences!\n");
  }
}

template<typename Real, typename Id, typename EvalMap>
void readAndCompareTapes(std::string const& tapeDirectory, std::string const& fileName, std::vector<Id> const& x_grad,
                         std::vector<Id> const& y_grad, size_t const& n, FILE*& resultHandle,
                         EvalMap const& evalHandlesTxt, EvalMap const& evalHandlesBin) {
  std::vector<Id> x_grad_txt(n);
  std::vector<Id> x_grad_bin(n);

  if constexpr (codi::TapeTraits::isPrimalValueTape<typename Real::Tape>) {
    auto textRead = codi::readTapeFile<Real>(tapeDirectory + fileName + "Text.txt", evalHandlesTxt);
    evalTape(textRead->getTape(), textRead->getInputs(), textRead->getOutputs(), y_grad, x_grad_txt);
    fprintf(resultHandle, "Comparing base vs. txt\n");
    compare(x_grad, x_grad_txt, resultHandle);

    auto binRead = codi::readTapeFile<Real>(tapeDirectory + fileName + "Binary.dat", evalHandlesBin);
    evalTape(binRead->getTape(), binRead->getInputs(), binRead->getOutputs(), y_grad, x_grad_bin);
    fprintf(resultHandle, "Comparing base vs. binary\n");
    compare(x_grad, x_grad_bin, resultHandle);

  } else {
    auto textRead = codi::readTapeFile<Real>(tapeDirectory + fileName + "Text.txt");
    evalTape(textRead->getTape(), textRead->getInputs(), textRead->getOutputs(), y_grad, x_grad_txt);
    fprintf(resultHandle, "Comparing base vs. txt\n");
    compare(x_grad, x_grad_txt, resultHandle);

    auto binRead = codi::readTapeFile<Real>(tapeDirectory + fileName + "Binary.dat");
    evalTape(binRead->getTape(), binRead->getInputs(), binRead->getOutputs(), y_grad, x_grad_bin);
    fprintf(resultHandle, "Comparing base vs. binary\n");
    compare(x_grad, x_grad_bin, resultHandle);
  }
}

template<typename Real, typename EvalMap>
void checkResults(std::string const& tapeDirectory, std::string const& name, FILE*& resultHandle,
                  EvalMap const& evalHandlesTxt, EvalMap const& evalHandlesBin) {
  fprintf(resultHandle, "Running %s:\n", name.c_str());

  using Id = typename Real::Identifier;
  using Tape = typename Real::Tape;

  size_t const n = 5;

  std::vector<Real> x(n);
  std::vector<Real> y(n);

  std::vector<Id> x_id(n);
  std::vector<Id> y_id(n);

  std::vector<typename Real::Real> x_grad(n);
  std::vector<typename Real::Real> y_grad(n);

  Tape& tape = Real::getTape();
  tape.setActive();

  func(tape, x, y, x_id, y_id);

  tape.setPassive();

  for (size_t i = 0; i < n; i += 1) {
    y_grad[i] = 1.0;
  }

  evalTape(tape, x_id, y_id, y_grad, x_grad);

  tape.resetHard();

  readAndCompareTapes<Real>(tapeDirectory, name, x_grad, y_grad, n, resultHandle, evalHandlesTxt, evalHandlesBin);
  fprintf(resultHandle, "\n");
}

int main(int nargs, char** args) {
  std::string tapeDirectory = "../../../results/tapeWrite/";

  FILE* resultHandle;
  openFile(resultHandle, "tapeRead.out", "w");

  std::vector<typename codi::RealReverse::Tape::EvalHandle> emptyMap;

  //----Jacobian Readers----
  checkResults<codi::RealReverse>(tapeDirectory, "jacobian_linear", resultHandle, emptyMap, emptyMap);
  checkResults<codi::RealReverseIndex>(tapeDirectory, "jacobian_multiuse", resultHandle, emptyMap, emptyMap);
  checkResults<codi::RealReverseIndexGen<double, double, codi::ReuseIndexManager<int>>>(
      tapeDirectory, "jacobian_multiuse", resultHandle, emptyMap, emptyMap);

  //----Primal Value Readers----
  checkResults<codi::RealReversePrimal>(tapeDirectory, "primal_linear", resultHandle,
                                        primal_linearTextCreateEvalHandles<typename codi::RealReversePrimal::Tape>(),
                                        primal_linearBinaryCreateEvalHandles<typename codi::RealReversePrimal::Tape>());

  checkResults<codi::RealReversePrimalIndex>(
      tapeDirectory, "primal_multiuse", resultHandle,
      primal_multiuseTextCreateEvalHandles<typename codi::RealReversePrimalIndex::Tape>(),
      primal_multiuseBinaryCreateEvalHandles<typename codi::RealReversePrimalIndex::Tape>());

  checkResults<codi::RealReversePrimalIndexGen<double, double, codi::ReuseIndexManager<int>>>(
      tapeDirectory, "primal_reuse", resultHandle,
      primal_reuseTextCreateEvalHandles<
          typename codi::RealReversePrimalIndexGen<double, double, codi::ReuseIndexManager<int>>::Tape>(),
      primal_reuseBinaryCreateEvalHandles<
          typename codi::RealReversePrimalIndexGen<double, double, codi::ReuseIndexManager<int>>::Tape>());
}
