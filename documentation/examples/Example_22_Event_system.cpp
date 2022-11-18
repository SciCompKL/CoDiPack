//! [Example 22 - Event system]

#include <codi.hpp>
#include <iostream>

using ActiveType = codi::RealReverse;
using Tape = ActiveType::Tape;
using Real = Tape::Real;
using Identifier = Tape::Identifier;
using Position = Tape::Position;
using VectorAccess = codi::VectorAccessInterface<Real, Identifier>;

//! [AD Workflow callback definitions]

void onTapeStartRecording(Tape&, void*) {
  std::cout << "TapeStartRecording" << std::endl;
}

void onTapeStopRecording(Tape&, void*) {
  std::cout << "TapeStopRecording" << std::endl;
}

void onTapeRegisterInput(Tape&, Real& value, Identifier& identifier, void* customData) {
  std::cout << "TapeRegisterInput value " << value << " identifier " << identifier << std::endl;
  if (customData != nullptr) {
    std::cout << "\tcustom data " << *(int*)customData << std::endl;
  }
}

void onTapeRegisterOutput(Tape&, Real& value, Identifier& identifier, void*) {
  std::cout << "TapeRegisterOutput value " << value << " identifier " << identifier << std::endl;
}

void onTapeEvaluate(Tape&, Position const& start, Position const& end, VectorAccess*,
                    codi::EventHints::EvaluationKind, codi::EventHints::Endpoint endpoint, void*) {

  auto to_string = [](codi::EventHints::Endpoint endpoint) -> std::string {
    switch (endpoint) {
      case codi::EventHints::Endpoint::Begin:
        return "begin";
      case codi::EventHints::Endpoint::End:
        return "end";
      default:
        return "unknown";
    }
  };

  std::cout << "TapeEvaluate " << to_string(endpoint) << " from " << start << " to "
            << end << std::endl;
}

void onTapeReset(Tape&, Position const& position, codi::EventHints::Reset, bool clearAdjoints, void*) {
  std::cout << "TapeReset position " << position << " clear adjoints " << clearAdjoints
            << std::endl;
}

//! [AD Workflow callback definitions]

//! [Statement callback definitions]

void onStatementStoreOnTape(Tape&, Identifier const& lhsIdentifier, Real const& newValue, size_t numActiveVariables,
                            Identifier const* rhsIdentifiers, Real const* jacobians, void*) {
  std::cout << "StatementStoreOnTape lhsIdentifier " << lhsIdentifier << " newValue " << newValue
            << " numActiveVariables " << numActiveVariables << std::endl
            << "\t";
  for (size_t i = 0; i < numActiveVariables; ++i) {
    if (i != 0) {
      std::cout << " ";
    }
    std::cout << rhsIdentifiers[i] << " " << jacobians[i] << ";";
  }
  std::cout << std::endl;
}

void onStatementEvaluate(Tape&, Identifier const& lhsIdentifier, size_t numAdjoints, Real const* adjoints, void*) {
  std::cout << "StatementEvaluate lhsIdentifier " << lhsIdentifier << " numAdjoints " << numAdjoints << std::endl
            << "\t";
  for (size_t i = 0; i < numAdjoints; ++i) {
    if (i != 0) {
      std::cout << " ";
    }
    std::cout << adjoints[i];
  }
  std::cout << std::endl;
}

//! [Statement callback definitions]

//! [Function]
void func(const ActiveType* x, size_t l, ActiveType* y) {
  y[0] = 0.0;
  y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}
//! [Function]

int main(int nargs, char** args) {

//! [Callback registration]
  codi::EventSystem<Tape>::registerTapeStartRecordingListener(onTapeStartRecording);
  codi::EventSystem<Tape>::registerTapeStopRecordingListener(onTapeStopRecording);
  codi::EventSystem<Tape>::registerTapeRegisterInputListener(onTapeRegisterInput);

  int myCustomData = 42;
  codi::EventSystem<Tape>::registerTapeRegisterInputListener(onTapeRegisterInput, &myCustomData);

  codi::EventSystem<Tape>::registerTapeRegisterOutputListener(onTapeRegisterOutput);
  codi::EventSystem<Tape>::registerTapeEvaluateListener(onTapeEvaluate);
  codi::EventSystem<Tape>::registerTapeResetListener(onTapeReset);

  codi::EventSystem<Tape>::registerStatementStoreOnTapeListener(onStatementStoreOnTape);
  codi::EventSystem<Tape>::registerStatementEvaluateListener(onStatementEvaluate);

//! [Callback registration]

  ActiveType x[5];
  ActiveType y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  Tape& tape = ActiveType::getTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }

  func(x, 5, y);

  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();

  y[0].setGradient(1.0);
  y[1].setGradient(2.0);

  tape.evaluate();

  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) [1 2]^T = (" << x[0].getGradient() << " "
                                            << x[1].getGradient() << " "
                                            << x[2].getGradient() << " "
                                            << x[3].getGradient() << " "
                                            << x[4].getGradient() << ")" << std::endl;

  tape.reset();

  return 0;
}

//! [Example 22 - Event system]
