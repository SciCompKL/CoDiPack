#define CODI_IDE 0

#include <codi.hpp>
#include <fstream>
#include <iostream>

using Real = codi::RealReverseTag;
using Tape = typename Real::Tape;

Real func(const Real& x, const Real& y) {
  return x * y;
}

static void tagLhsChangeErrorCallback(double const& currentValue, double const& newValue, void* userData) {
  std::ofstream* out = (std::ofstream*)userData;

  *out << "Wrong tag use detected '" << currentValue << "' is set to '" << newValue << "'." << std::endl;
}

static void tagErrorCallback(int const& correctTag, int const& wrongTag, bool tagError, bool useError,
                                    void* userData) {
  std::ofstream* out = (std::ofstream*)userData;

  // output default warning if no handle is defined.
  if (useError) {
    *out << "Wrong variable use detected." << std::endl;
  }
  if (tagError) {
    *out << "Wrong tag detected '" << wrongTag << "' should be '" << correctTag << "'." << std::endl;
  }
}

int main(int nargs, char** args) {

  std::ofstream out("run.out");

  Real x = 4.0;
  Real y = 3.0;
  Real z = 1.0;

  codi::PreaccumulationHelper<Real> ph;
  Tape& tape = Real::getTape();
  tape.setTagErrorCallback(tagErrorCallback, &out);
  tape.setTagLhsChangeErrorCallback(tagLhsChangeErrorCallback, &out);
  tape.setCurTag(42);
  tape.setActive();

  tape.registerInput(x);
  tape.registerInput(y);

  out << "Default test:" << std::endl;
  ph.start(x, y);
  Real w = func(x, y);
  ph.finish(false, w);
  w = w * z;

  out << "Input error test:" << std::endl;
  ph.start(x);
  w = func(x, y);
  ph.finish(false, w);
  w = w * z;

  out << "Output error test:" << std::endl;
  ph.start(x, y);
  w = func(x, y);
  ph.finish(false);
  w = w * z;

  out << "Do not use error:" << std::endl;
  tape.setTagPropertyOnVariable(x, codi::TagFlags::DoNotUse);
  w = func(x, y);
  tape.clearTagPropertiesOnVariable(x);

  out << "Do not change with same value:" << std::endl;
  tape.setTagPropertyOnVariable(w, codi::TagFlags::DoNotChange);
  w = func(x, y);

  out << "Do not change error test:" << std::endl;
  tape.setTagPropertyOnVariable(w, codi::TagFlags::DoNotChange);
  w = func(x, z);

  tape.registerOutput(y);

  tape.setPassive();
  tape.reset();

  out.close();

  return 0;
}
