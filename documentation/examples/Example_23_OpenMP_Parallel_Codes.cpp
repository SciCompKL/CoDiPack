#if CODI_EnableOpDiLib
//! [Example 23 - OpenMP Parallel Codes]
#include <codi.hpp>
#include <opdi.hpp>
#include <iostream>

#include <opdi/backend/macro/macroBackend.hpp>
#include <opdi.hpp>

using Real = codi::RealReverseIndexOpenMP;  // use a CoDiPack type suitable for OpenMP parallel applications
using Tape = typename Real::Tape;

int main(int nargs, char** args) {

  // initialize OpDiLib

  opdi::backend = new opdi::MacroBackend();
  opdi::backend->init();
  opdi::logic = new opdi::OmpLogic;
  opdi::logic->init();
  opdi::tool = new CoDiOpDiLibTool<Real>;

  // usual AD workflow in the serial parts of the code

  Real x = 4.0;

  Tape& tape = Real::getTape();
  tape.setActive();
  tape.registerInput(x);

  // parallel computation

  Real a[1000];
  Real y = 0.0;

  OPDI_PARALLEL()  // this example uses OpDiLib's macro backend, where OpenMP pragmas are replaced by such macros
  {
    OPDI_FOR()
    for (int i = 0; i < 1000; ++i)
    {
      a[i] = sin(x * i);
    }
    OPDI_END_FOR
  }
  OPDI_END_PARALLEL

  for (int i = 0; i < 1000; ++i) {
    y += a[i];
  }

  // usual AD workflow

  tape.registerOutput(y);
  tape.setPassive();
  y.setGradient(1.0);

  opdi::logic->prepareEvaluate();  // prepare OpDiLib for evaluation

  tape.evaluate();

  std::cout << "f(" << x << ") = " << y << std::endl;
  std::cout << "df/dx(" << x << ") = " << x.getGradient() << std::endl;

  // finalize OpDiLib

  opdi::backend->finalize();
  delete opdi::backend;
  delete opdi::logic;
  delete opdi::tool;

  return 0;
}

// don't forget to include the OpDiLib source file
#include <opdi.cpp>

//! [Example 23 - OpenMP Parallel Codes]
#else
#include <iostream>

int main(int nargs, char** args) {
  std::cout << "Please compile with 'make OPENMP=yes OPDILIB=yes OPDI_DIR=<path to OpDiLib>'." << std::endl;
  return 0;
}
#endif
