#if CODI_EnableMPI
//! [Example 13 - MPI communication]
#include <codi.hpp>
#include <medi/medi.hpp>
#include <iostream>
using namespace medi;

#include "codi/tools/mpi/codiMpiTypes.hpp"

using Real = codi::RealReverse;
using Tape = typename Real::Tape;
using MpiTypes = codi::CoDiMpiTypes<Real>;
MpiTypes* mpiTypes;

int main(int nargs, char** args) {
  AMPI_Init(&nargs, &args);  // Step 1: Replace all MPI_* functions and types with AMPI_*

  mpiTypes = new MpiTypes(); // Step 2: Create the CoDiPack MPI types.

  int rank;
  int size;

  AMPI_Comm_rank(AMPI_COMM_WORLD, &rank);
  AMPI_Comm_size(AMPI_COMM_WORLD, &size);

  if(size != 2) {
    std::cout << "Please start the tutorial with two processes." << std::endl;
  } else {

    Tape& tape = Real::getTape();
    tape.setActive();

    Real a = 3.0;
    if( 0 == rank ) {
      tape.registerInput(a);

      AMPI_Send(&a, 1, mpiTypes->MPI_TYPE, 1, 42, AMPI_COMM_WORLD);   // Step 3: Use the CoDiPack MPI type as the data type.
    } else {
      AMPI_Recv(&a, 1, mpiTypes->MPI_TYPE, 0, 42, AMPI_COMM_WORLD, AMPI_STATUS_IGNORE);

      tape.registerOutput(a);

      a.setGradient(100.0);
    }

    tape.setPassive();

    tape.evaluate();

    if(0 == rank) {
      std::cout << "Adjoint of 'a' on rank 0 is: " << a.getGradient() << std::endl;
    }
  }

  delete mpiTypes;           // Step 4: Cleanup the created CoDiPack MPI types.

  AMPI_Finalize();
}

#include <medi/medi.cpp>

//! [Example 13 - MPI communication]
#else
#include <iostream>

int main(int nargs, char** args) {
  std::cout << "Please compile with 'make MPI=yes MEDI_DIR=<path to medipack>' (You have to install MeDiPack, too)" << std::endl;
  return 0;
}
#endif
