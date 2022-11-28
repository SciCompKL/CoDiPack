//! [Example 21 - Special handling of linear system solvers]

#include <codi.hpp>
#include <iostream>

#if CODI_EnableEigen

using Real = codi::RealReverse;
using Tape = typename Real::Tape;
using Identifier = typename Real::Identifier;
using RealBase = typename Real::Real;

//! [Specialization of Eigen solver]
template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

//! [Function]
template<typename Type>
void func(Matrix<Type> const& A, Vector<Type> const& rhs, Vector<Type>& sol) {
  sol = A.colPivHouseholderQr().solve(rhs);
}
//! [Function]

template<typename Number>
struct EigenSolver : public codi::EigenLinearSystem<Number, Matrix, Vector> {
  public:

    using Base = codi::EigenLinearSystem<Number, Matrix, Vector>;
    using MatrixReal = typename Base::MatrixReal;
    using VectorReal = typename Base::VectorReal;

    void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
      std::cout << "Solve system says hello!!!" << std::endl;
      func(*A, *b, *x);
    }
};
//! [Specialization of Eigen solver]
#else
  #warning EIGEN_DIR not set. Skipping Eigen example.
#endif


int main(int nargs, char** args) {
#if CODI_EnableEigen
  int size = 10;

  Matrix<Real> A(size, size);
  Vector<Real> rhs(size);
  Vector<Real> sol(size);

  Tape& tape = Real::getTape();
  tape.setActive();

  Real matrixEntry = 1.0;
  Real rhsEntry = 1.0;

  tape.registerInput(matrixEntry);
  tape.registerInput(rhsEntry);

  for(int i = 0; i < size; i += 1) {
    A(i,i) = matrixEntry;
    if(i + 1 != size) {
      A(i, i + 1) = matrixEntry;
    }
    rhs(i) = rhsEntry;
  }

  std::cout << "Solving primal system:" << std::endl;
  codi::solveLinearSystem(EigenSolver<Real>(), A, rhs, sol);


  Real y = 0.0;
  for(int i = 0; i < size; i += 1) {
    y += sol(i);
  }

  tape.registerOutput(y);

  tape.setPassive();

  y.setGradient(1.0);
  std::cout << "Running reverse evaluation:" << std::endl;
  tape.evaluate();

  std::cout << "y = " << y << std::endl;
  std::cout << "dy/d matrixEntry = " << matrixEntry.getGradient() << std::endl;
  std::cout << "dy/d rhsEntry = " << rhsEntry.getGradient() << std::endl;

  tape.reset();
#else
  std::cerr << "EIGEN_DIR not set. Skipping Eigen example." << std::endl;
#endif

  return 0;
}
//! [Example 21 - Special handling of linear system solvers]
