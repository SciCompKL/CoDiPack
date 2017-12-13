#include <codi.hpp>

#include <iostream>

template<typename Real>
void solve2(const Real* A, const Real* b, Real* x) {

  // A = a[0] a[1]  A^-1 = 1/det *  a[3] -a[1]
  //     a[2] a[3]                 -a[2]  a[0]
  Real det = A[0] * A[3] - A[1] * A[2];

  x[0] = (A[3] * b[0] - A[1] * b[1]) / det;
  x[1] = (-A[2] * b[0] + A[0] * b[1]) / det;
}

void solve2_primal(const codi::RealReverse::Real* x, size_t m, codi::RealReverse::Real* y, size_t n, codi::DataStore* d) {
  CODI_UNUSED(m);
  CODI_UNUSED(n);
  CODI_UNUSED(d);

  solve2(&x[0], &x[4], y);
}

void solve2_rev(const codi::RealReverse::Real* x, codi::RealReverse::Real* x_b, size_t m, const codi::RealReverse::Real* y, const codi::RealReverse::Real* y_b, size_t n, codi::DataStore* d) {
  CODI_UNUSED(m);
  CODI_UNUSED(n);
  CODI_UNUSED(y);
  CODI_UNUSED(d);

  codi::RealReverse::Real ATrans[4] = {x[0], x[2], x[1], x[3]};

  codi::RealReverse::Real s[2];
  solve2(ATrans, y_b, s);

  // Adjoint of A (\bar A = -s*x^T) (In local terms x[0-3] = -s*y^T)
  x_b[0] = -s[0] * y[0];
  x_b[1] = -s[0] * y[1];
  x_b[2] = -s[1] * y[0];
  x_b[3] = -s[1] * y[1];

  // Adjoint of b (\bar b = s) (In local terms x[4-5] = s)
  x_b[4] = s[0];
  x_b[5] = s[1];
}

void primal() {
  std::cout << "double:" << std::endl;

  double u = 3.0;

  double A[4] = {u * 1.0, 0.5,  0.25, u * -1.0};
  double b[2] = {u * 10.0, u * 20.0};

  double x[2];

  solve2(A, b, x);

  double w = sqrt(x[0] * x[0] + x[1] * x[1]);

  std::cout << "Solution w: " << w << std::endl;
}

void derivative() {
  std::cout << "codi::RealReverse:" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  codi::RealReverse b[2] = {u * 10.0, u * 20.0};

  codi::RealReverse x[2];

  solve2(A, b, x);

  codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void externalFunction() {
  std::cout << "codi::RealReverse(External function):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  codi::RealReverse b[2] = {u * 10.0, u * 20.0};

  codi::RealReverse x[2];

  // external function helper start
  codi::ExternalFunctionHelper<codi::RealReverse> eh;
  for(int i = 0; i < 4; ++i) {
    eh.addInput(A[i]);
  }
  for(int i = 0; i < 2; ++i) {
    eh.addInput(b[i]);
  }

  for(int i = 0; i < 2; ++i) {
    eh.addOutput(x[i]);
  }

  eh.callPrimalFunc(solve2_primal);
  eh.addToTape(solve2_rev);
  // external function helper end

  codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void externalFunctionPassive() {
  std::cout << "codi::RealReverse(External function passive):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  codi::RealReverse b[2] = {u * 10.0, u * 20.0};

  codi::RealReverse x[2];

  // external function helper start
  codi::ExternalFunctionHelper<codi::RealReverse> eh;
  for(int i = 0; i < 4; ++i) {
    eh.addInput(A[i]);
  }
  for(int i = 0; i < 2; ++i) {
    eh.addInput(b[i]);
  }

  eh.callPassiveFunc(solve2<codi::RealReverse>, A, b, x);

  for(int i = 0; i < 2; ++i) {
    eh.addOutput(x[i]);
  }

  eh.addToTape(solve2_rev);
  // external function helper end

  codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

int main(int nargs, char** args) {

  primal();
  derivative();
  externalFunction();
  externalFunctionPassive();

  return 0;
}
