Tutorial A2: Preaccumulation {#TutorialA2}
============

The preaccumulation technique tries to reduce the memory of an AD tape
by reducing the data that is stored.

The function we want to optimize for memory is defined as
\f[
    y = f(x)
\f]
with \f$x \in \R^n\f$ and \f$y \in \R^m\f$.The reverse AD mode for this
function is
\f[
    \begin{aligned}
        \bar x \aeq & \frac{\d f}{\d x}(x) \bar y \\
        \bar y = & 0
    \end{aligned}
\f]
The CoDiPack tape stores an abstract definition of \f$\frac{\d f}{\d x}\f$
which is usually much more efficient then storing the Jacobian directly.
If \f$n\f$ and \f$m\f$ are quite large e.g. \f$10^6\f$ the memory for the
Jacobian would be 7450 Gb. If \f$n\f$ and \f$m\f$ are quite small e.g. 2
the Jacobian requires only 32 bytes. Depending on the computational complexity
of \f$f\f$ it can be more efficient to store the Jacobian on the tape
instead of the abstract AD representation.

Lets assume \f$f\f$ is implemented as:
~~~~{.cpp}
// time step ode in a explicit euler sheme
// x'(t) = Ax(t)
// x_n = x_c + dt * Ax(t)
void ode(const codi::RealReverse* start, codi::RealReverse* end, int steps, codi::RealReverse* A, double dt, size_t n) {
  codi::RealReverse* cur = new codi::RealReverse[n];

  for(size_t i = 0; i < n; ++i) {
    end[i] = start[i];
  }

  for(int curStep = 0; curStep < steps; ++curStep) {

    std::swap(end, cur);

    for(size_t i = 0; i < n; ++i) {
      end[i] = 0.0;
      for(size_t j = 0; j < n; ++j) {
        end[i] += A[i + j * n] * cur[j];
      }

      end[i] = cur[i] + dt * end[i];
    }
  }

  // we need to copy the result again if the number of steps is uneven
  if(steps % 2 == 1) {
    std::swap(end, cur);

    for(size_t i = 0; i < n; ++i) {
      end[i] = cur[i];
    }
  }

  delete[] cur;
}
~~~~
The inputs for this function are `start` and `A` and the outputs are `end`.
Depending on the number of steps the number of operations can be quite
large. In contrast, the number of inputs and outputs does not depend on the
number of steps and remains constant.

The aim is now to store the Jacobi of the call to `ode` instead of the
tape representation. The technique which we use for the is called preaccumulation.
The basic idea is to record the code section and immediately perform a reverse
evaluation of that code section. This reverse evaluation will yield the Jacobi
matrix and this matrix is stored instead of the tape representation.

The function `ode` is called in this small example:
~~~~{.cpp}
    codi::RealReverse u = 3.0;

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();
    tape.registerInput(u);

    codi::RealReverse A[4] = {u * 1.0, 0.5,  0.0, u * -1.0};
    codi::RealReverse start[2] = {u * 10.0, u * 20.0};

    codi::RealReverse end[2];

    ode(start, end, 1000, A, 1.0 / 1000.0, 2);

    codi::RealReverse w = sqrt(end[0] * end[0] + end[1] * end[1]);

    tape.registerOutput(w);

    tape.setPassive();
    w.setGradient(1);

    tape.evaluate();

    tape.printStatistics();

    std::cout << "Solution w: " << w << std::endl;
    std::cout << "Adjoint u: " << u.getGradient() << std::endl;
~~~~
This example requires 185 Kb of memory.

In order to perform the preaccumulation for the ode solver the following
steps are required:
 - Start the preaccumulation
 - Add the inputs
 - Call the function
 - Add the outputs
 - Finish the preaccumulation

The first step is quite simple. Inputs can be directly given in the
start of the preaccumulation region or added via a separate function:
~~~~{.cpp}
    codi::PreaccumulationHelper<codi::RealReverse> ph;

    ph.start(start[0], start[1]);
    for(size_t i = 0; i < 4; ++i) {
      ph.addInput(a[i]);
    }
~~~~

After all the inputs have been added the function can be called without any
modifications:
~~~~{.cpp}
    ode(start, end, 1000, A, 1.0 / 1000.0, 2);
~~~~

Finally, the outputs can be added and the preaccumulation needs to be
finished:
~~~~{.cpp}
    ph.finish(false, end[0], end[1]);
~~~~

The boolean argument to the finish function indicates if the preaccumulation
can use the adjoint vector as it is (false) or if the adjoint vector needs
to be stored and restored (true). The second case is only required if
the user has already started evaluating the tape.

The full code for the preaccumulation is now:
~~~~{.cpp}
    codi::RealReverse u = 3.0;

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();
    tape.registerInput(u);

    codi::RealReverse A[4] = {u * 1.0, 0.5,  0.0, u * -1.0};
    codi::RealReverse start[2] = {u * 10.0, u * 20.0};

    codi::RealReverse end[2];

    codi::PreaccumulationHelper<codi::RealReverse> ph;

    ph.start(start[0], start[1]);
    for(size_t i = 0; i < 4; ++i) {
      ph.addInput(A[i]);
    }

    ode(start, end, 1000, A, 1.0 / 1000.0, 2);

    ph.finish(false, end[0], end[1]);

    codi::RealReverse w = sqrt(end[0] * end[0] + end[1] * end[1]);

    tape.registerOutput(w);

    tape.setPassive();
    w.setGradient(1);

    tape.evaluate();

    tape.printStatistics();

    std::cout << "Solution w: " << w << std::endl;
    std::cout << "Adjoint u: " << u.getGradient() << std::endl;
~~~~
The preaccumulated example uses now only 269 bytes of memory.

The preaccumulation of a code section yields a reduction in memory only if
the number of input values and the number of output values is smaller than
the number of operations evaluated in that section.

The codi::PreaccumulationHelper is programmed such that it can be used
multiple times. After the finish call, the helper structure is in a state
such that start can be called again.

The full code for the example is:
~~~~{.cpp}
#include <codi.hpp>

#include <iostream>
#include <algorithm>

// time step ode in a explicit euler sheme
// x'(t) = Ax(t)
// x_n = x_c + dt * Ax(t)
void ode(const codi::RealReverse* start, codi::RealReverse* end, int steps, codi::RealReverse* A, double dt, size_t n) {
  codi::RealReverse* cur = new codi::RealReverse[n];

  for(size_t i = 0; i < n; ++i) {
    end[i] = start[i];
  }

  for(int curStep = 0; curStep < steps; ++curStep) {

    std::swap(end, cur);

    for(size_t i = 0; i < n; ++i) {
      end[i] = 0.0;
      for(size_t j = 0; j < n; ++j) {
        end[i] += A[i + j * n] * cur[j];
      }

      end[i] = cur[i] + dt * end[i];
    }
  }

  // we need to copy the result again if the number of steps is uneven
  if(steps % 2 == 1) {
    std::swap(end, cur);

    for(size_t i = 0; i < n; ++i) {
      end[i] = cur[i];
    }
  }

  delete[] cur;
}

void derivative() {
  std::cout << "codi::RealReverse:" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = {u * 1.0, 0.5,  0.0, u * -1.0};
  codi::RealReverse start[2] = {u * 10.0, u * 20.0};

  codi::RealReverse end[2];

  ode(start, end, 1000, A, 1.0 / 1000.0, 2);

  codi::RealReverse w = sqrt(end[0] * end[0] + end[1] * end[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  tape.printStatistics();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void preaccumulation() {
  std::cout << "codi::RealReverse(preaccumulation):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = {u * 1.0, 0.5,  0.0, u * -1.0};
  codi::RealReverse start[2] = {u * 10.0, u * 20.0};

  codi::RealReverse end[2];

  codi::PreaccumulationHelper<codi::RealReverse> ph;

  ph.start(start[0], start[1]);
  for(size_t i = 0; i < 4; ++i) {
    ph.addInput(A[i]);
  }

  ode(start, end, 1000, A, 1.0 / 1000.0, 2);

  ph.finish(false, end[0], end[1]);

  codi::RealReverse w = sqrt(end[0] * end[0] + end[1] * end[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  tape.printStatistics();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

int main(int nargs, char** args) {

  derivative();
  codi::RealReverse::getGlobalTape().reset();
  preaccumulation();

  return 0;
}
~~~~
