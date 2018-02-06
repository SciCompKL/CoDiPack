/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2018 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */
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
        end[i] += A[j + i * n] * cur[j];
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

  codi::RealReverse A[4] = {u * 1.0, 0.5,
                                0.0, u * -1.0};
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

  codi::RealReverse A[4] = {u * 1.0, 0.5,
                                0.0, u * -1.0};
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
