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

#include <toolDefines.h>

#include <iostream>

IN(2)
OUT(1)
POINTS(1) = {{2.0, 3.0}};

void func_forward(NUMBER& z, const NUMBER& w, const NUMBER& v) {
  z = w*v;
}

const int ITER = 5;

#if REVERSE_TAPE
void func_reverse(const NUMBER::Real* x, NUMBER::Real* x_b, size_t m, const NUMBER::Real* y, const NUMBER::Real* y_b, size_t n, codi::DataStore* d) {
  CODI_UNUSED(m);
  CODI_UNUSED(n);
  CODI_UNUSED(y);
  CODI_UNUSED(d);

  x_b[0] = x[1] * y_b[0];
  x_b[1] = x[0] * y_b[0];
}

void func(NUMBER* x, NUMBER* y) {
  NUMBER w[ITER];

  w[0] = x[0];
  for(int i = 1; i < ITER; ++i) {

    codi::ExternalFunctionHelper<NUMBER> eh(true);

    eh.addInput(x[1]);
    eh.addInput(w[i-1]);

    eh.callPassiveFunc(func_forward, w[i], w[i-1], x[1]);
    eh.addOutput(w[i]);

    eh.addToTape(func_reverse);
  }

  y[0] = w[ITER - 1]*w[ITER - 1];
}
#else
void func(NUMBER* x, NUMBER* y) {
  NUMBER w[ITER];

  w[0] = x[0];
  for(int i = 1; i < ITER; ++i) {

    func_forward(w[i], w[i-1], x[1]);
  }

  y[0] = w[ITER - 1]*w[ITER - 1];
}
#endif
