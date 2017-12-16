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

#include <vector>

IN(2)
OUT(2)
POINTS(1) = {
  {  1.0,     0.5}
};


void evalFunc(NUMBER* x, NUMBER* y, size_t size) {
  y[0] = x[0];
  y[1] = x[0];
  for(size_t i = 1; i < size; ++i) {
    y[0] += x[i];
    y[1] = max(y[1], x[i]);
  }
}

void func(NUMBER* x, NUMBER* y) {

#ifdef REVERSE_TAPE
  codi::PreaccumulationHelper<NUMBER> ph;
#else
  codi::ForwardPreaccumulationHelper<NUMBER> ph;
#endif

  const size_t size = codi::MaxStatementIntSize * 3;
  NUMBER intermediate[size];

  for(size_t i = 0; i < size; ++i) {
    intermediate[i] = x[0] * (double)i + x[1];
  }

  ph.start();
  for(size_t i = 0; i < size; ++i) {
    ph.addInput(intermediate[i]);
  }

  evalFunc(intermediate, y, size);

  ph.addOutput(y[0]);
  ph.addOutput(y[1]);
  ph.finish(false);
}
