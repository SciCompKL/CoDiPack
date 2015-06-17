/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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

#include <tools/dataStore.hpp>

#include <iostream>

IN(2)
OUT(1)
POINTS(1) = {{2.0, 3.0}};

void func_forward(NUMBER& z, const NUMBER& w, const NUMBER& v){
  z = w*v;
}

const int ITER = 5;

#ifdef REVERSE_TAPE
static void extFunc(void* checkpoint){
  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);
  NUMBER *x, w0, w1;
  check->getData(x);
  check->getData(w0);
  check->getData(w1);
  w0.gradient() += x[1].getValue()*w1.getGradient();
  x[1].gradient() += w0.getValue()*w1.getGradient();
}

static void delFunc(void* checkpoint){
  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);
  delete check;
  std::cout << "Delete" << std::endl;
}

void func(NUMBER* x, NUMBER* y) {
  NUMBER::TapeType& tape = NUMBER::getGlobalTape();
  NUMBER w[ITER];

  w[0] = x[0];
  for(int i = 1; i < ITER; ++i) {
    tape.setPassive();
    func_forward(w[i],w[i - 1],x[1]);
    tape.setActive();

    codi::DataStore *checkpoint = new codi::DataStore();
    tape.registerInput(w[i]);
    checkpoint->addData(x);
    checkpoint->addData(w[i-1]);
    checkpoint->addData(w[i]);
    tape.pushExternalFunctionHandle(&extFunc, checkpoint, delFunc);
  }


  y[0] = w[ITER - 1]*w[ITER - 1];
}
#else
void func(NUMBER*x, NUMBER* y){
  NUMBER w[ITER];

  w[0] = x[0];
  for(int i = 1; i < ITER; ++i) {
    func_forward(w[i],w[i - 1],x[1]);
  }


  y[0] = w[ITER - 1]*w[ITER - 1];
}
#endif

