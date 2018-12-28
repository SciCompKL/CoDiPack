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

void func_forward(NUMBER& z, const NUMBER& w, const NUMBER& v){
  z = w*v;
}

static void extFunc(void* t, void* checkpoint, void* i){
  CODI_UNUSED(t);

  codi::AdjointInterface<typename NUMBER::Real, typename NUMBER::GradientData>* ra = (codi::AdjointInterface<typename NUMBER::Real, typename NUMBER::GradientData>*)i;

  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);

  typename NUMBER::Real x1_v, x2_v;
  typename NUMBER::GradientData x1_i, x2_i, w_i;
  check->getData(x1_v);
  check->getData(x1_i);
  check->getData(x2_v);
  check->getData(x2_i);
  check->getData(w_i);

  size_t dim = ra->getVectorSize();

  for(size_t i = 0; i < dim; ++i) {

    typename NUMBER::Real w_b = ra->getAdjoint(w_i, i);
    ra->resetAdjoint(w_i, i);

    ra->updateAdjoint(x1_i, i, x2_v*w_b);
    ra->updateAdjoint(x2_i, i, x1_v*w_b);
  }
}

static void delFunc(void* tape, void* checkpoint){
  (void) tape;

  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);
  delete check;
  std::cout << "Delete" << std::endl;
}

void func(NUMBER* x, NUMBER* y) {
  NUMBER::TapeType& tape = NUMBER::getGlobalTape();
  codi::DataStore *checkpoint = new codi::DataStore;
  NUMBER w;
  tape.setPassive();
  func_forward(w,x[0],x[1]);
  tape.setActive();

  tape.registerInput(w);
  checkpoint->addData(x[0].getValue());
  checkpoint->addData(x[0].getGradientData());
  checkpoint->addData(x[1].getValue());
  checkpoint->addData(x[1].getGradientData());
  checkpoint->addData(w.getGradientData());
  tape.pushExternalFunctionHandle(&extFunc, checkpoint, delFunc);

  y[0] = w*w;
}
