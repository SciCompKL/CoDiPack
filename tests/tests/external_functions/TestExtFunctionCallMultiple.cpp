/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

const int ITER = 5;

#if REVERSE_TAPE
static void extFunc(void* t, void* checkpoint, void* i){
  CODI_UNUSED(t);

  codi::AdjointInterface<typename NUMBER::Real, typename NUMBER::GradientData>* ra = (codi::AdjointInterface<typename NUMBER::Real, typename NUMBER::GradientData>*)i;

  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);

  typename NUMBER::Real x_v, w0_v;
  typename NUMBER::GradientData x_i, w0_i, w1_i;
  check->getData(x_v);
  check->getData(x_i);
  check->getData(w0_v);
  check->getData(w0_i);
  check->getData(w1_i);

  size_t dim = ra->getVectorSize();

  for(size_t i = 0; i < dim; ++i) {

    typename NUMBER::Real w1_b = ra->getAdjoint(w1_i, i);
    ra->resetAdjoint(w1_i, i);

    ra->updateAdjoint(w0_i, i, x_v*w1_b);
    ra->updateAdjoint(x_i, i, w0_v*w1_b);
  }
}

static void extFuncPrimal(void* t, void* checkpoint, void* i) {
  CODI_UNUSED(t);

  codi::AdjointInterface<typename NUMBER::Real, typename NUMBER::GradientData>* ra = (codi::AdjointInterface<typename NUMBER::Real, typename NUMBER::GradientData>*)i;

  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);

  typename NUMBER::GradientData x1_i, x2_i, w_i;
  Real& x1_v = check->getDataRef<Real>();
  check->getData(x1_i);
  Real& x2_v = check->getDataRef<Real>();
  check->getData(x2_i);
  check->getData(w_i);

  x1_v = ra->getPrimal(x1_i); // Data is overwritten here
  x2_v = ra->getPrimal(x2_i); // Data is overwritten here

  typename NUMBER::Real z = x1_v * x2_v;
  ra->setPrimal( w_i, z);
}


static void delFunc(void* tape, void* checkpoint){
  (void) tape;

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
    checkpoint->addData(x[1].getValue());
    checkpoint->addData(x[1].getGradientData());
    checkpoint->addData(w[i-1].getValue());
    checkpoint->addData(w[i-1].getGradientData());
    checkpoint->addData(w[i].getGradientData());
    tape.pushExternalFunctionHandle(&extFunc, checkpoint, delFunc, nullptr, extFuncPrimal);
  }


  y[0] = w[ITER - 1]*w[ITER - 1];
}
#else
void func(NUMBER* x, NUMBER* y) {
  NUMBER w[ITER];

  w[0] = x[0];
  for(int i = 1; i < ITER; ++i) {
    func_forward(w[i],w[i - 1],x[1]);
  }


  y[0] = w[ITER - 1]*w[ITER - 1];
}
#endif
