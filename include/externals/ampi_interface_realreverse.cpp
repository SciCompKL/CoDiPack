/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include "ampi_tape.hpp"
#include "../tools/dataStore.hpp"
#include "../codi.hpp"

//#define INT64 int

#ifndef AD_TYPE
  #error Please specify the ad type with AD_TYPE.
#endif

typedef AD_TYPE type;

extern "C" {
  //forward declare von AMPI

  void ampi_get_val(void *buf, int *i, double *x) {
    *x=static_cast<type*>(buf)[*i].getValue();
  }
  void ampi_set_val(void* buf, int *i, double *v) {
    type &dummy= static_cast<type*>(buf)[*i];
    dummy.setValue(*v);
  }

  void ampi_get_idx(void *buf, int *i, INT64 *idx) {
    type &var = static_cast<type*>(buf)[*i];
    *idx = var.getGradientData();
  }


  void ampi_get_adj(INT64 *idx, double *x) {
    int index = *idx;
    if(index!=0) {
        *x = type::getGlobalTape().getGradient(index);
        type::getGlobalTape().setGradient(index, 0.0);
    } else {
        *x = 0.0;
    }
  }
  void ampi_set_adj(INT64 *idx, double *x) {
   int index = *idx;
    if(*idx!=0 && *x != 0.0){
        type::getGlobalTape().gradient(index) += *x;
    }
  }

  void delFunc(void* tape, void* cp) {
    CODI_UNUSED(tape);
    ampi_reset_entry(cp);
  }


  void execFunc(void*tape, void* handle) {
    CODI_UNUSED(tape);
    ampi_interpret_tape(handle);
  }

  void ampi_create_tape_entry(void* handle) {
     if (type::getGlobalTape().isActive()){
         type::getGlobalTape().pushExternalFunctionHandle(execFunc, handle, delFunc);
      }
  }

  void ampi_create_dummies_displ(void *buf, int* displ, int *size) {
      if (type::getGlobalTape().isActive()){

        type *values=static_cast<type*>(buf);
        for(int i=0;i<*size;++i) {
          //type &dummy=values[i];
          values[*displ + i]=0;
          type::getGlobalTape().registerInput(values[*displ + i]);
        }
     }
  }

  int ampi_is_tape_active () {
     return type::getGlobalTape().isActive();
  }
}
