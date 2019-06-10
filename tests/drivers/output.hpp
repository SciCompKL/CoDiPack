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

#pragma once

#include <codi/tools/data/hessian.hpp>

char const *const HEADER_FORMAT = "%6s_%03zd";
char const *const VALUE_FORMAT  = "%10g";
char const *const COL_SEPERATOR = " ";
char const *const LINE_END      = "\n";

template<typename Vec>
void writeOutputHessian(codi::Hessian<Vec> const& hes) {

  for(size_t curOut = 0; curOut < hes.m; curOut += 1) {

    // print header
    printf(HEADER_FORMAT, "out", curOut);
    for(size_t curIn = 0; curIn < hes.n; curIn += 1) {
      printf(COL_SEPERATOR);
      printf(HEADER_FORMAT, "in", curIn);
    }
    printf(LINE_END);

    for(size_t curIn1st = 0; curIn1st < hes.n; curIn1st += 1) {
      printf(HEADER_FORMAT, "in", curIn1st);

      for(size_t curIn2nd = 0; curIn2nd < hes.n; curIn2nd += 1) {
        printf(COL_SEPERATOR);
        printf(VALUE_FORMAT, hes(curOut, curIn1st, curIn2nd));
      }

      printf(LINE_END);
    }

    printf(LINE_END);
  }
}
