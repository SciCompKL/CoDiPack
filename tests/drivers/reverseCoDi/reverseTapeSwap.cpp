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
#include <codi.hpp>

#include "reverseDriverBase.hpp"

struct ReverseDriverSwapTape : public ReverseDriverBase {

    NUMBER::TapeType swapTape;

    NUMBER::GradientValue& getGradient(NUMBER &number) {
      return swapTape.gradient(number.getGradientData());
    }

    void evaluate() {
      swapTape.evaluate();
    }

    void doPreEvaluate() {
      swapTape.swap(NUMBER::getGlobalTape());
    }

    void doLoopCleanup() {
      swapTape.reset();
    }
};

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  ReverseDriverSwapTape driver;

  driver.run();
}
