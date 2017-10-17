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


#pragma once

#include <vector>

#include "../configure.h"
#include "../exceptions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename CoDiType, typename GradientValue>
  struct TapeVectorHelper {

      typedef typename CoDiType::Real Real;
      typedef typename CoDiType::GradientData GradientData;
      typedef typename CoDiType::TapeType Tape;
      typedef typename Tape::Position Position;

      std::vector<GradientValue> adjointVector;
      Tape& tape;

      GradientValue zeroValue;
      const GradientValue constZeroValue;

      TapeVectorHelper() :
        adjointVector(0),
        tape(CoDiType::getGlobalTape()),
        zeroValue(),
        constZeroValue() {
      }

      void setTape(Tape& tape) {
        this->tape = tape;
      }

      void setGradient(GradientData& value, const GradientValue& gradientValue) {
        gradient(value) = gradientValue;
      }

      GradientValue getGradient(const GradientData& value) {
        return gradient(value);
      }

      GradientValue& gradient(GradientData& value) {
        checkAdjointVectorSize();

        if(0 != value) {
          return adjointVector[value];
        } else {
          zeroValue = GradientValue();
          return zeroValue;
        }
      }

      const GradientValue& gradient(const GradientData& value) const {
        if(0 != value && value < adjointVector.size()) {
          return adjointVector[value];
        } else {
          return constZeroValue;
        }
      }

      void evaluate(const Position& start, const Position& end) {
        checkAdjointVectorSize();

        tape.evaluate(start, end, adjointVector.data());
      }

      void evaluate() {
        evaluate(tape.getPosition(), tape.getZeroPosition());
      }

      void clearAdjoints() {
        for(size_t i = 0; i < adjointVector.size(); i += 1) {
          adjointVector[i] = GradientValue();
        }
      }

    private:

      void checkAdjointVectorSize() {
        if(adjointVector.size() <= tape.getAdjointSize()) {
          adjointVector.resize(tape.getAdjointSize() + 1);
        }
      }
  };
}
