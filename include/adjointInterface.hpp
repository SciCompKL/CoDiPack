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

namespace codi {

  struct AdjointInterface {
      virtual ~AdjointInterface() {}

      virtual void setLhsAdjoint(const int index) = 0;
      virtual void resetAdjoint(const int index) = 0;
      virtual void updateAdjoint(const int index, double jacobi) = 0;
  };

  template<typename GradientValue>
  struct AdjointHandler final : public AdjointInterface {
      GradientValue* adjointVector;

      GradientValue lhsSeed;

      AdjointHandler(GradientValue* adjointVector) :
        adjointVector(adjointVector),
        lhsSeed() {}

      void setLhsAdjoint(const int index) {
        lhsSeed = adjointVector[index];
      }

      void resetAdjoint(const int index) {
        adjointVector[index] = GradientValue();
      }

      void updateAdjoint(const int index, const double jacobi) {
        adjointVector[index] += jacobi * lhsSeed;
      }
  };
}
