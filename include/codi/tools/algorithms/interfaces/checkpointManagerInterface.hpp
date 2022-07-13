/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <string>
#include <vector>

#include "../../../config.h"
#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../misc/enumBitset.hpp"
#include "../../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    struct CheckpointHandle {
        int getIteration() const;

        void setAlgorithmData(void* data);
        void* getAlgorithmData() const;
    };

    struct CheckpointBase : public CheckpointHandle {
      private:
        int iteration;

        void* appData;

      public:

        CheckpointBase(int iteration) : iteration(iteration), appData(NULL) {}

        int getIteration() const {
          return iteration;
        }

        void setAlgorithmData(void* data) {
          this->appData = data;
        }

        void* getAlgorithmData() const {
          return appData;
        }
    };

    struct CheckpointManagerInterface {
        virtual ~CheckpointManagerInterface() {}

        virtual CheckpointHandle* create() = 0;
        virtual std::vector<CheckpointHandle*> list() = 0;
        virtual void load(CheckpointHandle* cp) = 0;
        virtual void remove(CheckpointHandle* cp) = 0;

        virtual void write(CheckpointHandle* cp) = 0;
        virtual void read(CheckpointHandle* cp) = 0;
    };
  }
}
