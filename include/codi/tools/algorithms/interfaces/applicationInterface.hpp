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
#include "../../../misc/macros.hpp"
#include "../enums/applicationHints.hpp"
#include "../tools/residuum.hpp"
#include "applicationIOInterface.hpp"
#include "checkpointManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    template<typename T_Type>
    struct ApplicationInterface {
      public:

        using Type = CODI_DD(T_Type, CODI_LHS_EXPRESSION_PROXY);


        using Real = RealTraits::Real<Type>;
        using Res = Residuum<Real>;

        using InitFunc = void (*)(Type& value, void* data);

        InitFunc initFunc;
        void* initData;

        ApplicationInterface() : initFunc(nullptr), initData(nullptr) {}

        virtual ~ApplicationInterface() {}

        template<typename Func>
        void iterateY(Func&& func);

        template<typename Func>
        void iterateX(Func&& func);

        template<typename Func>
        void iterateP(Func&& func);

        template<typename Func>
        void iterateZ(Func&& func);

        size_t getSizeY();
        size_t getSizeX();
        size_t getSizeP();
        size_t getSizeZ();

        virtual int getNumberOfFunctionals() = 0;

        void evaluateG();
        void evaluateF();
        void evaluateP();

        Res residuumY(std::vector<Real> const& v1, std::vector<Real> const& v2);
        Res residuumX(std::vector<Real> const& v1, std::vector<Real> const& v2);
        Res residuumP(std::vector<Real> const& v1, std::vector<Real> const& v2);

        CheckpointManagerInterface* getCheckpointInterface();
        ApplicationIOInterface<Type>* getIOInterface();

        void initialize();
        ApplicationHints getHints();
        int getIteration();

        void print(std::string const& line);
        bool isConverged();  ///< Check if primal application has converged.
        bool isStop();       ///< External stop to abort algorithm

        // Init variable handling
        void setInitializationHandlingFunction(InitFunc const& func, void* data = nullptr) {
          initFunc = func;
          initData = data;
        }
        void handleInitializationVariable(Type& value) {
          if(nullptr != initFunc) {
            initFunc(value, initData);
          }
        }
    };
  }
}
