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

#include <map>


#include "../config.h"
#include "../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  namespace SparseEvaluation {

    enum class ElemeniationMissingOutput {
      Ignore = 0,
      Add,
      Throw
    };

    template<typename Real, typename Identifier>
    using NodeDependencies = std::map<Identifier, Real>;

    template<typename Real, typename Identifier>
    using DependencyMap = std::map<Identifier, NodeDependencies<Real, Identifier>>;


    template<typename Real, typename Identifier>
    CODI_INLINE static bool getIncomingDependencies(
        DependencyMap<Real, Identifier>& dependencies,
        Identifier const& lhsIdentifier,
        NodeDependencies<Real, Identifier>& result,
        ElemeniationMissingOutput missingOutputHandling) {

      auto iter = dependencies.find(lhsIdentifier);
      if(iter != dependencies.end()) {
        result = std::move(iter->second);
        dependencies.erase(iter);

        return true;
      } else {
        switch (missingOutputHandling) {
          case ElemeniationMissingOutput::Ignore:
            return false;
            break;
          case ElemeniationMissingOutput::Add:
            result.clear();
            result[lhsIdentifier] = Real(1.0); // Add a self reference.
            return true;
            break;
          case ElemeniationMissingOutput::Throw:
            CODI_EXCEPTION("Node for '%d' not in depency map. It needs to be declared as an output.", (int)lhsIdentifier);
            return false;
            break;
          default:
            CODI_EXCEPTION("Missing case '%d'.", (int)missingOutputHandling);
            return false;
            break;
        }
      }
    }
  }
}
