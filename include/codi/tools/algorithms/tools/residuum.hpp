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

#include <vector>
#include <cmath>

#include "../../../misc/macros.hpp"
#include "../../../misc/stringUtil.hpp"
#include "../../../config.h"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    template<typename T_Real>
    struct Residuum {
        using Real = CODI_DD(T_Real, double);

        Real l2;
        Real l1;
        Real lMax;
        size_t lMaxPos;

        Residuum() = default;

        static Residuum<Real> vectorBasedResiduum(std::vector<Real> const& v1, std::vector<Real> const& v2) {
          Residuum<Real> res{};
          res.lMax = -1e300;

          codiAssert(v1.size() == v2.size());

          for (size_t i = 0; i < v1.size(); i += 1) {
            Real diff = std::abs(v1[i] - v2[i]);
            res.l1 += diff;
            res.l2 += diff * diff;
            if (res.lMax < diff) {
              res.lMax = diff;
              res.lMaxPos = i;
            }
          }

          res.l2 = std::sqrt(res.l2);

          return res;
        }

        std::string formatHeader(std::string prefix) const {

          return StringUtil::format("%sY_L1 %sY_L2 %sY_LMax %sY_LMaxPos", prefix.c_str(),
                                    prefix.c_str(), prefix.c_str(), prefix.c_str());
        }

        std::string formatEntry(int width = 6) const {
          return StringUtil::format("%0.*e %0.*e %0.*e %d", width, l1, width, l2, width, lMax, lMaxPos);
        }
    };
  }
}
