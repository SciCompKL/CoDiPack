/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "jacobianInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper class for the delayed write access of a reference.
  ///
  /// This class can be returned instead of a reference when the owner of the reference wants to be informed about write
  /// actions to the reference. Each assign call is forwarded to `data.setLogic(i,j, v)`.
  ///
  /// @tparam T_Issuer  The issuing class of the delay accessor.
  template<typename T_Issuer>
  struct JacobianDelayAccessor {
    public:

      using Issuer = CODI_DD(T_Issuer, CODI_T(JacobianInterface<CODI_ANY>));  ///< See DelayAccessor.

    private:

      size_t i;
      size_t j;

      Issuer& data;

    public:

      /// Constructor
      JacobianDelayAccessor(size_t const i, size_t const j, Issuer& data) : i(i), j(j), data(data) {}

      /// Forwards to `data.setLogic(i, j, v)`.
      template<typename T>
      JacobianDelayAccessor& operator=(T const& v) {
        data.setLogic(i, j, v);

        return *this;
      }

      /// Convert to the underlying type.
      operator typename Issuer::T() const {
        return const_cast<Issuer const&>(data).operator()(i, j);
      }
  };
}
