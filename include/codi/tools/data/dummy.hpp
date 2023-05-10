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

#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Dummy value that can be assigned.
  struct DummyValue {
    public:

      /// Ignore the assignment.
      template<typename T>
      CODI_INLINE void operator=(T const& v) {
        CODI_UNUSED(v);
      }
  };

  /// Dummy vector that provides a dummy element access and size function.
  struct DummyVector {
    public:

      /// Return a dummy value.
      CODI_INLINE DummyValue operator[](size_t const i) {
        CODI_UNUSED(i);
        return DummyValue();
      }

      /// Return a dummy value.
      CODI_INLINE DummyValue operator[](size_t const i) const {
        CODI_UNUSED(i);
        return DummyValue();
      }

      /// Vector is always zero size.
      CODI_INLINE size_t size() const {
        return (size_t)0;
      }
  };
}
