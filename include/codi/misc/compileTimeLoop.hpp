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

#include <utility>

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Compile time loop evaluation.
   *
   * pos is counted backwards until zero excluding zero.
   *
   * Called range is: (0,pos]
   *
   * @tparam T_pos  Starting value for the loop. Counted downwards.
   */
  template<size_t T_pos>
  struct CompileTimeLoop {
    public:

      static size_t constexpr pos = T_pos;  ///< See CompileTimeLoop.

      /// Func is evaluated with args as func(pos, args...)
      template<typename Func, typename... Args>
      static CODI_INLINE void eval(Func&& func, Args&&... args) {
        func(std::integral_constant<size_t, pos>{}, std::forward<Args>(args)...);

        CompileTimeLoop<pos - 1>::eval(std::forward<Func>(func), std::forward<Args>(args)...);
      }
  };

  /// Termination of loop evaluation.
  template<>
  struct CompileTimeLoop<0> {
    public:

      static size_t constexpr pos = 0;  ///< See CompileTimeLoop.

      /// Nothing is evaluated.
      template<typename... Args>
      static CODI_INLINE void eval(Args&&... args) {
        CODI_UNUSED(args...);
      }
  };
}
