/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
   * step is added to pos until end is reached.
   *
   * Be careful if step not -1 or 1, the end position needs to be reached directly.
   *
   * Example:
   * \code
   * int a[10];
   *
   * // Via lambda:
   * CompileTimeLoop<0,10,1>::eval([&] (auto i) { a[i.value] = i.value; });
   *
   * // Via functor:
   * struct SetEntry {
   *
   *   int* vec;
   *   template<size_t pos>
   *   void operator()(std::integral_constant<size_t, pos>, Type& v, Derivative const& d) {
   *     vec[pos] = pos;
   *   }
   * };
   * CompileTimeLoop<0,10,1>::eval(SetEntry{a});
   * \endcode
   *
   * Called range is: [pos, end)
   *
   * @tparam T_pos   Start value for the loop.
   * @tparam T_end   End value for the loop.
   * @tparam T_step  Step value for increment or decrement.
   */
  template<size_t T_pos, size_t T_end, int T_step>
  struct CompileTimeLoop {
    public:

      static size_t constexpr pos = T_pos;  ///< See CompileTimeLoop.
      static size_t constexpr end = T_end;  ///< See CompileTimeLoop.
      static int constexpr step = T_step;   ///< See CompileTimeLoop.

      /// Func is evaluated with args as func(pos, args...)
      template<typename Func, typename... Args>
      static CODI_INLINE void eval(Func&& func, Args&&... args) {
        func(std::integral_constant<size_t, pos>{}, std::forward<Args>(args)...);
        CompileTimeLoop<pos + step, end, step>::eval(std::forward<Func>(func), std::forward<Args>(args)...);
      }
  };

  /// Termination of loop evaluation.
  template<size_t T_pos, int T_step>
  struct CompileTimeLoop<T_pos, T_pos, T_step> {
    public:

      static size_t constexpr pos = T_pos;  ///< See CompileTimeLoop.
      static size_t constexpr end = T_pos;  ///< See CompileTimeLoop.
      static int constexpr step = T_step;   ///< See CompileTimeLoop.

      /// Nothing is evaluated.
      template<typename... Args>
      static CODI_INLINE void eval(Args&&... args) {
        CODI_UNUSED(args...);
      }
  };

  /// Static for with i = 0 .. (N - 1). See CompileTimeLoop for details.
  template<std::size_t N, typename F, typename... Args>
  CODI_INLINE void static_for(F func, Args&&... args) {
    CompileTimeLoop<0, N, 1>::eval(func, std::forward<Args>(args)...);
  }
}
