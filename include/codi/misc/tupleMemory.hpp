/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <tuple>

#include "../config.h"
#include "../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Leaf for tuple implementation see TupleMemory for details.
   */
  template<size_t i, typename T>
  struct TupleMemoryLeaf {
    public:
      T value;  ///< Memory value.

      /// Constructor.
      template<typename TT>
      CODI_INLINE TupleMemoryLeaf(TT&& v) : value(std::forward<TT>(v)) {}

      /// Constructor.
      CODI_INLINE TupleMemoryLeaf(TupleMemoryLeaf const& v) = default;
  };

  /**
   * @brief Base for tuple implementation see TupleMemory for details.
   */
  template<typename Ids, typename... Ts>
  struct TupleMemoryBase;

  /**
   * @brief Base for tuple implementation see TupleMemory for details.
   */
  template<size_t... Ids, typename... Ts>
  struct TupleMemoryBase<std::integer_sequence<size_t, Ids...>, Ts...> : public TupleMemoryLeaf<Ids, Ts>... {
    public:

      /// Constructor.
      template<typename... TTs>
      CODI_INLINE TupleMemoryBase(TTs&&... v) : TupleMemoryLeaf<Ids, Ts>(std::forward<TTs>(v))... {}

      /// Constructor.
      CODI_INLINE TupleMemoryBase(TupleMemoryBase const& v) = default;
  };

  /**
   * @brief Tuple implementation which allows to force inline of the construction of the tuple.
   *
   * This is just a minimal implementation which we need for ComputeExpression.
   *  - TupleMemoryLeaf stores one entry from the specified tuple types.
   *  - TupleMemoryBase helper class unrolling the integer_sequence into a variadic template argument.
   *
   * @tparam Ts  Type of the stored elements. References are kept.
   */
  template<typename... Ts>
  struct TupleMemory : public TupleMemoryBase<std::make_index_sequence<sizeof...(Ts)>, Ts...> {
    public:

      using Args = std::tuple<Ts...>;  ///< Helper for getting the i-th element.

      using Base = TupleMemoryBase<std::make_index_sequence<sizeof...(Ts)>, Ts...>;  ///< Abbreviation for base class.

      /// Constructor.
      template<typename... TTs>
      CODI_INLINE TupleMemory(TTs&&... v) : Base(std::forward<TTs>(v)...) {}

      /// Constructor.
      CODI_INLINE TupleMemory(TupleMemory const& v) = default;

      /// Get specific element.
      template<size_t i>
      CODI_INLINE auto& get() const {
        return static_cast<TupleMemoryLeaf<i, typename std::tuple_element<i, Args>::type> const*>(this)->value;
      }
  };
}
