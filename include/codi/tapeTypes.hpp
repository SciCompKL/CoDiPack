/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /*
   * Helper macro that creates definition of all types that are defined in the ForwardTapeTypes structure.
   */
#define CODI_INLINE_FORWARD_TAPE_TYPES(name) \
  typedef typename name::Real Real; /**< The floating point calculation type in the CoDiPack types. */ \
  typedef typename name::GradientValue GradientValue; /**< The type for the gradient computation */ \
  typedef typename name::PassiveReal PassiveReal; /**< The most inner floating point type if CoDiPack types are nested. */

  /**
   * @brief Defines all the basic types that forward tapes use.
   *
   * @tparam          RealType  The floating point computation type for the CoDiPack type. This type needs to implement all
   *                            mathematical functions and operators (e.g. sin, cos, +, *)
   * @tparam GradientValueType  The type for the computation of the gradient value. This type needs to implement the
   *                            addition operator and a left hand side scalar multiplication.
   */
  template<typename RealType, typename GradientValueType>
  struct ForwardTapeTypes {

      typedef RealType Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef GradientValueType GradientValue; /**< The type for the gradient computation */

      /**
       * The most inner floating point type if CoDiPack types are nested.
       */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
  };

  /*
   * Helper macro that creates definition of all types that are defined in the ReverseTapeTypes structure.
   */
#define CODI_INLINE_REVERSE_TAPE_TYPES(name) \
  typedef typename name::Real Real; /**< The floating point calculation type in the CoDiPack types. */ \
  typedef typename name::GradientValue GradientValue; /**< The type for the gradient computation */ \
  typedef typename name::PassiveReal PassiveReal; /**< The most inner floating point type if CoDiPack types are nested. */ \
  typedef typename name::IndexHandler IndexHandler; /**< The type of the index handler */ \
  typedef typename name::Index Index; /**< The actual type for the adjoint identification. */


  /**
   * @brief Defines all the basic types that reverse tapes use.
   *
   * @tparam          RealType  The floating point computation type for the CoDiPack type. This type needs to implement all
   *                            mathematical functions and operators (e.g. sin, cos, +, *)
   * @tparam GradientValueType  The type for the computation of the gradient value. This type needs to implement the
   *                            addition operator and a left hand side scalar multiplication.
   * @tparam  IndexHandlerType  The index handler for the identification of the adjoint values. It needs to implement the
   *                            common interface from the index handlers in include/tapes/indices
   */
  template<typename RealType, typename GradientValueType, typename IndexHandlerType>
  struct ReverseTapeTypes {

      typedef RealType Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef GradientValueType GradientValue; /**< The type for the gradient computation */
      typedef IndexHandlerType IndexHandler; /**< The type of the index handler */
      typedef typename IndexHandlerType::Index Index; /**< The actual type for the adjoint identification. */

      /**
       * The most inner floating point type if CoDiPack types are nested.
       */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
  };
}
