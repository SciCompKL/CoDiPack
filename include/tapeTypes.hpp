/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#define CODI_INLINE_FORWARD_TAPE_TYPES(name) \
  typedef typename name::Real Real; \
  typedef typename name::GradientValue GradientValue; \
  typedef typename name::PassiveReal PassiveReal;

  /**
   * @brief Defines all the basic types that all tapes use.
   */
  template<typename RealType, typename GradientValueType>
  struct ForwardTapeTypes {

      typedef RealType Real;
      typedef GradientValueType GradientValue;

      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
  };

#define CODI_INLINE_REVERSE_TAPE_TYPES(name) \
  typedef typename name::Real Real; \
  typedef typename name::GradientValue GradientValue; \
  typedef typename name::PassiveReal PassiveReal; \
  typedef typename name::IndexHandler IndexHandler; \
  typedef typename name::IndexType IndexType;


  /**
   * @brief Defines all the basic types that all reverse tapes use.
   */
  template<typename RealType, typename GradientValueType, typename IndexHandlerType>
  struct ReverseTapeTypes {

      typedef RealType Real;
      typedef GradientValueType GradientValue;
      typedef IndexHandlerType IndexHandler;
      typedef typename IndexHandlerType::IndexType IndexType;

      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
  };
}
