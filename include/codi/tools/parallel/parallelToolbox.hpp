/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "atomicInterface.hpp"
#include "readWriteMutex.hpp"
#include "staticThreadLocalPointerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<template<typename> T_Atomic, template<typename, typename> T_StaticThreadLocalPointer>
  struct ParallelToolbox {
    public:
      template<typename Type>
      using Atomic = CODI_DD(T_Atomic<Type>, CODI_T(AtomicInterface<Type, CODI_ANY>));

      using ReadWriteMutex = codi::ReadWriteMutex<Atomic<int>>;

      template<typename Type, typename Owner>
      using StaticThreadLocalPointer = CODI_DD(CODI_T(T_StaticThreadLocalPointer<Type, Owner>),
                                               CODI_T(StaticThreadLocalPointerInterface<Type, Owner, CODI_ANY>));
  };
}
