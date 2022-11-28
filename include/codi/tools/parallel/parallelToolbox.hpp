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

#include "atomicInterface.hpp"
#include "mutexInterface.hpp"
#include "readWriteMutex.hpp"
#include "staticThreadLocalPointerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Constructs required to make CoDiPack suitable for a shared memory parallel environment.
   *
   * Hides details of a specific shared memory parallelization approach, such as OpenMP.
   *
   * Please refer to the individual constructs for documentation.
   *
   * @tparam T_Atomic                    Atomic implementation used by the toolbox. See codi::Atomic.
   * @tparam T_Mutex                     Mutex implementation used by the toolbox. See codi::Mutex.
   * @tparam T_StaticThreadLocalPointer  Static thread-local pointer implementation used by the toolbox.
   *                                     See codi::StaticThreadLocalPointer.
   */
  template<template<typename> class T_Atomic, typename T_Mutex,
           template<typename, typename> class T_StaticThreadLocalPointer>
  struct ParallelToolbox {
    public:
      template<typename Type>
      using Atomic = CODI_DD(T_Atomic<Type>, CODI_T(AtomicInterface<Type, CODI_ANY>));  ///< See codi::AtomicInterface.

      using Mutex = CODI_DD(T_Mutex, MutexInterface);            ///< See codi::MutexInterface.
      using Lock = codi::Lock<Mutex>;                            ///< See codi::Lock.
      using ReadWriteMutex = codi::ReadWriteMutex<Atomic<int>>;  ///< See codi::ReadWriteMutex.
      using LockForRead = codi::LockForRead<ReadWriteMutex>;     ///< See codi::LockForRead.
      using LockForWrite = codi::LockForWrite<ReadWriteMutex>;   ///< See codi::LockForWrite.

      /// See codi::StaticThreadLocalPointerInterface.
      template<typename Type, typename Owner>
      using StaticThreadLocalPointer = CODI_DD(CODI_T(T_StaticThreadLocalPointer<Type, Owner>),
                                               CODI_T(StaticThreadLocalPointerInterface<Type, Owner, CODI_ANY>));
  };
}
