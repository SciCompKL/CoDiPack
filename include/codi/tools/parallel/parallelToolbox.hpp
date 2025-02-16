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

#include "atomicInterface.hpp"
#include "mutexInterface.hpp"
#include "readWriteMutex.hpp"
#include "reverseAtomicInterface.hpp"
#include "staticThreadLocalPointerInterface.hpp"
#include "synchronizationInterface.hpp"
#include "threadInformationInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Collects parallel programming facilties required to make CoDiPack applicable in a shared memory parallel
   *        environment.
   *
   * The parallel programming facilities provided as template parameters have to be implemented against a specific
   * shared memory parallelism API, e.g., OpenMP. This class redeclares them and provides further derived types.
   * CoDiPack uses a ParallelToolbox for defining thread-safe active types and tapes while abstracting away the details
   * of the specific underlying shared memory parallelism API.
   *
   * Please refer to the individual using declarations for documentation.
   *
   * @tparam T_ThreadInformation         Thread information used by the toolbox. See codi::ThreadInformationInterface.
   * @tparam T_Atomic                    Atomic implementation used by the toolbox. See codi::AtomicInterface.
   * @tparam T_ReverseAtomic             Reverse atomic implementation used by the toolbox. See
   *                                     codi::ReverseAtomicInterface.
   * @tparam T_Mutex                     Mutex implementation used by the toolbox. See codi::MutexInterface.
   * @tparam T_StaticThreadLocalPointer  Static thread-local pointer implementation used by the toolbox.
   *                                     See codi::StaticThreadLocalPointerInterface.
   * @tparam T_Synchronization           Synchronization facalities that comply with the toolbox.
   *                                     See codi::SynchronizationInterface.
   */
  template<typename T_ThreadInformation, template<typename> class T_Atomic, template<typename> class T_ReverseAtomic,
           typename T_Mutex, template<typename, typename> class T_StaticThreadLocalPointer, typename T_Synchronization>
  struct ParallelToolbox {
    public:
      /// See codi::ParallelToolbox.
      using ThreadInformation = CODI_DD(T_ThreadInformation, ThreadInformationInterface);
      template<typename Type>
      using Atomic = CODI_DD(T_Atomic<Type>, CODI_DEFAULT_ATOMIC<Type>);  ///< See codi::ParallelToolbox.
      /// See codi::ParallelToolbox.
      template<typename Type>
      using ReverseAtomic = CODI_DD(T_ReverseAtomic<Type>, CODI_DEFAULT_REVERSE_ATOMIC<Type>);
      using Mutex = CODI_DD(T_Mutex, MutexInterface);  ///< See codi::ParallelToolbox.

      /// See codi::StaticThreadLocalPointerInterface.
      template<typename Type, typename Owner>
      using StaticThreadLocalPointer = CODI_DD(CODI_T(T_StaticThreadLocalPointer<Type, Owner>),
                                               CODI_T(StaticThreadLocalPointerInterface<Type, Owner, CODI_ANY>));

      using Synchronization = CODI_DD(T_Synchronization, DefaultSynchronization);  ///< See codi::ParallelToolbox.

      using Lock = codi::Lock<Mutex>;                                               ///< See codi::Lock.
      using ReadWriteMutex = codi::ReadWriteMutex<ThreadInformation, Atomic<int>>;  ///< See codi::ReadWriteMutex.
      using LockForRead = codi::LockForRead<ReadWriteMutex>;                        ///< See codi::LockForRead.
      using LockForWrite = codi::LockForWrite<ReadWriteMutex>;                      ///< See codi::LockForWrite.
  };

#if CODI_IDE
  /// Helper for IDE code completion.
  using CODI_DEFAULT_PARALLEL_TOOLBOX =
      ParallelToolbox<DefaultThreadInformation, CODI_DEFAULT_ATOMIC, MutexInterface,
                      CODI_DEFAULT_STATIC_THREAD_LOCAL_POINTER, DefaultSynchronization>;
#endif
}
