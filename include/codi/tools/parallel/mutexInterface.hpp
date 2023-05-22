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

  /**
   * @brief Abstracts a mutex.
   *
   * Simple mutex with lock and unlock operations. Lock provides RAII locking.
   */
  struct MutexInterface {
    public:
      CODI_INLINE MutexInterface() {}  ///< Constructor
      ~MutexInterface() {}             ///< Destructor

      CODI_INLINE void initialize();  ///< Initialize the mutex.
      CODI_INLINE void finalize();    ///< Finalize the mutex.

      CODI_INLINE void lock();    ///< Lock the mutex.
      CODI_INLINE void unlock();  ///< Unlock the mutex.
  };

  /**
   * @brief RAII mutex locking.
   *
   * @tparam T_Mutex  The underlying mutex type.
   */
  template<typename T_Mutex>
  struct Lock {
    public:
      using Mutex = CODI_DD(T_Mutex, MutexInterface);  ///< See Lock.

    private:
      Mutex& mutex;

    public:
      /// Constructor. Locks the mutex.
      /// @param mutex  The mutex to be locked.
      Lock(Mutex& mutex) : mutex(mutex) {
        mutex.lock();
      }

      /// Destructor. Unlocks the mutex.
      ~Lock() {
        mutex.unlock();
      }
  };
}
