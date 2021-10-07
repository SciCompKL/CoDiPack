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

#include "../../../misc/macros.hpp"

#include "atomicInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_AtomicInt>
  struct ReadWriteMutex {
    public:
      using AtomicInt = CODI_DD(T_AtomicInt, CODI_T(AtomicInterface<int, CODI_ANY>));

    private:
      AtomicInt numReaders;
      AtomicInt numWriters;

    public:
      CODI_INLINE ReadWriteMutex() : numReaders(0), numWriters(0) {}
      ~ReadWriteMutex() {}

      void lockRead() {
        int currentWriters;
        while (true) {
          // wait until there are no writers
          do { currentWriters = numWriters; } while (currentWriters > 0);
          // register reader
          ++numReaders;
          // success if there are still no writers
          currentWriters = numWriters;
          if (currentWriters == 0) {
            break;
          }
          // otherwise let writers go first and try again
          --numReaders;
        }
      }

      void unlockRead() {
        --numReaders;
      }

      void lockWrite() {
        int currentWriters;
        while (true) {
          // register writer
          currentWriters = ++numWriters;
          // success if we are the first/only writer
          if (currentWriters == 1) {
            break;
          }
          // otherwise try again
          --numWriters;
        }

        int currentReaders;
        // wait until there are no readers
        do { currentReaders = numReaders; } while (currentReaders != 0);
      }

      void unlockWrite() {
        --numWriters;
      }
  };

  template<typename T_ReadWriteMutex>
  struct LockForRead {
    public:
      using ReadWriteMutex = CODI_DD(T_ReadWriteMutex, codi::ReadWriteMutex<CODI_ANY>);

    private:
      ReadWriteMutex& mutex;

      LockForRead(ReadWriteMutex& mutex) : mutex(mutex) {
        mutex.lockRead();
      }

      ~LockForRead() {
        mutex.unlockRead();
      }
  };

  template<typename T_ReadWriteMutex>
  struct LockForWrite {
    public:
      using ReadWriteMutex = CODI_DD(T_ReadWriteMutex, codi::ReadWriteMutex<CODI_ANY>);

    private:
      ReadWriteMutex& mutex;

      LockForWrite(ReadWriteMutex& mutex) : mutex(mutex) {
        mutex.lockWrite();
      }

      ~LockForWrite() {
        mutex.unlockWrite();
      }
  };

}
