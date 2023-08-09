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
#include "atomicInterface.hpp"
#include "threadInformationInterface.hpp"

#ifdef __SANITIZE_THREAD__
  #define ANNOTATE_RWLOCK_CREATE(lock) AnnotateRWLockCreate(__FILE__, __LINE__, (void*)lock)
  #define ANNOTATE_RWLOCK_DESTROY(lock) AnnotateRWLockDestroy(__FILE__, __LINE__, (void*)lock)
  #define ANNOTATE_RWLOCK_ACQUIRED(lock, isWrite) AnnotateRWLockAcquired(__FILE__, __LINE__, (void*)lock, isWrite)
  #define ANNOTATE_RWLOCK_RELEASED(lock, isWrite) AnnotateRWLockReleased(__FILE__, __LINE__, (void*)lock, isWrite)

extern "C" void AnnotateRWLockCreate(const char* f, int l, void* addr);
extern "C" void AnnotateRWLockDestroy(const char* f, int l, void* addr);
extern "C" void AnnotateRWLockAcquired(const char* f, int l, void* addr, size_t isWrite);
extern "C" void AnnotateRWLockReleased(const char* f, int l, void* addr, size_t isWrite);
#endif

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Mutex construct that distinguishes between lock for read and lock for write.
   *
   * Since not all shared memory parallel APIs provide such mutexes, this is a custom implementation based on the atomic
   * type.
   *
   * The custom locking mechanism is annotated for the thread sanitizer so that the synchronization due to this mutex
   * is captured correctly when checking for data races.
   *
   * The user is responsible for correct pairing of lockForRead, unlockForRead and lockForWrite, unlockForWrite,
   * respectively. Use of the RAII locks LockForRead and LockForWrite is advised.
   *
   * Recursive locking for read is supported.
   *
   * @tparam T_ThreadInformation Implementation of ThreadInformationInterface.
   * @tparam T_AtomicInt  Implementation of AtomicInterface, instantiated with an underlying integer type.
   */
  template<typename T_ThreadInformation, typename T_AtomicInt>
  struct ReadWriteMutex {
    public:
      using ThreadInformation = CODI_DD(T_ThreadInformation, ThreadInformationInterface);  ///< See ReadWriteMutex.
      using AtomicInt = CODI_DD(T_AtomicInt, CODI_DEFAULT_ATOMIC<int>);                    ///< See ReadWriteMutex.

    private:
      AtomicInt numReaders;
      AtomicInt numWriters;
      int* nestingDepth;

#ifdef __SANITIZE_THREAD__
      int dummy;
#endif

    public:
      /// Constructor
      CODI_INLINE ReadWriteMutex()
          : numReaders(0),
            numWriters(0),
            nestingDepth(new int[ThreadInformation::getMaxThreads()]{})
#ifdef __SANITIZE_THREAD__
            ,
            dummy(0)
#endif
      {
#ifdef __SANITIZE_THREAD__
        ANNOTATE_RWLOCK_CREATE(&dummy);
#endif
      }

      /// Destructor
      ~ReadWriteMutex() {
        delete[] nestingDepth;
#ifdef __SANITIZE_THREAD__
        ANNOTATE_RWLOCK_DESTROY(&dummy);
#endif
      }

      /**
       * @brief Acquire mutex for read access.
       *
       * Waits until there are no writers. Multiple simultaneous acquisitions for reading are allowed.
       */
      void lockRead() {
        int currentWriters;
        while (true) {
          // nested lock for read
          int const threadId = ThreadInformation::getThreadId();
          if (nestingDepth[threadId] > 0) {
            ++nestingDepth[threadId];
            break;
          }

          // wait until there are no writers
          do {
            currentWriters = numWriters;
          } while (currentWriters > 0);
          // register reader
          ++numReaders;
          // success if there are still no writers
          currentWriters = numWriters;
          if (currentWriters == 0) {
            ++nestingDepth[threadId];
            break;
          }
          // otherwise let writers go first and try again
          --numReaders;
        }

#ifdef __SANITIZE_THREAD__
        ANNOTATE_RWLOCK_ACQUIRED(&dummy, false);
#endif
      }

      /// Release mutex that was acquired for read access.
      void unlockRead() {
#ifdef __SANITIZE_THREAD__
        ANNOTATE_RWLOCK_RELEASED(&dummy, false);
#endif
        int const threadId = ThreadInformation::getThreadId();
        --nestingDepth[threadId];
        if (nestingDepth[threadId] == 0) {
          --numReaders;
        }
      }

      /**
       * @brief Acquire mutex for write access.
       *
       * First writer comes first, as soon as there are no readers. Other writers, if any, wait until the first one is
       * done.
       */
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
        do {
          currentReaders = numReaders;
        } while (currentReaders != 0);

#ifdef __SANITIZE_THREAD__
        ANNOTATE_RWLOCK_ACQUIRED(&dummy, true);
#endif
      }

      /// Release mutex that was acquired for write access.
      void unlockWrite() {
#ifdef __SANITIZE_THREAD__
        ANNOTATE_RWLOCK_RELEASED(&dummy, true);
#endif

        --numWriters;
      }
  };

  /**
   * @brief RAII lock for read.
   *´
   * @tparam T_ReadWriteMutex  The underlying ReadWriteMutex.
   */
  template<typename T_ReadWriteMutex>
  struct LockForRead {
    public:
      using ReadWriteMutex =
          CODI_DD(T_ReadWriteMutex,
                  CODI_T(ReadWriteMutex<DefaultThreadInformation, CODI_DEFAULT_ATOMIC<int>>));  ///< See LockForRead.

    private:
      ReadWriteMutex& mutex;

    public:
      /// Constructor. Acquires lock for read access.
      LockForRead(ReadWriteMutex& mutex) : mutex(mutex) {
        mutex.lockRead();
      }

      /// Destructor. Releases lock for read access.
      ~LockForRead() {
        mutex.unlockRead();
      }
  };

  /**
   * @brief RAII lock for write
   *
   * @tparam T_ReadWriteMutex  The underlying ReadWriteMutex.
   */
  template<typename T_ReadWriteMutex>
  struct LockForWrite {
    public:
      using ReadWriteMutex =
          CODI_DD(T_ReadWriteMutex,
                  CODI_T(ReadWriteMutex<DefaultThreadInformation, CODI_DEFAULT_ATOMIC<int>>));  ///< See LockForWrite.

    private:
      ReadWriteMutex& mutex;

    public:
      /// Constructor. Acquires lock for write access.
      LockForWrite(ReadWriteMutex& mutex) : mutex(mutex) {
        mutex.lockWrite();
      }

      /// Destructor. Releases lock for write access.
      ~LockForWrite() {
        mutex.unlockWrite();
      }
  };
}
