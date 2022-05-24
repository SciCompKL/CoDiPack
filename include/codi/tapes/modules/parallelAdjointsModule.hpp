/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include "../../configure.h"
#include "../tapeTraits.hpp"
#include "../reverseTapeInterface.hpp"

#include <omp.h>

#ifdef __SANITIZE_THREAD__
  #define ANNOTATE_RWLOCK_CREATE(lock) \
    AnnotateRWLockCreate(__FILE__, __LINE__, (void*)lock)
  #define ANNOTATE_RWLOCK_DESTROY(lock) \
    AnnotateRWLockDestroy(__FILE__, __LINE__, (void*)lock)
  #define ANNOTATE_RWLOCK_ACQUIRED(lock, isWrite) \
    AnnotateRWLockAcquired(__FILE__, __LINE__, (void*)lock, isWrite)
  #define ANNOTATE_RWLOCK_RELEASED(lock, isWrite) \
    AnnotateRWLockReleased(__FILE__, __LINE__, (void*)lock, isWrite)

  extern "C" void AnnotateRWLockCreate(const char* f, int l, void* addr);
  extern "C" void AnnotateRWLockDestroy(const char* f, int l, void* addr);
  extern "C" void AnnotateRWLockAcquired(const char* f, int l, void* addr, size_t isWrite);
  extern "C" void AnnotateRWLockReleased(const char* f, int l, void* addr, size_t isWrite);
#endif

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename TapeTypes, typename Tape>
  class ParallelAdjointsModule : public virtual ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, Tape, typename TapeTypes::Position >{
    private:
      // typedefs
      typedef typename TapeTypes::Real Real;
      typedef typename TapeTypes::Index Index;
      typedef typename TapeTypes::GradientValue GradientValue;
      typedef typename TapeTypes::GradientData GradientData;
      typedef typename TapeTypes::Position Position;

      /**
       * @brief Cast this class to the full.
       *
       * The full type is able to access all functions from the tape and other modules.
       *
       * @return  The full tape implemenation.
       */
      Tape& cast() {
        return *static_cast<Tape*>(this);
      }

      /**
       * @brief Cast this class to the full.
       *
       * The full type is able to access all functions from the tape and other modules.
       *
       * @return  The full tape implemenation.
       */
      const Tape& cast() const {
        return *static_cast<const Tape*>(this);
      }

    private:

      #ifdef __SANITIZE_THREAD__
        static int lockDummy;
      #endif

      /**
       * @brief The adjoint vector.
       *
       * The size of the adjoint vector is set according to the requested positions.
       * But the positions should not be greater than the current expression counter.
       */
      struct AdjointsWrapper {
        GradientValue* adjoints;
        int lockForUse;
        int lockForRealloc;

        AdjointsWrapper() : adjoints(NULL) {
          #ifdef __SANITIZE_THREAD__
            ANNOTATE_RWLOCK_CREATE(&lockDummy);
          #endif
        }

        ~AdjointsWrapper() {
          if (adjoints != NULL) {
            free(adjoints);
            adjoints = NULL;
          }

          #ifdef __SANITIZE_THREAD__
            ANNOTATE_RWLOCK_DESTROY(&lockDummy);
          #endif
        }
      };

      static AdjointsWrapper adjointsWrapper;

      /** @brief The current size of the adjoint vector. */
      static Index adjointsSize;

      static void setAdjointsSize(const Index& newSize) {
        adjointsSize = newSize;
      }

    public:

      static void lockForUse() {
        int numReallocators;
        while (true) {
          // wait until there are no reallocators
          do {
            #pragma omp atomic read
            numReallocators = adjointsWrapper.lockForRealloc;
          } while (numReallocators > 0);

          // increment lock for use
          #pragma omp atomic update
          ++adjointsWrapper.lockForUse;

          // check if there are still no reallocators
          #pragma omp atomic read
          numReallocators = adjointsWrapper.lockForRealloc;

          // if so, let reallocators go first and try again
          if (numReallocators > 0) {
            #pragma omp atomic update
            --adjointsWrapper.lockForUse;
            continue;
          }
          break;
        }

        #ifdef __SANITIZE_THREAD__
          ANNOTATE_RWLOCK_ACQUIRED(&lockDummy, false);
        #endif
      }

      static void unlockAfterUse() {
        #ifdef __SANITIZE_THREAD__
          ANNOTATE_RWLOCK_RELEASED(&lockDummy, false);
        #endif

        #pragma omp atomic update
        --adjointsWrapper.lockForUse;
      }

      static void lockForRealloc() {
        // wait until there is exactly one lock for realloc
        int numReallocators;
        while (true) {
          #pragma omp atomic capture
          numReallocators = ++adjointsWrapper.lockForRealloc;

          if (numReallocators != 1) {
            #pragma omp atomic update
            --adjointsWrapper.lockForRealloc;
            continue;
          }
          break;
        }

        // wait until there are no users
        int users;
        do {
          #pragma omp atomic read
          users = adjointsWrapper.lockForUse;
        } while (users != 0);

        #ifdef __SANITIZE_THREAD__
          ANNOTATE_RWLOCK_ACQUIRED(&lockDummy, true);
        #endif
      }

      static void unlockAfterRealloc() {
        #ifdef __SANITIZE_THREAD__
          ANNOTATE_RWLOCK_RELEASED(&lockDummy, true);
        #endif

        #pragma omp atomic update
        --adjointsWrapper.lockForRealloc;
      }

      struct LockUse {
          LockUse() {
            lockForUse();
          }

          ~LockUse() {
            unlockAfterUse();
          }
      };

      struct LockRealloc {
          LockRealloc() {
            lockForRealloc();
          }

          ~LockRealloc() {
            unlockAfterRealloc();
          }
      };

    protected:

      ParallelAdjointsModule() {
      }

      ~ParallelAdjointsModule() {
      }

      void initAdjointsModule() {
        // Nothing to do
      }

      void initializeAdjointsModule() {
      }

      void finalizeAdjointsModule() {
      }

    protected:

      CODI_INLINE GradientValue* getAdjoints() const {
        LockUse lock;
        return adjointsWrapper.adjoints;
      }

    // ----------------------------------------------------------------------
    // Protected functions for the communication with the including class
    // ----------------------------------------------------------------------

      /**
      * @brief Adds information about adjoint vector.
      *
      * Adds the number of adjoint vector entries and the size of the adjoint vector.
      *
      * @param[in,out] values  The information is added to the values
      */
      void CODI_INLINE addAdjointValues(TapeValues& values) const {

        size_t nAdjoints      = cast().indexHandler.getMaximumGlobalIndex() + 1;
        double memoryAdjoints = static_cast<double>(nAdjoints) * static_cast<double>(sizeof(GradientValue)) * BYTE_TO_MB;

        values.addSection("Adjoint vector");
        values.addData("Number of adjoints", nAdjoints);
        values.addData("Memory allocated", memoryAdjoints, true, true);

        cast().indexHandler.addValues(values);
      }

    private:
      /**
       * @brief Helper function: Sets the adjoint vector to a new size.
       *
       * @param[in] size The new size for the adjoint vector.
       */
      static CODI_NO_INLINE void resizeAdjoints(const Index& size) {

        LockRealloc lock;
        Index oldSize = adjointsSize;

        if (size > oldSize) {
          setAdjointsSize(size);

          for(Index i = size; i < oldSize; ++i) {
            adjointsWrapper.adjoints[i].~GradientValue();
          }

          adjointsWrapper.adjoints = (GradientValue*)realloc((void*)adjointsWrapper.adjoints, sizeof(GradientValue) * static_cast<size_t>(size));

          if(NULL == adjointsWrapper.adjoints) {
            throw std::bad_alloc();
          }

          for(Index i = oldSize; i < size; ++i) {
            new (adjointsWrapper.adjoints + i) GradientValue();
          }
        }
      }

    protected:
      /**
       * @brief Resize the adjoint vector such that it fits the number of indices.
       */
      void CODI_INLINE resizeAdjointsToIndexSize() {
        if(getAdjointsSize() <= cast().indexHandler.getMaximumGlobalIndex()) {
          resizeAdjoints(cast().indexHandler.getMaximumGlobalIndex() + 1);
        }
      }

      /**
       * @brief Helper function: Deletes all arrays
       */
      static CODI_NO_INLINE void cleanAdjoints() {
        if(adjointsValid()) {
          LockRealloc lock;
          free(adjointsWrapper.adjoints);
          adjointsWrapper.adjoints = NULL;
          setAdjointsSize(0);
        }
      }

      static CODI_INLINE bool adjointsValid() {
        LockUse lock;
        return adjointsWrapper.adjoints != NULL;
      }

      /**
       * @brief Swap the data of the tape base module with the data of the other tape base module.
       *
       * @param[in] other  The object with the other tape base module.
       */
      void swapAdjointsModule(Tape& other) {
        // makes no longer sense as the adjoint vector is a static member

        CODI_UNUSED(other);
      }

    public:
      // access to adjoints
      static CODI_INLINE Index getAdjointsSize() {
        LockUse lock;
        return adjointsSize;
      }

    public:
      // no boundary check access for derived class

      template<typename AdjointData>
      static CODI_INLINE void setAdjoint(const Index& index, const GradientValue& value, AdjointData* data) {
        data[index] = value;
      }

      static CODI_INLINE void setAdjoint(const Index& index, const GradientValue& value) {
        LockUse lock;
        setAdjoint(index, value, adjointsWrapper.adjoints);
      }

      template<typename AdjointData>
      static CODI_INLINE void incrementAdjoint(const Index& index, const AdjointData& adj, const Real& jacobi, AdjointData* data) {
        data[index] += adj * jacobi;
      }

      static CODI_INLINE void incrementAdjoint(const Index& index, const GradientValue& adj, const GradientValue& jacobi) {
        incrementAdjoint(index, adj, jacobi, adjointsWrapper.adjoints);
      }

      template<typename AdjointData>
      static CODI_INLINE void incrementTangent(AdjointData& adj, const AdjointData* data, const Index& index, const Real& jacobi) {
        adj += data[index] * jacobi;
      }

      template<typename AdjointData>
      static CODI_INLINE void clearAdjoint(const Index& index, AdjointData* data) {
        setAdjoint(index, GradientValue(), data);
      }

      static CODI_INLINE void clearAdjoint(const Index& index) {
        clearAdjoint(index, adjointsWrapper.adjoints);
      }

    public:
      /**
       * @brief Get the gradient value of the corresponding index.
       *
       * @param[in] index The index of the active type.
       * @return The gradient value corresponding to the given index.
       */
      CODI_INLINE GradientValue getGradient(const Index& index) const {
        if(0 == index || getAdjointsSize() <= index) {
          return GradientValue();
        } else {
          LockUse lock;
          return adjointsWrapper.adjoints[index];
        }
      }

      CODI_INLINE void setGradient(const Index& index, const GradientValue& gradientValue) {
        codiAssert(index <= cast().indexHandler.getMaximumGlobalIndex());

        if (index == 0) {
          resizeAdjointsToIndexSize();
        }
        else {

          //TODO: Add error when index is bigger than expression count
          if(getAdjointsSize() <= index) {
            resizeAdjoints(cast().indexHandler.getMaximumGlobalIndex() + 1);
          }

          setAdjoint(index, gradientValue);
        }
      }

      CODI_INLINE void setGradient(Index& index, const GradientValue& gradientValue) {
        Index const& indexConst = index;
        setGradient(indexConst, gradientValue);
      }

      /**
       * @brief Get a reference to the gradient value of the corresponding index.
       *
       * An index of 0 will raise an codiAssert exception.
       *
       * @param[in] index The index of the active type.
       * @return The reference to the gradient data.
       */
      CODI_INLINE GradientValue& gradient(const Index& index) {
        codiAssert(0 != index);
        codiAssert(index <= cast().indexHandler.getMaximumGlobalIndex());

        if(getAdjointsSize() <= index) {
          return adjointsWrapper.adjoints[0];
        } else {
          return adjointsWrapper.adjoints[index];
        }
      }

      /**
       * @brief Get a reference to the gradient value of the corresponding index.
       *
       * An index of 0 will raise an codiAssert exception.
       *
       * @param[in] index The index of the active type.
       * @return The reference to the gradient data.
       */
      CODI_INLINE GradientValue& gradient(Index& index) {
        Index const& indexConst = index;
        return gradient(indexConst);
      }

      /**
       * @brief Get a constant reference to the gradient value of the corresponding index.
       *
       * @param[in] index The index of the active type.
       * @return The constant reference to the gradient data.
       */
      CODI_INLINE const GradientValue& gradient(const Index& index) const {
        if(getAdjointsSize() <= index) {
          return adjointsWrapper.adjoints[0];
        } else {
          return adjointsWrapper.adjoints[index];
        }
      }

      /**
       * @brief Sets all adjoint/gradients to zero.
       */
      CODI_INLINE void clearAdjoints(){

        LockUse lock;

        if(adjointsValid()) {
          Index adjointsSize = getAdjointsSize();
          for(Index i = 0; i < adjointsSize; ++i) {
            clearAdjoint(i);
          }
        }
      }

      /**
       * @brief Clear the adjoint vector and delete it.
       */
      CODI_INLINE void deleteAdjointVector() {
        cleanAdjoints();
      }
  };

  #if __SANITIZE_THREAD__
    template<typename TapeTypes, typename Tape>
    int ParallelAdjointsModule<TapeTypes, Tape>::lockDummy;
  #endif

  template<typename TapeTypes, typename Tape>
  typename ParallelAdjointsModule<TapeTypes, Tape>::AdjointsWrapper ParallelAdjointsModule<TapeTypes, Tape>::adjointsWrapper;

  template<typename TapeTypes, typename Tape>
  typename ParallelAdjointsModule<TapeTypes, Tape>::Index ParallelAdjointsModule<TapeTypes, Tape>::adjointsSize = 0;
}
