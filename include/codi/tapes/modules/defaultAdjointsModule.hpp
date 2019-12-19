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

#include "../../configure.h"
#include "../reverseTapeInterface.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {
  template<typename TapeTypes, typename Tape>
  class DefaultAdjointsModule : public virtual ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, Tape, typename TapeTypes::Position >{
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
      GradientValue* adjoints;

      /** @brief The current size of the adjoint vector. */
      Index adjointsSize;

      void setAdjointsSize(const Index& newSize) {
        adjointsSize = newSize;
      }

    protected:

      DefaultAdjointsModule() : adjoints(NULL), adjointsSize(0) {
      }

      ~DefaultAdjointsModule() {
        cleanAdjoints();
      }

      void initAdjointsModule() {
        // Nothing to do
      }

    protected:

      CODI_INLINE GradientValue* getAdjoints() const {
        return adjoints;
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
      CODI_INLINE void addAdjointValues(TapeValues& values) const {

        size_t nAdjoints      = cast().indexHandler.getMaximumGlobalIndex() + 1;
        double memoryAdjoints = static_cast<double>(nAdjoints) * static_cast<double>(sizeof(GradientValue)) * BYTE_TO_MB;

        values.addSection("Adjoint vector");
        values.addData("Number of adjoints", nAdjoints);
        values.addData("Memory allocated", memoryAdjoints, true, true);

        cast().indexHandler.addValues(values);
      }

      /**
       * @brief Helper function: Sets the adjoint vector to a new size.
       *
       * @param[in] size The new size for the adjoint vector.
       */
      CODI_NO_INLINE void resizeAdjoints(const Index& size) {
        Index oldSize = getAdjointsSize();
        setAdjointsSize(size);

        for(Index i = size; i < oldSize; ++i) {
          adjoints[i].~GradientValue();
        }

        adjoints = (GradientValue*)realloc((void*)adjoints, sizeof(GradientValue) * static_cast<size_t>(size));

        if(NULL == adjoints) {
          throw std::bad_alloc();
        }

        for(Index i = oldSize; i < size; ++i) {
          new (adjoints + i) GradientValue();
        }
      }

      /**
       * @brief Resize the adjoint vector such that it fits the number of indices.
       */
      CODI_INLINE void resizeAdjointsToIndexSize() {
        if(getAdjointsSize() <= cast().indexHandler.getMaximumGlobalIndex()) {
          resizeAdjoints(cast().indexHandler.getMaximumGlobalIndex() + 1);
        }
      }

      /**
       * @brief Helper function: Deletes all arrays
       */
      CODI_NO_INLINE void cleanAdjoints() {
        if(adjointsValid()) {
          free(adjoints);
          adjoints = NULL;
          setAdjointsSize(0);
        }
      }

      CODI_INLINE bool adjointsValid() const {
        return adjoints != NULL;
      }

      /**
       * @brief Swap the data of the tape base module with the data of the other tape base module.
       *
       * @param[in] other  The object with the other tape base module.
       */
      void swapAdjointsModule(Tape& other) {
        std::swap(adjoints, other.adjoints);
        std::swap(adjointsSize, other.adjointsSize);

        // the index handler is not swaped because it is either swaped in the recursive call to of the data vectors
        // or it is handled by the including class
      }

    public:
      // access to adjoints
      CODI_INLINE Index getAdjointsSize() const {
        return adjointsSize;
      }

    public:
      // no boundary check access for derived class

      template<typename AdjointData>
      static CODI_INLINE void setAdjoint(const Index& index, const GradientValue& value, AdjointData* data) {
        data[index] = value;
      }

      CODI_INLINE void setAdjoint(const Index& index, const GradientValue& value) {
        setAdjoint(index, value, adjoints);
      }

      template<typename AdjointData>
      static CODI_INLINE void incrementAdjoint(const Index& index, const AdjointData& adj, const Real& jacobi, AdjointData* data) {
        data[index] += adj * jacobi;
      }

      CODI_INLINE void incrementAdjoint(const Index& index, const GradientValue& adj, const GradientValue& jacobi) {
        incrementAdjoint(index, adj, jacobi, adjoints);
      }

      template<typename AdjointData>
      static CODI_INLINE void incrementTangent(AdjointData& adj, const AdjointData* data, const Index& index, const Real& jacobi) {
        adj += data[index] * jacobi;
      }

      template<typename AdjointData>
      static CODI_INLINE void clearAdjoint(const Index& index, AdjointData* data) {
        setAdjoint(index, GradientValue(), data);
      }

      CODI_INLINE void clearAdjoint(const Index& index) {
        clearAdjoint(index, adjoints);
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
          return adjoints[index];
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
      CODI_INLINE GradientValue& gradient(const Index& index) {
        codiAssert(0 != index);
        codiAssert(index <= cast().indexHandler.getMaximumGlobalIndex());

        //TODO: Add error when index is bigger than expression count
        if(getAdjointsSize() <= index) {
          resizeAdjoints(cast().indexHandler.getMaximumGlobalIndex() + 1);
        }

        return adjoints[index];
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
          return adjoints[0];
        } else {
          return adjoints[index];
        }
      }

      /**
       * @brief Sets all adjoint/gradients to zero.
       */
      CODI_INLINE void clearAdjoints(){
        if(adjointsValid()) {
          for(Index i = 0; i < getAdjointsSize(); ++i) {
            clearAdjoint(i);
          }
        }
      }

      /**
       * @brief Clear the adjoint vector and delete it.
       */
      void CODI_INLINE deleteAdjointVector() {
        cleanAdjoints();
      }
  };
}
