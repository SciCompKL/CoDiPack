/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <tuple>

#include "../chunk.hpp"
#include "../chunkVector.hpp"
#include "../externalFunctions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template <typename Tape, typename ChildVector, typename Real, typename IndexType>
  class JacobiModule {
		public:

    typedef typename ChildVector::Position ChildPosition;
//    typedef typename Tape::IndexType IndexType;
//    typedef typename Tape::Real Real;

    typedef Chunk2< Real, IndexType> Chunk;
    /** @brief The chunk vector for the jacobi data. */
    typedef ChunkVector<Chunk, ChildVector > Vector;

    typedef typename Vector::Position Position;

    /** @brief The data for the jacobies of each statements. */
    Vector vector;

    JacobiModule(ChildVector& childVector) :
      vector(DefaultChunkSize, childVector) {
      }

    inline const Tape& cast() const {
      return static_cast<const Tape&>(*this);
    }

    inline Tape& cast() {
      return static_cast<Tape&>(*this);
    }


    /**
     * @brief Set the size of the jacobi data chunks.
     *
     * @param[in] dataChunkSize The new size for the jacobi data chunks.
     */
    void setDataChunkSize(const size_t& dataChunkSize) {
      vector.setChunkSize(dataChunkSize);
    }

    /**
     * @brief Return the number of used data entries.
     * @return The number of used data entries.
     */
    size_t getUsedDataEntriesSize() {
      return vector.getDataSize();
    }

    void resize(const size_t& dataSize) {
      vector.resize(dataSize);
    }

    /**
     * @brief Stores the jacobi with the value 1.0 on the tape if the index is active.
     *
     * @param[in]  data Not used in this implementation.
     * @param[in] value Not used in this implementation.
     * @param[in] index Used to check if the variable is active.
     *
     * @tparam Data  The type of the data for the tape.
     */
    template<typename Data>
    inline void pushJacobi(Data& data, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      ENABLE_CHECK(OptCheckZeroIndex, 0 != index) {
        this->vector.setDataAndMove(std::make_tuple(1.0, index));
      }
    }

    /**
     * @brief Stores the jacobi on the tape if the index is active.
     *
     * @param[in]   data Not used in this implementation.
     * @param[in] jacobi Stored on the tape if the variable is active.
     * @param[in]  value Not used in this implementation.
     * @param[in]  index Used to check if the variable is active.
     *
     * @tparam Data  The type of the data for the tape.
     */
    template<typename Data>
    inline void pushJacobi(Data& data, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      ENABLE_CHECK(OptCheckZeroIndex, 0 != index) {
        ENABLE_CHECK(OptIgnoreInvalidJacobies, isfinite(jacobi)) {
          ENABLE_CHECK(OptJacobiIsZero, 0.0 != jacobi) {
            this->vector.setDataAndMove(std::make_tuple(jacobi, index));
          }
        }
      }
    }

    /**
     * @brief Evaluate a part of the jacobi vector.
     *
     * It has to hold start >= end.
     *
     * It calls the evaluation method for the expression counter.
     *
     * @param[in]       start The starting point for the jacobi vector.
     * @param[in]         end The ending point for the jacobi vector.
     * @param[inout]  stmtPos The current position in the statement vector. This value is used in the next invocation of this method.
     * @param[in]  statements The pointer to the statement vector.
     */
    inline void evaluateJacobies(const Position& start, const Position& end, size_t& stmtPos, StatementInt* &statementData) {
      Real* jacobiData;
      IndexType* indexData;
      size_t dataPos = start.data;
      ChildPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(jacobiData, indexData) = vector.getDataAtPosition(curChunk, 0);

        ChildPosition endInnerPos = vector.getInnerPosition(curChunk);
        cast().evalJacobiesCallback(curInnerPos, endInnerPos, stmtPos, statementData, dataPos, jacobiData, indexData);

        curInnerPos = endInnerPos;

        dataPos = vector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(jacobiData, indexData) = vector.getDataAtPosition(end.chunk, 0);
      cast().evalJacobiesCallback(curInnerPos, end.inner, stmtPos, statementData, dataPos, jacobiData, indexData);
    }

    /**
     * @brief Prints statistics about the tape on the screen
     *
     * Prints information such as stored statements/adjoints and memory usage on screen.
     */
    void printStatistics(){
      size_t nChunksData  = vector.getNumChunks();
      size_t totalData    = (nChunksData-1)*vector.getChunkSize()
                             +vector.getChunkUsedData(nChunksData-1);
      double  memoryUsedData = (double)totalData*(double)(sizeof(Real)+sizeof(IndexType))* BYTE_TO_MB;
      double  memoryAllocData= (double)nChunksData*(double)vector.getChunkSize()
                                *(double)(sizeof(Real)+sizeof(IndexType))* BYTE_TO_MB;

      std::cout << "-------------------------------------" << std::endl
                << "Jacobi entries "                       << std::endl
                << "-------------------------------------" << std::endl
                << "  Number of Chunks: " << std::setw(10) << nChunksData << std::endl
                << "  Total Number:     " << std::setw(10) << totalData   << std::endl
                << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryAllocData << " MB" << std::endl;
    }
  };
}
