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

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template <typename Tape, typename ChildVector>
  class StatementModule {

		public:

    typedef typename ChildVector::Position ChildPosition;

    /** @brief The data for each statement. */
    typedef Chunk1<StatementInt> Chunk;
    /** @brief The chunk vector for the statement data. */
    typedef ChunkVector<Chunk, ChildVector> Vector;

    typedef typename Vector::Position Position;

    /** @brief The data for the statements. */
    Vector vector;

    StatementModule(ChildVector& childVector) :
      vector(DefaultChunkSize, childVector) {
      }

    inline const Tape& cast() const {
      return static_cast<const Tape&>(*this);
    }

    inline Tape& cast() {
      return static_cast<Tape&>(*this);
    }

    /**
     * @brief Set the size of the statement data chunks.
     *
     * @param[in] statementChunkSize The new size for the statement data chunks.
     */
    void setStatementChunkSize(const size_t& statementChunkSize) {
      vector.setChunkSize(statementChunkSize);
    }

    /**
     * @brief Return the number of used statements.
     * @return The number of used statements.
     */
    size_t getUsedStatementsSize() {
      return vector.getDataSize();
    }

    void resize(const size_t& statementSize) {
      vector.resize(statementSize);
    }

    /**
     * @brief Evaluate a part of the statement vector.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the jacobi vector.
     *
     * @param[in] start The starting point for the statement vector.
     * @param[in]   end The ending point for the statement vector.
     */
    inline void evaluateStmt(const Position& start, const Position& end) {
      StatementInt* statementData;
      size_t dataPos = start.data;
      ChildPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(statementData) = vector.getDataAtPosition(curChunk, 0);

        ChildPosition endInnerPos = vector.getInnerPosition(curChunk);
        cast().evalStmtCallback(curInnerPos, endInnerPos, dataPos, statementData);

        curInnerPos = endInnerPos;

        dataPos = vector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(statementData) = vector.getDataAtPosition(end.chunk, 0);
      cast().evalStmtCallback(curInnerPos, end.inner, dataPos, statementData);
    }

    /**
     * @brief Prints statistics about the tape on the screen
     *
     * Prints information such as stored statements/adjoints and memory usage on screen.
     */
    void printStatistics(){
      size_t nChunksStmts  = vector.getNumChunks();
      size_t totalStmts    = (nChunksStmts-1)*vector.getChunkSize()
                             +vector.getChunkUsedData(nChunksStmts-1);
      double  memoryUsedStmts = (double)totalStmts*(double)sizeof(StatementInt)* BYTE_TO_MB;
      double  memoryAllocStmts= (double)nChunksStmts*(double)vector.getChunkSize()
                                *(double)sizeof(StatementInt)* BYTE_TO_MB;
      std::cout << std::endl
                << "Statements " << std::endl
                << "-------------------------------------" << std::endl
                << "  Number of Chunks: " << std::setw(10) << nChunksStmts << std::endl
                << "  Total Number:     " << std::setw(10) << totalStmts   << std::endl
                << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryAllocStmts << " MB" << std::endl
                << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryUsedStmts << " MB" << std::endl;

    }
  };
}
