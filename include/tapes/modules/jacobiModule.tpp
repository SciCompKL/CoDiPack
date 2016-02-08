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

#ifndef CHILD_VECTOR_TYPE
  #error Please define the type of the child vector
#endif
#ifndef TAPE_NAME
  #error Please define the name of the tape.
#endif

		public:

    typedef CHILD_VECTOR_TYPE JacobiChildVector;
    typedef typename JacobiChildVector::Position JacobiChildPosition;

    typedef Chunk2< Real, IndexType> JacobiChunk;
    /** @brief The chunk vector for the jacobi data. */
    typedef ChunkVector<JacobiChunk, JacobiChildVector > JacobiVector;

    typedef typename JacobiVector::Position JacobiPosition;

    /** @brief The data for the jacobies of each statements. */
    JacobiVector jacobiVector;

    /**
     * @brief Set the size of the jacobi data chunks.
     *
     * @param[in] dataChunkSize The new size for the jacobi data chunks.
     */
    void setDataChunkSize(const size_t& dataChunkSize) {
      jacobiVector.setChunkSize(dataChunkSize);
    }

    /**
     * @brief Return the number of used data entries.
     * @return The number of used data entries.
     */
    size_t getUsedDataEntriesSize() {
      return jacobiVector.getDataSize();
    }

    void resizeJacobi(const size_t& dataSize) {
      jacobiVector.resize(dataSize);
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
        this->jacobiVector.setDataAndMove(std::make_tuple(1.0, index));
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
            this->jacobiVector.setDataAndMove(std::make_tuple(jacobi, index));
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
    template<typename ... Args>
    inline void evaluateJacobies(const JacobiPosition& start, const JacobiPosition& end, Args&&... args) {
      Real* jacobiData;
      IndexType* indexData;
      size_t dataPos = start.data;
      JacobiChildPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(jacobiData, indexData) = jacobiVector.getDataAtPosition(curChunk, 0);

        JacobiChildPosition endInnerPos = jacobiVector.getInnerPosition(curChunk);
        evalJacobiesCallback(curInnerPos, endInnerPos, dataPos, jacobiData, indexData, std::forward<Args>(args)...);

        curInnerPos = endInnerPos;

        dataPos = jacobiVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(jacobiData, indexData) = jacobiVector.getDataAtPosition(end.chunk, 0);
      evalJacobiesCallback(curInnerPos, end.inner, dataPos, jacobiData, indexData, std::forward<Args>(args)...);
    }

    /**
     * @brief Prints statistics about the tape on the screen
     *
     * Prints information such as stored statements/adjoints and memory usage on screen.
     */
    void printJacobiStatistics(){
      size_t nChunksData  = jacobiVector.getNumChunks();
      size_t totalData    = (nChunksData-1)*jacobiVector.getChunkSize()
                             +jacobiVector.getChunkUsedData(nChunksData-1);
      double  memoryUsedData = (double)totalData*(double)(sizeof(Real)+sizeof(IndexType))* BYTE_TO_MB;
      double  memoryAllocData= (double)nChunksData*(double)jacobiVector.getChunkSize()
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

#undef CHILD_VECTOR_TYPE
