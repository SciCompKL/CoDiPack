/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

/*
 * In order to include this file the user has to define the preprocessor macro CHILD_VECTOR_TYPE and
 * VECTOR_TYPE.
 *
 * CHILD_VECTOR_TYPE defines the type of the nested vector for the data vector.
 * VECTOR_TYPE defines the type of the data vector.
 *
 * All these macros are undefined at the end of the file.
 *
 * The module defines the structures jacobiVector.
 * The module defines the types JacobiChildVector, JacobiChildPosition, JacobiVector, JacobiChunk,
 * JacobiPosition.
 *
 * It defines the methods pushJacobi(1.0), pushJacobi(Mul) printJacobiStatistics from the TapeInterface and ReverseTapeInterface.
 *
 * It defines the methods evaluateJacobies, incrementAdjoints, setDataChunkSize, getUsedJacobiesSize, resizeJacobi as interface functions for the
 * including class.
 */

#ifndef CHILD_VECTOR_TYPE
  #error Please define the type of the child vector
#endif
#ifndef VECTOR_TYPE
  #error Please define the name of the chunk vector type.
#endif

		private:

  // ----------------------------------------------------------------------
  // All definitions of the module
  // ----------------------------------------------------------------------

    /** @brief The child vector for the Jacobi data vector. */
    typedef CHILD_VECTOR_TYPE JacobiChildVector;

    /** @brief The position type of the Jacobi child vector */
    typedef typename JacobiChildVector::Position JacobiChildPosition;

    /** @brief The vector for the jacobi data. */
    typedef VECTOR_TYPE JacobiVector;

    /** @brief The data for the jacobies */
    typedef typename JacobiVector::ChunkType JacobiChunk;

    /** @brief The position type of the jacobi module. */
    typedef typename JacobiVector::Position JacobiPosition;

    /** @brief The data for the jacobies of each statements. */
    JacobiVector jacobiVector;

  private:

  // ----------------------------------------------------------------------
  // Private function for the communication with the including class
  // ----------------------------------------------------------------------

    /**
     * @brief Evaluate a part of the jacobi vector.
     *
     * It has to hold start >= end.
     *
     * It calls the evaluation method for the expression counter.
     *
     * @param[in]    start  The starting point for the jacobi vector.
     * @param[in]      end  The ending point for the jacobi vector.
     * @param[in,out] args  The data from the other vectors.
     *
     * @tparam Args The types for the arguments from the other vectors.
     */
    template<typename ... Args>
    CODI_INLINE void evaluateJacobies(const JacobiPosition& start, const JacobiPosition& end, Args&&... args) {
      Real* jacobiData;
      IndexType* indexData;
      size_t dataPos = start.data;
      JacobiChildPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        jacobiVector.getDataAtPosition(curChunk, 0, jacobiData, indexData);

        JacobiChildPosition endInnerPos = jacobiVector.getInnerPosition(curChunk);
        evalJacobiesCallback(curInnerPos, endInnerPos, dataPos, jacobiData, indexData, std::forward<Args>(args)...);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        curInnerPos = endInnerPos;

        dataPos = jacobiVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      jacobiVector.getDataAtPosition(end.chunk, 0, jacobiData, indexData);
      evalJacobiesCallback(curInnerPos, end.inner, dataPos, jacobiData, indexData, std::forward<Args>(args)...);
    }

    /**
     * @brief Perform the adjoint update of the reverse AD sweep
     *
     * Evaluates the equation
     *
     * \f[ \bar v_i += \frac{\d \phi}{\d v_i} \bar w \f]
     *
     * The \f[ v_i \f] are the arguments of the statement and are taken from the input jacobi and indices.
     * The value \f[ \bar w \f] is taken from the input adj.
     *
     * @param[in]                  adj  The adjoint of the lhs of the statement.
     * @param[in,out]         adjoints  The adjoint vector containing the adjoints of all variables.
     * @param[in]      activeVariables  The number of active arguments on the rhs.
     * @param[int,out]         dataPos  The position inside the jacobi and indices vectors. It is decremented by the number of active variables.
     * @param[in]             jacobies  The jacobies from the arguments of the statement.
     * @param[in]              indices  The indices from the arguments of the statements.
     */
     CODI_INLINE void incrementAdjoints(const GradientValue& adj, GradientValue* adjoints, const StatementInt& activeVariables, size_t& dataPos, Real* &jacobies, IndexType* &indices) {
      ENABLE_CHECK(OptZeroAdjoint, !isTotalZero(adj)){
        for(StatementInt curVar = 0; curVar < activeVariables; ++curVar) {
          --dataPos;
          adjoints[indices[dataPos]] += adj * jacobies[dataPos];

        }
      } else {
        dataPos -= activeVariables;
      }
    }

    /**
     * @brief Resize the jacobi data.
     *
     * Ensure that enough size is allocated such that dataSize number of items
     * can be stored.
     *
     * @param[in] dataSize  The size that should be allocated for the jacobi data.
     */
    void resizeJacobi(const size_t& dataSize) {
      jacobiVector.resize(dataSize);
    }

  public:

  // ----------------------------------------------------------------------
  // Public function from the TapeInterface and ReverseTapeInterface
  // ----------------------------------------------------------------------

    /**
     * @brief Set the size of the jacobi data chunks.
     *
     * @param[in] dataChunkSize The new size for the jacobi data chunks.
     */
    void setDataChunkSize(const size_t& dataChunkSize) {
      jacobiVector.setChunkSize(dataChunkSize);
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
    CODI_INLINE void pushJacobi(Data& data, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      ENABLE_CHECK(OptCheckZeroIndex, 0 != index) {
        this->jacobiVector.setDataAndMove(1.0, index);
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
    CODI_INLINE void pushJacobi(Data& data, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      ENABLE_CHECK(OptCheckZeroIndex, 0 != index) {
        ENABLE_CHECK(OptIgnoreInvalidJacobies, isfinite(jacobi)) {
          ENABLE_CHECK(OptJacobiIsZero, !isTotalZero(jacobi)) {
            this->jacobiVector.setDataAndMove(jacobi, index);
          }
        }
      }
    }

    /**
     * @brief Prints statistics about the stored Jacobie entries.
     *
     * Displays the number of chunks, the total number of jacobies, the
     * allocated memory and the used memory.
     *
     * @param[in,out]   out  The information is written to the stream.
     * @param[in]     hLine  The horizontal line that separates the sections of the output.
     *
     * @tparam Stream The type of the stream.
     */
    template<typename Stream>
    void printJacobiStatistics(Stream& out, const std::string hLine) const {
      size_t nChunksData   = jacobiVector.getNumChunks();
      size_t totalData     = jacobiVector.getDataSize();
      size_t sizeDataEntry = JacobiChunk::EntrySize;

      double  memoryUsedData = (double)totalData*(double)(sizeDataEntry)* BYTE_TO_MB;
      double  memoryAllocData= (double)nChunksData*(double)jacobiVector.getChunkSize()
                                *(double)(sizeDataEntry)* BYTE_TO_MB;

      out << hLine
          << "Jacobi entries\n"
          << hLine
          << "  Total Number:     " << std::setw(10) << totalData   << "\n"
          << "  Number of Chunks: " << std::setw(10) << nChunksData << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedData << " MB" << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocData << " MB" << "\n";
    }

    /**
     * @brief Return the number of used data entries.
     * @return The number of used data entries.
     */
    size_t getUsedDataEntriesSize() const {
      return jacobiVector.getDataSize();
    }

#undef CHILD_VECTOR_TYPE
#undef VECTOR_TYPE
