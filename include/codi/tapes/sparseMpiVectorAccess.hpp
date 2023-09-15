/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include <cstddef>
#include <map>
#include <vector>
#include <mpi.h>

#include "../config.h"
#include "../misc/macros.hpp"
#include "sparseEvaluation.hpp"
#include "misc/vectorAccessInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Real, typename T_Identifier>
  struct SparseMPIVectorAccess : public VectorAccessInterface<T_Real, T_Identifier> {
    public:

      static int constexpr MASK = 0x42FF42FF;

      using Real = CODI_DD(T_Real, double);           ///< See SparseMPIVectorAccess.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See SparseMPIVectorAccess.

      using NodeDependencies = SparseEvaluation::NodeDependencies<Real, Identifier>;
      using DependencyMap = SparseEvaluation::DependencyMap<Real, Identifier>;

      DependencyMap& dependencies;
      SparseEvaluation::ElemeniationMissingOutput missingOutputHandling;

      DependencyMap outputDependencies;
      DependencyMap inputDependencies;

      MPI_Comm mpiComm;
      int mpiRank;
      int mpiSize;
      int identifierMPIRankOffset;
      int identifierMPIRankMask;

      Real temp;

      int mpiIndex;

      SparseMPIVectorAccess(DependencyMap& dependencies,
                            SparseEvaluation::ElemeniationMissingOutput missingOutputHandling,
                            MPI_Comm mpiComm) :
            dependencies(dependencies),
            missingOutputHandling(missingOutputHandling),
            outputDependencies(),
            inputDependencies(),
            mpiComm(mpiComm),
            mpiRank(),
            mpiSize(),
            identifierMPIRankOffset(),
            identifierMPIRankMask(),
            temp(),
            mpiIndex(1)
      {
        MPI_Comm_size(mpiComm, &mpiSize);
        MPI_Comm_rank(mpiComm, &mpiRank);

        identifierMPIRankOffset = (int)std::ceil(std::log2((double)mpiSize));
        identifierMPIRankMask = (1 << identifierMPIRankOffset) - 1;
      }

      /*******************************************************************************/
      /// @name Misc


      size_t getVectorSize() const {return 1; }
      bool isLhsZero() {return false; }
      SparseMPIVectorAccess* clone() const {
        CODI_EXCEPTION("Not Supported.");
        return nullptr;
      }

      /*******************************************************************************/
      /// @name Indirect adjoint access

      void setLhsAdjoint(Identifier const& index) {
        CODI_UNUSED(index);
        CODI_EXCEPTION("Not used.");
      }

      void updateAdjointWithLhs(Identifier const& index, Real const& jacobian) {
        CODI_UNUSED(index, jacobian);
        CODI_EXCEPTION("Not used.");
      }

      /*******************************************************************************/
      /// @name Indirect tangent access

      void setLhsTangent(Identifier const& index) {
        CODI_UNUSED(index);
        CODI_EXCEPTION("Not used.");
      }

      void updateTangentWithLhs(Identifier const& index, Real const& jacobian) {
        CODI_UNUSED(index, jacobian);
        CODI_EXCEPTION("Not used.");
      }

      /*******************************************************************************/
      /// @name Direct adjoint access

      void resetAdjoint(Identifier const& index, size_t dim) {
        CODI_UNUSED(index, dim);

        // Done in getAdjoint.
      }

      void resetAdjointVec(Identifier const& index) {
        resetAdjoint(index, 0);
      }

      Real getAdjoint(Identifier const& index, size_t dim) {
        CODI_UNUSED(dim);

        SparseEvaluation::NodeDependencies<Real, Identifier> incomingDependencies;

        if(SparseEvaluation::getIncomingDependencies(dependencies, index, incomingDependencies, missingOutputHandling)) {
          Identifier outgoingIndex = generateIdentifier();

          outputDependencies[outgoingIndex] = incomingDependencies;
          return maskIdentifier(outgoingIndex);
        } else {
          return Real();
        }
      }

      void getAdjointVec(Identifier const& index, Real* const vec) {
        vec[0] = getAdjoint(index, 0);
      }

      Real const* getAdjointVec(Identifier const& index) {
        getAdjointVec(index, &temp);

        return &temp;
      }

      void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        CODI_UNUSED(dim);

        Identifier incommingIndex = unmaskIdentifier(adjoint);

        inputDependencies[incommingIndex] = NodeDependencies();

        dependencies[index][incommingIndex] += 1.0; // Just note the dependency here.
      }

      void updateAdjointVec(Identifier const& index, Real const* const vec) {
        updateAdjoint(index, 0, vec[0]);
      }

      /*******************************************************************************/
      /// @name Primal access

      void setPrimal(Identifier const& index, Real const& primal) {
        CODI_UNUSED(index, primal);
      }
      Real getPrimal(Identifier const& index) {
        CODI_UNUSED(index);
        return Real();
      }

      bool hasPrimals() {
        return false;
      }

      std::vector<Identifier> communicateDependenciesStage(std::vector<Identifier> const& requested) {
        // Step 1: Communicate number of input dependencies.
        std::vector<int> inputDependenciesOnRanksCount(mpiSize);
        std::vector<int> inputDependenciesOnRanksDispl(mpiSize);
        inputDependenciesOnRanksCount[mpiRank] = requested.size();
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, inputDependenciesOnRanksCount.data(), 1, MPI_INT, mpiComm);

        size_t totalInputDependencies = computeTotalAndDisplacements(inputDependenciesOnRanksCount, inputDependenciesOnRanksDispl);

        // Step 2: Communicate the input dependencies.
        std::vector<int> inputDependenciesOnRanks(totalInputDependencies);

        // Set our dependencies into the vector.
        int pos = inputDependenciesOnRanksDispl[mpiRank];
        for(auto const& cur : requested) {
          inputDependenciesOnRanks[pos] = cur;
          pos += 1;
        }
        MPI_Allgatherv(MPI_IN_PLACE, (int)inputDependencies.size(), MPI_INT, inputDependenciesOnRanks.data(), inputDependenciesOnRanksCount.data(), inputDependenciesOnRanksDispl.data(), MPI_INT, mpiComm);

        // Step 3: Compute for each rank what we have to send and receive what each rank will send us.
        std::vector<int> outputDependenciesOnRanksSendCount(mpiSize);
        std::vector<int> outputDependenciesOnRanksSendDispl(mpiSize);
        std::vector<int> outputDependenciesOnRanksRecvCount(mpiSize);
        std::vector<int> outputDependenciesOnRanksRecvDispl(mpiSize);
        for(int curRank = 0; curRank < mpiSize; curRank += 1) {
          if(curRank == mpiRank) {
            continue;  // Do not count for own rank.
          }

          int curOffset = inputDependenciesOnRanksDispl[curRank];
          for(int curPos = 0; curPos < inputDependenciesOnRanksCount[curRank]; curPos += 1) {
            int curIdentifier = inputDependenciesOnRanks[curOffset + curPos];
            int curInputRank = rankFromIdentifier(curIdentifier);

            if(curInputRank == mpiRank) {
              // This process got an identifier from us, we will send him our dependencies.
              outputDependenciesOnRanksSendCount[curRank] += outputDependencies[curIdentifier].size();
            }
          }
        }
        MPI_Alltoall(outputDependenciesOnRanksSendCount.data(), 1, MPI_INT, outputDependenciesOnRanksRecvCount.data(), 1, MPI_INT, mpiComm);

        size_t totalOutputDependenciesSend = computeTotalAndDisplacements(outputDependenciesOnRanksSendCount, outputDependenciesOnRanksSendDispl);
        size_t totalOutputDependenciesRecv = computeTotalAndDisplacements(outputDependenciesOnRanksRecvCount, outputDependenciesOnRanksRecvDispl);

        // Step 4: Send output dependencies.
        std::vector<int> outputDependenciesOnRanksSendInputIdentifier(totalOutputDependenciesSend);
        std::vector<int> outputDependenciesOnRanksSendOutputIdentifier(totalOutputDependenciesSend);
        std::vector<double> outputDependenciesOnRanksSendOutputJacobian(totalOutputDependenciesSend);
        std::vector<int> outputDependenciesOnRanksRecvInputIdentifier(totalOutputDependenciesRecv);
        std::vector<int> outputDependenciesOnRanksRecvOutputIdentifier(totalOutputDependenciesRecv);
        std::vector<double> outputDependenciesOnRanksRecvOutputJacobian(totalOutputDependenciesRecv);

        // Set our own dependencies in the vector.
        for(int curRank = 0; curRank < mpiSize; curRank += 1) {
          if(curRank == mpiRank) {
            continue;  // Do set in own rank.
          }

          int sendPos = outputDependenciesOnRanksSendDispl[curRank];

          int curOffset = inputDependenciesOnRanksDispl[curRank];
          for(int curPos = 0; curPos < inputDependenciesOnRanksCount[curRank]; curPos += 1) {
            int curIdentifier = inputDependenciesOnRanks[curOffset + curPos];
            int curInputRank = rankFromIdentifier(curIdentifier);

            if(curInputRank == mpiRank) {
              // This process got an identifier from us, we will send him our dependencies.

              NodeDependencies const& nodeDependencies = outputDependencies[curIdentifier];
              for(auto curDep : nodeDependencies) {
                outputDependenciesOnRanksSendInputIdentifier[sendPos] = curIdentifier;
                outputDependenciesOnRanksSendOutputIdentifier[sendPos] = curDep.first;
                outputDependenciesOnRanksSendOutputJacobian[sendPos] = curDep.second;

                sendPos += 1;
              }
            }
          }
        }
        MPI_Alltoallv(outputDependenciesOnRanksSendInputIdentifier.data(), outputDependenciesOnRanksSendCount.data(), outputDependenciesOnRanksSendDispl.data(), MPI_INT,
                      outputDependenciesOnRanksRecvInputIdentifier.data(), outputDependenciesOnRanksRecvCount.data(), outputDependenciesOnRanksRecvDispl.data(), MPI_INT,
                      mpiComm);
        MPI_Alltoallv(outputDependenciesOnRanksSendOutputIdentifier.data(), outputDependenciesOnRanksSendCount.data(), outputDependenciesOnRanksSendDispl.data(), MPI_INT,
                      outputDependenciesOnRanksRecvOutputIdentifier.data(), outputDependenciesOnRanksRecvCount.data(), outputDependenciesOnRanksRecvDispl.data(), MPI_INT,
                      mpiComm);
        MPI_Alltoallv(outputDependenciesOnRanksSendOutputJacobian.data(), outputDependenciesOnRanksSendCount.data(), outputDependenciesOnRanksSendDispl.data(), MPI_DOUBLE,
                      outputDependenciesOnRanksRecvOutputJacobian.data(), outputDependenciesOnRanksRecvCount.data(), outputDependenciesOnRanksRecvDispl.data(), MPI_DOUBLE,
                      mpiComm);

        // Step 5: Resolve requested dependencies.
        std::vector<Identifier> newlyUnresolved = {};
        Identifier lastIdentifier = 0;
        NodeDependencies* curNodeDependencies = nullptr;
        for(size_t curPos = 0; curPos < totalOutputDependenciesRecv; curPos += 1) {

          // Update the pointer to the current node dependencies.
          if(lastIdentifier != outputDependenciesOnRanksRecvInputIdentifier[curPos]) {
            lastIdentifier = outputDependenciesOnRanksRecvInputIdentifier[curPos];
            curNodeDependencies = &inputDependencies[lastIdentifier];
          }

          // Inserted all dependencies into the mpi  input dependencies.


          Identifier curOutputIdentifier = outputDependenciesOnRanksRecvOutputIdentifier[curPos];
          (*curNodeDependencies)[curOutputIdentifier] += outputDependenciesOnRanksRecvOutputJacobian[curPos];

          if(curOutputIdentifier < 0) {
            if(mpiRank == rankFromIdentifier(curOutputIdentifier)) {
              // Depency from our rank. Just add it as a new mpi dependency.
              inputDependencies[curOutputIdentifier] = outputDependencies[curOutputIdentifier];
            } else {
              // Dependency from another rank. Add it to the unresolved dependencies and update the dependency map.
              newlyUnresolved.push_back(curOutputIdentifier);
            }
          }
        }

        return newlyUnresolved;
      }

      void communicateDependencies() {
        std::vector<Identifier> required = {};
        required.reserve(inputDependencies.size());

        for(auto const& cur : inputDependencies) {
          required.push_back(cur.first);
        }

        int maxRemaining = 0;
        do {
          std::vector<Identifier> remainingDependencies = communicateDependenciesStage(required);

          maxRemaining = remainingDependencies.size();
          MPI_Allreduce(MPI_IN_PLACE, &maxRemaining, 1, MPI_INT, MPI_MAX, mpiComm);

          required = std::move(remainingDependencies);
        } while(maxRemaining != 0);
      }


      int rankFromIdentifier(Identifier index) {
        int rank = index;
        rank = -rank;                        // Undo negative value.
        rank = rank & identifierMPIRankMask; // Mask away the index.

        return rank;
      }


      Identifier indexFromIdentifier(Identifier index) {
        int id = index;
        id = -id;                           // Undo negative value.
        id = id >> identifierMPIRankOffset; // Shift away the rank.

        return id;
      }

    private:

      size_t computeTotalAndDisplacements(std::vector<int> const& vec, std::vector<int>& displacements) {
        size_t total = 0;
        for(size_t i = 0; i < vec.size(); i += 1) {
          displacements[i] = total;
          total += vec[i];
        }

        return total;
      }

      Identifier generateIdentifier(Identifier identifier, int rank) {
        Identifier index = identifier;  // Basic index is the current counter.
        index = index << identifierMPIRankOffset;   // Shift the bytes for the mpi rank.
        index = index +  rank;     // Set the mpi rank in the lower ranks.
        index = -index;               // Identify mpi indices with negative values.

        return index;
      }

      Identifier generateIdentifier() {
        mpiIndex += 1;
        return generateIdentifier(mpiIndex, mpiRank);
      }

      Real maskIdentifier(Identifier index) {
        codiAssert(sizeof(double) == sizeof(long int));
        codiAssert(2 * sizeof(Identifier) == sizeof(long int));

        int masked[2] = {index, MASK};

        return *reinterpret_cast<Real*>(masked);
      }

      Identifier unmaskIdentifier(Real maskedIndex) {
        int* unmasked = reinterpret_cast<int*>(&maskedIndex);

        if(MASK != unmasked[1]) {
          CODI_EXCEPTION("Adjoint was modified in MPI communication.");
        }

        return unmasked[0];
      }

  };
}
