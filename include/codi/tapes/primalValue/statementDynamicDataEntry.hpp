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

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../data/dataInterface.hpp"
#include "../misc/statementSizes.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /** Helper structure for reading and writing the dynamic statement data.
   *
   *  A call to readForward or readReverse will create a data structure in a static context and populate the pointers
   *  with references from the byte array.
   *
   *  A call to first reserve and then to store will first reserve the data on the data stream and afterwards populate
   *  the pointers with references to the data stream.
   *
   *  @tparam T_Real  Primal computation type, e.g. double.
   *  @tparam T_Identifier  Identifier for the internal management, e.g. int.
   *  @tparam T_StatementCallDefaultArguments  Structure representing the default arguments for the reverse evaluations.
   *  @tparam T_DynamicSizeData  Data stream that stores the dynamic data.
   *  @tparam T_LinearIndexHandling  If linear index handling or reuse index handling is used.
   */
  template<typename T_Real, typename T_Identifier, typename T_StatementCallDefaultArguments, typename T_DynamicSizeData,
           bool T_LinearIndexHandling>
  struct StatementDynamicDataEntry {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See StatementDynamicDataEntry.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See StatementDynamicDataEntry.
      using StatementCallDefaultArguments =
          CODI_DD(T_StatementCallDefaultArguments,
                  CODI_T(StatementCallDefaultArgumentsBase<Identifier>));  ///< See StatementDynamicDataEntry.
      using DynamicSizeData = CODI_DD(T_DynamicSizeData, DataInterface);   ///< See StatementDynamicDataEntry.
      static bool constexpr LinearIndexHandling = T_LinearIndexHandling;   ///< See StatementDynamicDataEntry.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      Real* __restrict__ oldPrimalValues;       ///< Overwritten primal values. Only used for reuse index handling.
      Identifier* __restrict__ lhsIdentifiers;  ///< Lhs identifiers. Only used for reuse index handling.

      PassiveReal* __restrict__ constantValues;  ///< Constant values in the statement.
      Identifier* __restrict__ rhsIdentifiers;   ///< Rhs identifiers.
      Real* __restrict__ passiveValues;          ///< Passive values from active arguments of the statement.

      /// The dynamic data position in stmtArgs is updated.
      CODI_INLINE static StatementDynamicDataEntry readForward(StatementSizes const& stmtSizes,
                                                               StatementCallDefaultArguments& stmtArgs) {
        return readForward(stmtSizes, stmtArgs.dynamicSizeValues, stmtArgs.curDynamicSizePos,
                           stmtArgs.numberOfPassiveArguments);
      }

      /// The dynamic data position in stmtArgs is updated.
      CODI_INLINE static StatementDynamicDataEntry readReverse(StatementSizes const& stmtSizes,
                                                               StatementCallDefaultArguments& stmtArgs) {
        StatementDynamicDataEntry data;

        size_t pos = stmtArgs.curDynamicSizePos;

        if (!LinearIndexHandling) {
          pos -= stmtSizes.outputArgs * sizeof(Real);
          data.oldPrimalValues = ((Real*)(&stmtArgs.dynamicSizeValues[pos]));
          pos -= stmtSizes.outputArgs * sizeof(Identifier);
          data.lhsIdentifiers = ((Identifier*)(&stmtArgs.dynamicSizeValues[pos]));
        }

        pos -= stmtSizes.constantArgs * sizeof(PassiveReal);
        data.constantValues = (PassiveReal*)(&stmtArgs.dynamicSizeValues[pos]);
        pos -= stmtSizes.inputArgs * sizeof(Identifier);
        data.rhsIdentifiers = (Identifier*)(&stmtArgs.dynamicSizeValues[pos]);
        pos -= stmtArgs.numberOfPassiveArguments * sizeof(Real);
        data.passiveValues = (Real*)(&stmtArgs.dynamicSizeValues[pos]);

        stmtArgs.curDynamicSizePos = pos;

        return data;
      }

      /// Overwrite the old primal values with new ones. Usually called in a forward tape evaluation.
      CODI_INLINE void updateOldPrimalValues(StatementSizes const& stmtSizes, StatementCallDefaultArguments& revArgs,
                                             Real* __restrict__ primalVector) {
        for (size_t iLhs = 0; iLhs < stmtSizes.outputArgs; iLhs += 1) {
          Identifier const lhsIdentifier = revArgs.getLhsIdentifier(iLhs, lhsIdentifiers);
          oldPrimalValues[iLhs] = primalVector[lhsIdentifier];
        }
      }

      /// Copies the passive values into the first entries of the primal value vector.
      CODI_INLINE void copyPassiveValuesIntoPrimalVector(StatementCallDefaultArguments& revArgs,
                                                         Real* __restrict__ primalVector) {
        for (Config::ArgumentSize curPos = 0; curPos < revArgs.numberOfPassiveArguments; curPos += 1) {
          primalVector[curPos] = passiveValues[curPos];
        }
      }

      /// Reserve the data for the dynamic size data stream.
      CODI_INLINE static void reserve(DynamicSizeData& data, StatementSizes const& sizes,
                                      size_t const& activeArguments) {
        data.reserveItems(computeSize(sizes, activeArguments));
      }

      /// Store the data for the dynamic data stream.
      /// @return The pointers need to be populated with the data.
      CODI_INLINE static StatementDynamicDataEntry store(DynamicSizeData& data, StatementSizes const& sizes,
                                                         size_t const& activeArguments) {
        char* dataPointer = nullptr;
        data.getDataPointers(dataPointer);
        data.addDataSize(computeSize(sizes, activeArguments));

        size_t pos = 0;
        return readForward(sizes, dataPointer, pos, sizes.inputArgs - activeArguments);
      }

    private:

      /// Common read forward function for static context and data reservation.
      CODI_INLINE static StatementDynamicDataEntry readForward(StatementSizes const& stmtSizes,
                                                               char* __restrict__ dynamicSizeValues,
                                                               size_t& __restrict__ curDynamicSizePos,
                                                               size_t const& __restrict__ numberOfPassiveArguments) {
        StatementDynamicDataEntry data;

        size_t pos = curDynamicSizePos;

        data.passiveValues = (Real*)(&dynamicSizeValues[pos]);
        pos += numberOfPassiveArguments * sizeof(Real);
        data.rhsIdentifiers = (Identifier*)(&dynamicSizeValues[pos]);
        pos += stmtSizes.inputArgs * sizeof(Identifier);
        data.constantValues = (PassiveReal*)(&dynamicSizeValues[pos]);
        pos += stmtSizes.constantArgs * sizeof(PassiveReal);

        if (!LinearIndexHandling) {
          data.lhsIdentifiers = ((Identifier*)(&dynamicSizeValues[pos]));
          pos += stmtSizes.outputArgs * sizeof(Identifier);
          data.oldPrimalValues = ((Real*)(&dynamicSizeValues[pos]));
          pos += stmtSizes.outputArgs * sizeof(Real);
        }

        curDynamicSizePos = pos;

        return data;
      }

      /// Size in bytes of the dynamic statement data.
      CODI_INLINE static size_t computeSize(StatementSizes const& sizes, size_t const& activeArguments) {
        size_t dynamicSize = (sizes.inputArgs - activeArguments) * sizeof(Real) +
                             (sizes.inputArgs) * sizeof(Identifier) + sizes.constantArgs * sizeof(PassiveReal);
        if (!LinearIndexHandling) {
          dynamicSize += sizes.outputArgs * (sizeof(Identifier) + sizeof(Real));
        }

        return dynamicSize;
      }
  };
}
