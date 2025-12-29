/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#include "../../config.h"
#include "commonReaderWriterBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief  Used to register a statement for a Jacobian tape.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be restored.
   */
  template<typename T_Type>
  struct JacobianBaseTapeReader : public CommonBaseTapeReader<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeReaderInterface.
      using Tape = typename Type::Tape;                           ///< See TapeReaderInterface.
      using Identifier = typename Type::Identifier;               ///< See TapeReaderInterface.
      using Real = typename Type::Real;                           ///< See TapeReaderInterface.

      JacobianBaseTapeReader() : CommonBaseTapeReader<T_Type>() {};  ///< Constructor.

      /// Used to register the currently read statement form the tape file onto the new tape.
      void registerStatement(Identifier const& lhsIdentifier, Config::ArgumentSize const& nArgs,
                             std::vector<Identifier> const& rhsIdentifiers, std::vector<Real> const& rhsJacobians,
                             Identifier& lowestIndex, bool& isFirstIdentifier) {
        // Update the lowest identifier for the Linear case. The lowest identifier will be the first identifier.
        if (isFirstIdentifier && Tape::TapeTypes::IsLinearIndexHandler) {
          lowestIndex = lhsIdentifier - 1;
          isFirstIdentifier = false;
        }

        // Record the lhsIdentifier and the number of arguments of the current statement on the tape.
        Identifier lhsIdentifierWithOffset = lhsIdentifier - lowestIndex;

        // Store the rhsIdentifiers and rhsJacobians for the current statement.
        std::vector<Identifier> rhsIdentifierWithOffset;

        if (nArgs == Config::StatementLowLevelFunctionTag) CODI_Unlikely {
          // TODO.
        } else if (nArgs == Config::StatementInputTag) CODI_Unlikely {
          // Do nothing.
        } else CODI_Likely {
          for (size_t rhsCount = 0; rhsCount < nArgs; rhsCount++) {
            rhsIdentifierWithOffset.push_back(rhsIdentifiers[rhsCount] - lowestIndex);
          }
        }

        this->tape.createStatementManual(0, lhsIdentifierWithOffset, nArgs, rhsJacobians.data(),
                                         rhsIdentifierWithOffset.data());
      }
  };
}
