/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <fstream>
#include <iostream>
#include <vector>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../tapes/statementEvaluators/statementEvaluatorInterface.hpp"
#include "../identifierCacheOptimizer.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Write connectivity data for a tape.
   *
   * Creates two files
   *  - \<name\>_rhs.dat
   *  - \<name\>_lhs.dat
   *  the both contain two columns of data. The first column is the statement id. The second one contains either the
   *  rhs or lhs identifier.
   *
   * @tparam T_Tape  Tape tape on which the optimization is applied.
   */
  template<typename T_Tape>
  struct WriteConnectivityData : public ApplyIdentifierModification<T_Tape, WriteConnectivityData<T_Tape>> {
      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_TAPE);  ///< See WriteConnectivityData.

      using Base = ApplyIdentifierModification<T_Tape, WriteConnectivityData<T_Tape>>;  ///< Base class abbreviation.

      using Real = typename Tape::Real;              ///< See FullTapeInterface.
      using Identifier = typename Tape::Identifier;  ///< See FullTapeInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See FullTapeInterface.

      std::ofstream rhsConnectivity = {};  ///< Output stream for rhs connectivity.
      std::ofstream lhsConnectivity = {};  ///< Output stream for lhs connectivity.

      Identifier stmtId = 0;  ///< counter for statement data.

      Tape& tape;  ///< Tape that is modified.

      /// Constructor.
      WriteConnectivityData(Tape& t) : Base(t), tape(t) {}

      /// Write to rhs connectivity.
      CODI_INLINE void applyToInput(Identifier& id) {
        rhsConnectivity << stmtId << " " << id << "\n";
      }

      /// Write to lhs connectivity.
      CODI_INLINE void applyToOutput(Identifier& id) {
        lhsConnectivity << stmtId << " " << id << "\n";
      }

      /// Prepare for next statement.
      CODI_INLINE void applyPostOutputLogic() {
        stmtId += 1;
      }

      /// Create the two files. See the class description for details.
      void eval(std::string const& name) {
        rhsConnectivity.open(name + "_rhs.dat");
        lhsConnectivity.open(name + "_lhs.dat");

        tape.iterateForward(*this);

        rhsConnectivity.close();
        lhsConnectivity.close();
      }
  };
}
