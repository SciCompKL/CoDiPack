/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
   * @brief  Generates a file with the math representation for each of the statements.
   *
   * An example of a statement in the generated file:
   * <tt>T33 = X37*X6</tt>
   *
   * X, Y and T are used to indicate an input, output or a temporary variable respectively.
   *
   * See codi::CommonTextTapeWriter for the methods that are used in this class.
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct MathRepWriter : public CommonTextTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Base = CommonTextTapeWriter<T_Type>;                  ///< See CommonTextTapeWriter.

      using Tape = typename Type::Tape;              ///< See TapeWriterInterface.
      using Identifier = typename Type::Identifier;  ///< See TapeWriterInterface.
      using Real = typename Type::Real;              ///< See TapeWriterInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See TapeWriterInterface.

      /// Constructor.
      MathRepWriter(std::string const& name, std::vector<Identifier> const& in, std::vector<Identifier> const& out)
          : Base(false, name, in, out) {};

      /**
       * \copydoc codi::TapeWriterInterface::writeStatement(WriteInfo const&, Identifier const&, Real const&,
       *                                                    Config::ArgumentSize const&, size_t const&,
       *                                                    Identifier const* const, size_t const&,
       *                                                    Real const* const, size_t&, Real const* const,
       *                                                    EvalHandle)
       */
      void writeStatement(WriteInfo const& info, Identifier const& curLhsIdentifier, Real const& primalValue,
                          Config::ArgumentSize const& nPassiveValues, size_t const& curRhsIdentifiersPos,
                          Identifier const* const rhsIdentifiers, size_t const& curPassiveValuePos,
                          Real const* const passiveValues, size_t& curConstantPos, Real const* const constantValues,
                          EvalHandle stmtEvalHandle) {
        if (nPassiveValues == Config::StatementInputTag) CODI_Unlikely {
          // Do nothing.
        } else CODI_Likely {
          // The mathRep string is modified to include the identifier and value.
          std::string mathRep =
              Base::modifyMathRep(info.mathRepresentation, curLhsIdentifier, &rhsIdentifiers[curRhsIdentifiersPos],
                                  info.numberOfActiveArguments);
          // Print the statement string.
          fprintf(this->fileHandleGraph, "%s\n", mathRep.c_str());
        }
      }
  };
}
