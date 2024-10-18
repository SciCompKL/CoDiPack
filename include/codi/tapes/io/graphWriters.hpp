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
#include "../../traits/tapeTraits.hpp"
#include "commonReaderWriterBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * @brief  Generates a graphical .dot file for a Jacobian tape.
   *
   * An example of how Nodes and Edges are produced in the .dot file:
   *
   * Nodes: <tt>A35_1 [label = "T35"];</tt>
   * The  <tt>T35</tt> indicates that identifier 35 is a temporary variable.
   *
   * Edges: <tt>A35_1 -> A56_2 [label="0.909297"];</tt>
   * The extension is used to record multiple unique nodes for a identifier and the label represents the Jacobian.
   *
   * See codi::CommonTextTapeWriter for the methods that are used in this class.
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct JacobianGraphTapeWriter : public CommonTextTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Base = CommonTextTapeWriter<T_Type>;                  ///< See CommonTextTapeWriter.

      using Tape = typename Type::Tape;              ///< See TapeWriterInterface.
      using Identifier = typename Type::Identifier;  ///< See TapeWriterInterface.
      using Real = typename Type::Real;              ///< See TapeWriterInterface.

      bool printJacobians;  ///< Flag that enables the jacobians on the edges in the graph.

      /// Constructor.
      JacobianGraphTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                              std::vector<Identifier> const& out, bool ifJacobians)
          : Base(true, name, in, out), printJacobians(ifJacobians) {};

      /**
       * \copydoc codi::TapeWriterInterface::writeStatement(Identifier const&, size_t&, Real const* const,
       *                                                    Identifier const* const, Config::ArgumentSize const&)
       */
      void writeStatement(Identifier const& curLhsIdentifier, size_t& curJacobianPos, Real const* const rhsJacobians,
                          Identifier const* const rhsIdentifiers, Config::ArgumentSize const& nJacobians) {
        std::string node;
        if (nJacobians == Config::StatementInputTag) CODI_Unlikely {
          // Do nothing.
        } else CODI_Likely {
          //  Ensure that the all the rhsIdentifiers have been added before connecting edges to them.
          Base::placeUnusedRhsNodes(&rhsIdentifiers[curJacobianPos], nJacobians);

          // Add the curLhsIdentifier node.
          Base::createNode(curLhsIdentifier, this->formatNodeLabel(curLhsIdentifier));

          // Loop through rhsIdentifiers and create the edges. The this->identifierExtensions is used to find the
          // current extension for each of the rhsIdentifiers.
          for (size_t argCount = 0; argCount < nJacobians; argCount++) {
            std::string labelText = "";
            if (printJacobians) {
              labelText = std::to_string(rhsJacobians[curJacobianPos + argCount]);
            }
            Base::createEdge(rhsIdentifiers[curJacobianPos + argCount], curLhsIdentifier, labelText);
          }

          this->identifierExtensions[curLhsIdentifier] += 1;
        }
      }
  };

  /**
   * @brief  Generates a graphical .dot file for a primal value tape. The writer adds the math representation of
   * statements into the node labels.
   *
   * An example of how Nodes and Edges are produced in the .dot file:
   *
   * Nodes: <tt>A35_1 [label = <T33 = X37*X6>]];</tt>
   * The label contains the math rep of the current statement.
   *
   * Edges: <tt>A35_1 -> A56_2;</tt>
   * The extension is used to record multiple unique nodes for a identifier.
   *
   * See codi::CommonTextTapeWriter for the methods that are used in this class.
   * See codi::TapeWriterInterface for a general description on how to use tape writers.
   *
   * @tparam T_Type The CoDiPack type of the tape that is to be written out.
   */
  template<typename T_Type>
  struct PrimalGraphTapeWriter : public CommonTextTapeWriter<T_Type> {
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeWriterInterface.
      using Base = CommonTextTapeWriter<T_Type>;                  ///< See CommonTextTapeWriter.

      using Tape = typename Type::Tape;              ///< See TapeWriterInterface.
      using Identifier = typename Type::Identifier;  ///< See TapeWriterInterface.
      using Real = typename Type::Real;              ///< See TapeWriterInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See TapeWriterInterface.

      /// Constructor.
      PrimalGraphTapeWriter(std::string const& name, std::vector<Identifier> const& in,
                            std::vector<Identifier> const& out)
          : Base(true, name, in, out) {};

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
        CODI_UNUSED(primalValue, curPassiveValuePos, passiveValues, curConstantPos, constantValues, stmtEvalHandle);

        std::string node;
        if (nPassiveValues == Config::StatementInputTag) CODI_Unlikely {
          // Do nothing.
        } else CODI_Likely {
          //  Ensure that the all the rhsIdentifiers have been added before connecting edges to them. The mathRep string
          //  is also modified to include the identifier and value.
          Base::placeUnusedRhsNodes(&rhsIdentifiers[curRhsIdentifiersPos], info.numberOfActiveArguments);

          // The mathRep string is modified to include the identifier and value.
          std::string mathRep =
              Base::modifyMathRep(info.mathRepresentation, curLhsIdentifier, &rhsIdentifiers[curRhsIdentifiersPos],
                                  info.numberOfActiveArguments);
          // Add the curLhsIdentifier node.
          Base::createNode(curLhsIdentifier, mathRep);

          // Loop through rhsIdentifiers and create the edges. The identifierExtensions is used to find the current
          // extension for each of the rhsIdentifiers.
          for (size_t argCount = 0; argCount < info.numberOfActiveArguments; argCount++) {
            Base::createEdge(rhsIdentifiers[curRhsIdentifiersPos + argCount], curLhsIdentifier);
          }

          // Only record the increased extension after recording the edges. This ensure that if the lhs equals the rhs,
          // that it results in tow unique nodes.
          this->identifierExtensions[curLhsIdentifier] += 1;
        }
      }
  };
}
