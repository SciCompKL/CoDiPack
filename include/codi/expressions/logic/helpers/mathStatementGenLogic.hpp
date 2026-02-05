/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <complex>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../../config.h"
#include "../../../misc/macros.hpp"
#include "../../../traits/expressionTraits.hpp"
#include "../traversalLogic.hpp"
#include "forEachLeafLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /** @brief Creates a math representation from a given rhs in the form of a string.
   *
   *  @tparam Identifier The Codi identifier type for internal management, e.g. int.
   */
  template<typename Identifier>
  struct MathStatementGenLogic : public ForEachLeafLogic<MathStatementGenLogic<Identifier>> {
    public:
      Identifier passiveThreshold;  ///< The identifiers that are allocated for passive values.

      /// Constructor.
      MathStatementGenLogic(Identifier passiveThreshold = 0) : passiveThreshold(passiveThreshold) {}

      /// Produces a math representation string for a given statement.
      template<typename Node>
      CODI_INLINE void eval(NodeInterface<Node> const& node, std::string& mathRep) {
        std::vector<std::string> topNodeRep;
        this->toNode(node.cast(), topNodeRep);

        // Check for unnecessary brackets at the start of the statement.
        // If found, remove first and last bracket.
        if (topNodeRep[0].find("(") == 0) {
          topNodeRep[0] = topNodeRep[0].substr(1, topNodeRep[0].size() - 2);
        }

        mathRep = topNodeRep[0];
      }

      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// Links two nodes with a math operation.
      template<typename Node>
      CODI_INLINE void node(Node const& node, std::vector<std::string>& nodeRep) {
        std::vector<std::string> linkRep;
        std::string stmtOperator = node.getMathRep();

        this->toLinks(node, linkRep);

        if (linkRep.size() == 2) {
          // Determine if the operator comes before the two args.
          // This is indicated by the addition of "()" at the end of the stmtOperator.
          if (stmtOperator.find("()") != std::string::npos) {
            // Remove ')'.
            stmtOperator.pop_back();
            nodeRep.push_back(stmtOperator + linkRep[0] + ", " + linkRep[1] + ")");
          } else {
            nodeRep.push_back("(" + linkRep[0] + " " + stmtOperator + " " + linkRep[1] + ")");
          }
        } else {
          nodeRep.push_back(stmtOperator + "(" + linkRep[0] + ")");
        }
      }

      /// Called for leaf nodes which implement LhsExpressionInterface.
      template<typename Node>
      void handleActive(Node const& node, std::vector<std::string>& linkRep) {
        // Check for passive values
        if (node.getIdentifier() <= passiveThreshold) {
          linkRep.push_back("p(" + std::to_string(node.getValue()) + ")");
        } else {
          linkRep.push_back("x" + std::to_string(node.getIdentifier()));
        }
      }

      /// Called for leaf nodes which implement ConstantExpression.
      template<typename Node>
      void handleConstant(Node const& node, std::vector<std::string>& linkRep) {
        linkRep.push_back("c(" + convert_value(node.getValue()) + ")");
      }

      /// Called for leaf nodes which have an EmptyOperation
      template<typename Node>
      void handleEmpty(Node const& node, std::vector<std::string>& linkRep) {
        CODI_UNUSED(node);

        linkRep.push_back("");
      }

    private:

      template<typename T>
      std::string convert_value(T const& v) {
        return std::to_string(v);
      }

      template<typename T>
      std::string convert_value(std::complex<T> const& v) {
        return "(" + std::to_string(std::real(v)) + " + " + std::to_string(std::imag(v)) + ")";
      }

      template<typename T>
      std::string convert_value(T* const& v) {
        return "p" + std::to_string(reinterpret_cast<size_t>(v));
      }
  };
}
