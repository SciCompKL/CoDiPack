/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../traits/tapeTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base class for manual statement pushes on the tape.
   *
   * See StatementPushHelper for details.
   *
   * @tparam T_Type  The CoDiPack type on whose tape the statements are pushed.
   * @tparam T_Impl  The actual implementation.
   */
  template<typename T_Type, typename T_Impl>
  struct StatementPushHelperBase {
    public:

      /// See StatementPushHelperBase.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      /// See StatementPushHelperBase.
      using Impl = CODI_DD(T_Impl, CODI_T(StatementPushHelperBase<CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;  ///< See LhsExpressionInterface.

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Finish the push of a statement. Performs lhs = primal and cleans up all data.
      void endPushStatement(Type& lhs, Real const& primal);

      /// Add the Jacobian of an argument of the statement.
      void pushArgument(Type const& arg, Real const& jacobian);

      /// Initialize all data for the push of a statement.
      void startPushStatement();

      /// @}
      /*******************************************************************************/
      /// @name Start of general implementation
      /// @{

      /// Push a complete statement where the Jacobians and arguments are provided as iterator objects.
      template<typename ArgIter, typename JacobiIter>
      CODI_INLINE void pushStatement(Type& lhs, Real const& primal, ArgIter const startArg, ArgIter const endArg,
                                     JacobiIter const startJac) {
        cast().startPushStatement();

        JacobiIter jacPos = startJac;
        ArgIter argPos = startArg;
        while (argPos != endArg) {
          cast().pushArgument(*argPos, *jacPos);

          ++jacPos;
          ++argPos;
        }

        cast().endPushStatement(lhs, primal);
      }

      /// Push a complete statement where the Jacobians and arguments are provided as arrays.
      template<typename ArgVector, typename JacobiVector>
      CODI_INLINE void pushStatement(Type& lhs, Real const& primal, ArgVector const& arguments,
                                     JacobiVector const& jacobians, size_t const size) {
        cast().startPushStatement();

        for (size_t i = 0; i < size; ++i) {
          cast().pushArgument(arguments[i], jacobians[i]);
        }

        cast().endPushStatement(lhs, primal);
      }

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };

  /**
   * @brief Add statements to the tape where the Jacobians are computed manually.
   *
   * This helper class can be used to optimize the storage of a statement or to handle involved functions that cannot
   * be handled with CoDiPack.
   *
   * The three basic use cases are:
   * \snippet examples/statementPushHelper.cpp Statement push helper
   *
   * After a statement is pushed, the helper can be used again for the next statement.
   *
   * @tparam T_Type  The CoDiPack type on whose tape the statements are pushed.
   */
  template<typename T_Type, typename = void>
  struct StatementPushHelper : public StatementPushHelperBase<T_Type, StatementPushHelper<T_Type>> {
    public:

      ///< See StatementPushHelper.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.
      /// See LhsExpressionInterface.
      using Tape = typename Type::Tape;

    protected:
      Identifier indexData[Config::MaxArgumentSize];  ///< Storage for the identifiers of the arguments.
      Real jacobianData[Config::MaxArgumentSize];     ///< Storage for the Jacobians of the arguments.
      size_t dataPos;                                 ///< Current number of arguments.

    public:

      /*******************************************************************************/
      /// @name Implementation of StatementPushHelperBase
      /// @{

      /// \copydoc codi::StatementPushHelperBase::startPushStatement()
      void startPushStatement() {
        dataPos = 0;
      }

      /// \copydoc codi::StatementPushHelperBase::pushArgument()
      void pushArgument(Type const& arg, Real const& jacobian) {
        Tape& tape = Type::getTape();

        if (Config::MaxArgumentSize <= dataPos) {
          CODI_EXCEPTION("Adding more than %zu arguments to a statement.", Config::MaxArgumentSize);
        }

        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, tape.isActive())) {
          if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, 0 != arg.getIdentifier())) {
            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) {
              if (CODI_ENABLE_CHECK(Config::CheckJacobianIsZero, !RealTraits::isTotalZero(jacobian))) {
                indexData[dataPos] = arg.getIdentifier();
                jacobianData[dataPos] = jacobian;
                dataPos += 1;
              }
            }
          }
        }
      }

      /// \copydoc codi::StatementPushHelperBase::endPushStatement()
      void endPushStatement(Type& lhs, Real const& primal) {
        Tape& tape = Type::getTape();

        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, tape.isActive())) {
          if (0 != dataPos) {
            tape.storeManual(primal, lhs.getIdentifier(), dataPos);

            for (size_t i = 0; i < dataPos; ++i) {
              tape.pushJacobianManual(jacobianData[i], 0.0, indexData[i]);
            }
          }
        }

        lhs.value() = primal;
      }

      /// @}
  };

#ifndef DOXYGEN_DISABLE

  /// Specialization for forward tapes.
  template<typename T_Type>
  struct StatementPushHelper<T_Type, TapeTraits::EnableIfForwardTape<typename T_Type::Tape>>
      : public StatementPushHelperBase<T_Type, StatementPushHelper<T_Type>> {
    public:

      /// See StatementPushHelper.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      using Real = typename Type::Real;          ///< See LhsExpressionInterface.
      using Gradient = typename Type::Gradient;  ///< See LhsExpressionInterface.

    protected:

      Gradient lhsTangent;  ///< Tangent value for the left hand side.

    public:

      /*******************************************************************************/
      /// @name Implementation of StatementPushHelperBase
      /// @{

      /// \copydoc codi::StatementPushHelperBase::startPushStatement()
      void startPushStatement() {
        lhsTangent = Gradient();
      }

      /// \copydoc codi::StatementPushHelperBase::pushArgument()
      void pushArgument(Type const& arg, Real const& jacobian) {
        if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) {
          lhsTangent += jacobian * arg.getGradient();
        }
      }

      /// \copydoc codi::StatementPushHelperBase::endPushStatement()
      void endPushStatement(Type& lhs, Real const& primal) {
        lhs.gradient() = lhsTangent;
        lhs.value() = primal;
      }

      /// @}
  };

  /// Specialization for double.
  template<>
  struct StatementPushHelper<double, void> {
    public:

      using Type = double;  ///< See LhsExpressionInterface.
      using Real = double;  ///< See LhsExpressionInterface.

      /*******************************************************************************/
      /// @name Implementation of StatementPushHelperBase and overwrite of all methods.
      /// @{

      /// \copydoc codi::StatementPushHelperBase::startPushStatement()
      void startPushStatement() {}

      /// \copydoc codi::StatementPushHelperBase::pushArgument()
      void pushArgument(Type const& arg, Real const& jacobian) {
        CODI_UNUSED(arg, jacobian);
      }

      /// \copydoc codi::StatementPushHelperBase::endPushStatement()
      void endPushStatement(Type& lhs, Real const& primal) {
        lhs = primal;
      }

      /// \copydoc codi::StatementPushHelperBase::pushStatement(Type&, Real const&, ArgIter const, ArgIter const,
      /// JacobiIter const)
      template<typename ArgIter, typename JacobiIter>
      void pushStatement(Type& lhs, Real const& primal, ArgIter const startArg, ArgIter const endArg,
                         JacobiIter const startJac) {
        CODI_UNUSED(startArg, endArg, startJac);

        endPushStatement(lhs, primal);
      }

      /// \copydoc codi::StatementPushHelperBase::pushStatement(Type&, Real const&,ArgVector const&, JacobiVector
      /// const&, size_t const)
      template<typename ArgVector, typename JacobiVector>
      void pushStatement(Type& lhs, Real const& primal, ArgVector const& arguments, JacobiVector const& jacobians,
                         size_t const size) {
        CODI_UNUSED(arguments, jacobians, size);

        endPushStatement(lhs, primal);
      }

      /// @}
  };
#endif
}
