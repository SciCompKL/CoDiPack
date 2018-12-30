/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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


#pragma once

#include "../configure.h"
#include "../exceptions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Helper for the manual pushing of a statement on the CoDiPack tape.
   *
   * This helper class can be used to optimize the storage of a statement or to handle involved functions, that can not
   * be handled with CoDiPack.
   *
   * Lets assume that the the function is represent by
   * \f[ w = \phi(x) \quad, \f]
   * where \f$ w \in \R \f$. The reverse AD model of this function is
   * \f[ \bar{x} = \frac{d \phi}{dx}^T(x)\cdot \bar{w}. \f]
   * where \f$ \frac{d \phi}{dx} \f$ is just a vector. The statement push helper allows the user now to compute
   * \f$ \phi \f$ and \f$ \frac{d \phi}{dx} \f$ manually and store the information on the tape.
   *
   * Lets assume that we want to optimize the statement y = x * x.
   * The derivative of that statement is 2.0*x, so we can do this with the statement helper in three ways.
   *
   * The most basic way in which the helper can be used is:
   * \code{.cpp}
   * StatementPushHelper<CoDiType> sh;
   *
   * sh.startPushStatement();
   * sh.pushArgument(x, 2.0 * x.value());
   * sh.endPushStatement(y, x.value() * x.value());
   * \endcode
   *
   * The other two ways simplify the procedure is arrays are available. The iterator version:
   * \code{.cpp}
   * StatementPushHelper<CoDiType> sh;
   *
   * std::vector<CoDiType> values;
   * std::vector<double> jacobies;
   * values.push_back(x);
   * jacobies.push_back(2.0 * x.value());
   *
   * sh.pushStatement(y, x.value() * x.value(), values.begin(), values.end(), jacobies.begin());
   * \endcode
   *
   * The array version:
   * \code{.cpp}
   * StatementPushHelper<CoDiType> sh;
   *
   * CoDiType values[] = {x};
   * double jacobies[] = {2.0 * x.value()};
   *
   * sh.pushStatement(y, x.value() * x.value(), values, jacobies, 1);
   * \endcode
   *
   * @tparam CoDiType  This needs to be one of the CoDiPack types defined through an ActiveReal
   */
  template<typename CoDiType>
  struct StatementPushHelper {

      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientData GradientData; /**< The type for the identification of gradients. */

      /** The type of the tape implementation. */
      typedef typename CoDiType::TapeType Tape;

      GradientData indexVector[MaxStatementIntSize]; /**< Store the identification data for the inputs */
      Real jacobiVector[MaxStatementIntSize]; /**< Store the Jacobi for each input */
      size_t vectorPos; /**< Current position in the storage vectors */

      /**
       * @brief Resets the internal state such that a new expression can be pushed.
       */
      void startPushStatement() {
        vectorPos = 0;
      }

      /**
       * @brief Adds an argument to the expression.
       *
       * The value is tested if it should be pushed on the tape. Invalid values and values that would not change the
       * result are ignored.
       *
       * @param[in]    arg  The CoDiPack value that represents the argument.
       * @param[in] jacobi  The corresponding Jacobi value for the argument.
       */
      void pushArgument(const CoDiType& arg, const Real& jacobi) {
        Tape& tape = CoDiType::getGlobalTape();

        if(MaxStatementIntSize == vectorPos ) {
          CODI_EXCEPTION("Adding more than %zu arguments to a statement.", MaxStatementIntSize);
        }

        ENABLE_CHECK (OptTapeActivity, tape.isActive()) {
          ENABLE_CHECK(OptCheckZeroIndex, 0 != arg.getGradientData()) {
            ENABLE_CHECK(OptIgnoreInvalidJacobies, codi::isfinite(jacobi)) {
              ENABLE_CHECK(OptJacobiIsZero, !isTotalZero(jacobi)) {

                indexVector[vectorPos] = arg.getGradientData();
                jacobiVector[vectorPos] = jacobi;
                vectorPos += 1;
              }
            }
          }
        }
      }

      /**
       * @brief Finalize the statement and push it on the tape.
       *
       * The left hand side value is updated with the new primal value and marked as active.
       *
       * This method will also perform an activity analysis and my disable the left hand side.
       *
       * @param[in,out] lhs  The CoDiPack value on the left hand side of the statement.
       * @param[in]  primal  The new primal value that is assigned to the left hand side.
       */
      void endPushStatement(CoDiType& lhs, const Real& primal) {
        Tape& tape = CoDiType::getGlobalTape();

        ENABLE_CHECK (OptTapeActivity, tape.isActive()) {
          if(0 != vectorPos) {
            tape.storeManual(primal, lhs.getGradientData(), vectorPos);

            for(size_t i = 0; i < vectorPos; ++i) {
              tape.pushJacobiManual(jacobiVector[i], 0.0, indexVector[i]);
            }
          }
        }

        lhs.value() = primal;
      }

      /**
       * @brief Helper function if the Jacobies and arguments are stored in a list which can be iterated.
       *
       * See startPushStatement, pushArgument and endPushStatement for details.
       *
       * @param[in,out]  lhs  The CoDiPack value on the left hand side of the statement.
       * @param[in]   primal  The new primal value that is assigned to the left hand side.
       * @param[in] startArg  The start iterator for the arguments.
       * @param[in]   endArg  The end iterator for the arguments.
       * @param[in] startJac  The start iterator for the Jacobies.
       *
       * @tparam    ArgIter  Iterator type for the arguments.
       * @tparam JacobiIter  Iterator type for the Jacobies.
       */
      template<typename ArgIter, typename JacobiIter>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgIter startArg, const ArgIter endArg, const JacobiIter startJac) {

        startPushStatement();

        JacobiIter jacPos = startJac;
        ArgIter argPos = startArg;
        while(argPos != endArg) {
          pushArgument(*argPos, *jacPos);

          ++jacPos;
          ++argPos;
        }

        endPushStatement(lhs, primal);

      }

      /**
       * @brief Helper function if the Jacobies and arguments are stored in arrays.
       *
       * See startPushStatement, pushArgument and endPushStatement for details.
       *
       * @param[in,out]   lhs  The CoDiPack value on the left hand side of the statement.
       * @param[in]    primal  The new primal value that is assigned to the left hand side.
       * @param[in] arguments  The array with the arguments.
       * @param[in]  jacobies  The array with the Jacobies.
       * @param[in]      size  The size of the array.
       *
       * @tparam    ArgVector  Array type for the arguments.
       * @tparam JacobiVector  Array type for the Jacobies.
       */
      template<typename ArgVector, typename JacobiVector>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgVector& arguments, const JacobiVector& jacobies, const size_t size) {

        startPushStatement();

        for(size_t i = 0; i < size; ++i) {
          pushArgument(arguments[i], jacobies[i]);
        }

        endPushStatement(lhs, primal);
      }
  };

  /**
   * @brief Helper for the manual pushing of a statement on the CoDiPack tape.
   *
   * For a general explanation see StatementPushHelper.
   *
   * Instead of evaluating the reverse mode AD equation:
   * \f[ \bar{x} = \frac{d \phi}{dx}^T(x)\cdot \bar{w} \f]
   *
   * This helper evaluates the the forward mode AD equation:
   * \f[ \dot{w} = \frac{d \phi}{dx}(x)\cdot \dot{x} \f]
   * @tparam CoDiType
   */
  template<typename CoDiType>
  struct ForwardStatementPushHelper {

      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientValue GradientValue; /**< The type for the values of gradients. */

      /** The gradient value for the left hand side */
      GradientValue lhsTangent;

      /**
       * @brief Resets the internal state such that a new expression can be pushed.
       */
      void startPushStatement() {
        lhsTangent = GradientValue();
      }

      /**
       * @brief Adds an argument to the expression.
       *
       * @param[in]    arg  The CoDiPack value that represents the argument.
       * @param[in] jacobi  The corresponding Jacobi value for the argument.
       */
      void pushArgument(const CoDiType& arg, const Real& jacobi) {
        ENABLE_CHECK(OptIgnoreInvalidJacobies, codi::isfinite(jacobi)) {
          lhsTangent += jacobi * arg.getGradient();
        }
      }

      /**
       * @brief Finalize the statement and update the gradient of the left hand side value.
       *
       * @param[in,out] lhs  The CoDiPack value on the left hand side of the statement.
       * @param[in]  primal  The new primal value that is assigned to the left hand side.
       */
      void endPushStatement(CoDiType& lhs, const Real& primal) {
        lhs.gradient() = lhsTangent;
        lhs.value() = primal;
      }

      /**
       * @brief Helper function if the Jacobies and arguments are stored in a list which can be iterated.
       *
       * See startPushStatement, pushArgument and endPushStatement for details.
       *
       * @param[in,out]  lhs  The CoDiPack value on the left hand side of the statement.
       * @param[in]   primal  The new primal value that is assigned to the left hand side.
       * @param[in] startArg  The start iterator for the arguments.
       * @param[in]   endArg  The end iterator for the arguments.
       * @param[in] startJac  The start iterator for the Jacobies.
       *
       * @tparam    ArgIter  Iterator type for the arguments.
       * @tparam JacobiIter  Iterator type for the Jacobies.
       */
      template<typename ArgIter, typename JacobiIter>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgIter startArg, const ArgIter endArg, const JacobiIter startJac) {

        startPushStatement();

        JacobiIter jacPos = startJac;
        ArgIter argPos = startArg;
        while(argPos != endArg) {
          pushArgument(*argPos, *jacPos);

          ++jacPos;
          ++argPos;
        }

        endPushStatement(lhs, primal);
      }

      /**
       * @brief Helper function if the Jacobies and arguments are stored in arrays.
       *
       * See startPushStatement, pushArgument and endPushStatement for details.
       *
       * @param[in,out]   lhs  The CoDiPack value on the left hand side of the statement.
       * @param[in]    primal  The new primal value that is assigned to the left hand side.
       * @param[in] arguments  The array with the arguments.
       * @param[in]  jacobies  The array with the Jacobies.
       * @param[in]      size  The size of the array.
       *
       * @tparam    ArgVector  Array type for the arguments.
       * @tparam JacobiVector  Array type for the Jacobies.
       */
      template<typename ArgVector, typename JacobiVector>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgVector& arguments, const JacobiVector& jacobies, const size_t size) {

        startPushStatement();

        for(size_t i = 0; i < size; ++i) {
          pushArgument(arguments[i], jacobies[i]);
        }

        endPushStatement(lhs, primal);
      }
  };
}
