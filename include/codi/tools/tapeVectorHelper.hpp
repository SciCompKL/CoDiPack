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

#include <vector>

#include "../adjointInterface.hpp"
#include "../adjointInterfaceImpl.hpp"
#include "../configure.h"
#include "../exceptions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Allows for an arbitrary adjoint evaluation of a recorded tape.
   *
   * For the full documentation see TapeVectorHelper.
   *
   * This class can be used in a generalized context. All modifications of the adjoint vector need
   * to be performed via the AdjointInterface obtained via the getAdjointInterface method.
   *
   * This interface needs to be renewed everytime the tape changes or is changed.
   *
   * @tparam CoDiType  A CoDiPack type that which is defined via an ActiveReal.
   */
  template<typename CoDiType>
  struct TapeVectorHelperInterface {

      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientData GradientData; /**< The type for the identification of gradients. */
      typedef typename CoDiType::TapeType Tape; /**< The type of the tape implementation. */
      typedef typename Tape::Position Position; /**< The position for the tape. */

      /**
       * @brief The reference to the tape that is used in the evaluation.
       */
      Tape& tape;

      /**
       * @brief Create a new instance which uses the global tape as the default tape in the background.
       */
      TapeVectorHelperInterface() :
        tape(CoDiType::getGlobalTape()) {
      }

      /**
       * @brief Virtual destructor.
       */
      virtual ~TapeVectorHelperInterface() {}

      /**
       * @brief Set the tape which should be used in the evaluation.
       *
       * @param[in] tape  The tape which is used in the evaluation.
       */
      void setTape(Tape& tape) {
        this->tape = tape;
      }

      /**
       * @brief Delete the adjoint vector.
       */
      virtual void deleteAdjointVector() = 0;

      /**
       * @brief Evaluate the tape from start to end with the adjoint vector of this helper.
       *
       * It has to hold start >= end.
       *
       * @param[in] start  The starting position for the adjoint evaluation.
       * @param[in]   end  The ending position for the adjoint evaluation.
       */
      virtual void evaluate(const Position& start, const Position& end) = 0;

      /**
       * @brief Evaluate the full tape with the adjoint vector of this helper.
       */
      void evaluate() {
        evaluate(tape.getPosition(), tape.getZeroPosition());
      }

      /**
       * @brief Reset all adjoint to there default value.
       */
      virtual void clearAdjoints() = 0;

      /**
       * @brief Obtain a general interface to the adjoint vector in order to modify it.
       *
       * Everytime the tape is modified, the interface needs to be renewed.
       *
       * @return A general interface to the adjoint vector.
       */
      virtual AdjointInterface<Real, GradientData>* getAdjointInterface() = 0;
  };

  /**
   * @brief Allows for an arbitrary adjoint evaluation of a recorded tape.
   *
   * The evaluation of a reverse AD tape is independent of the recording of the tape.
   * That is, the reverse evaluation can be performed simultaneously multiple times or with a different vector mode.
   * Lets assumed that the recorded tape is represented by the function
   * \f[ y = F(x) \eqdot \f]
   * The reverse AD mode evaluates then the directional adjoint derivative in the form
   * \f[ \bar{x} = \frac{d F}{d x}^T(x)\cdot \bar{y}. \f]
   * \f$ \bar{x} \f$ and \f$ \bar{y} \f$ are real valued vectors. It is now possible to extend the AD theory such that
   * multiple directions are evaluated simultaneously which yields
   * \f[ \bar{X} = \frac{d F}{d x}^T(x)\cdot \bar{Y}. \f]
   * \f$ \bar{X} \f$ and \f$ \bar{Y} \f$ are now real valued matrices. E.g. \f$ \bar{Y} \in \R^{n\times d} \f$ where
   * \f$ d \in \N \f$ represents the number of adjoint directions. Since the tape representation is independent of
   * the number of directions, the compile time selection of the appropriate vector with e.g. RealReverseVec is not
   * necessary.
   *
   * The vector helper class allows the user to evaluate a tape with an arbitrary vector mode. In order to record and
   * evaluate the function \f$ f(x) = {a * b, a + b, a - b, a / b} \f$ the following steps are necessary:
   * \code{.cpp}
   *  void func(RealReverse& a, RealReverse& b, RealReverse* y) {
   *    y[0] = a * b;
   *    y[1] = a + b;
   *    y[2] = a - b;
   *    y[3] = a / b;
   *  }
   *
   *  ...
   *  RealReverse::TapeType& tape = RealReverse::globalTape;
   *  tape.setActive();
   *  RealReverse a = 3.0;
   *  RealReverse b = 2.0;
   *  tape.registerInput(a);
   *  tape.registerInput(b);
   *
   *  RealReverse y[4];
   *  func(a, b, y);
   *
   *  for(int i = 0; i < 4; ++i) {
   *    tape.registerOutput(y[i]);
   *  }
   *  tape.setPassive();
   *
   *  // First evaluation
   *  y[0].setGradient(1.0);
   *  tape.evaluate();
   *  assert(a.getGradient() == 2.0);
   *
   *  // Second evaluation
   *  tape.clearAdjoints();
   *  y[1].setGradient(1.0);
   *  tape.evaluate();
   *  assert(a.getGradient() == 1.0);
   *
   *  ...
   *  \endcode
   *
   *  The same evaluation can be achieved with the helper structure:
   *  \code{.cpp}
   *  void func(RealReverse& a, RealReverse& b, RealReverse* y) {
   *    y[0] = a * b;
   *    y[1] = a + b;
   *    y[2] = a - b;
   *    y[3] = a / b;
   *  }
   *
   *  ...
   *  RealReverse::TapeType& tape = RealReverse::globalTape;
   *  tape.setActive();
   *  RealReverse a = 3.0;
   *  RealReverse b = 2.0;
   *  tape.registerInput(a);
   *  tape.registerInput(b);
   *
   *  RealReverse y[4];
   *  func(a, b, y);
   *
   *  for(int i = 0; i < 4; ++i) {
   *    tape.registerOutput(y[i]);
   *  }
   *  tape.setPassive();
   *
   *  TapeVectorHelper<RealReverse, Direction<double, 4>> vh;
   *  for(int i = 0; i < 4; ++i) {
   *    vh.gradient(y[i].getGradientData())[i] = 1.0;
   *  }
   *  vh.evaluate();
   *  assert(vh.gradient(a.getGradientData())[0] == 2.0);
   *  assert(vh.gradient(a.getGradientData())[1] == 1.0);
   *  assert(vh.gradient(a.getGradientData())[2] == 1.0);
   *  assert(vh.gradient(a.getGradientData())[3] == 0.5);
   *  \endcode
   *
   *  The major difference in using the vector helper is that the adjoints need to be set on the helper instead of the
   *  tape. Therefore, the #codi::ActiveReal::getGradientData() functions need to be used. It returns the identifier for the
   *  tape and this identifier can be used to set the adjoint values in the vector helper.
   *
   *  In the default configuration of CoDiPack the TapeVectorHelper works only with Jacobi tapes. For primal value tapes
   *  the preprocessor option CODI_EnableVariableAdjointInterfaceInPrimalTapes needs to be set.
   *
   * @tparam      CoDiType  A CoDiPack type that which is defined via an ActiveReal.
   * @tparam GradientValue  The gradient data for each value.
   */
  template<typename CoDiType, typename GradientValue>
  struct TapeVectorHelper : public TapeVectorHelperInterface<CoDiType> {

      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientData GradientData; /**< The type for the identification of gradients. */
      typedef typename CoDiType::TapeType Tape; /**< The type of the tape implementation. */
      typedef typename Tape::Position Position; /**< The position for the tape. */

      /**
       * @brief Storage for the adjoint values.
       */
      std::vector<GradientValue> adjointVector;

      GradientValue zeroValue; /**< Helper value for out of bounds access */
      const GradientValue constZeroValue; /**< Helper value for out of bounds access */

      AdjointInterfaceImpl<Real, GradientData, GradientValue>* adjointInterface; /**< General access to the adjoint vector for the generalized interface. */

      /**
       * @brief Create a new instance which uses the global tape as the default tape in the background.
       */
      TapeVectorHelper() :
        TapeVectorHelperInterface<CoDiType>(),
        adjointVector(0),
        zeroValue(),
        constZeroValue(),
        adjointInterface(nullptr) {
      }

      /**
       * @brief Destructor
       */
      ~TapeVectorHelper() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }
      }

      /**
       * @brief Delete the adjoint vector.
       */
      void deleteAdjointVector() {
        adjointVector.resize(0);
        adjointVector.shrink_to_fit();
      }

      /**
       * @brief Set the gradient value in the internal adjoint vector.
       *
       * @param[in]         value  The identifier for the corresponding primal value.
       * @param[in] gradientValue  The gradient value which is set in the internal adjoint vector.
       */
      void setGradient(GradientData& value, const GradientValue& gradientValue) {
        gradient(value) = gradientValue;
      }

      /**
       * @brief Get the gradient value from the internal adjoint vector.
       *
       * @param[in] value  The identifier for the corresponding primal value.
       *
       * @return The gradient value from the internal adjoint vector.
       */
      const GradientValue& getGradient(const GradientData& value) const {
        return gradient(value);
      }

      /**
       * @brief Get the gradient value from the internal adjoint vector.
       *
       * This is the unchecked version of the gradient function.
       *
       * @param[in] value  The identifier for the corresponding primal value.
       *
       * @return The gradient value from the internal adjoint vector.
       */
      GradientValue& gradientAt(const GradientData& value) {
        return adjointVector[arrayAccess(value)];
      }

      /**
       * @brief Get the gradient value from the internal adjoint vector.
       *
       * This is the unchecked version of the gradient function.
       *
       * @param[in] value  The identifier for the corresponding primal value.
       *
       * @return The gradient value from the internal adjoint vector.
       */
      const GradientValue& gradientAt(const GradientData& value) const {
        return adjointVector[arrayAccess(value)];
      }

      /**
       * @brief Get the reference to the gradient value from the internal adjoint vector.
       *
       * @param[in] value  The identifier for the corresponding primal value.
       *
       * @return The reference to the gradient value from the internal adjoint vector.
       */
      GradientValue& gradient(const GradientData& value) {
        checkAdjointVectorSize();

        if(0 != value && value < (GradientData)adjointVector.size()) {
          return adjointVector[arrayAccess(value)];
        } else {
          zeroValue = GradientValue();
          return zeroValue;
        }
      }

      /**
       * @brief Get the reference to the gradient value from the internal adjoint vector.
       *
       * @param[in] value  The identifier for the corresponding primal value.
       *
       * @return The reference to the gradient value from the internal adjoint vector.
       */
      const GradientValue& gradient(const GradientData& value) const {
        if(0 != value && value < (GradientData)adjointVector.size()) {
          return adjointVector[arrayAccess(value)];
        } else {
          return constZeroValue;
        }
      }

      /**
       * @brief Evaluate the tape from start to end with the adjoint vector of this helper.
       *
       * It has to hold start >= end.
       *
       * @param[in] start  The starting position for the adjoint evaluation.
       * @param[in]   end  The ending position for the adjoint evaluation.
       */
      void evaluate(const Position& start, const Position& end) {
        checkAdjointVectorSize();

        this->tape.evaluate(start, end, adjointVector.data());
      }

      using TapeVectorHelperInterface<CoDiType>::evaluate;

      /**
       * @brief Reset all adjoint to there default value.
       */
      void clearAdjoints() {
        for(size_t i = 0; i < adjointVector.size(); i += 1) {
          adjointVector[arrayAccess(i)] = GradientValue();
        }
      }

      /**
       * @brief Obtain a general interface to the adjoint vector in order to modify it.
       *
       * Everytime the tape is modified, the interface needs to be renewed.
       *
       * @return A general interface to the adjoint vector.
       */
      AdjointInterface<Real, GradientData>* getAdjointInterface() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }

        checkAdjointVectorSize();
        adjointInterface = new AdjointInterfaceImpl<Real, GradientData, GradientValue> (adjointVector.data());
        return adjointInterface;
      }

    private:

      /**
       * @brief Update the vector size of the internal adjoint vector such that it can hold the adjoint values of the
       * tape evaluation.
       */
      void checkAdjointVectorSize() {
        if(adjointVector.size() <= this->tape.getAdjointSize()) {
          adjointVector.resize(this->tape.getAdjointSize() + 1);
        }
      }
  };
}
