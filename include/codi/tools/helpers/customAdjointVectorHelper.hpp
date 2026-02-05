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

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../tapes/misc/adjointVectorAccess.hpp"
#include "../../tapes/misc/vectorAccessInterface.hpp"
#include "../../traits/tapeTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief General interface for an arbitrary adjoint evaluation.
   *
   * In contrast to CustomAdjointVectorHelper this interface provides a generalization to the former class. This
   * interface can be used if at compile time the required number of vector entries is not known.
   *
   * See CustomAdjointVectorHelper for details.
   *
   * Access to the adjoint vector is gained through the VectorAccessInterface.
   *
   * @tparam T_Type  The underlying CoDiPack type.
   */
  template<typename T_Type>
  struct CustomAdjointVectorInterface {
    public:

      /// See CustomAdjointVectorInterface
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      using Real = typename Type::Real;              ///< See LhsExpressionInterface
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface

      /// See LhsExpressionInterface
      using Tape = typename Type::Tape;
      using Position = typename Tape::Position;  ///< See PositionalEvaluationTapeInterface

    private:

      Tape* tape;  ///< Current tape for evaluations. Default: the Type's current tape.

    public:

      /// Constructor
      CustomAdjointVectorInterface() : tape(&Type::getTape()) {}

      /// Destructor
      virtual ~CustomAdjointVectorInterface() {}

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Set all adjoints to zero
      virtual void clearAdjoints() = 0;

      /// Delete the adjoint vector
      virtual void deleteAdjointVector() = 0;

      /// True if the size for the adjoint vector is fixed. It will not be resized according to the largest identifier
      /// of the tape.
      virtual bool isAdjointSizeFixed() = 0;

      /// Set a fixed size for the adjoint vector. It will no longer be resized to the largest identifier of the tape. A
      /// size of 0 will disable the behavior.
      virtual void resizeAdjointToFixedSize(size_t size) = 0;

      /// Resize the adjoint vector to the largest identifier of the tape. This will remove any fixed size that is set.
      virtual void resizeAdjointToTapeSize() = 0;

      /// \copydoc codi::PositionalEvaluationTapeInterface::evaluate()
      virtual void evaluate(Position const& start, Position const& end) = 0;

      // clang-format off
      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward(T_Position const&, T_Position const&, AdjointsManagement)
      // clang-format on
      virtual void evaluateForward(Position const& start, Position const& end) = 0;

      /// Get a new general interface to the adjoint vector.
      virtual VectorAccessInterface<Real, Identifier>* getVectorInterface() = 0;

      /// @}
      /*******************************************************************************/
      /// @name Common methods
      /// @{

      /// \copydoc codi::ReverseTapeInterface::evaluate()
      void evaluate() {
        evaluate(tape->getPosition(), tape->getZeroPosition());
      }

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      void evaluateForward() {
        evaluate(tape->getPosition(), tape->getZeroPosition());
      }

      /// Set the tape for the evaluations.
      void setTape(Tape& tape) {
        this->tape = &tape;
      }

      /// Get the tape for the evaluations.
      Tape& getTape() {
        return *this->tape;
      }

      /// @}
  };

  /**
   * @brief Allows for an arbitrary adjoint evaluation of a recorded tape.
   *
   * The evaluation of a reverse AD tape is independent of the recording of the tape. For example, the reverse
   * evaluation can be performed simultaneously on differtent adjoint vectors or with a vector mode that is not the same
   * as the vector mode implied by the CoDiPack type.
   *
   * An example for a different vector mode and the general usage is:
   * \snippet documentation/examples/Example_02_Custom_adjoint_vector_evaluation.cpp Custom Adjoint Vector Helper
   *
   * The major difference in using the vector helper is that you cannot interact with the adjoints directly, e.g. via
   * the #codi::LhsExpressionInterface::getGradient or #codi::LhsExpressionInterface::setGradient functions.
   * Instead, you extract the identifier from the AD variable via the #codi::LhsExpressionInterface::getIdentifier()
   * functions. The identifier is then used to access the adjoint values via the vector helper, e.g. for seeding or
   * extraction.
   *
   * For details on custom vector evaluations see CustomAdjointVectorEvaluationTapeInterface.
   *
   * @tparam T_Type  The underlying CoDiPack type.
   * @tparam T_Gradient  Type of the entries of the custom adjoint vector.
   */
  template<typename T_Type, typename T_Gradient>
  struct CustomAdjointVectorHelper : public CustomAdjointVectorInterface<T_Type> {
    public:

      ///< See CustomAdjointVectorHelper.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using Gradient = CODI_DD(T_Gradient, double);  ///< See CustomAdjointVectorHelper.

      using Base = CustomAdjointVectorInterface<Type>;  ///< Abbreviation for the base class.

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.

      /// See LhsExpressionInterface
      using Tape = CODI_DD(typename Type::Tape, CODI_DEFAULT_TAPE);
      using Position = typename Tape::Position;  ///< See PositionalEvaluationTapeInterface.

    protected:
      std::vector<Gradient> adjointVector;  ///< Custom adjoint vector.

      Gradient zeroValue;             ///< Temporary zero value.
      Gradient const constZeroValue;  ///< Temporary constant zero value.

      /// Last created adjoint interface.
      AdjointVectorAccess<Real, Identifier, Gradient*>* adjointInterface;

      /// Fixed size for adjoint.
      ///  - 0:            Used largest identifier for size of adjoint vector.
      ///  - other values: Use this size for the adjoint vector.
      size_t fixedAdjointSize = {};

    public:

      /// Constructor
      CustomAdjointVectorHelper()
          : Base(), adjointVector(0), zeroValue(), constZeroValue(), adjointInterface(nullptr) {}

      /// Destructor
      ~CustomAdjointVectorHelper() {
        if (nullptr != adjointInterface) {
          delete adjointInterface;
        }
      }

      /*******************************************************************************/
      /// @name Implementation of CustomAdjointVectorInterface interface
      /// @{

      /// \copydoc codi::CustomAdjointVectorInterface::clearAdjoints()
      void clearAdjoints() {
        for (size_t i = 0; i < adjointVector.size(); i += 1) {
          adjointVector[i] = Gradient();
        }
      }

      /// \copydoc codi::CustomAdjointVectorInterface::deleteAdjointVector()
      void deleteAdjointVector() {
        adjointVector.resize(0);
        adjointVector.shrink_to_fit();
      }

      /// \copydoc codi::CustomAdjointVectorInterface::isAdjointSizeFixed()
      bool isAdjointSizeFixed() {
        return fixedAdjointSize != 0;
      }

      /// \copydoc codi::CustomAdjointVectorInterface::resizeAdjointToFixedSize()
      void resizeAdjointToFixedSize(size_t size) {
        fixedAdjointSize = size;
        adjointVector.resize(size);
        adjointVector.shrink_to_fit();
      }

      /// \copydoc codi::CustomAdjointVectorInterface::resizeAdjointToTapeSize();
      void resizeAdjointToTapeSize() {
        fixedAdjointSize = 0;
        checkAdjointVectorSize();
      }

      /// \copydoc codi::CustomAdjointVectorInterface::evaluate()
      void evaluate(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::getTape().evaluate(start, end, adjointVector.data());
      }
      using Base::evaluate;

      /// \copydoc codi::CustomAdjointVectorInterface::evaluateForward()
      void evaluateForward(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::getTape().evaluateForward(start, end, adjointVector.data());
      }
      using Base::evaluateForward;

      /// \copydoc codi::CustomAdjointVectorInterface::getVectorInterface()
      VectorAccessInterface<Real, Identifier>* getVectorInterface() {
        if (nullptr != adjointInterface) {
          delete adjointInterface;
        }

        checkAdjointVectorSize();
        adjointInterface = new AdjointVectorAccess<Real, Identifier, Gradient*>(adjointVector.data());
        return adjointInterface;
      }

      /// @}
      /*******************************************************************************/
      /// @name Gradient access methods
      /// @{

      /// Get a constant reference to the gradient.
      Gradient const& getGradient(Identifier const& identifier) const {
        return gradient(identifier);
      }

      /// Get a reference to the gradient. Unchecked access.
      Gradient& gradientUnchecked(Identifier const& identifier) {
        return adjointVector[identifier];
      }

      /// Get a constant reference to the gradient. Unchecked access.
      Gradient const& gradientUnchecked(Identifier const& identifier) const {
        return adjointVector[identifier];
      }

      /// Get a reference to the gradient. Checked access.
      Gradient& gradient(Identifier const& identifier) {
        checkAdjointVectorSize();

        if (0 != identifier && identifier < (Identifier)adjointVector.size()) {
          return adjointVector[identifier];
        } else {
          zeroValue = Gradient();
          return zeroValue;
        }
      }

      /// Get a constant reference to the gradient. Checked access.
      Gradient const& gradient(Identifier const& identifier) const {
        if (0 != identifier && identifier < (Identifier)adjointVector.size()) {
          return adjointVector[identifier];
        } else {
          return constZeroValue;
        }
      }

      /// Set the gradient. Checked access.
      void setGradient(Identifier& identifier, Gradient const& gradientValue) {
        gradient(identifier) = gradientValue;
      }

      /// @}

    private:

      void checkAdjointVectorSize() {
        if (0 == fixedAdjointSize) {  // Skip if fixed size is set
          if (adjointVector.size() <= Base::getTape().getParameter(TapeParameters::LargestIdentifier)) {
            adjointVector.resize(Base::getTape().getParameter(TapeParameters::LargestIdentifier) + 1);
          }
        }
      }
  };
}
