#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../tapes/aux/adjointVectorAccess.hpp"
#include "../../tapes/aux/vectorAccessInterface.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
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
   * @tparam _Type  A CoDiPack on which the evaluations take place.
   */
  template<typename _Type>
  struct CustomAdjointVectorInterface {
    public:

      /// See CustomAdjointVectorInterface
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real; ///< See LhsExpressionInterface
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface

      /// See LhsExpressionInterface
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position; ///< See PositionalEvaluationTapeInterface

    protected:

      Tape& tape; ///< Current tape for evaluations. Default: globalTape

      /// Constructor
      CustomAdjointVectorInterface() :
        tape(Type::getGlobalTape()) {
      }

      /// Destructor
      virtual ~CustomAdjointVectorInterface() {}

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Set all adjoints to zero
      virtual void clearAdjoints() = 0;

      /// Delete the adjoint vector
      virtual void deleteAdjointVector() = 0;

      /// \copydoc codi::PositionalEvaluationTapeInterface::evaluate()
      virtual void evaluate(Position const& start, Position const& end) = 0;

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward(Position const&, Position const&)
      virtual void evaluateForward(Position const& start, Position const& end) = 0;

      /// Get a new general interface to the adjoint vector.
      virtual VectorAccessInterface<Real, Identifier>* getVectorInterface() = 0;

      /// @}
      /*******************************************************************************/
      /// @name Common methods
      /// @{

      /// \copydoc codi::ReverseTapeInterface::evaluate()
      void evaluate() {
        evaluate(tape.getPosition(), tape.getZeroPosition());
      }


      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      void evaluateForward() {
        evaluate(tape.getPosition(), tape.getZeroPosition());
      }

      /// Set the tape for the evaluations.
      void setTape(Tape& tape) {
        this->tape = tape;
      }

      /// @}
  };

  /**
   * @brief Allows for an arbitrary adjoint evaluation of a recorded tape.
   *
   * The evaluation of a reverse AD tape is independent of the recording of the tape.
   * That is, the reverse evaluation can be performed simultaneously multiple times or with a different vector mode.
   *
   * An example for a different vector mode and the general usage is:
   * \snippet examples/customAdjointVectorHelper.cpp Custom Adjoint Vector Helper
   *
   * The major difference in using the vector helper is that the adjoints need to be set on the helper instead of the
   * tape. Therefore, the #codi::LhsExpressionInterface::getIdentifier() functions need to be used. It returns the
   * identifier for the tape and this identifier can be used to set the adjoint values in the vector helper.
   *
   * For details on custom vector evaluations see CustomAdjointVectorEvaluationTapeInterface.
   *
   * @tparam _Type  A CoDiPack on which the evaluations take place.
   * @tparam _Gradient  The gradient data for each value.
   */
  template<typename _Type, typename _Gradient>
  struct CustomAdjointVectorHelper : public CustomAdjointVectorInterface<_Type> {
    public:

      ///< See CustomAdjointVectorHelper
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double); ///< See CustomAdjointVectorHelper

      using Base = CustomAdjointVectorInterface<Type>; ///< Abbreviation for the base class

      using Real = typename Type::Real; ///< See LhsExpressionInterface
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface

      /// See LhsExpressionInterface
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position; ///< See PositionalEvaluationTapeInterface

    protected:
      std::vector<Gradient> adjointVector;  ///< Custom adjoint vector

      Gradient zeroValue; ///< Temporary zero value
      Gradient const constZeroValue; ///< Temporary constant zero value.

      AdjointVectorAccess<Real, Identifier, Gradient>* adjointInterface; ///< Last provided adjoint interface.

    public:

      /// Constructor
      CustomAdjointVectorHelper() :
        Base(),
        adjointVector(0),
        zeroValue(),
        constZeroValue(),
        adjointInterface(nullptr) {
      }

      /// Destructor
      ~CustomAdjointVectorHelper() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }
      }

      /*******************************************************************************/
      /// @name Implementation of CustomAdjointVectorInterface interface
      /// @{

      /// \copydoc codi::CustomAdjointVectorInterface::clearAdjoints()
      void clearAdjoints() {
        for(size_t i = 0; i < adjointVector.size(); i += 1) {
          adjointVector[i] = Gradient();
        }
      }

      /// \copydoc codi::CustomAdjointVectorInterface::deleteAdjointVector()
      void deleteAdjointVector() {
        adjointVector.resize(0);
        adjointVector.shrink_to_fit();
      }

      /// \copydoc codi::CustomAdjointVectorInterface::evaluate()
      void evaluate(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::tape.evaluate(start, end, adjointVector.data());
      }
      using Base::evaluate;

      /// \copydoc codi::CustomAdjointVectorInterface::evaluateForward()
      void evaluateForward(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::tape.evaluateForward(start, end, adjointVector.data());
      }
      using Base::evaluateForward;

      /// \copydoc codi::CustomAdjointVectorInterface::getVectorInterface()
      VectorAccessInterface<Real, Identifier>* getVectorInterface() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }

        checkAdjointVectorSize();
        adjointInterface = new AdjointVectorAccess<Real, Identifier, Gradient>(adjointVector.data());
        return adjointInterface;
      }

      /// @}
      /*******************************************************************************/
      /// @name Gradient access methods
      /// @{

      /// Get gradient element as value.
      Gradient const& getGradient(Identifier const& identifier) const {
        return gradient(identifier);
      }


      /// Get reference to the gradient value. Unchecked access.
      Gradient& gradientAt(Identifier const& identifier) {
        return adjointVector[identifier];
      }

      /// Get value to the gradient data. Unchecked access.
      Gradient const& gradientAt(Identifier const& identifier) const {
        return adjointVector[identifier];
      }

      /// Get reference to the gradient value. Checked access.
      Gradient& gradient(Identifier const& identifier) {
        checkAdjointVectorSize();

        if(0 != identifier && identifier < (Identifier)adjointVector.size()) {
          return adjointVector[identifier];
        } else {
          zeroValue = Gradient();
          return zeroValue;
        }
      }

      /// Get value to the gradient data. Checked access.
      Gradient const& gradient(Identifier const& identifier) const {
        if(0 != identifier && identifier < (Identifier)adjointVector.size()) {
          return adjointVector[identifier];
        } else {
          return constZeroValue;
        }
      }

      /// Set the gradient value. Checked access.
      void setGradient(Identifier& identifier, Gradient const& gradientValue) {
        gradient(identifier) = gradientValue;
      }

      /// @}

    private:

      void checkAdjointVectorSize() {
        if(adjointVector.size() <= Base::tape.getParameter(TapeParameters::LargestIdentifier)) {
          adjointVector.resize(Base::tape.getParameter(TapeParameters::LargestIdentifier) + 1);
        }
      }
  };
}
