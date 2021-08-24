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
   * @tparam _Type  The underlying CoDiPack type.
   */
  template<typename _Type>
  struct CustomAdjointVectorInterface {
    public:

      /// See CustomAdjointVectorInterface
      using Type = CODI_DD(_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;              ///< See LhsExpressionInterface
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface

      /// See LhsExpressionInterface
      using Tape = CODI_DD(typename Type::Tape, CODI_T(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position;  ///< See PositionalEvaluationTapeInterface

    protected:

      Tape& tape;  ///< Current tape for evaluations. Default: globalTape

    public:

      /// Constructor
      CustomAdjointVectorInterface() : tape(Type::getGlobalTape()) {}

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
   * @tparam _Type  The underlying CoDiPack type.
   * @tparam _Gradient  Type of the entries of the custom adjoint vector.
   */
  template<typename _Type, typename _Gradient>
  struct CustomAdjointVectorHelper : public CustomAdjointVectorInterface<_Type> {
    public:

      ///< See CustomAdjointVectorHelper.
      using Type = CODI_DD(_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Gradient = CODI_DD(_Gradient, double);  ///< See CustomAdjointVectorHelper.

      using Base = CustomAdjointVectorInterface<Type>;  ///< Abbreviation for the base class.

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.

      /// See LhsExpressionInterface
      using Tape = CODI_DD(typename Type::Tape, CODI_T(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position;  ///< See PositionalEvaluationTapeInterface.

    protected:
      std::vector<Gradient> adjointVector;  ///< Custom adjoint vector.

      Gradient zeroValue;             ///< Temporary zero value.
      Gradient const constZeroValue;  ///< Temporary constant zero value.

      AdjointVectorAccess<Real, Identifier, Gradient>* adjointInterface;  ///< Last created adjoint interface.

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
        if (nullptr != adjointInterface) {
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
        if (adjointVector.size() <= Base::tape.getParameter(TapeParameters::LargestIdentifier)) {
          adjointVector.resize(Base::tape.getParameter(TapeParameters::LargestIdentifier) + 1);
        }
      }
  };
}
