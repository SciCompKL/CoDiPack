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

  template<typename _Type>
  struct CustomAdjointVectorInterface {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position;

      Tape& tape;

      CustomAdjointVectorInterface() :
        tape(Type::getGlobalTape()) {
      }

      virtual ~CustomAdjointVectorInterface() {}

      /*******************************************************************************
       * Section: Start of interface definition
       *
       * Description: TODO
       *
       */

      virtual void clearAdjoints() = 0;
      virtual void deleteAdjointVector() = 0;
      virtual void evaluate(Position const& start, Position const& end) = 0;
      virtual void evaluateForward(Position const& start, Position const& end) = 0;
      virtual VectorAccessInterface<Real, Identifier>* getVectorInterface() = 0;

      /*******************************************************************************
       * Section: Common methods
       *
       * Description: TODO
       *
       */

      void evaluate() {
        evaluate(tape.getPosition(), tape.getZeroPosition());
      }


      void evaluateForward() {
        evaluate(tape.getPosition(), tape.getZeroPosition());
      }

      void setTape(Tape& tape) {
        this->tape = tape;
      }

  };

  template<typename _Type, typename _Gradient>
  struct CustomAdjointVectorHelper : public CustomAdjointVectorInterface<_Type> {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double);

      using Base = CustomAdjointVectorInterface<Type>;

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position;

      std::vector<Gradient> adjointVector;

      Gradient zeroValue;
      Gradient const constZeroValue;

      AdjointVectorAccess<Real, Identifier, Gradient>* adjointInterface;

      CustomAdjointVectorHelper() :
        Base(),
        adjointVector(0),
        zeroValue(),
        constZeroValue(),
        adjointInterface(nullptr) {
      }

      /**
       * @brief Destructor
       */
      ~CustomAdjointVectorHelper() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }
      }

      void clearAdjoints() {
        for(size_t i = 0; i < adjointVector.size(); i += 1) {
          adjointVector[i] = Gradient();
        }
      }

      void deleteAdjointVector() {
        adjointVector.resize(0);
        adjointVector.shrink_to_fit();
      }

      void evaluate(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::tape.evaluate(start, end, adjointVector.data());
      }
      using Base::evaluate;

      void evaluateForward(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::tape.evaluateForward(start, end, adjointVector.data());
      }
      using Base::evaluateForward;

      Gradient const& getGradient(Identifier const& identifier) const {
        return gradient(identifier);
      }

      VectorAccessInterface<Real, Identifier>* getVectorInterface() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }

        checkAdjointVectorSize();
        adjointInterface = new AdjointVectorAccess<Real, Identifier, Gradient>(adjointVector.data());
        return adjointInterface;
      }

      Gradient& gradientAt(Identifier const& identifier) {
        return adjointVector[identifier];
      }

      Gradient const& gradientAt(Identifier const& identifier) const {
        return adjointVector[identifier];
      }

      Gradient& gradient(Identifier const& identifier) {
        checkAdjointVectorSize();

        if(0 != identifier && identifier < (Identifier)adjointVector.size()) {
          return adjointVector[identifier];
        } else {
          zeroValue = Gradient();
          return zeroValue;
        }
      }

      Gradient const& gradient(Identifier const& identifier) const {
        if(0 != identifier && identifier < (Identifier)adjointVector.size()) {
          return adjointVector[identifier];
        } else {
          return constZeroValue;
        }
      }

      void setGradient(Identifier& identifier, Gradient const& gradientValue) {
        gradient(identifier) = gradientValue;
      }

    private:

      void checkAdjointVectorSize() {
        if(adjointVector.size() <= Base::tape.getParameter(TapeParameters::LargestIdentifier)) {
          adjointVector.resize(Base::tape.getParameter(TapeParameters::LargestIdentifier) + 1);
        }
      }
  };
}
