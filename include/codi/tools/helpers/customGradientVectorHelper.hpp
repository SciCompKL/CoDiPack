#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../../tapes/aux/adjointVectorAccess.hpp"
#include "../../tapes/aux/vectorAccessInterface.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"


/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type>
  struct CustomGradientVectorInterface {
    public:

      using Type = DECLARE_DEFAULT(_Type, TEMPLATE(LhsExpressionInterface<double, double, ANY, ANY>));

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Tape = DECLARE_DEFAULT(typename Type::Tape,TEMPLATE(FullTapeInterface<double, double, int, ANY>));
      using Position = typename Tape::Position;

      Tape& tape;

      CustomGradientVectorInterface() :
        tape(Type::getGlobalTape()) {
      }

      virtual ~CustomGradientVectorInterface() {}

      /*******************************************************************************
       * Section: Start of interface defintion
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
  struct CustomGradientVectorHelper : public CustomGradientVectorInterface<_Type> {

      using Type = DECLARE_DEFAULT(_Type, TEMPLATE(LhsExpressionInterface<double, double, ANY, ANY>));
      using Gradient = DECLARE_DEFAULT(_Gradient, double);

      using Base = CustomGradientVectorInterface<Type>;

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Tape = DECLARE_DEFAULT(typename Type::Tape,TEMPLATE(FullTapeInterface<double, double, int, ANY>));
      using Position = typename Tape::Position;

      std::vector<Gradient> gradientVector;

      Gradient zeroValue;
      Gradient const constZeroValue;

      AdjointVectorAccess<Real, Identifier, Gradient>* adjointInterface;

      CustomGradientVectorHelper() :
        Base(),
        gradientVector(0),
        zeroValue(),
        constZeroValue(),
        adjointInterface(nullptr) {
      }

      /**
       * @brief Destructor
       */
      ~CustomGradientVectorHelper() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }
      }

      void clearAdjoints() {
        for(size_t i = 0; i < gradientVector.size(); i += 1) {
          gradientVector[i] = Gradient();
        }
      }

      void deleteAdjointVector() {
        gradientVector.resize(0);
        gradientVector.shrink_to_fit();
      }

      void evaluate(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::tape.evaluate(start, end, gradientVector.data());
      }
      using Base::evaluate;

      void evaluateForward(Position const& start, Position const& end) {
        checkAdjointVectorSize();

        Base::tape.evaluateForward(start, end, gradientVector.data());
      }
      using Base::evaluateForward;

      Gradient const& getGradient(Identifier const& value) const {
        return gradient(value);
      }

      VectorAccessInterface<Real, Identifier>* getVectorInterface() {
        if(nullptr != adjointInterface) {
          delete adjointInterface;
        }

        checkAdjointVectorSize();
        adjointInterface = new AdjointVectorAccess<Real, Identifier, Gradient> (gradientVector.data());
        return adjointInterface;
      }

      Gradient& gradientAt(Identifier const& value) {
        return gradientVector[value];
      }

      Gradient const& gradientAt(Identifier const& value) const {
        return gradientVector[value];
      }

      Gradient& gradient(Identifier const& value) {
        checkAdjointVectorSize();

        if(0 != value && value < (Identifier)gradientVector.size()) {
          return gradientVector[value];
        } else {
          zeroValue = Gradient();
          return zeroValue;
        }
      }

      Gradient const& gradient(Identifier const& value) const {
        if(0 != value && value < (Identifier)gradientVector.size()) {
          return gradientVector[value];
        } else {
          return constZeroValue;
        }
      }

      void setGradient(Identifier& value, Gradient const& gradientValue) {
        gradient(value) = gradientValue;
      }

    private:

      void checkAdjointVectorSize() {
        if(gradientVector.size() <= Base::tape.getOption(ConfigurationOption::LargestIdentifier)) {
          gradientVector.resize(Base::tape.getOption(ConfigurationOption::LargestIdentifier) + 1);
        }
      }
  };
}
