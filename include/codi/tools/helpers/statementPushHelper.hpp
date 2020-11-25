#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../traits/tapeTraits.hpp"


/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename _Impl>
  struct StatementPushHelperBase {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(StatementPushHelperBase<CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;

      /*******************************************************************************
       * Section: Interface definition
       *
       */

      void endPushStatement(Type& lhs, Real const& primal);
      void pushArgument(Type const& arg, Real const& jacobian);
      void startPushStatement();

      /*******************************************************************************
       * Section: Start of general implementation
       *
       */

      template<typename ArgIter, typename JacobiIter>
      CODI_INLINE void pushStatement(Type& lhs, Real const& primal,
                                     ArgIter const startArg, ArgIter const endArg,
                                     JacobiIter const startJac) {

        cast().startPushStatement();

        JacobiIter jacPos = startJac;
        ArgIter argPos = startArg;
        while(argPos != endArg) {
          cast().pushArgument(*argPos, *jacPos);

          ++jacPos;
          ++argPos;
        }

        cast().endPushStatement(lhs, primal);

      }

      template<typename ArgVector, typename JacobiVector>
      CODI_INLINE void pushStatement(Type& lhs, Real const& primal,
                                     ArgVector const& arguments, JacobiVector const& jacobians,
                                     size_t const size) {

        cast().startPushStatement();

        for(size_t i = 0; i < size; ++i) {
          cast().pushArgument(arguments[i], jacobians[i]);
        }

        cast().endPushStatement(lhs, primal);
      }

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }


  };

  template<typename _Type, typename = void>
  struct StatementPushHelper : public StatementPushHelperBase<_Type, StatementPushHelper<_Type>> {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));

      Identifier indexData[Config::MaxArgumentSize];
      Real jacobianData[Config::MaxArgumentSize];
      size_t dataPos;

      void startPushStatement() {
        dataPos = 0;
      }

      void pushArgument(Type const& arg, Real const& jacobian) {
        Tape& tape = Type::getGlobalTape();

        if(Config::MaxArgumentSize < dataPos ) {
          CODI_EXCEPTION("Adding more than %zu arguments to a statement.", Config::MaxArgumentSize);
        }

        CODI_ENABLE_CHECK (Config::CheckTapeActivity, tape.isActive()) {
          CODI_ENABLE_CHECK(Config::CheckZeroIndex, 0 != arg.getIdentifier()) {
            CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobies, RealTraits::isTotalFinite(jacobian)) {
              CODI_ENABLE_CHECK(Config::CheckJacobiIsZero, !RealTraits::isTotalZero(jacobian)) {

                indexData[dataPos] = arg.getIdentifier();
                jacobianData[dataPos] = jacobian;
                dataPos += 1;
              }
            }
          }
        }
      }

      void endPushStatement(Type& lhs, Real const& primal) {
        Tape& tape = Type::getGlobalTape();

        CODI_ENABLE_CHECK (Config::CheckTapeActivity, tape.isActive()) {
          if(0 != dataPos) {
            tape.storeManual(primal, lhs.getIdentifier(), dataPos);

            for(size_t i = 0; i < dataPos; ++i) {
              tape.pushJacobiManual(jacobianData[i], 0.0, indexData[i]);
            }
          }
        }

        lhs.value() = primal;
      }
  };

  template<typename _Type>
  struct StatementPushHelper<_Type, TapeTraits::enableIfForwardTape<typename _Type::Tape> >
      : public StatementPushHelperBase<_Type, StatementPushHelper<_Type>>
  {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;
      using Gradient = typename Type::Gradient;

      Gradient lhsTangent;

      void startPushStatement() {
        lhsTangent = Gradient();
      }

      void pushArgument(Type const& arg, Real const& jacobian) {
        CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobies, RealTraits::isTotalFinite(jacobian)) {
          lhsTangent += jacobian * arg.getGradient();
        }
      }

      void endPushStatement(Type& lhs, Real const& primal) {
        lhs.gradient() = lhsTangent;
        lhs.value() = primal;
      }
  };

  template<>
  struct StatementPushHelper<double, void> {
    public:

      using Type = double;
      using Real = double;


      void startPushStatement() {}

      void pushArgument(Type const& arg, Real const& jacobian) {
        CODI_UNUSED(arg, jacobian);
      }

      void endPushStatement(Type& lhs, Real const& primal) {
        lhs = primal;
      }

      template<typename ArgIter, typename JacobiIter>
      void pushStatement(Type& lhs, Real const& primal,
                         ArgIter const startArg, ArgIter const endArg,
                         JacobiIter const startJac) {
        CODI_UNUSED(startArg, endArg, startJac);

        endPushStatement(lhs, primal);
      }

      template<typename ArgVector, typename JacobiVector>
      void pushStatement(Type& lhs, Real const& primal,
                         ArgVector const& arguments, JacobiVector const& jacobians,
                         size_t const size) {
        CODI_UNUSED(arguments, jacobians, size);

        endPushStatement(lhs, primal);
      }
  };
}
