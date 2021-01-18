#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../traits/tapeTraits.hpp"


/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base class for manual statement pushed on the tape.
   *
   * See StatementPushHelper for details.
   *
   * @tparam _Type  The CoDiPack type on which tape the statements are pushed.
   * @tparam _Impl  The actual implementation
   */
  template<typename _Type, typename _Impl>
  struct StatementPushHelperBase {
    public:

      /// See StatementPushHelperBase
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      /// See StatementPushHelperBase
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(StatementPushHelperBase<CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real; ///< See LhsExpressionInterface

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Finish the push of a statement. Performs lhs = primal and cleans up all data.
      void endPushStatement(Type& lhs, Real const& primal);

      /// Add the Jacobian of an argument for the statement.
      void pushArgument(Type const& arg, Real const& jacobian);

      /// Initialize all data for the push of a statement.
      void startPushStatement();

      /// @}
      /*******************************************************************************/
      /// @name Start of general implementation
      /// @{

      /// Push a complete statement where the Jacobians and arguments are provided as iterator objects.
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

      /// Push a complete statement where the Jacobians and arguments are provided as arrays.
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

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }


  };

  /**
   * @brief Add statements to the tape where the Jacobians are computed manually.
   *
   * This helper class can be used to optimize the storage of a statement or to handle involved functions, that can not
   * be handled with CoDiPack.
   *
   * The three basic use cases are:
   * \snippet examples/statementPushHelper.cpp Statement push helper
   *
   * After a statement is pushed the helper can be used again for the next statement.
   *
   * @tparam _Type  The CoDiPack type on which tape the statements are pushed.
   */
  template<typename _Type, typename = void>
  struct StatementPushHelper : public StatementPushHelperBase<_Type, StatementPushHelper<_Type>> {
    public:

      ///< See StatementPushHelper
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real; ///< See LhsExpressionInterface
      using Identifier = typename Type::Identifier; ///< See LhsExpressionInterface
      /// See LhsExpressionInterface
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));

    protected:
      Identifier indexData[Config::MaxArgumentSize];  ///< Storage for the identifiers of the arguments
      Real jacobianData[Config::MaxArgumentSize];     ///< Storage for the Jacobians of the arguments
      size_t dataPos; ///< Current number of arguments.

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
        Tape& tape = Type::getGlobalTape();

        if(Config::MaxArgumentSize <= dataPos ) {
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

      /// \copydoc codi::StatementPushHelperBase::endPushStatement()
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

      /// @}
  };

#ifndef DOXYGEN_DISABLE

  /// Specialization for forward tapes.
  template<typename _Type>
  struct StatementPushHelper<_Type, TapeTraits::EnableIfForwardTape<typename _Type::Tape> >
      : public StatementPushHelperBase<_Type, StatementPushHelper<_Type>>
  {
    public:

      /// See StatementPushHelper
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real; ///< See LhsExpressionInterface
      using Gradient = typename Type::Gradient; ///< See LhsExpressionInterface

    protected:

      Gradient lhsTangent; ///< Tangent value for the left hand side

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
        CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobies, RealTraits::isTotalFinite(jacobian)) {
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

  /// Specialization for double
  template<>
  struct StatementPushHelper<double, void> {
    public:

      using Type = double; ///< See LhsExpressionInterface
      using Real = double; ///< See LhsExpressionInterface

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

      /// \copydoc codi::StatementPushHelperBase::pushStatement(Type&, Real const&, ArgIter const, ArgIter const, JacobiIter const)
      template<typename ArgIter, typename JacobiIter>
      void pushStatement(Type& lhs, Real const& primal,
                         ArgIter const startArg, ArgIter const endArg,
                         JacobiIter const startJac) {
        CODI_UNUSED(startArg, endArg, startJac);

        endPushStatement(lhs, primal);
      }

      /// \copydoc codi::StatementPushHelperBase::pushStatement(Type&, Real const&,ArgVector const&, JacobiVector const&, size_t const)
      template<typename ArgVector, typename JacobiVector>
      void pushStatement(Type& lhs, Real const& primal,
                         ArgVector const& arguments, JacobiVector const& jacobians,
                         size_t const size) {
        CODI_UNUSED(arguments, jacobians, size);

        endPushStatement(lhs, primal);
      }

      /// @}
  };
#endif
}
