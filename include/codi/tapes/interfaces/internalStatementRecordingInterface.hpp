#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Internal tape interface that is used by active types to trigger the storing of an expression.
   *
   * The interface is the callback interface of values in the program to the tape implementation. Each value in the
   * program needs to allocate an identifier and this identifier needs to be initialized with a call to initIdentifier().
   * When the value in the program is destroyed the identifier freed by the tape before it is destroyed by the program.
   * That is destroyIdentifier() needs to be called before the value is destructed.
   *
   * The compile time switch AllowJacobianOptimization signals the values that the underlying tape is a Jacobian tape
   * that certain operations can be hidden from the tape recording process.
   *
   * store() needs to be called every time by the value when it is assigned. The left hand side value (lhs) needs to
   * implement the LhsExpressionInterface, the right hand side value (rhs) needs to implement the ExpressionInterface.
   *
   * ActiveType is the default implementation in CoDiPack which uses this interface and implements the behavior
   * described above.
   *
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   */
  template<typename _Identifier>
  struct InternalStatementRecordingInterface {
    public:

      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int); ///< See InternalStatementRecordingInterface

      /*******************************************************************************/
      /// @name Interface definition

      static bool constexpr AllowJacobianOptimization = CODI_UNDEFINED_VALUE; ///< If certain operations can be hidden from the tape.

      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier); ///< Needs to be called for each identifier, after it is allocated.
      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier); ///< Needs to be called for each identifier, before it is deallocated.

      /**
       * @brief Needs to called by a value every time it is assigned.
       *
       * Update of the value is performed by the tape. The tape will additionally store information for the e.g.
       * reversal of the statement.
       *
       * @tparam Lhs  Needs to implement LhsExpressionInterface.
       * @tparam Rhs  Needs to implement ExpressionInterface.
       */
      template<typename Lhs, typename Rhs>
      void store(Lhs& lhs, Rhs const& rhs);
  };
}
