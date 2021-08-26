#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Internal tape interface that is used by active types to trigger the storing of an expression.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * This interface contains callbacks used by AD variables to access the tape implementation. Each AD variable in the
   * program allocates an identifier and this identifier has to be initialized with a call to initIdentifier(). When an
   * AD variable in the program is destroyed, its identifier has to be freed by the tape by a call to
   * destroyIdentifier() before it is deallocated.
   *
   * The compile time switch AllowJacobianOptimization signals the AD variables that the underlying tape is a Jacobian
   * tape, indicating that certain operations can be hidden from the tape recording process.
   *
   * store() has to be called by the AD variable every time it is assigned. The left hand side value (lhs) has to
   * implement LhsExpressionInterface, the right hand side value (rhs) has to implement ExpressionInterface.
   *
   * ActiveType is the default implementation in CoDiPack which uses this interface and implements the behavior
   * described above.
   *
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Identifier>
  struct InternalStatementRecordingTapeInterface {
    public:

      using Identifier = CODI_DD(T_Identifier, int);  ///< See InternalStatementRecordingTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      static bool constexpr AllowJacobianOptimization =
          CODI_UNDEFINED_VALUE;  ///< If certain operations can be hidden from the tape.

      template<typename Real>
      void initIdentifier(Real& value,
                          Identifier& identifier);  ///< Has to be called for each identifier, after it is allocated.
      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier);  ///< Has to be called for each identifier, before
                                                                    ///< it is deallocated.

      /**
       * @brief Has to be called by an AD variable every time it is assigned.
       *
       * Update of the value is performed by the tape. The tape will additionally store information, e.g., for the
       * reversal of the statement.
       *
       * @tparam Lhs  Has to implement LhsExpressionInterface.
       * @tparam Rhs  Has to implement ExpressionInterface.
       */
      template<typename Lhs, typename Rhs>
      void store(Lhs& lhs, Rhs const& rhs);
  };
}
