#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "../traits/realTraits.hpp"
#include "activeType.hpp"
#include "assignmentOperators.hpp"
#include "incrementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Creates a pseudo active type from data references. Can be used to overlay existing data with active types.
   *
   * The class stores references to the value and identifier. The identifier is taken as it is and not initialized or
   * destroyed. The class only wraps the data in a CoDiPack expression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _ActiveType  The type of the active type which is wrapped.
   */
  template<typename _ActiveType>
  struct WritableActiveTypeWrapper
      : public LhsExpressionInterface<typename _ActiveType::Real, typename _ActiveType::Gradient,
                                      typename _ActiveType::Tape, WritableActiveTypeWrapper<_ActiveType>>,
        public AssignmentOperators<typename _ActiveType::Tape, WritableActiveTypeWrapper<_ActiveType>>,
        public IncrementOperators<typename _ActiveType::Tape, WritableActiveTypeWrapper<_ActiveType>> {
    public:

      using ActiveType = CODI_DD(_ActiveType, CODI_T(ActiveType<CODI_ANY>));  ///< See WritableActiveTypeWrapper.
      using Tape = typename ActiveType::Tape;                                 ///< See ActiveType.

      using Real = typename Tape::Real;                   ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.
      using Identifier = typename Tape::Identifier;       ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;           ///< See LhsExpressionInterface.

      using Base =
          LhsExpressionInterface<Real, Gradient, Tape, WritableActiveTypeWrapper>;  ///< Base class abbreviation.

    private:

      Real& primalValue;
      Identifier& identifier;

    public:

      /// The identifier is not initialized. It is assumed to be a valid identifier (either default or assigned by an
      /// expression).
      CODI_INLINE WritableActiveTypeWrapper(Real& value, Identifier& identifier)
          : primalValue(value), identifier(identifier) {
        // deliberately left empty
      }

      /// The identifier is not destroyed. It is assumed to be still valid, since this is only a reference to the
      /// actual value.
      CODI_INLINE ~WritableActiveTypeWrapper() {
        // deliberately left empty
      }

      /// See LhsExpressionInterface::operator =(ExpressionInterface const&)
      CODI_INLINE WritableActiveTypeWrapper<ActiveType>& operator=(WritableActiveTypeWrapper<ActiveType> const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, WritableActiveTypeWrapper>&>(*this) = v;
        return *this;
      }
      using LhsExpressionInterface<Real, Gradient, Tape, WritableActiveTypeWrapper>::operator=;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = WritableActiveTypeWrapper const&;  ///< \copydoc codi::ExpressionInterface::StoreAs

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier() {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::value()
      CODI_INLINE Real& value() {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::value() const
      CODI_INLINE Real const& value() const {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::getGlobalTape()
      static CODI_INLINE Tape& getGlobalTape() {
        return ActiveType::getGlobalTape();
      }

      /// @}
  };
}  // end of namespace codi
