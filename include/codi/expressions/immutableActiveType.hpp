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
   * @brief Creates a pseudo active type from a data value. Can be used to overlay existing data with immutable active
   * types.
   *
   * The class stores copies to the value and identifier. The identifier is taken as it is and not initialized or
   * destroyed. The class only wraps the data in a CoDiPack expression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _ActiveType  The type of the active type which is wrapped.
   */
  template<typename _ActiveType>
  struct ImmutableActiveType
      : public LhsExpressionInterface<typename _ActiveType::Real, typename _ActiveType::Gradient,
                                      typename _ActiveType::Tape, ImmutableActiveType<_ActiveType>>,
        public AssignmentOperators<typename _ActiveType::Tape, ImmutableActiveType<_ActiveType>>,
        public IncrementOperators<typename _ActiveType::Tape, ImmutableActiveType<_ActiveType>> {
    public:

      using ActiveType = CODI_DD(_ActiveType, CODI_T(ActiveType<CODI_ANY>));  ///< See ImmutableActiveType.
      using Tape = typename ActiveType::Tape;                                 ///< See ActiveType.

      using Real = typename Tape::Real;                   ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.
      using Identifier = typename Tape::Identifier;       ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;           ///< See LhsExpressionInterface.

      using Base = LhsExpressionInterface<Real, Gradient, Tape, ImmutableActiveType>;  ///< Base class abbreviation.

    private:

      Real const primalValue;
      Identifier const identifier;

    public:

      /// The identifier is not initialized. It is assumed to be a valid identifier (either default or assigned by an
      /// expression) and has to be valid throughout the lifespan of this object.
      CODI_INLINE ImmutableActiveType(Real value, Identifier identifier) : primalValue(value), identifier(identifier) {
        // deliberately left empty
      }

      /// Create an immutable copy of an active type. It is assumed that the identifier is valid throughout the lifespan
      /// of this object.
      CODI_INLINE ImmutableActiveType(ActiveType const& value)
          : primalValue(value.getValue()), identifier(value.getIdentifier()) {
        // deliberately left empty
      }

      /// The identifier is not destroyed. It is assumed to be still valid, since this is only an immutable copy of
      /// the actual value.
      CODI_INLINE ~ImmutableActiveType() {
        // deliberately left empty
      }

      /// This class is immutable, delete all assignment operators.
      CODI_INLINE ImmutableActiveType<ActiveType>& operator=(ImmutableActiveType<ActiveType> const& v) = delete;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ImmutableActiveType const&;              ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult = typename ActiveType::ActiveResult;  ///< \copydoc codi::ExpressionInterface::ActiveResult

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      ///
      /// Only the const access functions are implemented.
      ///
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
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
