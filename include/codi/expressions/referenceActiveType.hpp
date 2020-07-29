#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../traits/realTraits.hpp"
#include "assignmentOperators.hpp"
#include "incerementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type>
  struct ReferenceActiveType : public LhsExpressionInterface<typename _Type::Real, typename _Type::Gradient, typename _Type::Tape, ReferenceActiveType<_Type> >,
                               public AssignmentOperators<_Type, ReferenceActiveType<_Type>>,
                               public IncrementOperators<_Type, ReferenceActiveType<_Type>> {
    public:

      using Type = DECLARE_DEFAULT(_Type, TEMPLATE(LhsExpressionInterface<double, double, ANY, ANY>));
      using Tape = typename Type::Tape;

      using StoreAs = ReferenceActiveType const&;

      using Real = typename Tape::Real;
      using PassiveReal = PassiveRealType<Real>;
      using Identifier = typename Tape::Identifier;
      using Gradient = typename Tape::Gradient;

  private:

      Type& reference;

  public:

      // Additional data structures used by the tape implementation for optimizations.

      mutable Real jacobian;


      /*******************************************************************************
       * Section: Interface implementation
       *
       */

      CODI_INLINE ReferenceActiveType(Type & v) : reference(v), jacobian() {}

      CODI_INLINE ReferenceActiveType<Tape>& operator=(ReferenceActiveType<Tape > const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, ReferenceActiveType>&>(*this) = v;
        return *this;
      }
      using LhsExpressionInterface<Real, Gradient, Tape, ReferenceActiveType>::operator=;

      CODI_INLINE Identifier& getIdentifier() {
        return reference.getIdentifier();
      }

      CODI_INLINE Identifier const& getIdentifier() const {
        return reference.getIdentifier();
      }

      CODI_INLINE Real& value() {
        return reference.value();
      }

      CODI_INLINE Real const& value() const {
        return reference.value();
      }

      static CODI_INLINE Tape& getGlobalTape() {
        return Type::getGlobalTape();
      }
  };
}
