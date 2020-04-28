#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../data/position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Identifier>
  struct IdentifierInformationTapeInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      static bool constexpr LinearIndexHandling = CODI_UNDEFINED_VALUE;

      Identifier getPassiveIndex() const;
      Identifier getInvalidIndex() const;
      bool isIdentifierActive(Identifier const& index) const;

      template<typename Lhs>
      void deactivateValue(LhsExpressionInterface<Real, Gradient, IdentifierInformationTapeInterface, Lhs>& value);
  };
}
