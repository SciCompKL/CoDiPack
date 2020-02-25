#pragma once

#include "../../config.h"
#include "../../aux/macros.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Identifier>
  struct InternalExpressionTapeInterface {
    public:

      using Identifier = DECLARE_DEFAULT(_Identifier, int);

      static bool constexpr AllowJacobianOptimization = UNDEFINED_VALUE;

      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier);
      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier);

      template<typename Lhs, typename Rhs>
      void store(Lhs& lhs, Rhs const& rhs);
  };
}
