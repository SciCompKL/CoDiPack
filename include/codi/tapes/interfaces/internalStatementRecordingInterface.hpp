#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Identifier>
  struct InternalStatementRecordingInterface {
    public:

      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      static bool constexpr AllowJacobianOptimization = CODI_UNDEFINED_VALUE;

      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier);
      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier);

      template<typename Lhs, typename Rhs>
      void store(Lhs& lhs, Rhs const& rhs);
  };
}
