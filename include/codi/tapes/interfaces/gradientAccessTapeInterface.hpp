#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Gradient, typename _Identifier>
  struct GradientAccessTapeInterface {
    public:

      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double);
      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void setGradient(Identifier const& identifier, Gradient const& gradient);
      Gradient const& getGradient(Identifier const& identifier) const;

      Gradient& gradient(Identifier const& identifier);
      Gradient const& gradient(Identifier const& identifier) const;
  };
}
