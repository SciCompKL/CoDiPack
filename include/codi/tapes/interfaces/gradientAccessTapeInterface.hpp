#pragma once

#include "../../config.h"
#include "../../aux/macros.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Gradient, typename _Identifier>
  struct GradientAccessTapeInterface {
    public:

      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);

      void setGradient(Identifier const& identifier, Gradient const& gradient);
      Gradient const& getGradient(Identifier const& identifier) const;

      Gradient& gradient(Identifier const& identifier);
      Gradient const& gradient(Identifier const& identifier) const;
  };
}
