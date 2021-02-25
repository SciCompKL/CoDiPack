#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Allow for a direct access to the gradient information computed by the tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The gradient information is usually accessed via the helper functions of the ActiveType. For example
   * \code{.cpp}
   *   ActiveType<Tape> w = 1.0;
   *
   *   w.gradient() = 100.0;
   *
   *   std::cout << "Gradient of w: " << w.getGradient() << std::endl;
   * \endcode
   *
   * These helper function are shortcuts to the functions provided in this interface. But the functions here can also be
   * used to gain the sensitivity information of a variable that is no longer present. For example
   * (documentation/examples/gradientAccessTapeInterface.cpp):
   * \snippet examples/gradientAccessTapeInterface.cpp Gradient Access
   *
   * Implementation hint: Size of the gradient vector should be checked before access.
   *
   * @tparam _Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   */
  template<typename _Gradient, typename _Identifier>
  struct GradientAccessTapeInterface {
    public:

      using Gradient = CODI_DD(_Gradient, double); ///< See GradientAccessTapeInterface
      using Identifier = CODI_DD(_Identifier, int); ///< See GradientAccessTapeInterface

      /*******************************************************************************/
      /// @name Interface definition

      void setGradient(Identifier const& identifier, Gradient const& gradient); ///< Set the gradient
      Gradient const& getGradient(Identifier const& identifier) const; ///< Get the gradient

      Gradient& gradient(Identifier const& identifier); ///< Reference access to gradient
      Gradient const& gradient(Identifier const& identifier) const; ///< Constant reference access to gradient
  };
}
