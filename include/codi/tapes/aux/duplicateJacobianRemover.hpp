#pragma once

#include <array>

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Combines entries of Jacobians with the same identifier.
   *
   * This class is used in the storing process of the Jacobians for an expression. For each pushData, it checks if a
   * Jacobian with the same identifier has already been pushed. If so, then it combines these Jacobians.
   *
   * This behavior can be enabled with `-DCODI_RemoveDuplicateJacobianArguments=1`. See JacobianBaseTape::pushJacobians
   * for details.
   *
   * @tparam _Real  The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam _Identifier  The adjoint/tangent identifier type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename _Real, typename _Identifier>
  struct DuplicateJacobianRemover {
    public:

      using Real = CODI_DD(_Real, double);           ///< See DuplicateJacobianRemover.
      using Identifier = CODI_DD(_Identifier, int);  ///< See DuplicateJacobianRemover.
      using ArgumentSize = Config::ArgumentSize;     ///< Definition of ArgumentSize type.

    private:
      std::array<Identifier, Config::MaxArgumentSize> indices;
      std::array<Real, Config::MaxArgumentSize> jacobies;
      ArgumentSize size;

    public:

      /// Constructor
      DuplicateJacobianRemover() = default;

      /// For all added items, check if one matches the identifier. If yes combine, if no append.
      CODI_INLINE void pushData(Real const& jacobi, Identifier const& index) {
        bool found = false;
        ArgumentSize pos;
        for (pos = 0; pos < size; pos += 1) {
          if (indices[pos] == index) {
            found = true;
            break;
          }
        }

        if (!found) {
          size += 1;
          indices[pos] = index;
          jacobies[pos] = jacobi;
        } else {
          jacobies[pos] += jacobi;
        }
      }

      /// Add the data to the provided vector. Resets the internal data for a new statement push.
      /// @tparam Vec  DataInterface with Chunk2<Real, Identifier> as data.
      template<typename Vec>
      CODI_INLINE void storeData(Vec& vec) {
        for (ArgumentSize pos = 0; pos < size; pos += 1) {
          vec.pushData(jacobies[pos], indices[pos]);
        }

        // Reset the data for the next statement.
        size = 0;
      }
  };
}
