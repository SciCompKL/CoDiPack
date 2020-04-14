#pragma once

#include <array>


#include "../../aux/macros.h"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Identifier>
  struct JacobianSorter {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);
      using ArgumentSize = Config::ArgumentSize;

      std::array<Identifier, Config::MaxArgumentSize> indices;
      std::array<Real, Config::MaxArgumentSize> jacobies;
      ArgumentSize size;

      JacobianSorter() = default;

      CODI_INLINE void pushData(Real const& jacobi, Identifier const& index) {
        bool found = false;
        ArgumentSize pos;
        for(pos = 0; pos < size; pos += 1) {
          if(indices[pos] == index) {
            found = true;
            break;
          }
        }

        if(!found) {
          size += 1;
          indices[pos] = index;
          jacobies[pos] = jacobi;
        } else {
          jacobies[pos] += jacobi;
        }

      }

      template<typename Vec>
      CODI_INLINE void storeData(Vec& vec) {
        for(ArgumentSize pos = 0; pos < size; pos += 1) {
          vec.pushData(jacobies[pos], indices[pos]);
        }

        // Reset the data for the next statement
        size = 0;
      }
  };
}
