#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../aux/optionalArg.hpp"
#include "../../config.h"
#include "../aux/tapeValues.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Index>
  struct IndexManagerInterface {
    public:

      using Index = CODI_DECLARE_DEFAULT(_Index, int);

      static bool constexpr AssignNeedsStatement = CODI_UNDEFINED_VALUE;
      static bool constexpr IsLinear = CODI_UNDEFINED_VALUE;

      static Index constexpr UnusedIndex = 0;
      static Index constexpr InvalidIndex = -1;

      void addToTapeValue(TapeValues& values) const;

      bool assignIndex(Index& index);
      bool assignUnusedIndex(Index& index);
      void copyIndex(Index& lhs, Index const& rhs);
      void freeIndex(Index& index);

      Index getLargestAssignedIndex() const;

      void reset();
  };
}
