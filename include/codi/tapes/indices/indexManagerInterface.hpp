#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../config.h"
#include "../aux/tapeValues.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Index>
  struct IndexManagerInterface {
    public:

      using Index = DECLARE_DEFAULT(_Index, int);

      static bool constexpr AssignNeedsStatement = UNDEFINED_VALUE;
      static bool constexpr IsLinear = UNDEFINED_VALUE;

      static Index constexpr UnusedIndex = 0;
      static Index constexpr InvalidIndex = -1;

      void addToTapeValues(TapeValues& values) const;

      bool assignIndex(Index& index);
      bool assignUnusedIndex(Index& index);
      void copyIndex(Index& lhs, Index const& rhs);
      void freeIndex(Index& index);

      Index getLargestAssignedIndex() const;

      void reset();
  };
}
