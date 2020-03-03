#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../aux/optionalArg.hpp"
#include "../../config.h"
#include "../aux/tapeValues.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Index>
  struct IndexManagerInterface {
    public:

      using Index = DECLARE_DEFAULT(_Index, int);

      static bool const AssignNeedsStatement = UNDEFINED_VALUE;
      static bool const IsLinear = UNDEFINED_VALUE;

      static Index constexpr UnusedIndex = 0;
      static Index constexpr InvalidIndex = -1;

      void addToTapeValue(TapeValues& values) const;

      void assignIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value);
      void assignUnusedIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value);
      void copyIndex(Index& lhs, Index const& rhs);
      void freeIndex(Index& index);

      Index getLargestAssignedIndex() const;

      void reset();
  };
}
