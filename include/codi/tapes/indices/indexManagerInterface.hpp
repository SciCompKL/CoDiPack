#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../aux/optionalArg.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Index>
  struct IndexManagerInterface {
    public:

      using Index = DECLARE_DEFAULT(_Index, int);

      static bool const AssignNeedsStatement = UNDEFINED_VALUE;
      static bool const IsLinear = UNDEFINED_VALUE;

      void assignIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value);
      void assignUnusedIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value);
      void copyIndex(Index& lhs, Index const& rhs);
      void freeIndex(Index& index);

      Index getLargestAssignedIndexIndex() const;

      void reset();
  };
}
