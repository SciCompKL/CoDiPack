#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../aux/tapeValues.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Indices provide the mapping of primal values to their adjoint counterparts.
   *
   * In operator overloading AD each primal variable (e.g. \f$ w \f$) needs to be mapped to the adjoint counterpart
   * (e.g. \f$ \bar w \f$. Since adjoint can not be stored in the primal an identifier (usually an index) is associated
   * with each primal. This identifier is then used to access the adjoint variable.
   *
   * The interface defines the three basic operations which can be applied to a variable: assign, copy and free.
   * For each of these operations on the primal variable, the corresponding function on the identifier needs to be
   * called.
   *
   * freeIndex() only needs to be called in destructors. If a variable is overwritten only assign needs to be called on
   * the left hand side identifier. The index manager has to decide how the old identifier is handled.
   *
   * assignUnusedIndex() provides identifiers, that have not been used after the last reset. These identifiers can be
   * used for input values of the program because the adjoint will not be overwritten by intermediate variables.
   *
   * CopyNeedsStatement is a static check if the index manager implements a copy optimization, that is if it creates a
   * new identifier for the left hand side or copies the right hand side identifier over.
   *
   * IsLinear tells if the indices are coupled to the statements of a program. The tape needs to be managed accordingly.
   *
   * Mathematical and implementational details are explained in \ref SBG2021Index.
   *
   * @tparam _Index  Type for the identifier usually an integer type.
   */
  template<typename _Index>
  struct IndexManagerInterface {
    public:

      using Index = CODI_DD(_Index, int);  ///< See IndexManagerInterface

      /*******************************************************************************/
      /// @name Global constants

      static Index constexpr UnusedIndex = 0;  ///< Default inactive index for all index managers.
      static Index constexpr InvalidIndex = -1;  ///< Default invalid index for all index mangers. (Max value for unsigned types.)

      /*******************************************************************************/
      /// @name Identifier handling

      static bool constexpr CopyNeedsStatement = CODI_UNDEFINED_VALUE;  ///< If no copy optimization is implemented. See IndexManagerInterface.
      static bool constexpr IsLinear = CODI_UNDEFINED_VALUE;  ///< If identifiers are coupled to the statements. See IndexManagerInterface.

      bool assignIndex(Index& index);  ///< Call on assignment on a primal value e.g. On `w` for  `w = a + b`.
                                      ///< @return true if new indices have been generated internally.
      bool assignUnusedIndex(Index& index);  ///< Call on registering input values.
                                             ///< @return true if new indices have been generated internally.
      void copyIndex(Index& lhs, Index const& rhs);  ///< Call on copy on a primal value e.g. `w = a`.
      void freeIndex(Index& index);  ///< Call on destruction of primal value. Usually in the destructor.

      void reset();  ///< Reset for a new recording

      /*******************************************************************************/
      /// @name Misc functions

      /**
       * @brief Add storage and other information to the tape values.
       * @param[inout] values  Will only create new data entries and no new section.
       */
      void addToTapeValues(TapeValues& values) const;

      Index getLargestAssignedIndex() const;  ///< @return The largest assigned index.
  };
}
