/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../misc/macros.hpp"
#include "../tapes/interfaces/editingTapeInterface.hpp"
#include "../tools/parallel/parallelToolbox.hpp"
#include "activeTypeBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Represents a concrete lvalue in the CoDiPack expression tree.
   *
   * This active type implements a static thread-local tape, as suitable for parallel taping.
   *
   * @tparam T_Tape  The tape that manages all expressions created with this type.
   * @tparam T_ParallelToolbox  Toolbox used to parallelize this type.
   */
  template<typename T_Tape, typename T_ParallelToolbox>
  struct ParallelActiveType : public ActiveTypeBase<T_Tape, ParallelActiveType<T_Tape, T_ParallelToolbox>> {
    public:

      /// See ParallelActiveType.
      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_PARALLEL_TAPE);
      /// See ParallelActiveType.
      using ParallelToolbox = CODI_DD(T_ParallelToolbox, CODI_DEFAULT_PARALLEL_TOOLBOX);

      using Base = ActiveTypeBase<Tape, ParallelActiveType>;  ///< Base class abbreviation.

      using typename Base::Gradient;     ///< See ActiveTypeBase.
      using typename Base::Identifier;   ///< See ActiveTypeBase.
      using typename Base::PassiveReal;  ///< See ActiveTypeBase.
      using typename Base::Real;         ///< See ActiveTypeBase.

      using typename Base::ActiveResult;  ///< See ActiveTypeBase.
      using typename Base::StoreAs;       ///< See ActiveTypeBase.

      /// See ParallelToolbox.
      using ThreadLocalTapePointer =
          typename ParallelToolbox::template StaticThreadLocalPointer<Tape, ParallelActiveType>;

    private:

      static ThreadLocalTapePointer tape;

    public:

      /// Constructor
      CODI_INLINE ParallelActiveType(ParallelActiveType const& v) : Base(static_cast<Base const&>(v)) {}

      using Base::Base;  // Use constructors from base class.

      /// Destructor
      CODI_INLINE ~ParallelActiveType() {}

      /*******************************************************************************/
      /// @name Assignment operators (all forwarding to the base class)
      /// @{

      /// See ActiveTypeBase::operator=(ActiveTypeBase const&).
      CODI_INLINE ParallelActiveType& operator=(ParallelActiveType const& v) {
        static_cast<Base&>(*this) = static_cast<Base const&>(v);
        return *this;
      }

      using Base::operator=;

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getTape()
      static CODI_INLINE Tape& getTape() {
        return *(tape.get());
      }

      /// @}
      /*******************************************************************************/
      /// @name Additional functions used for parallel taping.
      /// @{

      /// Get the thread-local tape pointer.
      static CODI_INLINE Tape* getTapePtr() {
        return tape.get();
      }

      /// Set the thread-local tape pointer.
      static CODI_INLINE void setTapePtr(Tape* other) {
        tape.set(other);
      }

      /// @}
  };

  template<typename Tape, typename ParallelToolbox>
  typename ParallelActiveType<Tape, ParallelToolbox>::ThreadLocalTapePointer
      ParallelActiveType<Tape, ParallelToolbox>::tape;

#if CODI_IDE
  /// Helper for IDE code completion.
  using CODI_DEFAULT_PARALLEL_ACTIVE_TYPE =
      ParallelActiveType<CODI_DEFAULT_PARALLEL_TAPE, CODI_DEFAULT_PARALLEL_TOOLBOX>;
#endif
}
