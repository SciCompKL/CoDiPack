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

#include "activeTypeBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Represents a concrete lvalue in the CoDiPack expression tree.
   *
   * See ActiveTypeBase.
   *
   * This active type implements a tape with no state. getTape() constructs a new tape on every call.
   *
   * @tparam T_Tape  The tape that manages all expressions created with this type.
   */
  template<typename T_Tape>
  struct ActiveTypeNoTape : public ActiveTypeBase<T_Tape, ActiveTypeNoTape<T_Tape>> {
    public:

      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_TAPE);  ///< See ActiveTypeNoTape.

      using Base = ActiveTypeBase<T_Tape, ActiveTypeNoTape>;  ///< Base class abbreviation.

      using typename Base::Gradient;     ///< See ActiveTypeBase.
      using typename Base::Identifier;   ///< See ActiveTypeBase.
      using typename Base::PassiveReal;  ///< See ActiveTypeBase.
      using typename Base::Real;         ///< See ActiveTypeBase.

      using typename Base::ActiveResult;  ///< See ActiveTypeBase.
      using typename Base::StoreAs;       ///< See ActiveTypeBase.

    public:

      /// Constructor
      constexpr CODI_INLINE_NO_FA ActiveTypeNoTape() = default;

      /// Constructor
      constexpr CODI_INLINE ActiveTypeNoTape(PassiveReal const& value) : Base(value) {}

      /// Constructor
      CODI_INLINE ActiveTypeNoTape(ActiveTypeNoTape<Tape> const& v) : Base(static_cast<Base const&>(v)) {}

      using Base::Base;  // Use constructors from base class.

      /*******************************************************************************/
      /// @name Assignment operators (all forwarding to the base class)
      /// @{

      /// See ActiveTypeBase::operator=(ActiveTypeBase const&).
      CODI_INLINE ActiveTypeNoTape& operator=(ActiveTypeNoTape const& v) {
        static_cast<Base&>(*this) = static_cast<Base const&>(v);
        return *this;
      }

      using Base::operator=;

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getTape()
      static CODI_INLINE Tape getTape() {
        return Tape();
      }

      /// @}
  };
}
