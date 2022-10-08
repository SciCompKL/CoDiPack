/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
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
   * This active type implements a static tape.
   *
   * @tparam T_Tape  The tape that manages all expressions created with this type.
   */
  template<typename T_Tape>
  struct ActiveType : public ActiveTypeBase<T_Tape, ActiveType<T_Tape>> {
    public:

      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));  ///< See ActiveType.

      using Base = ActiveTypeBase<T_Tape, ActiveType>;  ///< Base class abbreviation.

      using typename Base::Real;          ///< See ActiveTypeBase.
      using typename Base::PassiveReal;   ///< See ActiveTypeBase.
      using typename Base::Identifier;    ///< See ActiveTypeBase.
      using typename Base::Gradient;      ///< See ActiveTypeBase.

      using typename Base::StoreAs;       ///< See ActiveTypeBase.
      using typename Base::ActiveResult;  ///< See ActiveTypeBase.

    private:

      static Tape tape;

    public:

      /*******************************************************************************/
      /// @name Constructors (all forwarding to the base class)
      /// @{

      /// Constructor
      CODI_INLINE ActiveType() : Base() {}

      /// Constructor
      CODI_INLINE ActiveType(ActiveType<Tape> const& v) : Base(static_cast<Base const&>(v)) {}

      /// Constructor
      CODI_INLINE ActiveType(Real const& value) : Base(value) {}

      /// Constructor
      template<typename U = Real, typename = RealTraits::EnableIfNotPassiveReal<U>>
      CODI_INLINE ActiveType(PassiveReal const& value) : Base(value) {}

      /// Constructor
      template<typename Rhs>
      CODI_INLINE ActiveType(ExpressionInterface<Real, Rhs> const& rhs) : Base(rhs) {}

      /// Constructor
      template<typename Rhs, typename U = Real, typename = RealTraits::EnableIfNotPassiveReal<U>>
      CODI_INLINE ActiveType(ExpressionInterface<typename U::Real, Rhs> const& rhs) : Base(rhs) {}

      /// @}

      /// Destructor
      CODI_INLINE ~ActiveType() {}

      /*******************************************************************************/
      /// @name Assignment operators (all forwarding to the base class)
      /// @{

      /// See ActiveTypeBase::operator=(ActiveTypeBase const&).
      CODI_INLINE ActiveType& operator=(ActiveType const& v) {
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
        return tape;
      }

      /// @}
  };

  template<typename Tape>
  Tape ActiveType<Tape>::tape{};
}
