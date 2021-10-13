/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "../misc/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "../tools/parallel/parallelToolbox.hpp"
#include "../traits/realTraits.hpp"
#include "assignmentOperators.hpp"
#include "incrementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Represents a concrete lvalue in the CoDiPack expression tree.
   *
   * This class provides a parallel alternative to ActiveType.
   *
   * @tparam T_Tape  The tape that manages all expressions created with this type.
   * @tparam T_ParallelToolbox  Toolbox used to parallelize this type.
   */
  template<typename T_Tape, typename T_ParallelToolbox>
  struct ParallelActiveType
      : public LhsExpressionInterface<typename T_Tape::Real, typename T_Tape::Gradient, T_Tape,
                                      ParallelActiveType<T_Tape, T_ParallelToolbox>>,
        public AssignmentOperators<T_Tape, ParallelActiveType<T_Tape, T_ParallelToolbox>>,
        public IncrementOperators<T_Tape, ParallelActiveType<T_Tape, T_ParallelToolbox>> {
    public:

      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));  /// See ActiveType.
      using ParallelToolbox = CODI_DD(T_ParallelToolbox, CODI_T());  ///< See ParallelActiveType.

      using Real = typename Tape::Real;                   ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.
      using Identifier = typename Tape::Identifier;       ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;           ///< See LhsExpressionInterface.

      using Base = LhsExpressionInterface<Real, Gradient, Tape, ParallelActiveType>;  ///< Base class abbreviation.

    private:

      Real primalValue;
      Identifier identifier;

      static typename ParallelToolbox::template StaticThreadLocalPointer<Tape, ParallelActiveType> tape;

    public:

      /// Constructor
      CODI_INLINE ParallelActiveType() : primalValue(), identifier() {
        Base::init();
      }

      /// Constructor
      CODI_INLINE ParallelActiveType(ParallelActiveType const& v) : primalValue(), identifier() {
        Base::init();
        this->getTape().store(*this, v);
      }

      /// Constructor
      CODI_INLINE ParallelActiveType(PassiveReal const& value) : primalValue(value), identifier() {
        Base::init();
      }

      /// Constructor
      template<class Rhs>
      CODI_INLINE ParallelActiveType(ExpressionInterface<Real, Rhs> const& rhs) : primalValue(), identifier() {
        Base::init();
        this->getTape().store(*this, rhs.cast());
      }

      /// Destructor
      CODI_INLINE ~ParallelActiveType() {
        Base::destroy();
      }

      /// See LhsExpressionInterface::operator =(ExpressionInterface const&).
      CODI_INLINE ParallelActiveType& operator=(ParallelActiveType const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, ParallelActiveType>&>(*this) = v;
        return *this;
      }
      using LhsExpressionInterface<Real, Gradient, Tape, ParallelActiveType>::operator=;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ParallelActiveType const&;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult = ParallelActiveType;    ///< \copydoc codi::ExpressionInterface::ActiveResult

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier() {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::value()
      CODI_INLINE Real& value() {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::value() const
      CODI_INLINE Real const& value() const {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::getTape()
      static CODI_INLINE Tape& getTape() {
        return *(tape.get());
      }

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
  typename ParallelToolbox::template StaticThreadLocalPointer<Tape, ParallelActiveType<Tape, ParallelToolbox>>
      ParallelActiveType<Tape, ParallelToolbox>::tape;
}
