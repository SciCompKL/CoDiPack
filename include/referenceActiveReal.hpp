/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "activeReal.hpp"
#include "expressions.hpp"
#include "typeTraits.hpp"
#include "expressionTraits.hpp"
#include <iostream>

namespace codi {

  template<typename Real, typename Tape>
  class ReferenceActiveReal : public Expression<Real, ReferenceActiveReal<Real, Tape> > {
  public:

    static const bool storeAsReference = true;

    static Tape globalTape;


    typedef Real RealType;

    typedef Tape TapeType;

    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    typedef typename Tape::GradientData GradientData;

  private:
    const ActiveReal<Tape>& reference;
    mutable Real jacobi;

  public:

    inline ReferenceActiveReal(const ActiveReal<Tape>& reference) :
      reference(reference),
      jacobi() {}


    template<typename Data>
    inline void calcGradient(Data& data) const {
      this->jacobi += 1.0;
    }

    template<typename Data>
    inline void calcGradient(Data& data, const Real& jacobi) const {
      this->jacobi += jacobi;
    }

    template<typename Data>
    inline void pushLazyJacobies(Data& data) const {
      if(0.0 != jacobi) {
        reference.calcGradient(data, jacobi);
        jacobi = 0.0; // reset jacobi for the next statement or the next call for this statement
      }
    }

    inline const GradientData& getGradientData() const {
      return reference.getGradientData();
    }

    inline Real getGradient() const {
      return reference.getGradient();
    }

    inline const Real& getValue() const {
      return reference.getValue();
    }

  private:
    inline ReferenceActiveReal<Real, Tape>& operator=(const ReferenceActiveReal& rhs){
    }
  };

  /**
   * @brief Specialization of the TypeTraits for the ActiveReal type.
   *
   * @tparam Real The floating point value of the active real.
   * @tparam Tape The tape of the active real.
   */
  template<typename Real, typename Tape>
  class TypeTraits<ReferenceActiveReal<Real, Tape> > {
    public:
      /**
       * @brief The passive type is the passive type of Real.
       */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

      /**
       * @brief Get the primal value of the origin of this type.
       * @param[in] t The value from which the primal is extracted.
       * @return The primal value of the origin of this type..
       */
      static const typename TypeTraits<Real>::PassiveReal getBaseValue(const ReferenceActiveReal<Real, Tape>& t) {
        return TypeTraits<Real>::getBaseValue(t.getValue());
      }
  };

  template<typename Real, typename Tape>
  struct ExpressionTraits<ReferenceActiveReal<Real, Tape> >  {
    /**
     * @brief The maximum number of active values for an ActiveReal is one.
     */
    static const size_t maxActiveVariables = 1;
  };
}
