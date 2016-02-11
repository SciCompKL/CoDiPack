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

namespace codi {
  template<typename Real, size_t dim>
  class AdVector {
    private:
      Real vector[dim];

    public:

      inline AdVector() :
        vector() {}

      inline Real& operator[] (const size_t& i) {
        return vector[i];
      }

      inline const Real& operator[] (const size_t& i) const {
        return vector[i];
      }

      inline AdVector<Real, dim>& operator = (const AdVector<Real, dim>& v) {
        for(size_t i = 0; i < dim; ++i) {
          this->vector[i] = v.vector[i];
        }

        return *this;
      }

      inline AdVector<Real, dim>& operator += (const AdVector<Real, dim>& v) {
        for(size_t i = 0; i < dim; ++i) {
          this->vector[i] += v.vector[i];
        }

        return *this;
      }
  };

  template<typename Real, size_t dim>
  inline AdVector<Real, dim> operator * (const Real& s, const AdVector<Real, dim>& v) {
    AdVector<Real, dim> r;
    for(size_t i = 0; i < dim; ++i) {
      r[i] = s * v[i];
    }

    return r;
  }

  template<typename Real, size_t dim>
  inline AdVector<Real, dim> operator * (const AdVector<Real, dim>& v, const Real& s) {
    return s * v;
  }

  template<typename A, typename Real, size_t dim>
  inline bool operator != (const A& s, const AdVector<Real, dim>& v) {
    for(size_t i = 0; i < dim; ++i) {
      if( s != v[i] ) {
        return true;
      }
    }

    return false;
  }

  template<typename A, typename Real, size_t dim>
  inline bool operator != (const AdVector<Real, dim>& v, const A& s) {
    return s != v;
  }
}
