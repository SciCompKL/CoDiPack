/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */

#include <codi.hpp>

template<typename Real, typename Id, typename Tape>
void func(Tape& tape, std::vector<Real>& x, std::vector<Real>& y, std::vector<Id>& x_id, std::vector<Id>& y_id) {
  Real sum = 0.0;
  Real mul = 1.0;

  for (size_t i = 0; i < x.size(); i += 1) {
    x[i] = (i + 1);
    tape.registerInput(x[i]);
    x_id[i] = x[i].getIdentifier();

    y[i] = sin(x[i]);

    // t1 and t2 will always have the same index but with different meanings.
    {  // Force delete of t1 after scope.
      Real t1 = sum + y[i];
      sum += t1;
    }
    {  // Force delete of t1 after scope.
      Real t2 = mul * y[i];
      mul *= t2;
    }
  }

  for (size_t i = 0; i < x.size(); i += 1) {
    if (0 == (i % 2)) {
      y[i] += sum;
    } else {
      y[i] *= mul;
    }
  }

  std::complex<Real> a(x[0], x[1]);
  std::complex<Real> b(x[2], x[3]);

  std::complex<Real> c = a * b;
  y[0] += norm(c);

  for (size_t i = 0; i < x.size(); i += 1) {
    y_id[i] = y[i].getIdentifier();
    tape.registerOutput(y[i]);
  }
}
