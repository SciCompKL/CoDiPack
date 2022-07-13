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

#include <string>
#include <vector>

#include "../../../config.h"
#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../misc/enumBitset.hpp"
#include "../../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    /// Flags should be one out of each category, that is
    /// {status} + {function} + {kind} + (optional: {version})
    enum class OutputFlags {
      // Category: status
      Intermediate,
      Final,
      // Category: function
      F,
      G,
      P,
      // Category: kind
      Primal,
      Derivative,
      // Category: version (optional)
      V1,
      V2,
      MaxElement
    };
    using OutputHints = EnumBitset<OutputFlags>;

#define ENUM OutputFlags
#include "../../../misc/enumOperations.tpp"

    template<typename T_Type>
    struct ApplicationIOInterface {
        using Type = CODI_DD(T_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
        using Real = typename Type::Real;

        virtual ~ApplicationIOInterface() {}

        virtual void writeRestartY(std::string const& fileName, std::vector<Real> const& v) = 0;
        virtual void writeRestartX(std::string const& fileName, std::vector<Real> const& v) = 0;
        virtual void writeRestartP(std::string const& fileName, std::vector<Real> const& v) = 0;

        virtual void writeRestartData(std::string const& filename, char* data, size_t length) = 0;

        virtual void readRestartY(std::string const& fileName, std::vector<Real>& v) = 0;
        virtual void readRestartX(std::string const& fileName, std::vector<Real>& v) = 0;
        virtual void readRestartP(std::string const& fileName, std::vector<Real>& v) = 0;

        virtual void readRestartData(std::string const& filename, char*& data, size_t& length) = 0;

        virtual void writeY(int iteration, std::vector<Real> const& v, OutputHints flags) = 0;
        virtual void writeX(int iteration, std::vector<Real> const& v, OutputHints flags) = 0;
        virtual void writeP(int iteration, std::vector<Real> const& v, OutputHints flags) = 0;
        virtual void writeZ(int iteration, std::vector<Real> const& v, OutputHints flags) = 0;
    };
  }
}
