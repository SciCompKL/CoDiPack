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
      // Category: hints (optional)
      Vector, // Force vector output
      MaxElement
    };
    using OutputHints = EnumBitset<OutputFlags>;

#define ENUM OutputFlags
#include "../../../misc/enumOperations.tpp"

    enum class OutputType {
      P,
      X,
      Y,
      Z,
      MaxElement
    };

    template<typename T_Type>
    struct ApplicationIOInterface {
        using Type = CODI_DD(T_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
        using Real = RealTraits::Real<Type>;

        using RealVector =  std::vector<Real>;

        virtual ~ApplicationIOInterface() {}

        virtual void writeRestartY(std::string const& fileName, RealVector const& v) = 0;
        virtual void writeRestartX(std::string const& fileName, RealVector const& v) = 0;
        virtual void writeRestartP(std::string const& fileName, RealVector const& v) = 0;

        virtual void writeRestartData(std::string const& filename, char* data, size_t length) = 0;

        virtual void readRestartY(std::string const& fileName, RealVector& v) = 0;
        virtual void readRestartX(std::string const& fileName, RealVector& v) = 0;
        virtual void readRestartP(std::string const& fileName, RealVector& v) = 0;

        virtual void readRestartData(std::string const& filename, char*& data, size_t& length) = 0;

        virtual void writeY(int iteration, RealVector const& v, OutputHints flags, int vec = 0) = 0;
        virtual void writeX(int iteration, RealVector const& v, OutputHints flags, int vec = 0) = 0;
        virtual void writeP(int iteration, RealVector const& v, OutputHints flags, int vec = 0) = 0;
        virtual void writeZ(int iteration, RealVector const& v, OutputHints flags, int vec = 0) = 0;


        // Utility functions

        /// Always relative to the base path. Only used in write*. Empty for reset.
        virtual void changeFolder(std::string const& path) = 0;

        virtual void createFolder(std::string const& path) = 0;

        // Helper functions

        void writeY(int iteration, std::vector<RealVector> const& v, OutputHints flags, int vecOffset = 0) {
          for(size_t i = 0; i < v.size(); i += 1) {
            writeY(iteration, v[i], flags, i + vecOffset);
          }
        }

        void writeX(int iteration, std::vector<RealVector> const& v, OutputHints flags, int vecOffset = 0) {
          for(size_t i = 0; i < v.size(); i += 1) {
            writeX(iteration, v[i], flags, i + vecOffset);
          }
        }

        void writeP(int iteration, std::vector<RealVector> const& v, OutputHints flags, int vecOffset = 0) {
          for(size_t i = 0; i < v.size(); i += 1) {
            writeP(iteration, v[i], flags, i + vecOffset);
          }
        }

        void writeZ(int iteration, std::vector<RealVector> const& v, OutputHints flags, int vecOffset = 0) {
          for(size_t i = 0; i < v.size(); i += 1) {
            writeZ(iteration, v[i], flags, i + vecOffset);
          }
        }

        virtual void write(OutputType type, int iteration, RealVector const& v, OutputHints flags, int vec = 0) {
          switch (type) {
            case OutputType::P: writeP(iteration, v, flags, vec); break;
            case OutputType::X: writeX(iteration, v, flags, vec); break;
            case OutputType::Y: writeY(iteration, v, flags, vec); break;
            case OutputType::Z: writeZ(iteration, v, flags, vec); break;
          default:
            CODI_EXCEPTION("Unimplemented switch case.");
          }
        }

        virtual void write(OutputType type, int iteration, std::vector<RealVector> const& v, OutputHints flags, int vecOffset = 0) {
          switch (type) {
            case OutputType::P: writeP(iteration, v, flags, vecOffset); break;
            case OutputType::X: writeX(iteration, v, flags, vecOffset); break;
            case OutputType::Y: writeY(iteration, v, flags, vecOffset); break;
            case OutputType::Z: writeZ(iteration, v, flags, vecOffset); break;
          default:
            CODI_EXCEPTION("Unimplemented switch case.");
          }
        }
    };
  }
}
