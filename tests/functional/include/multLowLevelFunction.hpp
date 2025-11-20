/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */

#include <codi.hpp>

/// w = a * b
template<typename Type>
struct MultLowLevelFunction {
    using Real = typename Type::Real;
    using Identifier = typename Type::Identifier;
    using AdjointVectorAccess = codi::VectorAccessInterface<Real, Identifier>*;
    using Tape = typename Type::Tape;

    using IterCallback = typename codi::LowLevelFunctionEntry<Tape, Real, Identifier>::IterCallback;

    /// Id for this function.
    static codi::Config::LowLevelFunctionToken ID;

    struct Data {

        static size_t constexpr size = sizeof(Real) * 3 + sizeof(Identifier) * 3;

        Real* a_v;
        Real* b_v;
        Real* w_v_old;

        Identifier* a_i;
        Identifier* b_i;
        Identifier* w_i;


        void read(codi::ByteDataView& dataStore) {
          a_v = dataStore.read<Real>(1);
          b_v = dataStore.read<Real>(1);
          w_v_old = dataStore.read<Real>(1);
          a_i = dataStore.read<Identifier>(1);
          b_i = dataStore.read<Identifier>(1);
          w_i = dataStore.read<Identifier>(1);
        }
    };

    /// Store on tape.
    CODI_INLINE static void evalAndStore(Type const& a, Type const& b, Type& w) {
      Tape& tape = Type::getTape();

      registerOnTape();

              // Perform passive primal operation.
      w.value() = a.getValue() * b.getValue();
      // Register on tape
      Real oldPrimal = tape.registerExternalFunctionOutput(w);

              // Store the data.
      size_t dataSize = Data::size;
      codi::ByteDataView dataStore = {};
      tape.pushLowLevelFunction(ID, dataSize, dataStore);

      dataStore.write(a.getValue());
      dataStore.write(b.getValue());
      dataStore.write(oldPrimal);

      dataStore.write(a.getIdentifier());
      dataStore.write(b.getIdentifier());
      dataStore.write(w.getIdentifier());
    }

    /// Function for forward interpretation.
    CODI_INLINE static void forward(Tape* tape, codi::ByteDataView& dataStore, AdjointVectorAccess adjoints) {
      Data data = {};

      data.read(dataStore);

      if (Tape::HasPrimalValues) {
        *data.a_v = adjoints->getPrimal(*data.a_i);
        *data.b_v = adjoints->getPrimal(*data.b_i);

        Real w = *data.a_v * *data.b_v;

        *data.w_v_old = adjoints->getPrimal(*data.w_i);
        adjoints->setPrimal(*data.w_i, w);
      }

      size_t vecDim = adjoints->getVectorSize();
      for (size_t curDim = 0; curDim < vecDim; curDim += 1) {
        Real w_d = *data.b_v * adjoints->getAdjoint(*data.a_i, curDim)
                   + *data.a_v * adjoints->getAdjoint(*data.b_i, curDim);

        adjoints->resetAdjoint(*data.w_i, curDim);
        adjoints->updateAdjoint(*data.w_i, curDim, w_d);
      }
    }

    /// Function for reverse interpretation.
    CODI_INLINE static void reverse(Tape* tape, codi::ByteDataView& dataStore, AdjointVectorAccess adjoints) {
      Data data = {};

      data.read(dataStore);

      if (Tape::HasPrimalValues) {
        adjoints->setPrimal(*data.w_i, *data.w_v_old);
      }

      size_t vecDim = adjoints->getVectorSize();
      for (size_t curDim = 0; curDim < vecDim; curDim += 1) {
        Real w_b = adjoints->getAdjoint(*data.w_i, curDim);
        adjoints->resetAdjoint(*data.w_i, curDim);

        adjoints->updateAdjoint(*data.a_i, curDim, *data.b_v * w_b);
        adjoints->updateAdjoint(*data.b_i, curDim, *data.a_v * w_b);
      }
    }

    /// Function for iteration over inputs.
    CODI_INLINE static void iterateInputs(Tape* tape, codi::ByteDataView& dataStore, IterCallback func, void* userData) {
      Data data = {};

      data.read(dataStore);

      func(data.a_i, userData);
      func(data.b_i, userData);
    }

    /// Function for iteration over inputs.
    CODI_INLINE static void iterateOutputs(Tape* tape, codi::ByteDataView& dataStore, IterCallback func, void* userData) {
      Data data = {};

      data.read(dataStore);

      func(data.w_i, userData);
    }


            /// Register function on tape.
    CODI_INLINE static void registerOnTape() {
      if (codi::Config::LowLevelFunctionTokenInvalid == ID) {
        using Entry = codi::LowLevelFunctionEntry<Tape, typename Type::Real, typename Type::Identifier>;
        ID = Type::getTape().registerLowLevelFunction(Entry(reverse, forward, nullptr, nullptr, iterateInputs, iterateOutputs));
      }
    }
};

template<typename Type>
codi::Config::LowLevelFunctionToken MultLowLevelFunction<Type>::ID = codi::Config::LowLevelFunctionTokenInvalid;
