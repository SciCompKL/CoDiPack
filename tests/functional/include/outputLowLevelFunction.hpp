/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
struct OutputLowLevelFunction {
    using Real = typename Type::Real;
    using Identifier = typename Type::Identifier;
    using AdjointVectorAccess = codi::VectorAccessInterface<Real, Identifier>*;
    using Tape = typename Type::Tape;

    using IterCallback = typename codi::LowLevelFunctionEntry<Tape, Real, Identifier>::IterCallback;

    /// Id for this function.
    static codi::Config::LowLevelFunctionToken ID;

    struct Data {

        static size_t constexpr size_per_entry = sizeof(Real) * 1 + sizeof(Identifier) * 1;

        static size_t size(int n) {
          return sizeof(int) + n * size_per_entry + sizeof(Real*);
        }

        int n;

        Real* w_v_old;

        Identifier* w_i;

        Real* buffer;

        void read(codi::ByteDataView& dataStore) {
          n = dataStore.read<int>();
          w_v_old = dataStore.read<Real>(n);
          w_i = dataStore.read<Identifier>(n);
          buffer = dataStore.read<Real*>();
        }
    };

    /// Store on tape.
    CODI_INLINE static void evalAndStore(Type* w, int n, Real* buffer) {
      std::vector<Real> oldPrimals = {};
      Tape& tape = Type::getTape();

      registerOnTape();

      // Perform passive primal operation.
      for(int i = 0; i < n; i += 1) {
        w[i].value() = buffer[i];
      }
      // Register on tape
      for(int i = 0; i < n; i += 1) {
        if(buffer[i + n] != 0) {
          oldPrimals.push_back(tape.registerExternalFunctionOutput(w[i]));
        } else {
          w[i] = buffer[i]; // Delete old index.
          oldPrimals.push_back(0.0);
        }
      }


      // Store the data.
      size_t dataSize = Data::size(n);
      codi::ByteDataView dataStore = {};
      tape.pushLowLevelFunction(ID, dataSize, dataStore);

      dataStore.write(n);
      for(int i = 0; i < n; i += 1) {
        dataStore.write(oldPrimals[i]);
      }

      for(int i = 0; i < n; i += 1) {
        dataStore.write(w[i].getIdentifier());
      }

      dataStore.write(buffer);
    }

    CODI_INLINE static void evalAndStore(Type& w, Real* buffer) {
      evalAndStore(&w, 1, buffer);
    }

    /// Function for forward interpretation.
    CODI_INLINE static void forward(Tape* tape, codi::ByteDataView& dataStore, AdjointVectorAccess adjoints) {
      Data data = {};

      data.read(dataStore);

      if (Tape::HasPrimalValues) {
        for(int i = 0; i < data.n; i += 1) {
          data.w_v_old[i] = adjoints->getPrimal(data.w_i[i]);
          adjoints->setPrimal(data.w_i[i], data.buffer[i]);
        }
      }

      size_t vecDim = adjoints->getVectorSize();
      for(int i = 0; i < data.n; i += 1) {
        for (size_t curDim = 0; curDim < vecDim; curDim += 1) {
          adjoints->resetAdjoint(data.w_i[i], curDim);
          adjoints->updateAdjoint(data.w_i[i], curDim, data.buffer[i + (vecDim + 1) * data.n]);
        }
      }
    }

    /// Function for reverse interpretation.
    CODI_INLINE static void reverse(Tape* tape, codi::ByteDataView& dataStore, AdjointVectorAccess adjoints) {
      Data data = {};

      data.read(dataStore);

      if (Tape::HasPrimalValues) {
        for(int i = 0; i < data.n; i += 1) {
          adjoints->setPrimal(data.w_i[i], data.w_v_old[i]);
        }
      }

      size_t vecDim = adjoints->getVectorSize();
      for(int i = 0; i < data.n; i += 1) {
        for (size_t curDim = 0; curDim < vecDim; curDim += 1) {
          data.buffer[i + (curDim + 1) * data.n] = adjoints->getAdjoint(data.w_i[i], curDim);
          adjoints->resetAdjoint(data.w_i[i], curDim);
        }
      }
    }

    /// Function for iteration over inputs.
    CODI_INLINE static void iterateInputs(Tape* tape, codi::ByteDataView& dataStore, IterCallback func, void* userData) {
      Data data = {};

      data.read(dataStore);
    }

    /// Function for iteration over inputs.
    CODI_INLINE static void iterateOutputs(Tape* tape, codi::ByteDataView& dataStore, IterCallback func, void* userData) {
      Data data = {};

      data.read(dataStore);

      for(int i = 0; i < data.n; i += 1) {
        func(&data.w_i[i], userData);
      }
    }

    ///Register function on tape.
    CODI_INLINE static void registerOnTape() {
      if (codi::Config::LowLevelFunctionTokenInvalid == ID) {
        using Entry = codi::LowLevelFunctionEntry<Tape, typename Type::Real, typename Type::Identifier>;
        ID = Type::getTape().registerLowLevelFunction(Entry(reverse, forward, nullptr, nullptr, iterateInputs, iterateOutputs));
      }
    }
};

template<typename Type>
codi::Config::LowLevelFunctionToken OutputLowLevelFunction<Type>::ID = codi::Config::LowLevelFunctionTokenInvalid;
