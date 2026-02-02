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
struct MultLowLevelFunction {
    using Real = typename Type::Real;
    using Identifier = typename Type::Identifier;
    using AdjointVectorAccess = codi::VectorAccessInterface<Real, Identifier>*;
    using Tape = typename Type::Tape;

    using IterCallback = typename codi::LowLevelFunctionEntry<Tape, Real, Identifier>::IterCallback;

    /// Id for this function.
    static codi::Config::LowLevelFunctionToken ID;

    struct Data {

        static size_t constexpr size_per_entry = sizeof(Real) * 3 + sizeof(Identifier) * 3;

        static size_t size(int n) {
          return sizeof(int) + n * size_per_entry;
        }

        int n;

        Real* a_v;
        Real* b_v;
        Real* w_v_old;

        Identifier* a_i;
        Identifier* b_i;
        Identifier* w_i;


        void read(codi::ByteDataView& dataStore) {
          n = dataStore.read<int>();
          a_v = dataStore.read<Real>(n);
          b_v = dataStore.read<Real>(n);
          w_v_old = dataStore.read<Real>(n);
          a_i = dataStore.read<Identifier>(n);
          b_i = dataStore.read<Identifier>(n);
          w_i = dataStore.read<Identifier>(n);
        }
    };

    /// Store on tape.
    CODI_INLINE static void evalAndStore(Type const* a, Type const* b, Type* w, int n) {
      std::vector<Real> oldPrimals = {};
      Tape& tape = Type::getTape();

      registerOnTape();

      // Perform passive primal operation.
      for(int i = 0; i < n; i += 1) {
        w[i].value() = a[i].getValue() * b[i].getValue();
      }
      // Register on tape
      for(int i = 0; i < n; i += 1) {
        oldPrimals.push_back(tape.registerExternalFunctionOutput(w[i]));
      }


      // Store the data.
      size_t dataSize = Data::size(n);
      codi::ByteDataView dataStore = {};
      tape.pushLowLevelFunction(ID, dataSize, dataStore);

      dataStore.write(n);
      for(int i = 0; i < n; i += 1) {
        dataStore.write(a[i].getValue());
      }
      for(int i = 0; i < n; i += 1) {
        dataStore.write(b[i].getValue());
      }
      for(int i = 0; i < n; i += 1) {
        dataStore.write(oldPrimals[i]);
      }

      for(int i = 0; i < n; i += 1) {
        dataStore.write(a[i].getIdentifier());
      }
      for(int i = 0; i < n; i += 1) {
        dataStore.write(b[i].getIdentifier());
      }
      for(int i = 0; i < n; i += 1) {
        dataStore.write(w[i].getIdentifier());
      }
    }

    CODI_INLINE static void evalAndStore(Type const& a, Type const& b, Type& w) {
      evalAndStore(&a, &b, &w, 1);
    }

    /// Function for forward interpretation.
    CODI_INLINE static void forward(Tape* tape, codi::ByteDataView& dataStore, AdjointVectorAccess adjoints) {
      Data data = {};

      data.read(dataStore);

      if (Tape::HasPrimalValues) {
        for(int i = 0; i < data.n; i += 1) {
          data.a_v[i] = adjoints->getPrimal(data.a_i[i]);
          data.b_v[i] = adjoints->getPrimal(data.b_i[i]);

          Real w = data.a_v[i] * data.b_v[i];

          data.w_v_old[i] = adjoints->getPrimal(data.w_i[i]);
          adjoints->setPrimal(data.w_i[i], w);
        }
      }

      size_t vecDim = adjoints->getVectorSize();
      for(int i = 0; i < data.n; i += 1) {
        for (size_t curDim = 0; curDim < vecDim; curDim += 1) {
          Real w_d = data.b_v[i] * adjoints->getAdjoint(data.a_i[i], curDim)
                     + data.a_v[i] * adjoints->getAdjoint(data.b_i[i], curDim);

          adjoints->resetAdjoint(data.w_i[i], curDim);
          adjoints->updateAdjoint(data.w_i[i], curDim, w_d);
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
          Real w_b = adjoints->getAdjoint(data.w_i[i], curDim);
          adjoints->resetAdjoint(data.w_i[i], curDim);

          adjoints->updateAdjoint(data.a_i[i], curDim, data.b_v[i] * w_b);
          adjoints->updateAdjoint(data.b_i[i], curDim, data.a_v[i] * w_b);
        }
      }
    }

    /// Function for iteration over inputs.
    CODI_INLINE static void iterateInputs(Tape* tape, codi::ByteDataView& dataStore, IterCallback func, void* userData) {
      Data data = {};

      data.read(dataStore);

      for(int i = 0; i < data.n; i += 1) {
        func(&data.a_i[i], userData);
      }
      for(int i = 0; i < data.n; i += 1) {
        func(&data.b_i[i], userData);
      }
    }

    /// Function for iteration over inputs.
    CODI_INLINE static void iterateOutputs(Tape* tape, codi::ByteDataView& dataStore, IterCallback func, void* userData) {
      Data data = {};

      data.read(dataStore);

      for(int i = 0; i < data.n; i += 1) {
        func(&data.w_i[i], userData);
      }
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
