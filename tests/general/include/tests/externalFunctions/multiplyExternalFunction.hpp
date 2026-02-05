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
#pragma once

#include <codi.hpp>

template<typename T_Number, typename T_Tape, typename = void>
struct MultiplyExternalFunction {
  public:
    using Number = CODI_DECLARE_DEFAULT(T_Number, codi::ActiveType<CODI_ANY>);
    using Tape = CODI_DECLARE_DEFAULT(T_Tape, CODI_ANY);

    static Number create(Number const& x1, Number const& x2, Tape const& tape) {
      codi::CODI_UNUSED(tape);
      return x1 * x2;
    }
};

template<typename T_Number, typename T_Tape>
struct MultiplyExternalFunction<T_Number, T_Tape, codi::TapeTraits::EnableIfReverseTape<T_Tape>> {
  public:
    using Number = CODI_DECLARE_DEFAULT(T_Number, codi::ActiveType<CODI_ANY>);

    using Tape = CODI_DECLARE_DEFAULT(T_Tape,
                                      CODI_TEMPLATE(codi::FullTapeInterface<double, double, int, codi::EmptyPosition>));

    using DataExtraction = codi::RealTraits::DataExtraction<Number>;
    using Real = typename DataExtraction::Real;
    using Identifier = typename DataExtraction::Identifier;

    using VAI = codi::VectorAccessInterface<typename Tape::Real, typename Tape::Identifier>;
    using NumberVAI = typename codi::AggregatedTypeVectorAccessWrapperFactory<Number>::RType;

    static void extFuncReverse(Tape* t, void* d, VAI* vaiReal) {
      codi::CODI_UNUSED(t);

      codi::ExternalFunctionUserData* data = static_cast<codi::ExternalFunctionUserData*>(d);

      NumberVAI* vai = codi::AggregatedTypeVectorAccessWrapperFactory<Number>::create(vaiReal);

      Real x1_v, x2_v;
      Identifier x1_i, x2_i, w_i;
      data->getData(x1_v);
      data->getData(x1_i);
      data->getData(x2_v);
      data->getData(x2_i);
      data->getData(w_i);

      size_t dim = vai->getVectorSize();

      for (size_t i = 0; i < dim; ++i) {
        Real w_b = vai->getAdjoint(w_i, i);
        vai->resetAdjoint(w_i, i);

        vai->updateAdjoint(x1_i, i, codi::ComputationTraits::transpose(x2_v) * w_b);
        vai->updateAdjoint(x2_i, i, codi::ComputationTraits::transpose(x1_v) * w_b);
      }

      codi::AggregatedTypeVectorAccessWrapperFactory<Number>::destroy(vai);
    }

    static void extFuncPrimal(Tape* t, void* d, VAI* vaiReal) {
      codi::CODI_UNUSED(t);

      codi::ExternalFunctionUserData* data = static_cast<codi::ExternalFunctionUserData*>(d);

      NumberVAI* vai = codi::AggregatedTypeVectorAccessWrapperFactory<Number>::create(vaiReal);

      Identifier x1_i, x2_i, w_i;
      Real& x1_v = data->getDataRef<Real>();
      data->getData(x1_i);
      Real& x2_v = data->getDataRef<Real>();
      data->getData(x2_i);
      data->getData(w_i);

      x1_v = vai->getPrimal(x1_i);  // Data is overwritten here
      x2_v = vai->getPrimal(x2_i);  // Data is overwritten here

      Real z = x1_v * x2_v;
      vai->setPrimal(w_i, z);

      codi::AggregatedTypeVectorAccessWrapperFactory<Number>::destroy(vai);
    }

    static void extFuncForward(Tape* t, void* d, VAI* vaiReal) {
      codi::CODI_UNUSED(t);

      codi::ExternalFunctionUserData* data = static_cast<codi::ExternalFunctionUserData*>(d);

      NumberVAI* vai = codi::AggregatedTypeVectorAccessWrapperFactory<Number>::create(vaiReal);

      Identifier x1_i, x2_i, w_i;
      Real& x1_v = data->getDataRef<Real>();
      data->getData(x1_i);
      Real& x2_v = data->getDataRef<Real>();
      data->getData(x2_i);
      data->getData(w_i);

      if (vai->hasPrimals()) {
        x1_v = vai->getPrimal(x1_i);  // Data is overwritten here
        x2_v = vai->getPrimal(x2_i);  // Data is overwritten here
      }

      size_t dim = vai->getVectorSize();

      for (size_t i = 0; i < dim; ++i) {
        Real x1_d = vai->getAdjoint(x1_i, i);
        Real x2_d = vai->getAdjoint(x2_i, i);

        Real w_d = x1_d * x2_v + x1_v * x2_d;
        vai->resetAdjoint(w_i, i);
        vai->updateAdjoint(w_i, i, w_d);
      }

      Real z = x1_v * x2_v;
      vai->setPrimal(w_i, z);

      codi::AggregatedTypeVectorAccessWrapperFactory<Number>::destroy(vai);
    }

    static void delFunc(Tape* tape, void* d) {
      codi::CODI_UNUSED(tape);

      codi::ExternalFunctionUserData* data = static_cast<codi::ExternalFunctionUserData*>(d);
      delete data;
    }

    static void iterInFunc(Tape* t, void* d, typename codi::ExternalFunction<Tape>::IterCallback func, void* userData) {
      codi::CODI_UNUSED(t);

      codi::ExternalFunctionUserData* data = static_cast<codi::ExternalFunctionUserData*>(d);

      Real x1_v, x2_v;
      data->getData(x1_v);
      Identifier& x1_i = data->template getDataRef<Identifier>();
      data->getData(x2_v);
      Identifier& x2_i = data->template getDataRef<Identifier>();
      Identifier& w_i = data->template getDataRef<Identifier>();

      if constexpr (std::is_integral_v<Identifier>) {
        func(&x1_i, userData);
        func(&x2_i, userData);
      } else {
        for (auto& i : x1_i) {
          func(&i, userData);
        }
        for (auto& i : x2_i) {
          func(&i, userData);
        }
      }

      codi::CODI_UNUSED(w_i);
    }

    static void iterOutFunc(Tape* t, void* d, typename codi::ExternalFunction<Tape>::IterCallback func,
                            void* userData) {
      codi::CODI_UNUSED(t);

      codi::ExternalFunctionUserData* data = static_cast<codi::ExternalFunctionUserData*>(d);

      Real x1_v, x2_v;
      data->getData(x1_v);
      Identifier& x1_i = data->template getDataRef<Identifier>();
      data->getData(x2_v);
      Identifier& x2_i = data->template getDataRef<Identifier>();
      Identifier& w_i = data->template getDataRef<Identifier>();

      if constexpr (std::is_integral_v<Identifier>) {
        func(&w_i, userData);
      } else {
        for (auto& i : w_i) {
          func(&i, userData);
        }
      }

      codi::CODI_UNUSED(x1_i, x2_i);
    }

    static Number create(Number const& x1, Number const& x2, Tape& tape) {
      codi::ExternalFunctionUserData* data = new codi::ExternalFunctionUserData;
      Number w;
      tape.setPassive();
      w = x1 * x2;
      tape.setActive();

      codi::RealTraits::registerExternalFunctionOutput(w);
      data->addData(codi::RealTraits::getValue(x1));
      data->addData(codi::RealTraits::getIdentifier(x1));
      data->addData(codi::RealTraits::getValue(x2));
      data->addData(codi::RealTraits::getIdentifier(x2));
      data->addData(codi::RealTraits::getIdentifier(w));
      tape.pushExternalFunction(codi::ExternalFunction<Tape>(extFuncReverse, extFuncForward, extFuncPrimal, data,
                                                             delFunc, iterInFunc, iterOutFunc));

      return w;
    }
};
