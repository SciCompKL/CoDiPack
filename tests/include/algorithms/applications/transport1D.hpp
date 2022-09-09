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

#include <codi.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "applicationBase.hpp"

template<typename T_Type>
struct Transport1DSettings {
  public:
    using Type = CODI_DD(T_Type, CODI_T(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

    Type rho = 1.0;
    Type gamma = 0.1;
    Type L = 1.0;

    int N = 101;
    int maxIterSolve = 100;
    int maxT = 10000;

    int functionalNumber = 1;

    std::vector<Type> control = std::vector<Type>(101, 1.0);
};

template<typename T_Type>
struct Transport1D : public TestApplicationBase<T_Type, Transport1D<T_Type>> {
  public:
    using Type = CODI_DD(T_Type, CODI_T(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
    using Base = TestApplicationBase<Type, Transport1D>;
    using Res = typename Base::Res;
    using Settings = Transport1DSettings<Type>;

    static size_t constexpr FunctionalMax = 16;

    Settings settings;

    Type dx;
    Type dt;

    Type coeff;

    std::vector<Type> A_d;
    std::vector<Type> A_m;
    std::vector<Type> A_p;
    std::vector<Type> b;

    std::vector<Type> phi;
    std::vector<double> phi_old;
    Type z[FunctionalMax];

    double res;

    Transport1D(Settings settings) : Base(this), settings(settings) {}

    Transport1D() : Base(this), settings() {}

  private:

    void solve(std::vector<Type>& phi) {
      for (int k = 0; k < settings.maxIterSolve; k++) {
        for (int i = 1; i < settings.N - 1; i++) {
          phi[i] = (b[i] - A_m[i] * phi[i - 1] - A_p[i] * phi[i + 1]) / A_d[i];
        }
      }
    }

    Type computeExplicitTerm(int i) {
      Type temp = -settings.control[i] * (phi[i + 1] - phi[i - 1]) / (2.0 * dx) +
                  coeff * ((phi[i + 1] - 2.0 * phi[i] + phi[i - 1]) / (dx * dx));
      return temp * 0.5 + phi[i] / dt;
    }

  public:

    void initialize() {
      dx = settings.L / (settings.N - 1);
      dt = dx;
      coeff = settings.gamma / settings.rho;

      A_d.resize(settings.N);
      A_m.resize(settings.N);
      A_p.resize(settings.N);
      b.resize(settings.N);
      phi.resize(settings.N);
      phi_old.resize(settings.N);

      for (int i = 0; i < settings.N; i++) {
        settings.control[i] = 1.0;
        if(Base::getHints() & codi::algorithms::ApplicationFlags::InitializationComputesP) {
          this->handleInitializationVariable(settings.control[i]);
        }
      }

      evaluateP();

      for (int i = 0; i < settings.N; i++) {
        phi[i] = 0.0;
      }
      phi[0] = 0.0;
      phi[settings.N - 1] = 1.0;
      Base::setIteration(0);
    }

    void evaluateG() {
      for (int i = 0; i < settings.N; i++) {
        phi_old[i] = codi::RealTraits::getPassiveValue(phi[i]);
      }

      for (int i = 1; i < settings.N - 1; i++) {
        b[i] = computeExplicitTerm(i);
      }
      b[0] = 0.0;
      b[settings.N - 1] = 1.0;

      solve(phi);
      phi[0] = 0.0;
      phi[settings.N - 1] = 1.0;

      res = 0.0;
      for (int i = 0; i < settings.N; i++) {
        double diff = codi::RealTraits::getPassiveValue(phi[i]) - phi_old[i];
        res += diff * diff;
      }
      res = sqrt(res);

      if(!(Base::hints & codi::algorithms::ApplicationFlags::FComputationIsAvailable)) {
        evaluateF();
      }

      Base::setIteration(Base::getIteration() + 1);
    }

    void evaluateF() {
      z[0] = 0.0;
      z[1] = 0.0;
      z[0] += 0.5 * sin(0) * phi[0];
      z[1] += 0.5 * cos(0) * phi[0];
      for (int i = 1; i < settings.N - 1; i += 1) {
        z[0] += sin(i * dx) * phi[i];
        z[1] += cos(i * dx) * phi[i];
      }
      z[0] += 0.5 * sin((settings.N - 1) * dx) * phi[settings.N - 1];
      z[1] += 0.5 * sin((settings.N - 1) * dx) * phi[settings.N - 1];

      z[0] *= dx;
      z[1] *= dx;

      for(size_t i = 2; i < FunctionalMax; i += 1) {
        z[i] = phi[settings.N - 1 - i];
      }
    }

    void evaluateP() {
      for (int i = 1; i < settings.N - 1; i += 1) {
        A_d[i] = 1.0 / dt + coeff / (dx * dx);
        A_m[i] = -settings.control[i] / (4.0 * dx) - (0.5 * coeff) / (dx * dx);
        A_p[i] = settings.control[i] / (4.0 * dx) - (0.5 * coeff) / (dx * dx);
      }
      // Initialize boundary conditions
      A_d[0] = 1.0;
      A_d[settings.N - 1] = 1.0;
    }

    template<typename Func>
    void iterateY(Func&& func) {
      for (int i = 0; i < settings.N; i += 1) {
        func(phi[i], i);
      }
    }

    template<typename Func>
    void iterateX(Func&& func) {
      for (int i = 0; i < settings.N; i += 1) {
        func(settings.control[i], i);
      }
    }

    template<typename Func>
    void iterateP(Func&& func) {
      for (int i = 0; i < settings.N; i += 1) {
        size_t pos = 3 * i;
        func(A_d[i], pos);
        func(A_p[i], pos + 1);
        func(A_m[i], pos + 2);
      }
    }

    template<typename Func>
    void iterateZ(Func&& func) {
      for (int i = 0; i < settings.functionalNumber; i += 1) {
        func(z[i], i);
      }
    }

    size_t getSizeY() {
      return settings.N;
    }

    size_t getSizeX() {
      return settings.N;
    }

    size_t getSizeP() {
      return settings.N * 3;
    }

    size_t getSizeZ() {
      return settings.functionalNumber;
    }

    int getNumberOfFunctionals() { return settings.functionalNumber; }

    void runPrimal() {
      initialize();

      Base::print(codi::StringUtil::format("Iter Res\n"));

      for (int T = 0; T < settings.maxT; T++) {
        evaluateG();

        Base::print(codi::StringUtil::format("%d %0.6e\n", Base::getIteration(), res));

        if (res < 0.00000001) {
          break;
        }
      }

      evaluateF();
    }
};
