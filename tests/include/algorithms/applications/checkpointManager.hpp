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

#include <cstdlib>
#include <iostream>
#include <fstream>

#include <codi.hpp>

struct TestCheckpoint : public codi::algorithms::CheckpointBase {
  private:

    std::vector<double> data;

  public:

    TestCheckpoint(int iteration, std::vector<double> const& data) :
      codi::algorithms::CheckpointBase(iteration),
      data(data)
    {}

    std::vector<double> const& getData() {
      return data;
    }
};


template<typename T_Type, typename T_Application>
struct TestCheckpointManager : public codi::algorithms::CheckpointManagerInterface {
  public:
    using Application = CODI_DD(T_Application, CODI_T(codi::algorithms::ApplicationInterface<CODI_ANY>));
    using Type = CODI_DD(T_Type, CODI_T(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

    using Checkpoint = TestCheckpoint;

  private:

    Application* app;

  public:

    TestCheckpointManager(Application* app) :
      app(app)
    {}

    Checkpoint* create() {
      std::vector<double> data(app->getSizeY());

      app->iterateY([&](Type& value, size_t pos){
        data[pos] = codi::RealTraits::getPassiveValue(value);
      });

      return new Checkpoint(app->getIteration(), data);
    }

    std::vector<codi::algorithms::CheckpointHandle*> list() {
      return std::vector<codi::algorithms::CheckpointHandle*>(0);
    }

    void load(codi::algorithms::CheckpointHandle* cp) {
      Checkpoint* check = cast(cp);

      app->iterateY([=](Type& value, size_t pos){
        value = check->getData()[pos];
      });

      app->setIteration(check->getIeration());
    }

    void remove(codi::algorithms::CheckpointHandle* cp) {
      Checkpoint* check = cast(cp);

      delete check;
    }

    void write(codi::algorithms::CheckpointHandle* cp) { (void)cp; } // Do nothing
    void read(codi::algorithms::CheckpointHandle* cp) { (void)cp; } // Do nothing

  private:

    Checkpoint* cast(codi::algorithms::CheckpointHandle* cp) {
      return static_cast<Checkpoint*>(cp);
    }
};
