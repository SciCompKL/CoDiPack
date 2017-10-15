/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */


#pragma once

#include <vector>

#include "../configure.h"
#include "../exceptions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   *
   *  Currently not working for vector types.
   *  Currently not tested for hihger order derivatives.
   */
  template<typename CoDiType>
  struct PreaccumulationHelper {

      typedef typename CoDiType::Real Real;
      typedef typename CoDiType::GradientData GradientData;
      typedef typename CoDiType::GradientValue GradientValue;

      typedef typename CoDiType::TapeType Tape;
      typedef typename Tape::Position Position;

      std::vector<GradientData> inputData;
      std::vector<CoDiType*> outputData;
      Position startPos;

      std::vector<GradientValue> storedAdjoints;

      std::vector<Real> jacobie;
      std::vector<int> nonZeros;

      template<typename ... Inputs>
      void addInput(const Inputs& ... inputs) {
        Tape& tape = CoDiType::getGlobalTape();

        if(tape.isActive()) {
          addInputRec(inputs...);
        }
      }

      template<typename ... Inputs>
      void start(const Inputs& ... inputs) {
        Tape& tape = CoDiType::getGlobalTape();

        if(tape.isActive()) {
          inputData.clear();
          outputData.clear();

          startPos = tape.getPosition();

          addInputRec(inputs...);
        }
      }

      template<typename ... Outputs>
      void addOutput(Outputs& ... outputs) {
        Tape& tape = CoDiType::getGlobalTape();

        if(tape.isActive()) {
          addOutputRec(outputs...);
        }
      }

      template<typename ... Outputs>
      void finish(bool storeAdjoints, Outputs& ... outputs) {
        Tape& tape = CoDiType::getGlobalTape();

        if(tape.isActive()) {

          addOutputRec(outputs...);

          if(storeAdjoints) {
            storeInputAdjoints();
          }

          doPreaccumulation();

          if(storeAdjoints) {
            restoreInputAdjoints();
          }
        }
      }

    private:

      void addInputLogic(const CoDiType& input) {
        const GradientData& gradData = input.getGradientData();
        if(0 != gradData) {
          inputData.push_back(gradData);
        }
      }

      void addInputRec() {
        // terminator implementation
      }

      template<typename ... Inputs>
      void addInputRec(const CoDiType& input, const Inputs& ... r) {
        addInputLogic(input);
        addInputRec(r...);
      }

      void addOutputLogic(CoDiType& output) {
        const GradientData& gradData = output.getGradientData();
        if(0 != gradData) {
          outputData.push_back(&output);
        }
      }

      void addOutputRec() {
        // terminator implementation
      }

      template<typename ... Outputs>
      void addOutputRec(CoDiType& output, Outputs& ... r) {
        addOutputLogic(output);
        addOutputRec(r...);
      }

      void storeInputAdjoints() {
        Tape& tape = CoDiType::getGlobalTape();

        if(storedAdjoints.size() < inputData.size()) {
          storedAdjoints.resize(inputData.size());
        }

        for(size_t i = 0; i < inputData.size(); ++i) {
          GradientData index = inputData[i];
          GradientValue& adjoint = tape.gradient(index);
          storedAdjoints[i] = adjoint;
          adjoint = GradientValue();
        }
      }

      void restoreInputAdjoints() {
        Tape& tape = CoDiType::getGlobalTape();

        for(size_t i = 0; i < inputData.size(); ++i) {
          GradientData index = inputData[i];
          tape.gradient(index) = storedAdjoints[i];
        }
      }

      void doPreaccumulation() {
        // perform the accumulation of the tape part
        Tape& tape = CoDiType::getGlobalTape();

        Position endPos = tape.getPosition();
        size_t jacobiSize = inputData.size() * outputData.size();
        if(jacobie.size() < jacobiSize) {
          jacobie.resize(jacobiSize);
        }
        if(nonZeros.size() < outputData.size()) {
          nonZeros.resize(outputData.size());
        }

        for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {

          nonZeros[curOut] = 0;
          size_t jacobiOffset = curOut * inputData.size();

          GradientData indexOut = outputData[curOut]->getGradientData();
          tape.setGradient(indexOut, 1.0);
          tape.evaluatePreacc(endPos, startPos);

          for (size_t curIn = 0; curIn < inputData.size(); ++curIn) {

            GradientData indexIn = inputData[curIn];
            GradientValue& adj = tape.gradient(indexIn);
            jacobie[curIn + jacobiOffset] = adj;

            if(0.0 != adj) {
              nonZeros[curOut] += 1;
            }

            adj = 0.0;
          }

          tape.clearAdjoints(endPos, startPos);
        }

        // store the Jacobi matrix
        tape.reset(startPos);

        for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {

          CoDiType& value = *outputData[curOut];
          if(0 != nonZeros[curOut]) {

            int nonZerosLeft = nonZeros[curOut];
            GradientData lastGradientData = value.getGradientData();
            bool staggeringActive = false;
            int curIn = 0;
            size_t jacobiOffset = curOut * inputData.size();

            // push statements as long as there are nonzeros left
            // if there are more than MaxStatementIntValue nonzeros, then we neeed to stager the
            // statement pushes
            while(nonZerosLeft > 0) {

              // calucluate the number of Jacobies for this statement
              int jacobiesForStatement = nonZerosLeft;
              if(jacobiesForStatement >= (int)MaxStatementIntValue) {

                jacobiesForStatement = MaxStatementIntValue;
                if(staggeringActive) { /* Space is used up but we need one Jacobi for the staggering */
                  jacobiesForStatement -= 1;
                }
              }
              nonZerosLeft -= jacobiesForStatement; /* update nonzeros so that we now if it is the last round */

              GradientData storedGradientData = lastGradientData;
              tape.storeManual(value.getValue(), lastGradientData, jacobiesForStatement + (int)staggeringActive);
              if(staggeringActive) { /* Not the first stagering so push the last output */
                tape.pushJacobiManual(1.0, 0.0, storedGradientData);
              }

              // push the rest of the Jacobies for the statement
              while(jacobiesForStatement > 0) {
                if(0.0 != jacobie[curIn + jacobiOffset]) {
                  tape.pushJacobiManual(jacobie[curIn + jacobiOffset], 0.0, inputData[curIn]);
                  jacobiesForStatement -= 1;
                }
                curIn += 1;
              }

              staggeringActive = true;
            }

            value.getGradientData() = lastGradientData; /* now set gradient data for the real output value */
          } else {
            // disable tape index since there is no dependency
            tape.destroyGradientData(value.value(), value.getGradientData());
          }
        }
      }
  };

  template<typename CoDiType>
  struct ForwardPreaccumulationHelper {

      template<typename ... Inputs>
      void addInput(const Inputs& ... inputs) {
        // do nothing
      }

      template<typename ... Inputs>
      void start(const Inputs& ... inputs) {
        // do nothing
      }

      template<typename ... Outputs>
      void addOutput(Outputs& ... outputs) {
        // do nothing
      }

      template<typename ... Outputs>
      void finish(bool storeAdjoints, Outputs& ... outputs) {
        CODI_UNUSED(storeAdjoints);

        // do nothing
      }
  };
}
