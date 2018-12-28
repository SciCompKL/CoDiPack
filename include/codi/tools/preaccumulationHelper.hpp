/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
   * @brief Stores the Jacobi matrix for a code section.
   *
   * The preaccumulation of a code section describes the process of replacing the recorded tape entries with the Jacobi
   * matrix of that section. If the code part is defined by the function
   * \f[ y = f(x) \eqdot \f]
   * then the reverse AD mode needs to compute the reverse AD equation
   * \f[ \bar{x} = \frac{df}{dx}^T(x)\cdot \bar{y}. \f]
   * for this section. If nothing is changed then for \f$ f \f$ several statements and arguments are recorded.
   * If the computation of requires 200 statements with a total of 600 arguments the storage for this is on a Jacobi tape
   * would be 7400 byte. If the function has two input arguments and two output arguments, the storage for the Jacobi
   * matrix of this function would require 50 byte.
   *
   * The procedure for the preaccumulation of a code section is:
   *
   * \code{.cpp}
   * PreaccumulationHelper<CoDiType> ph;
   *
   * ph.start(<list of input arguments>); // can also be empty
   * ph.addInput(<list of input arguments>); // optional, can be called multiple times
   *
   * ... < code section that will be preaccumulated >
   *
   * ph.addOutput(<list of output arguments>); // optional, can be called multiple times
   * ph.finish(false, <list of output arguments>); // can also be empty
   * \endcode
   *
   * The first argument to #finish is most of the times false. If it is set to true then the current adjoint values
   * in the tape will be stored before the preaccumulation is done and restored afterwards.
   *
   * The preaccumulation helper can be used multiple times, the start routine resets the state such that multiple
   * evaluations are possible. This improves the performance of the helper since stack allocations are only performed
   * once.
   *
   * Restrictions:
   *  - Currently not working for vector types.
   *  - Currently not tested for higher order derivatives.
   *
   * @tparam CoDiType  This needs to be one of the CoDiPack types defined through an ActiveReal
   */
  template<typename CoDiType>
  struct PreaccumulationHelper {

      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientData GradientData;  /**< The type for the gradient identification */
      typedef typename CoDiType::GradientValue GradientValue;  /**< The type for the gradient computation */

      typedef typename CoDiType::TapeType Tape; /**< The type for the tape */
      typedef typename Tape::Position Position; /**< The type for the position in the tape */

      std::vector<GradientData> inputData; /**< The identifiers for the input data of the preaccumulation section. */
      std::vector<CoDiType*> outputData; /**< The pointers to the output values of the preaccumulation section. */
      Position startPos; /**< The starting point for the preaccumulation. */

      std::vector<GradientValue> storedAdjoints; /**< The old values of the adjoints for the values of the preaccumulation */

      std::vector<Real> jacobie; /**< The Jacobi matrix used to hold the result of the preaccumulation. */
      std::vector<int> nonZeros; /**< The number of nonzero values for each output value in the Jacobi matrix. */

      /**
       * @brief Add extra inputs to the preaccumulated section.
       *
       * This function needs to be called after start and before any computations are performed.
       *
       * @param[in] inputs  Inputs for the preaccumulated section.
       *
       * @tparam Inputs  The data type for the inputs. This needs to be CoDiType given as a template parameter to the
       *                 helper.
       */
      template<typename ... Inputs>
      void addInput(const Inputs& ... inputs) {
        Tape& tape = CoDiType::getGlobalTape();

        if(tape.isActive()) {
          addInputRec(inputs...);
        }
      }

      /**
       * @brief Starts the section for the preaccumulation.
       *
       * This is the first function that needs to be called.
       *
       * @param[in] inputs  Inputs for the preaccumulated section.
       *
       * @tparam Inputs  The data type for the inputs. This needs to be CoDiType given as a template parameter to the
       *                 helper.
       */
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

      /**
       * @brief Adds extra outputs to the section of the preaccumulation.
       *
       * This function needs to be called after all computations in the preaccumulated section have been evaluated.
       *
       * @param[in,out] outputs  Outputs for the preaccumulated section.
       *
       * @tparam Outputs  The data type for the outputs. This needs to be CoDiType given as a template parameter to the
       *                  helper.
       */
      template<typename ... Outputs>
      void addOutput(Outputs& ... outputs) {
        Tape& tape = CoDiType::getGlobalTape();

        if(tape.isActive()) {
          addOutputRec(outputs...);
        }
      }

      /**
       * @brief Performs the preaccumulation of the code region defined between the start call and this function call.
       *
       * @param[in] storeAdjoints  Usually false, only required to be true if the preaccumulation takes place during an
       *                           tape evaluation.
       * @param[in,out]   outputs  Outputs for the preaccumulated section.
       *
       * @tparam Outputs  The data type for the outputs. This needs to be CoDiType given as a template parameter to the
       *                  helper.
       */
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

      /**
       * @brief Helper function for adding an input value.
       *
       * @param[in] input  The input value.
       */
      void addInputLogic(const CoDiType& input) {
        const GradientData& gradData = input.getGradientData();
        if(0 != gradData) {
          inputData.push_back(gradData);
        }
      }

      /** @brief terminator for the recursive implementation */
      void addInputRec() {
        // terminator implementation
      }

      /**
       * @brief Add all inputs in a recursive manner.
       *
       * @param[in] input  The input that is added during this recursive call.
       * @param[in]     r  The reminder that still needs to be added.
       *
       * @tparam Inputs  Input types that still need to be handled.
       */
      template<typename ... Inputs>
      void addInputRec(const CoDiType& input, const Inputs& ... r) {
        addInputLogic(input);
        addInputRec(r...);
      }

      /**
       * @brief Helper function for adding an output value.
       *
       * @param[in,out] output  The output value.
       */
      void addOutputLogic(CoDiType& output) {
        const GradientData& gradData = output.getGradientData();
        if(0 != gradData) {
          outputData.push_back(&output);
        }
      }

      /** @brief terminator for the recursive implementation */
      void addOutputRec() {
        // terminator implementation
      }

      /**
      * @brief Add all outputs in a recursive manner.
      *
      * @param[in,out] output  The output that is added during this recursive call.
      * @param[in]          r  The reminder that still needs to be added.
      *
      * @tparam Outputs  Output types that still need to be handled.
      */
      template<typename ... Outputs>
      void addOutputRec(CoDiType& output, Outputs& ... r) {
        addOutputLogic(output);
        addOutputRec(r...);
      }

      /**
       * @brief Helper function for storing the adjoints of all input values.
       *
       * The storing of input adjoints is required if the tape already contains meaning full values in the adjoint vector.
       */
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

      /**
       * @brief Helper function for restoring the adjoints of all input values.
       *
       * See also storeInputAdjoints.
       */
      void restoreInputAdjoints() {
        Tape& tape = CoDiType::getGlobalTape();

        for(size_t i = 0; i < inputData.size(); ++i) {
          GradientData index = inputData[i];
          tape.gradient(index) = storedAdjoints[i];
        }
      }

      /**
       * @brief Performs the actual preaccumulation.
       *
       * The first part of the function perform the preaccumulation of the code section that is enclosed between the
       * start and finish call. The sections creates the Jacobi matrix with respect to all input and output values.
       *
       * The second part of the function stores the Jacobi matrix on the tape. Since CoDiPack has a limit how many input
       * arguments a statement can have, the storing of the matrix may need to be staggered. This is done such that each
       * statement contains the previous left hand side value. Only the last statement will have the output value as the
       * left hand side.
       */
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


        if(inputData.size() < outputData.size()) {
          // forward accumulation of Jacobi

          for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {
            nonZeros[curOut] = 0;
          }

          for (size_t curIn = 0; curIn < inputData.size(); ++curIn) {

            GradientData indexIn = inputData[curIn];
            tape.setGradient(indexIn, 1.0);
            tape.evaluateForwardPreacc(startPos, endPos);

            for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {

              GradientData indexOut = outputData[curOut]->getGradientData();
              GradientValue& adj = tape.gradient(indexOut);
              jacobie[curIn + curOut * inputData.size()] = adj;

              if(0.0 != adj) {
                nonZeros[curOut] += 1;
              }
            }

            tape.setGradient(indexIn, 0.0);

            tape.clearAdjoints(endPos, startPos);
          }
        } else {
          // reverse accumulation of Jacobi

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
        }

        // store the Jacobi matrix
        tape.reset(startPos);

        for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {

          CoDiType& value = *outputData[curOut];
          if(0 != nonZeros[curOut]) {

            int nonZerosLeft = nonZeros[curOut];
            // we need to use here the value of the gradient data such that it is correctly deleted.
            GradientData lastGradientData = value.getGradientData();
            bool staggeringActive = false;
            int curIn = 0;
            size_t jacobiOffset = curOut * inputData.size();

            // push statements as long as there are non zeros left
            // if there are more than MaxStatementIntValue non zeros, then we need to stagger the
            // statement pushes
            while(nonZerosLeft > 0) {

              // calculate the number of Jacobies for this statement
              int jacobiesForStatement = nonZerosLeft;
              if(jacobiesForStatement >= (int)MaxStatementIntValue) {

                jacobiesForStatement = MaxStatementIntValue;
                if(staggeringActive) { /* Space is used up but we need one Jacobi for the staggering */
                  jacobiesForStatement -= 1;
                }
              }
              nonZerosLeft -= jacobiesForStatement; /* update non zeros so that we now if it is the last round */

              GradientData storedGradientData = lastGradientData;
              tape.storeManual(value.getValue(), lastGradientData, jacobiesForStatement + (int)staggeringActive);
              if(staggeringActive) { /* Not the first staggering so push the last output */
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

  /**
   * Helper implementation of the same interface as the PreaccumulationHelper for forward AD tapes.
   *
   * This implementation does nothing in all methods.
   * @tparam CoDiType
   */
  template<typename CoDiType>
  struct ForwardPreaccumulationHelper {

      /**
       * @brief Does nothing
       * @param[in] inputs  Not used
       * @tparam Inputs  Not used.
       */
      template<typename ... Inputs>
      void addInput(const Inputs& ... inputs) {
        CODI_UNUSED_VAR(inputs...);
        // do nothing
      }

      /**
       * @brief Does nothing
       * @param[in] inputs  Not used
       * @tparam Inputs  Not used.
       */
      template<typename ... Inputs>
      void start(const Inputs& ... inputs) {
        CODI_UNUSED_VAR(inputs...);
        // do nothing
      }

      /**
       * @brief Does nothing
       * @param[in,out] outputs  Not used
       * @tparam Outputs  Not used.
       */
      template<typename ... Outputs>
      void addOutput(Outputs& ... outputs) {
        CODI_UNUSED_VAR(outputs...);
        // do nothing
      }

      /**
       * @brief Does nothing
       * @param[in] storeAdjoints  Not used
       * @param[in,out]   outputs  Not used
       * @tparam Outputs  Not used.
       */
      template<typename ... Outputs>
      void finish(bool storeAdjoints, Outputs& ... outputs) {
        CODI_UNUSED(storeAdjoints);
        CODI_UNUSED_VAR(outputs...);

        // do nothing
      }
  };
}
