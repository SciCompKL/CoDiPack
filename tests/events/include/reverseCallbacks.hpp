/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <codi.hpp>

#include "string_conversions.hpp"

namespace ReverseCallbacks {

  /*******************************************************************************/
  /// @name AD workflow
  /// @{

  template<typename Tape>
  void onTapeStartRecording(Tape&, void*) {
    std::cout << "TapeStartRecording" << std::endl;
  }

  template<typename Tape>
  void onTapeStopRecording(Tape&, void*) {
    std::cout << "TapeStopRecording" << std::endl;
  }

  template<typename Tape>
  void onTapeRegisterInput(Tape&, typename Tape::Real& value, typename Tape::Identifier& identifier, void*) {
    std::cout << "TapeRegisterInput value " << value << " identifier " << identifier << std::endl;
  }

  template<typename Tape>
  void onTapeRegisterOutput(Tape&, typename Tape::Real& value, typename Tape::Identifier& identifier, void*) {
    std::cout << "TapeRegisterOutput value " << value << " identifier " << identifier << std::endl;
  }

  template<typename Tape>
  void onTapeEvaluate(Tape&, typename Tape::Position const& start, typename Tape::Position const& end,
                      codi::VectorAccessInterface<typename Tape::Real, typename Tape::Identifier>*,
                      codi::EventHints::EvaluationKind direction, codi::EventHints::Endpoint endpoint, void*) {
    std::cout << "TapeEvaluate " << to_string(direction) << " " << to_string(endpoint) << " from " << start << " to "
              << end << std::endl;
  }

  template<typename Tape>
  void onTapeReset(Tape&, typename Tape::Position const& position, codi::EventHints::Reset kind, bool clearAdjoints,
                   void*) {
    std::cout << "TapeReset " << to_string(kind) << " position " << position << " clear adjoints " << clearAdjoints
              << std::endl;
  }

  /// @}
  /*******************************************************************************/
  /// @name Preaccumulation
  /// @{

  template<typename Tape>
  void onPreaccStart(Tape&, void*) {
    std::cout << "PreaccStart" << std::endl;
  }

  template<typename Tape>
  void onPreaccFinish(Tape&, void* customData) {
    std::cout << "PreaccFinish" << std::endl;
  }

  template<typename Tape>
  void onPreaccAddInput(Tape&, typename Tape::Real const& value, typename Tape::Identifier const& identifier, void*) {
    std::cout << "PreaccAddInput value " << value << " identifier " << identifier << std::endl;
  }

  template<typename Tape>
  void onPreaccAddOutput(Tape&, typename Tape::Real& value, typename Tape::Identifier& identifier, void*) {
    std::cout << "PreaccAddOutput value " << value << " identifier " << identifier << std::endl;
  }

  /// @}
  /*******************************************************************************/
  /// @name Statements
  /// @{

  template<typename Tape>
  void onStatementPrimal(Tape&, typename Tape::Real const& lhsValue, typename Tape::Identifier const& lhsIdentifier,
                         typename Tape::Real const& newValue, codi::EventHints::Statement statement, void*) {
    std::cout << "StatementPrimal " << to_string(statement) << " lhsValue " << lhsValue << " lhsIdentifier "
              << lhsIdentifier << " newValue " << newValue << std::endl;
  }

  template<typename Tape>
  struct GlobalStatementCounters {
    public:
      static size_t storeOnTape;
      static size_t evaluate;
      static size_t evaluatePrimal;

      static void assertEqual() {
        if (storeOnTape != evaluate) {
          std::cerr << "StatementStoreOnTape count (" << storeOnTape << ") does not match StatementEvaluate count ("
                    << evaluate << ")" << std::endl;
          abort();
        }
        if (evaluatePrimal != 0 && storeOnTape != evaluate) {
          std::cerr << "StatementStoreOnTape count (" << storeOnTape
                    << ") does not match StatementEvaluatePrimal count (" << evaluatePrimal << ")" << std::endl;
          abort();
        }
      }
  };

  template<typename Tape>
  size_t GlobalStatementCounters<Tape>::storeOnTape = 0;

  template<typename Tape>
  size_t GlobalStatementCounters<Tape>::evaluate = 0;

  template<typename Tape>
  size_t GlobalStatementCounters<Tape>::evaluatePrimal = 0;

  template<typename Tape>
  void onStatementStoreOnTape(Tape&, typename Tape::Identifier const& lhsIdentifier,
                              typename Tape::Real const& newValue, size_t numActiveVariables,
                              typename Tape::Identifier const* rhsIdentifiers, typename Tape::Real const* jacobians,
                              void*) {
    std::cout << "StatementStoreOnTape lhsIdentifier " << lhsIdentifier << " newValue " << newValue
              << " numActiveVariables " << numActiveVariables << std::endl
              << "\t";
    for (size_t i = 0; i < numActiveVariables; ++i) {
      if (i != 0) {
        std::cout << " ";
      }
      std::cout << rhsIdentifiers[i] << " " << jacobians[i] << ";";
    }
    std::cout << std::endl;

    ++GlobalStatementCounters<Tape>::storeOnTape;
  }

  template<typename Tape>
  void onStatementEvaluate(Tape&, typename Tape::Identifier const& lhsIdentifier, size_t numAdjoints,
                           typename Tape::Real const* adjoints, void*) {
    std::cout << "StatementEvaluate lhsIdentifier " << lhsIdentifier << " numAdjoints " << numAdjoints << std::endl
              << "\t";
    for (size_t i = 0; i < numAdjoints; ++i) {
      if (i != 0) {
        std::cout << " ";
      }
      std::cout << adjoints[i];
    }
    std::cout << std::endl;

    ++GlobalStatementCounters<Tape>::evaluate;
  }

  template<typename Tape>
  void onStatementEvaluatePrimal(Tape&, typename Tape::Identifier const& lhsIdentifier,
                                 typename Tape::Real const& lhsValue, void*) {
    std::cout << "StatementEvaluatePrimal lhsIdentifier " << lhsIdentifier << " lhsValue " << lhsValue << std::endl;

    ++GlobalStatementCounters<Tape>::evaluatePrimal;
  }

  /// @}
  /*******************************************************************************/
  /// @name Index management
  /// @{

  template<typename Tape>
  void onIndexAssign(typename Tape::Identifier const& index, void*) {
    std::cout << "IndexAssign index " << index << std::endl;
  }

  template<typename Tape>
  void onIndexFree(typename Tape::Identifier const& index, void*) {
    std::cout << "IndexFree index " << index << std::endl;
  }

  template<typename Tape>
  void onIndexCopy(typename Tape::Identifier const& index, void*) {
    std::cout << "IndexCopy index " << index << std::endl;
  }

  /// @}

  template<typename Tape>
  std::list<typename codi::EventSystem<Tape>::Handle> registerAll() {
    std::list<typename codi::EventSystem<Tape>::Handle> handles;
    handles.push_back(codi::EventSystem<Tape>::registerTapeStartRecordingListener(onTapeStartRecording<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerTapeStopRecordingListener(onTapeStopRecording<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerTapeRegisterInputListener(onTapeRegisterInput<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerTapeRegisterOutputListener(onTapeRegisterOutput<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerTapeEvaluateListener(onTapeEvaluate<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerTapeResetListener(onTapeReset<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerPreaccStartListener(onPreaccStart<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerPreaccFinishListener(onPreaccFinish<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerPreaccAddInputListener(onPreaccAddInput<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerPreaccAddOutputListener(onPreaccAddOutput<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerStatementPrimalListener(onStatementPrimal<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerStatementStoreOnTapeListener(onStatementStoreOnTape<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerStatementEvaluateListener(onStatementEvaluate<Tape>));
    handles.push_back(
        codi::EventSystem<Tape>::registerStatementEvaluatePrimalListener(onStatementEvaluatePrimal<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerIndexAssignListener(onIndexAssign<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerIndexFreeListener(onIndexFree<Tape>));
    handles.push_back(codi::EventSystem<Tape>::registerIndexCopyListener(onIndexCopy<Tape>));
    return handles;
  }

}
