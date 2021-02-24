/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include <sstream>

#include <opdi/tool/toolInterface.hpp>

template<typename _CoDiType>
struct CoDiOpDiTool : public opdi::ToolInterface {
  public:
    using CoDiType = _CoDiType;
    using Tape = typename CoDiType::TapeType;
    using Position = typename Tape::Position;

  private:
    static void callHandleReverse(void*, void* handlePtr, void*) {
      opdi::Handle* handle = (opdi::Handle*) handlePtr;
      handle->reverseFunc(handle->data);
    }

    static void callHandleDelete(void*, void* handlePtr) {
      opdi::Handle* handle = (opdi::Handle*) handlePtr;
      if (handle->deleteFunc != nullptr) {
        handle->deleteFunc(handle->data);
      }
      delete handle;
    }

  public:
    void* createTape() {
      return (void*) new Tape;
    }

    void deleteTape(void* tapePtr) {
      Tape* tape = (Tape*) tapePtr;
      delete tape;
    }

    void* allocPosition() {
      return new Position();
    }

    void freePosition(void* positionPtr) {
      Position* position = (Position*) positionPtr;
      delete position;
    }

    size_t getPositionSize() {
      return sizeof(Position);
    }

    std::string positionToString(void* positionPtr) {
      Position* position = (Position*) positionPtr;
      std::stringstream conv;
      conv << *position;
      return conv.str();
    }

    void getTapePosition(void* tapePtr, void* positionPtr) {
      Tape* tape = (Tape*) tapePtr;
      Position* position = (Position*) positionPtr;

      *position = tape->getPosition();
    }

    void getZeroPosition(void* tapePtr, void* positionPtr) {
      Tape* tape = (Tape*) tapePtr;
      Position* position = (Position*) positionPtr;

      *position = tape->getZeroPosition();
    }

    void copyPosition(void* dstPtr, void* srcPtr) {
      Position* dst = (Position*) dstPtr;
      Position* src = (Position*) srcPtr;

      *dst = *src;
    }

    int comparePosition(void* lhsPtr, void* rhsPtr) {
      Position* lhs = (Position*) lhsPtr;
      Position* rhs = (Position*) rhsPtr;

      if (*lhs <= *rhs) {
        if (*rhs <= *lhs) {
          return 0;
        }
        else {
          return -1;
        }
      }
      else {
        return 1;
      }
    }

    bool isActive(void* tapePtr) {
      Tape* tape = (Tape*) tapePtr;
      return tape->isActive();
    }

    void setActive(void* tapePtr, bool active) {
      Tape* tape = (Tape*) tapePtr;
      if (active) {
        tape->setActive();
      }
      else {
        tape->setPassive();
      }
    }

    void evaluate(void* tapePtr, void* startPtr, void* endPtr, bool useAtomics = true) {
      Tape* tape = (Tape*) tapePtr;
      Position* start = (Position*) startPtr;
      Position* end = (Position*) endPtr;

      typename Tape::GradientValue* adjoints = &tape->gradient(0);
      using NonAtomicGradientValue = codi::RemoveAtomic<typename Tape::GradientValue>;
      using AtomicGradientValue = codi::Atomic<NonAtomicGradientValue>;

      if (useAtomics) {
        AtomicGradientValue* safeAdjoints = (AtomicGradientValue*) adjoints;
        tape->evaluate(*start, *end, safeAdjoints);
      }
      else {
        NonAtomicGradientValue* unsafeAdjoints = (NonAtomicGradientValue*) adjoints;
        tape->evaluate(*start, *end, unsafeAdjoints);
      }
    }

    void reset(void* tapePtr, bool clearAdjoints) {
      Tape* tape = (Tape*) tapePtr;
      tape->reset(clearAdjoints);
    }

    void reset(void* tapePtr, void* positionPtr, bool clearAdjoints) {
      Tape* tape = (Tape*) tapePtr;
      Position* position = (Position*) positionPtr;
      tape->reset(*position, clearAdjoints);
    }

    void* getThreadLocalTape() {
      return (void*) CoDiType::getGlobalTapePtr();
    }

    void setThreadLocalTape(void* tapePtr) {
      Tape* tape = (Tape*) tapePtr;
      CoDiType::setGlobalTapePtr(tape);
    }

    void pushExternalFunction(void* tapePtr, opdi::Handle const* handle) {
      Tape* tape = (Tape*) tapePtr;
      tape->pushExternalFunctionHandle(CoDiOpDiTool::callHandleReverse, (void*) handle, CoDiOpDiTool::callHandleDelete);
    }

    void erase(void* tapePtr, void* startPtr, void* endPtr) {
      Tape* tape = (Tape*) tapePtr;
      Position* start = (Position*) startPtr;
      Position* end = (Position*) endPtr;

      tape->erase(*start, *end);
    }

    void append(void* dstTapePtr, void* srcTapePtr, void* startPtr, void* endPtr) {
      Tape* dstTape = (Tape*) dstTapePtr;
      Tape* srcTape = (Tape*) srcTapePtr;
      Position* start = (Position*) startPtr;
      Position* end = (Position*) endPtr;

      dstTape->append(*srcTape, *start, *end);
    }
};
