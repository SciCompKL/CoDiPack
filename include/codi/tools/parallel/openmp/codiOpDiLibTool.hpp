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

#include <opdi/tool/toolInterface.hpp>
#include <sstream>

#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../expressions/parallelActiveType.hpp"
#include "../../../misc/macros.hpp"
#include "../../../tapes/interfaces/editingTapeInterface.hpp"
#include "../../../tapes/misc/vectorAccessInterface.hpp"
#include "../../../traits/atomicTraits.hpp"
#include "openMPAtomic.hpp"

#ifndef DOXYGEN_DISABLE

template<typename T_CoDiType>
struct CoDiOpDiLibTool : public opdi::ToolInterface {
  public:
    using CoDiType = CODI_DD(T_CoDiType, codi::CODI_DEFAULT_PARALLEL_ACTIVE_TYPE);
    using Real = typename CoDiType::Real;
    using Identifier = typename CoDiType::Identifier;
    using Tape = typename CoDiType::Tape;
    using Position = typename Tape::Position;

    using VAI = codi::VectorAccessInterface<Real, Identifier>;

  private:
    static void callHandleReverse(Tape*, void* handlePtr, VAI*) {
      opdi::Handle* handle = (opdi::Handle*)handlePtr;
      handle->reverseFunc(handle->data);
    }

    static void callHandleDelete(Tape*, void* handlePtr) {
      opdi::Handle* handle = (opdi::Handle*)handlePtr;
      if (handle->deleteFunc != nullptr) {
        handle->deleteFunc(handle->data);
      }
      delete handle;
    }

  public:
    void init() {}

    void finalize() {}

    void* createTape() {
      return (void*)new Tape;
    }

    void deleteTape(void* tapePtr) {
      Tape* tape = (Tape*)tapePtr;
      delete tape;
    }

    void* allocPosition() {
      return new Position();
    }

    void freePosition(void* positionPtr) {
      Position* position = (Position*)positionPtr;
      delete position;
    }

    size_t getPositionSize() {
      return sizeof(Position);
    }

    std::string positionToString(void* positionPtr) {
      Position* position = (Position*)positionPtr;
      std::stringstream conv;
      conv << *position;
      return conv.str();
    }

    void getTapePosition(void* tapePtr, void* positionPtr) {
      Tape* tape = (Tape*)tapePtr;
      Position* position = (Position*)positionPtr;

      *position = tape->getPosition();
    }

    void getZeroPosition(void* tapePtr, void* positionPtr) {
      Tape* tape = (Tape*)tapePtr;
      Position* position = (Position*)positionPtr;

      *position = tape->getZeroPosition();
    }

    void copyPosition(void* dstPtr, void* srcPtr) {
      Position* dst = (Position*)dstPtr;
      Position* src = (Position*)srcPtr;

      *dst = *src;
    }

    int comparePosition(void* lhsPtr, void* rhsPtr) {
      Position* lhs = (Position*)lhsPtr;
      Position* rhs = (Position*)rhsPtr;

      if (*lhs <= *rhs) {
        if (*rhs <= *lhs) {
          return 0;
        } else {
          return -1;
        }
      } else {
        return 1;
      }
    }

    bool isActive(void* tapePtr) {
      Tape* tape = (Tape*)tapePtr;
      return tape->isActive();
    }

    void setActive(void* tapePtr, bool active) {
      Tape* tape = (Tape*)tapePtr;
      if (active) {
        tape->setActive();
      } else {
        tape->setPassive();
      }
    }

    void evaluate(void* tapePtr, void* startPtr, void* endPtr, bool useAtomics = true) {
      Tape* tape = (Tape*)tapePtr;
      Position* start = (Position*)startPtr;
      Position* end = (Position*)endPtr;

      if (tape->isActive()) {
        std::cerr << "Warning: OpDiLib evaluation of an active tape." << std::endl;
      }

      typename Tape::Gradient* adjoints = &tape->gradient(0);
      using NonAtomicGradient = codi::AtomicTraits::RemoveAtomic<typename Tape::Gradient>;
      using AtomicGradient = codi::OpenMPAtomic<NonAtomicGradient>;

      if (useAtomics) {
        AtomicGradient* safeAdjoints = (AtomicGradient*)adjoints;
        tape->evaluate(*start, *end, safeAdjoints);
      } else {
        NonAtomicGradient* unsafeAdjoints = (NonAtomicGradient*)adjoints;
        tape->evaluate(*start, *end, unsafeAdjoints);
      }
    }

    void reset(void* tapePtr, bool clearAdjoints) {
      Tape* tape = (Tape*)tapePtr;
      tape->reset(clearAdjoints);
    }

    void reset(void* tapePtr, void* positionPtr, bool clearAdjoints) {
      Tape* tape = (Tape*)tapePtr;
      Position* position = (Position*)positionPtr;
      tape->resetTo(*position, clearAdjoints);
    }

    void* getThreadLocalTape() {
      return (void*)CoDiType::getTapePtr();
    }

    void setThreadLocalTape(void* tapePtr) {
      Tape* tape = (Tape*)tapePtr;
      CoDiType::setTapePtr(tape);
    }

    void pushExternalFunction(void* tapePtr, opdi::Handle const* handle) {
      Tape* tape = (Tape*)tapePtr;
      tape->pushExternalFunction(codi::ExternalFunction<Tape>::create(CoDiOpDiLibTool::callHandleReverse, (void*)handle,
                                                                      CoDiOpDiLibTool::callHandleDelete));
    }

    void erase(void* tapePtr, void* startPtr, void* endPtr) {
      Tape* tape = (Tape*)tapePtr;
      Position* start = (Position*)startPtr;
      Position* end = (Position*)endPtr;

      tape->erase(*start, *end);
    }

    void append(void* dstTapePtr, void* srcTapePtr, void* startPtr, void* endPtr) {
      Tape* dstTape = (Tape*)dstTapePtr;
      Tape* srcTape = (Tape*)srcTapePtr;
      Position* start = (Position*)startPtr;
      Position* end = (Position*)endPtr;

      dstTape->append(*srcTape, *start, *end);
    }
};

#endif
