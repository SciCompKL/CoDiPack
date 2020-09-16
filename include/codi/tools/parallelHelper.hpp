/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include <vector>
#include <map>
#include <stack>
#include <thread>
#include <mutex>
#include <utility>
#include <set>
#include <queue>
#include <list>
#include <atomic>
#include <deque>
#include <string>
#include <fstream>

#include "../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Provides tools for handling multiple inter-dependent tapes.
   * @tparam Tape The type of tape used.
   *
   * In multithreaded applications, multiple tapes must be managed, and their evaluation
   * must respect the synchronization that likely occured during the forward pass. This
   * helper provides
   * - routines for creating and managing multiple tapes,
   * - broadcast calls to multiple tapes,
   * - management of thread-local tapes,
   * - recording of meta information that reflects the tapes' dependencies,
   * - helpers for custom reverse evaluation,
   * - a built-in evaluation routine as a convenient option for an automatically parallelized reverse pass.
   *
   * The public member functions can safely be used in a multithreaded application.
   *
   * With respect to tape management, care should be taken about correct pairing (register - forget, create - delete).
   * Also, subsequent calls of those functions with the same tape (e.g. register the same tape twice) can cause unexpected behaviour and should be avoided.
   *
   * Such misuse is (amongst others) detected by setting CODI_EnableAssert.
   *
   * Setting CODI_EnableParallelHelperDebugOutput provides additional information about tapes.
   * - 1: Report beginnings and ends of frames.
   * - 2: Report setting and clearing of threadlocal tapes.
   * - 4: Report beginnings and ends of frame evaluations.
   * - 8: Scheduler reports on submitted and finished frames.
   *
   * Multiple options can be enabled simultaneously by forming sums.
   */
  template<typename Tape>
  class ParallelHelper {
    public:

      /**
       * @brief Type used for tape identification.
       */
      typedef int TapeId;

      /**
       * @brief Invalid tape id.
       */
      static const TapeId InvalidTapeId = -1;

    private:
      /**
       * @brief Type used to order synchronization events.
       */
      typedef unsigned long SyncEvent;

      /**
       * @brief Convenience typedef. The tape's position type.
       */
      typedef typename Tape::Position Position;

      /**
       * @brief Mutex that supports read and write locking.
       *
       * Since there is no shared mutex support in C++11, this is a custom
       * mutex that uses two atomic ints to mimick the behaviour.
       *
       * Tape data is stored in a map which is not thread-safe. Access is synchronized by means of this read-write mutex.
       */
      struct SharedMutex {

        /**
         * @brief Indicates lock for write.
         */
        std::atomic<int> hasWriter;

        /**
         * @brief Counts locks for read.
         */
        std::atomic<int> numReaders;

        /**
         * @brief Constructor.
         */
        SharedMutex() : hasWriter(false), numReaders(0) {
        }

        /**
         * @brief Lock for writing.
         *
         * Busy-waits to acquire hasWriter, then busy-waits until there are no readers.
         */
        CODI_INLINE void lockWrite() {
          while (1 == this->hasWriter.fetch_or(1)) {} // wait until hasWriter is false, then set it to true
          while (this->numReaders.load() != 0) {} // wait until there are no readers
        }

        /**
         * @brief Lock for reading.
         *
         * After a successful lock, the readers counter remains incremented.
         */
        CODI_INLINE void lockRead() {
          while (true) {
            while (this->hasWriter.load()) {}; // wait until there is no writer
            this->numReaders++;
            if (this->hasWriter.load()) { // check if there is still no writer
              this->numReaders--; // otherwise delay reading
            }
            else {
               break;
            }
          }
        }

        /**
         * @brief Unlock from writing.
         *
         * Releases hasWriter.
         */
        CODI_INLINE void unlockWrite() {
          int previous = this->hasWriter.fetch_and(0); // set has Writer to false
          codiAssert(previous); // assert if it was false already
          CODI_UNUSED(previous);
        }

        /**
         * @brief Unlock from reading.
         *
         * Decrements the readers counter.
         */
        CODI_INLINE void unlockRead() {
          this->numReaders--;
        }
      };

      /**
       * @brief RAII lock for read.
       */
      struct ReadLock {
        SharedMutex& mutex;

        /**
         * @brief Constructor.
         * @param mutex The mutex to lock.
         *
         * Locks the mutex for reading.
         */
        ReadLock(SharedMutex& mutex) : mutex(mutex) {
          this->mutex.lockRead();
        }

        /**
         * @brief Destructor.
         *
         * Unlocks the mutex.
         */
        ~ReadLock() {
          this->mutex.unlockRead();
        }
      };

      /**
       * @brief RAII lock for write.
       */
      struct WriteLock {
        SharedMutex& mutex;

        /**
         * @brief Constructor.
         * @param mutex The mutex to lock.
         *
         * Locks the mutex for writing.
         */
        WriteLock(SharedMutex& mutex) : mutex(mutex) {
          this->mutex.lockWrite();
        }

        /**
         * @brief Destructor.
         *
         * Unlocks the mutex.
         */
        ~WriteLock() {
          this->mutex.unlockWrite();
        }
      };

      /**
       * @brief Stores a tape pointer together with meta information.
       */
      struct TapeData {
        /**
         * @brief Pointer to the tape.
         */
        Tape* tape;

        /**
         * @brief Application wide unique tape id.
         */
        TapeId tapeId;

        /**
         * @brief User-defined name of the tape.
         *
         * Makes debugging output more readable.
         */
        std::string name;

        /**
         * @brief Indicates ownership of tape pointer.
         *
         * As an example, the master tape might be externally managed wheres tapes of local workers are not.
         */
        bool externallyManaged;

        /**
         * @brief Meta information required for joint reverse evaluation of multiple tapes.
         *
         * The general idea is that a frame marks a part of the tape that arises from a continuous,
         * unsynchronized sequence of computations that does not depend on activity on other tapes.
         */
        struct Frame {
          /**
           * @brief Starting position of the frame.
           */
          Position start;

          /**
           * @brief Past-the-end position of the frame.
           */
          Position end;

          /**
           * @brief Time stamp that marks the beginning of the forward evaluation.
           */
          SyncEvent startEvent;

          /**
           * @brief Time stamp that marks the end of the forward evaluation.
           */
          SyncEvent endEvent;

          /**
           * @brief Constructor for the internally relevant use case.
           * @param start Start position of the frame.
           * @param startEvent Sync event at start of the frame.
           */
          Frame(Position start, SyncEvent startEvent) : start(start), end(), startEvent(startEvent), endEvent() {
          }
        };

        /**
         * @brief Subdivide the tape into multiple frames.
         *
         * New frames are inserted at the beginning so that the past-the-end iterator
         * corresponds to "all frames evaluated".
         */
        std::deque<Frame> frames;

        /**
         * @brief Frame iterator.
         *
         * Used to indicate the next frame to evaluate.
         * Evaluated frames are not popped so that the tape collection can be evaluated multiple times.
         */
        typename std::deque<Frame>::iterator frameIterator;

        /**
         * @brief Constructor.
         * @param tape Pointer to the tape.
         * @param tapeId Id of the tape.
         * @param externallyManaged Memory management flag.
         * @param name User-defined name of the tape.
         */
        TapeData(Tape* tape, TapeId tapeId, bool externallyManaged, const std::string& name) : tape(tape), tapeId(tapeId), name(name), externallyManaged(externallyManaged) {
        }

        /**
         * @brief Default constructor.
         */
        TapeData() : tape(nullptr), tapeId(InvalidTapeId), name(), externallyManaged(false) {
        }
      };

      /**
       * @brief Stores multiple tapes and corresponding data.
       */
      std::map<TapeId, TapeData> tapeData;

      /**
       * @brief Used to create application wide unique tape ids.
       */
      static std::atomic<TapeId> nextTapeId;

      /**
       * @brief Protects the tape data map.
       */
      SharedMutex tapeDataMutex;

      /**
       * @brief Thread safe generation of sync events.
       */
      std::atomic<SyncEvent> nextEvent;

      /**
       * @brief Stores the default tape of a thread when a specialized tape is set for the tape.
       *
       * In general this variable will be null if the default thread for the tape is set to the global tape pointer
       * of CoDiPack.
       */
      static thread_local Tape* threadDefaultTape;

      /*************** members related to debug output ****************/

      #if CODI_EnableParallelHelperDebugOutput
        /**
         * @brief Used to create unique thread ids.
         */
        static std::atomic<int> nextThreadId;

        /**
         * @brief Thread id.
         *
         * Used to distinguish threads in debug output.
         */
        static thread_local int threadId;

        /**
         * @brief Protects stdout.
         *
         * Used during debugging to synchronize output.
         */
        static std::mutex outputMutex;
      #endif

      /*************** internal helpers ****************/

      /*
       * Note that synchronization of map access takes place in the public member functions.
       * Likewise, no error checks are performed at this level.
       */

      /**
       * @brief Adds a new tape.
       * @param tape The new tape to register.
       * @param externallyManaged Ownership of the tape pointer.
       * @return The id assigned to the tape.
       */
      CODI_INLINE TapeId internalRegisterTape(Tape* tape, bool externallyManaged, const std::string& name) {
        TapeId id = nextTapeId++;
        this->tapeData[id] = TapeData(tape, id, externallyManaged, name);
        return id;
      }

      /**
       * @brief Forgets an externally managed or already deleted tape.
       * @param id Id of the tape.
       */
      CODI_INLINE void internalForgetTape(TapeId id) {
        this->tapeData.erase(id);
      }

      /**
       * @brief Delete a tape, then forget it.
       * @param id Id of the tape.
       */
      CODI_INLINE void internalDeleteTape(TapeId id) {
        delete this->tapeData.at(id).tape;
        this->internalForgetTape(id);
      }

      /**
       * @brief Delete or forget a tape, depending on ownership.
       * @param id Id of the tape.
       */
      CODI_INLINE void internalClearTape(TapeId id) {
        if (this->tapeData.at(id).externallyManaged) {
          this->internalForgetTape(id);
        }
        else {
          this->internalDeleteTape(id);
        }
      }

      /**
       * @brief Minimal validity check.
       * @param id Id of the tape to check.
       * @return Whether the last recorded frame of the given tape has admissible boundaries.
       *
       * Valid means: no frames so far, or: start position <= end position and start event <= end event.
       */
      CODI_INLINE bool validFrame(TapeId id) {
        return this->tapeData.at(id).frames.empty() || ((this->tapeData.at(id).frames.front().start <= this->tapeData.at(id).frames.front().end) && (this->tapeData.at(id).frames.front().startEvent <= this->tapeData.at(id).frames.front().endEvent));
      }

      /**
       * @brief Convenience accessor.
       * @param id Id of the tape.
       */
      CODI_INLINE bool isExternallyManaged(TapeId id) {
        return this->tapeData.at(id).externallyManaged;
      }

      /*************** internal helpers related to debug output ****************/

      #if CODI_EnableParallelHelperDebugOutput
        /**
         * @brief Terminator.
         */
        CODI_INLINE void dataOutput() {}

        /**
         * @brief Output of arbitrary data.
         * @tparam Head Type of the next datum.
         * @tparam Tail Types of remaining data.
         * @param head Next datum to print.
         * @param tail Remaining data.
         */
        template<typename Head, typename... Tail>
        CODI_INLINE void dataOutput(const Head& head, const Tail&... tail) {
          std::cout << head << " ";
          dataOutput(tail...);
        }

        /**
         * @brief Synchronized debug output.
         * @param data Data to print.
         *
         * Synchronizes output via outputMutex. Also appends newline at the end.
         */
        template<typename... Data>
        CODI_INLINE void debugOutput(const Data&... data) {

          std::lock_guard<std::mutex> lock(outputMutex);
          this->dataOutput(data...);
          std::cout << std::endl;
        }

      #endif

    public:
      /**
       * @brief Constructor.
       */
      ParallelHelper() : nextEvent(0) {
      }

      /**
       * @brief Destructor.
       */
      ~ParallelHelper() {
      }

      /*************** parallel helper management ****************/

      /**
       * @brief Initialize the parallel helper.
       *
       * Nothing to do right now.
       */
      CODI_INLINE void init() {}

      /**
       * @brief Clear all tapes known to the parallel helper.
       *
       * Forget all externally managed tapes and delete the others.
       */
      CODI_INLINE void clear() {
        WriteLock lock(this->tapeDataMutex);
        while (!this->tapeData.empty()) {
          this->internalClearTape(this->tapeData.begin()->first);
        }
      }

      /*************** tape management ****************/

      /**
       * @brief Register a tape that is memory managed elsewhere.
       * @return The id assigned to the tape.
       *
       * The counterpart is forgetTape.
       */
      CODI_INLINE TapeId registerTape(Tape* tape, std::string name = "") {
        codiAssert(!hasTape(tape));
        WriteLock lock(this->tapeDataMutex);
        return this->internalRegisterTape(tape, true, name);
      }

      /**
       * @brief Register the thread-local tape of the calling thread.
       * @return The id assigned to the tape.
       *
       * The counterpart is forgetTape.
       */
      CODI_INLINE TapeId registerTape(std::string name = "") {
        ActiveReal<Tape> dummy;
        Tape* tape = dummy.getGlobalTapePtr();
        codiAssert(tape != nullptr);
        return this->registerTape(tape, name);
      }

      /**
       * @brief Forget an externally managed tape.
       * @param id The id of the tape.
       *
       * The counterpart is registerTape.
       */
      CODI_INLINE void forgetTape(TapeId id) {
        codiAssert(this->hasTape(id));
        codiAssert(this->isExternallyManaged(id));

        WriteLock lock(this->tapeDataMutex);
        this->internalForgetTape(id);
      }

      /**
       * @brief Forget the thread-local tape of the calling thread.
       *
       * The counterpart is registerTape.
       */
      CODI_INLINE void forgetTape() {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        this->forgetTape(this->getTapeId(tape));
      }

      /**
       * @brief Create a new tape that is managed by the helper.
       * @return The id of the tape.
       *
       * The counterpart is deleteTape.
       */
      CODI_INLINE TapeId createTape(std::string name = "") {
        WriteLock lock(this->tapeDataMutex);
        Tape* tape = new Tape;
        return this->internalRegisterTape(tape, false, name);
      }

      /**
       * @brief Delete a tape that is managed by the helper.
       * @param id Id of the tape.
       *
       * The counterpart is createTape.
       */
      CODI_INLINE void deleteTape(TapeId id) {
        codiAssert(this->hasTape(id));
        codiAssert(!this->isExternallyManaged(id));

        WriteLock lock(this->tapeDataMutex);
        this->internalDeleteTape(id);
      }

      /**
       * @brief Check if a tape id is known to the parallel helper.
       * @param id The id of the tape.
       * @return True if there is such tape.
       */
      CODI_INLINE bool hasTape(TapeId id) {
        ReadLock lock(this->tapeDataMutex);
        return this->tapeData.find(id) != this->tapeData.end();
      }

      /**
       * @brief Check if a tape is known to the parallel helper.
       * @param id Address of the tape.
       * @return True if there is such tape.
       */
      CODI_INLINE bool hasTape(Tape* tape) {
        ReadLock lock(this->tapeDataMutex);
        for (auto& tapeDataPair : this->tapeData) {
          if (tapeDataPair.second.tape == tape) {
            return true;
          }
        }
        return false;
      }

      /**
       * @brief Convert a tape id into the corresponding tape pointer.
       * @param id Tape id.
       * @return Pointer to the associated tape.
       */
      CODI_INLINE Tape* getTape(TapeId id) {
        codiAssert(this->hasTape(id));
        ReadLock lock(this->tapeDataMutex);
        return this->tapeData.at(id).tape;
      }

      /**
       * @brief Convert a tape pointer into the corresponding tape id.
       * @param tape A tape pointer.
       * @return The corresponding tape id, InvalidTapeId if no such tape is known.
       */
      CODI_INLINE TapeId getTapeId(Tape* tape) {
        ReadLock lock(this->tapeDataMutex);
        for (auto& tapeDataPair : this->tapeData) {
          if (tapeDataPair.second.tape == tape) {
            return tapeDataPair.first;
          }
        }
        codiAssert(false && "tape not known to helper");
        return InvalidTapeId;
      }

      /**
       * @brief Tape name accessor.
       * @param id Tape id.
       * @return Name of the tape with the given id.
       */
      std::string getTapeName(TapeId id) {
        codiAssert(this->hasTape(id));
        ReadLock lock(this->tapeDataMutex);
        return this->tapeData.at(id).name;
      }

      /**
       * @brief Access tape name of the thread-local tape of the calling thread.
       * @return Name of the tape.
       *
       * Assumes that this tape is known to the parallel helper.
       */
      std::string getTapeName() {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        return this->getTapeName(this->getTapeId(tape));
      }

      /*************** broadcasts to all tapes ****************/

      /**
       * @brief Sets all tapes active.
       */
      CODI_INLINE void setActive() {
        ReadLock lock(this->tapeDataMutex);
        for (auto& tapeDataPair : this->tapeData) {
          tapeDataPair.second.tape->setActive();
        }
      }

      /**
       * @brief Sets all tapes passive.
       */
      CODI_INLINE void setPassive() {
        ReadLock lock(this->tapeDataMutex);
        for (auto& tapeDataPair : this->tapeData) {
          tapeDataPair.second.tape->setPassive();
        }
      }

      /**
       * @brief Reset all tapes.
       *
       * The shared adjoint vector is only reset once. Also, frame information is cleared.
       * No tapes are forgotten or deleted.
       */
      CODI_INLINE void reset() {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(!this->tapeData.empty());

        this->tapeData.begin()->second.tape->clearAdjoints(); // all tapes share the same adjoint vector, it suffices if it is cleared once...
        for (auto& tapeDataPair : this->tapeData) {
          tapeDataPair.second.tape->reset(false); // ...as indicated by the "false" here.
          tapeDataPair.second.frames.clear();
        }
      }

      /**
       * @brief Prints statistics of all tapes, together with tape id, name and address.
       */
      template<typename Stream = std::ostream>
      CODI_INLINE void printStatistics(Stream& out = std::cout) const {
        ReadLock lock(this->tapeDataMutex);
        for (auto& tapeDataPair : this->tapeData) {
          out << "------------- Statistics of tape " << tapeDataPair.first << " (" << tapeDataPair.second.name << ", " << tapeDataPair.second.tape << ") -------------" << std::endl;
          tapeDataPair.second.tape->printStatistics(out);
        }
      }

      /**
       * @brief Tabular representation of tape statistics including tape id, name and address.
       */
      template<typename Stream = std::ostream>
      CODI_INLINE void printTable(Stream& out = std::cout) {
        ReadLock lock(this->tapeDataMutex);

        // generate table header
        if (!this->tapeData.empty()) {
          out << "id; name; address; ";
          this->tapeData.begin()->second.tape->printTableHeader(out);
        }

        // generate table rows
        for (auto& tapeDataPair : this->tapeData) {
          out << tapeDataPair.first << "; " << tapeDataPair.second.name << "; " << tapeDataPair.second.tape << "; ";
          tapeDataPair.second.tape->printTableRow(out);
        }
      }

      /*************** tape meta information management ****************/

      /**
       * @brief Begins a frame on the tape with the given id.
       * @param id Tape id.
       *
       * CODI_EnableAssert detects:
       * - incomplete previous frame, that is, beginFrame has not been closed by endFrame
       * - gap between this frame and the previous one, that is, the tape has parts that are not contained in any frame
       *
       * The counterpart is endFrame.
       */
      CODI_INLINE void beginFrame(TapeId id) {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(this->hasTape(id));

        SyncEvent event = this->nextEvent++;

        #if CODI_EnableParallelHelperDebugOutput & 1
          std::string spaces;
          for (int i = 0; i < id; ++i)
            spaces.push_back(' ');
          debugOutput(spaces, this->tapeData.at(id).name, "begins frame at", this->tapeData.at(id).tape->getPosition(), event);
        #endif

        codiAssert(this->validFrame(id)); // incomplete frame detection
        codiAssert(this->tapeData.at(id).frames.empty() || this->tapeData.at(id).frames.front().end == this->tapeData.at(id).tape->getPosition()); // gap detection

        this->tapeData.at(id).frames.push_front(typename codi::ParallelHelper<Tape>::TapeData::Frame(this->tapeData.at(id).tape->getPosition(), event));
      }

      /**
       * @brief Begins a frame on the thread-local tape of the calling thread.
       *
       * Assumes that this tape is known to the parallel helper.
       *
       * The counterpart is endFrame.
       */
      CODI_INLINE void beginFrame() {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        this->beginFrame(this->getTapeId(tape));
      }

      /**
       * @brief Ends a frame on the tape with the given id.
       * @param id Tape id.
       * @param discardIfEmpty Control whether empty frames are kept or discarded.
       * @return true if the frame was discarded, false otherwise.
       *
       * Automatically skips empty frames, that is, start == end.
       *
       * CODI_EnableAssert detects:
       * - most recent frame already closed
       *
       * The counterpart is beginFrame.
       */
      CODI_INLINE bool endFrame(TapeId id, bool discardIfEmpty = true) {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(this->hasTape(id));
        codiAssert(this->tapeData.at(id).frames.front().end == Position()); // frame already closed

        SyncEvent event = this->nextEvent++;

        #if CODI_EnableParallelHelperDebugOutput & 1
          std::string spaces;
          for (int i = 0; i < id; ++i)
            spaces.push_back(' ');
          debugOutput(spaces, this->tapeData.at(id).name, "ends frame at", this->tapeData.at(id).tape->getPosition(), event);
        #endif

        Position end = this->tapeData.at(id).tape->getPosition();

        if (this->tapeData.at(id).frames.front().start == end && discardIfEmpty) { // automatically skip empty frames
          this->tapeData.at(id).frames.pop_front();
          return true;
        }
        else {
          this->tapeData.at(id).frames.front().end = end;
          this->tapeData.at(id).frames.front().endEvent = event;
          return false;
        }
      }

      /**
       * @brief Ends a frame on the thread-local tape of the calling thread.
       * @param discardIfEmpty Control whether empty frames are kept or discarded.
       * @return true if the frame was discarded, false otherwise.
       *
       * Assumes that this tape is known to the parallel helper.
       *
       * The counterpart is beginFrame.
       */
      CODI_INLINE bool endFrame(bool discardIfEmpty = true) {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        return this->endFrame(this->getTapeId(tape), discardIfEmpty);
      }

      /**
       * @brief Check if the most recent frame on the given tape is valid.
       * @param id Tape id.
       * @return Result of the validFrame check.
       *
       * Checks in particular if a beginFrame has been close by an endFrame.
       */
      CODI_INLINE bool lastFrameValid(TapeId id) {
        codiAssert(this->hasTape(id));
        ReadLock lock(this->tapeDataMutex);
        return this->validFrame(id);
      }

      /**
       * @brief Checks if the most recent frame on the thread-local tape of the calling thread is valid.
       * @return Result of the validFrame check.
       *
       * Assumes that this tape is known to the parallel helper.
       */
      CODI_INLINE bool lastFrameValid() {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        return this->lastFrameValid(this->getTapeId(tape));
      }

      /**
       * @brief Export the timeline indicated by the frame boundaries.
       * @param filename Name of the file to write to.
       *
       * Creates or truncates the specified file and exports the timeline. Each line contains
       *
       * tapeName startEvent endEvent startEvent endEvent ...
       */
      CODI_INLINE void exportTimeline(const std::string& filename) {
        ReadLock lock(this->tapeDataMutex);

        std::ofstream out(filename);
        codiAssert(out.is_open());

        for (auto& tapeDataPair : this->tapeData) {

          if (tapeDataPair.second.name.empty()) {
            out << "tape_" << tapeDataPair.first;
          }
          else {
            out << tapeDataPair.second.name;
          }

          for (auto& frame : tapeDataPair.second.frames) {
            out << " " << frame.startEvent << " " << frame.endEvent;
          }

          out << std::endl;
        }

        out.close();
      }

      /*************** management of thread-local tapes ****************/

      /**
       * @brief Set the thread-local tape of the calling thread to the one with the given id.
       *
       * The current tape of the thread is stored. A call to clearThisThreadsTape will restore the current tape.
       * @param id Tape id.
       */
      CODI_INLINE void setThisThreadsTape(TapeId id) {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(this->hasTape(id));
        codiAssert(nullptr == ParallelHelper::threadDefaultTape); // If not null, then clearThisThreadsTape has not been called.

        #if CODI_EnableParallelHelperDebugOutput & 2
          debugOutput("thread", threadId, "uses now tape", id, this->tapeData.at(id).name);
        #endif

        ParallelHelper::threadDefaultTape = ActiveReal<Tape>::getGlobalTapePtr();
        ActiveReal<Tape>::setGlobalTapePtr(this->tapeData.at(id).tape);
      }

      /**
       * @brief Set the thread-local tape of the calling thread to the old stored tape.
       */
      CODI_INLINE void clearThisThreadsTape() {
        codiAssert(nullptr != ParallelHelper::threadDefaultTape); // If null, then setThisThreadsTape has not been called.

        #if CODI_EnableParallelHelperDebugOutput & 2
          debugOutput("thread", threadId, "cleared its thread-local tape");
        #endif

        ActiveReal<Tape>::setGlobalTapePtr(ParallelHelper::threadDefaultTape);
        ParallelHelper::threadDefaultTape = nullptr;
      }

      /*************** reverse pass management ****************/

      /**
       * @brief Initialize the frame iterators.
       *
       * Must be called prior to custom reverse evaluation.
       * Internally, sets all frame iterators to the respective most recent frame.
       */
      CODI_INLINE void prepareEvaluation() {
        ReadLock lock(this->tapeDataMutex);
        for (auto& tapeDataPair : this->tapeData) {
          tapeDataPair.second.frameIterator = tapeDataPair.second.frames.begin();
        }
      }

      /**
       * @brief Returns the number of frames recorded for the tape with the given id.
       * @param id Tape id.
       * @return Number of frames.
       */
      CODI_INLINE size_t totalNumberOfFrames(TapeId id) {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(this->hasTape(id));
        return this->tapeData.at(id).frames.size();
      }

      /**
       * @brief For the given tape, returns the number of frames left to evaluate.
       * @param id Tape id.
       * @return Number of frames left to evaluate.
       */
      CODI_INLINE size_t numberOfFramesLeftToEvaluate(TapeId id) {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(this->hasTape(id));
        return std::distance(this->tapeData.at(id).frameIterator, this->tapeData.at(id).frames.end());
      }

      /**
       * @brief For the thread local tape of the calling thread, returns the number of frames left to evaluate.
       * @return Number of frames left to evaluate.
       */
      CODI_INLINE size_t numberOfFramesLeftToEvaluate() {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        return this->numberOfFramesLeftToEvaluate(this->getTapeId(tape));
      }

      /**
       * @brief Checks if all frames of the given tape have been evaluated.
       * @param id Tape id.
       * @return True if all frames are evaluated.
       */
      CODI_INLINE bool evaluationDone(TapeId id) {
        return this->numberOfFramesLeftToEvaluate(id) == 0;
      }

      /**
       * @brief Checks if all frames of the thread local tape of the calling thread have been evaluated.
       * @return True if all frames are evaluated.
       */
      CODI_INLINE bool evaluationDone() {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        return this->evaluationDone(this->getTapeId(tape));
      }

      /**
       * @brief Evaluates the next frame of the given tape with a custom adjoint vector
       * @param id Tape id.
       *
       * Internally, advances the frame iterator.
       */
      template<typename AdjointVec>
      CODI_INLINE void evaluateNextFrame(TapeId id, AdjointVec& adjoints) {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(this->hasTape(id));
        codiAssert(!this->evaluationDone(id));

        #if CODI_EnableParallelHelperDebugOutput & 4
          std::string spaces;
          for (int i = 0; i < id; ++i)
            spaces.push_back(' ');
          debugOutput(spaces, "thread", threadId, "evaluates", this->tapeData.at(id).name, "from", this->tapeData.at(id).frameIterator->end, "to", this->tapeData.at(id).frameIterator->start, "from", this->tapeData.at(id).frameIterator->endEvent, "to", this->tapeData.at(id).frameIterator->startEvent);
        #endif

        auto oldIter = this->tapeData.at(id).frameIterator++;
        this->tapeData.at(id).tape->evaluate(oldIter->end, oldIter->start, adjoints);

        #if CODI_EnableParallelHelperDebugOutput & 4
          debugOutput(spaces, "thread", threadId, "finished evaluating", this->tapeData.at(id).name, "from", oldIter->end, "to", oldIter->start, "from", oldIter->endEvent, "to", oldIter->startEvent);
        #endif
      }

      /**
       * @brief Evaluates the next frame of the given tape.
       * @param id Tape id.
       *
       * Internally, advances the frame iterator.
       */
      CODI_INLINE void evaluateNextFrame(TapeId id) {
        ReadLock lock(this->tapeDataMutex);
        codiAssert(this->hasTape(id));
        codiAssert(!this->evaluationDone(id));

        #if CODI_EnableParallelHelperDebugOutput & 4
          std::string spaces;
          for (int i = 0; i < id; ++i)
            spaces.push_back(' ');
          debugOutput(spaces, "thread", threadId, "evaluates", this->tapeData.at(id).name, "from", this->tapeData.at(id).frameIterator->end, "to", this->tapeData.at(id).frameIterator->start, "from", this->tapeData.at(id).frameIterator->endEvent, "to", this->tapeData.at(id).frameIterator->startEvent);
        #endif

        auto oldIter = this->tapeData.at(id).frameIterator++;
        this->tapeData.at(id).tape->evaluate(oldIter->end, oldIter->start);

        #if CODI_EnableParallelHelperDebugOutput & 4
          debugOutput(spaces, "thread", threadId, "finished evaluating", this->tapeData.at(id).name, "from", oldIter->end, "to", oldIter->start, "from", oldIter->endEvent, "to", oldIter->startEvent);
        #endif
      }

      /**
       * @brief Evaluates the next frame of the thread local tape of the calling thread.
       *
       * Retrieves the tape id and forwards the call to evaluateNextFrame(TapeId).
       */
      CODI_INLINE void evaluateNextFrame() {
        Tape* tape = ActiveReal<Tape>::getGlobalTapePtr();
        this->evaluateNextFrame(this->getTapeId(tape));
      }

      /**
       * @brief Built-in parallel reverse evaluation routine.
       * @param numThreads Number of worker threads spawned.
       *
       * Takes care of the call to prepareEvaluation.
       *
       * An admissible reverse evaluation schedule is recovered from the recorded frame information, in particular the timings.
       *
       * Spawns numThreads worker threads. The calling thread servers as scheduler thread and hands out frames for evaluation to the workers.
       *
       * It might not be as efficient as a custom reverse pass implementation since the user has most likely additional frame
       * dependency information available. Nonetheless, it can be used as a convenient fallback, or to verify a custom implementation.
       */
      CODI_INLINE void evaluate(int numThreads) {
        ReadLock lock(this->tapeDataMutex);
        this->prepareEvaluation();

        // The environment is evaluated by a local mini scheduler.

        // scheduler helpers
        std::map<TapeId, bool> queueBlocked; // marks whether a frame in this tape is currently being evaluated
        for (auto& tapeDataPair : this->tapeData) {
          queueBlocked[tapeDataPair.first] = false;
        }

        std::set<SyncEvent> runningStartEvents; // contains the start events of all frames that are currently being evaluated

        // queue of launched frames, ready for pick-up by threads
        std::queue<std::pair<TapeId, SyncEvent> > reverseEvalQueue;
        std::atomic_flag queueFlag = ATOMIC_FLAG_INIT;

        // queue of finished events, ready for pick-up by scheduler
        std::queue<std::pair<TapeId, SyncEvent> > finished;
        std::atomic_flag finishedFlag = ATOMIC_FLAG_INIT;

        // stop signal for threads
        std::atomic<bool> shouldStop(false);

        // thread body
        auto threadFun = [&](int threadId) {
          (void)threadId;

          // thread main loop
          while (true) {
            // find something to do
            while (queueFlag.test_and_set()) {}

            if (!reverseEvalQueue.empty()) {
              std::pair<TapeId, SyncEvent> job = reverseEvalQueue.front();
              reverseEvalQueue.pop();
              queueFlag.clear();

              this->evaluateNextFrame(job.first);

              while (finishedFlag.test_and_set()) {}
              finished.push(job);
              finishedFlag.clear();
            }
            else {
              queueFlag.clear();

              // guarantee that queue is empty before stop, therefore check only in else branch
              if (shouldStop.load()) {
                break;
              }
            }
          }
        };

        // create threads
        std::list<std::thread> threads;
        for (int i = 0; i < numThreads; ++i) {
          threads.emplace(threads.end(), std::thread(threadFun, i));
        }

        bool schedulerStop = false;

        // scheduler main loop
        while (!schedulerStop) {
          // identify frame with latest end event among unblocked queues as candidate
          TapeId candidate;
          SyncEvent endEvent;
          bool foundJob = false;

          schedulerStop = true;

          for (auto& tapeDataPair : this->tapeData) {
            if (this->evaluationDone(tapeDataPair.first)) {
              continue;
            }
            else {
              schedulerStop = false;
              if (queueBlocked[tapeDataPair.first]) {
                continue;
              }
            }
            if (!foundJob || tapeDataPair.second.frameIterator->endEvent > endEvent) { // first job found or additional job found with later end event
              foundJob = true;
              candidate = tapeDataPair.first;
              endEvent = tapeDataPair.second.frameIterator->endEvent;
            }
          }

          // launch it if its end event is ordered after the latest start event among running jobs
          if (foundJob && (runningStartEvents.empty() || endEvent > *(--runningStartEvents.end()))) {
            auto startEvent = this->tapeData[candidate].frameIterator->startEvent;

            while (queueFlag.test_and_set()) {}
            reverseEvalQueue.push(std::make_pair(candidate, startEvent));
            queueFlag.clear();

            #if CODI_EnableParallelHelperDebugOutput & 8
              debugOutput("scheduler submitted frame on tape", candidate, this->tapeData.at(candidate).name, "with events", startEvent, endEvent);
            #endif

            queueBlocked[candidate] = true;
            runningStartEvents.insert(startEvent);
          }

          // check for finished jobs
          while (finishedFlag.test_and_set()) {}
          if (!finished.empty()) {
            std::pair<TapeId, SyncEvent> job = finished.front();
            finished.pop();
            finishedFlag.clear();

            queueBlocked[job.first] = false;

            runningStartEvents.erase(runningStartEvents.find(job.second));

            #if CODI_EnableParallelHelperDebugOutput & 8
              debugOutput("evaluation on tape", job.first, this->tapeData.at(job.first).name, "of a frame with start event", job.second, "done");
            #endif
          }
          else {
            finishedFlag.clear();
          }
        }

        // stop all threads
        shouldStop.store(true);

        for (auto threadsIter = threads.begin(); threadsIter != threads.end(); ++threadsIter) {
          threadsIter->join();
        }
      }
  };

  /*************** static member initialization ****************/

  template<typename Tape>
  std::atomic<typename ParallelHelper<Tape>::TapeId> ParallelHelper<Tape>::nextTapeId(0);

  template<typename Tape>
  thread_local Tape* ParallelHelper<Tape>::threadDefaultTape = nullptr;

  #if CODI_EnableParallelHelperDebugOutput
    template<typename Tape>
    std::atomic_int ParallelHelper<Tape>::nextThreadId(0);

    template<typename Tape>
    thread_local int ParallelHelper<Tape>::threadId = ParallelHelper<Tape>::nextThreadId++;

    template<typename Tape>
    std::mutex ParallelHelper<Tape>::outputMutex;
  #endif

}
