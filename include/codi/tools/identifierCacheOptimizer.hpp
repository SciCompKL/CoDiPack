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

#include <fstream>
#include <iostream>
#include <vector>

#include "../config.h"
#include "../misc/macros.hpp"
#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "../tapes/statementEvaluators/statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper class for iterating or changing the identifiers of a tape.
   *
   * This class implements the CallbacksInterface and can be used in a custom tape evaluation.
   *
   * Implementations may overwrite these methods:
   *  - applyToInput
   *  - applyPostInputLogic
   *  - applyToOutput
   *  - applyPostOutputLogic
   *
   *  They are called in the same order as listed above. See the method documentation for further information.
   *
   *  @tparam T_Tape  Tape tape on which the modification is applied.
   *  @tparam T_Impl  Final implementation of the modification.
   */
  template<typename T_Tape, typename T_Impl>
  struct ApplyIdentifierModification : public CallbacksInterface<typename T_Tape::Real, typename T_Tape::Identifier> {
      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_TAPE);            ///< See ApplyIdentifierModification.
      using Impl = CODI_DD(T_Impl, ApplyIdentifierModification);  ///< See ApplyIdentifierModification.

      using Real = typename Tape::Real;              ///< See FullTapeInterface.
      using Identifier = typename Tape::Identifier;  ///< See FullTapeInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See FullTapeInterface.

    private:

      Tape& tape;
      Real* primals = nullptr;

    public:

      /// Constructor.
      ApplyIdentifierModification(Tape& tape) : tape(tape) {
        if constexpr (codi::TapeTraits::isPrimalValueTape<Tape>) {
          primals = tape.getPrimalVector();
        }
      }

      /// Called for each input of each statement or low level function.
      CODI_INLINE void applyToInput(Identifier& id) {
        CODI_UNUSED(id);
        // Empty
      }

      /// Called for each output of each statement or low level function.
      CODI_INLINE void applyToOutput(Identifier& id) {
        CODI_UNUSED(id);
        // Empty
      }

      /// Called after \c applyToInput has been called for all inputs.
      CODI_INLINE void applyPostInputLogic() {
        // Empty
      }

      /// Called after \c applyToOutput has been called for all outputs.
      CODI_INLINE void applyPostOutputLogic() {
        // Empty
      }

    private:

      /// Helper function for low level function inputs.
      static void applyToInput_func(Identifier* id, ApplyIdentifierModification* data) {
        data->cast().applyToInput(*id);
      }

      /// Helper function for low level function outputs.
      static void applyToOutput_func(Identifier* id, ApplyIdentifierModification* data) {
        data->cast().applyToOutput(*id);
      }

    public:

      /// Implementation of CallbacksInterface::handleStatement for Jacobian tapes.
      ///
      /// Calls applyToInput for all inputs, then applyPostInputLogic, afterwards applyToOutput, and finally
      ///  applyPostOutputLogic.
      CODI_INLINE void handleStatement(Identifier& lhsIndex, codi::Config::ArgumentSize const& size,
                                       Real const* jacobians, Identifier* rhsIdentifiers) {
        CODI_UNUSED(jacobians);

        Impl& impl = cast();

        for (codi::Config::ArgumentSize i = 0; i < size; i += 1) {
          impl.applyToInput(rhsIdentifiers[i]);
        }
        impl.applyPostInputLogic();

        impl.applyToOutput(lhsIndex);
        impl.applyPostOutputLogic();
      }

      /// Implementation of CallbacksInterface::handleStatement for primal value tapes.
      ///
      /// Calls applyToInput for all inputs, then applyPostInputLogic, afterwards applyToOutput, and finally
      ///  applyPostOutputLogic.
      CODI_INLINE void handleStatement(EvalHandle const& evalHandle, codi::Config::ArgumentSize const& nPassiveValues,
                                       size_t& linearAdjointPosition, char* stmtData) {
        using StatementEvaluator = typename Tape::StatementEvaluator;

        Impl& impl = cast();

        codi::WriteInfo writeInfo;
        StatementEvaluator::template call<codi::StatementCall::WriteInformation, Tape>(evalHandle, writeInfo, primals,
                                                                                       nPassiveValues, stmtData);

        StatementEvaluator::template call<codi::StatementCall::IterateInputs, Tape>(
            evalHandle, linearAdjointPosition, reinterpret_cast<void (*)(Identifier*, void*)>(applyToInput_func), this,
            nPassiveValues, stmtData);
        impl.applyPostInputLogic();

        StatementEvaluator::template call<codi::StatementCall::IterateOutputs, Tape>(
            evalHandle, linearAdjointPosition, reinterpret_cast<void (*)(Identifier*, void*)>(applyToOutput_func), this,
            nPassiveValues, stmtData);
        impl.applyPostOutputLogic();

        if (Tape::LinearIndexHandling) {
          linearAdjointPosition += writeInfo.numberOfOutputArguments;
        }
      }

      /// Implementation of CallbacksInterface::handleLowLevelFunction.
      ///
      /// Calls applyToInput for all inputs, then applyPostInputLogic, afterwards applyToOutput for all outputs, and
      /// finally applyPostOutputLogic.
      CODI_INLINE void handleLowLevelFunction(codi::LowLevelFunctionEntry<Tape, Real, Identifier> const& func,
                                              codi::ByteDataView& llfData) {
        Impl& impl = cast();

        func.template call<codi::LowLevelFunctionEntryCallKind::IterateInputs>(
            &tape, llfData, reinterpret_cast<void (*)(Identifier*, void*)>(applyToInput_func), this);
        impl.applyPostInputLogic();

        llfData.reset();
        func.template call<codi::LowLevelFunctionEntryCallKind::IterateOutputs>(
            &tape, llfData, reinterpret_cast<void (*)(Identifier*, void*)>(applyToOutput_func), this);
        impl.applyPostOutputLogic();
      }

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };

  /// Helper for generating new identifiers.
  ///
  /// Will generate start + direction * i identifiers with i running from 0 to infinity.
  ///
  /// @tparam T_Identifier  The identifier of the CoDiPack type.
  template<typename T_Identifier>
  struct IdentifierGenerator {
      using Identifier = CODI_DD(T_Identifier, int);
      using UnsignedIdentifier = std::make_signed_t<Identifier>;

      Identifier start = 0;          ///< Start of the identifier range.
      Identifier nextFree = 0;       ///< Next generated identifier.
      UnsignedIdentifier nextDirection = 0;  ///< Step for the next generated identifier.

      std::set<Identifier> stack = {};  ///< List of all unused identifiers.

      /// Initialize the range.
      CODI_INLINE void init(Identifier start, UnsignedIdentifier dir) {
        this->start = start;
        nextFree = start;
        nextDirection = dir;
      }

      /// Free an identifier.
      CODI_INLINE void free(Identifier id) {
        stack.insert(id);
      }

      /// Generate an unused identifier.
      CODI_INLINE Identifier generateFresh() {
        Identifier gen = nextFree;
        nextFree += nextDirection;

        return gen;
      }

      /// Generate an identifier.
      CODI_INLINE Identifier generate() {
        Identifier gen = {};
        if (stack.empty()) {
          gen = generateFresh();
        } else {
          gen = *stack.begin();
          stack.erase(stack.begin());
        }

        return gen;
      }

      /// Check if the identifier is handled by this generator.
      CODI_INLINE bool isHandledByThis(Identifier id) {
        if (id == 0) {
          return false;
        }

        if (nextDirection < 0) {
          // Reverse generation
          return id > nextFree;
        } else {
          // Regular generation
          return id < nextFree;
        }
      }

      /// Get number of generated identifiers.
      CODI_INLINE Identifier getGeneratedSize() {
        return (nextFree - start) / nextDirection;
      }
  };

  /**
   *  @brief Helper class for the lifetime management of the identifiers.
   *
   *  For each statement the lifetime of the output identifier is stored. Since low level functions
   *  or primal value tapes can have multiple outputs, an extra list is created for each tape entry that has more
   *  than one output. These can be looked up by statement id and the output id.
   *
   *  @tparam T_Identifier  The identifier of the CoDiPack type.
   *  @tparam T_Lifetime    The lifetime type. For large tapes this needs to be increased. Signed type required since
   *                        negative values are used for statements with multiple outputs.
   */
  template<typename T_Identifier, typename T_Lifetime = int>
  struct LifetimeManager {
    private:
      using Identifier = CODI_DD(T_Identifier, int);  ///< See LifetimeManager.
      using Lifetime = CODI_DD(T_Lifetime, int);      ///< See LifetimeManager.

      Identifier invalidId;

      /// Lifetime of the output of each statement. Negative for statements with more than one output.
      /// This is then the lookup index in llfLivetimesId and llfLivetimeOffsets.
      std::vector<Lifetime> stmtLivetime = {};

      std::vector<Lifetime> llfLivetimeOffsets = {};  ///< Offset into llfLivetimesId and llfLivetimes for each
                                                    ///< statement.
      std::vector<Identifier> llfLivetimesId = {};  ///< Output ids of the statements. Sorted when statement is
                                                    ///< finalized.
      std::vector<Lifetime> llfLivetimes = {};      ///< Lifetimes for each id. Sorted when statement is finalized.

      Lifetime curLLFOutputPos = 0;  ///< Helper for statement id lookup. Initialized once for each statement.
      Lifetime endLLFOutputPos = 0;  ///< Helper for statement id lookup. Initialized once for each statement.

      int outputSize = 0;  ///< Helper for writing statements. Counts the number of outputs for the current statement.

    public:

      /// Constructor
      CODI_INLINE LifetimeManager(Identifier invalidId) : invalidId(invalidId) {
        llfLivetimeOffsets.push_back(0);  // First starting range.
      }

      /*******************************************************************************/
      /// @name Functions for the translation phase. Lifetimes can only be read.
      /// @{

      /// Check if the statement has just one output.
      CODI_INLINE bool isLLFStatement(Lifetime const& stmtId) {
        return stmtLivetime[stmtId] < invalidId;
      }

      /// Call before lifetimes are read for a statement.
      CODI_INLINE void prepareStatementRead(Lifetime const& stmtId) {
        // Check if this is a low level function and perform the setup.
        if (stmtLivetime[stmtId] < invalidId) {
          auto [start, end] = getLLFRange(stmtLivetime[stmtId]);

          curLLFOutputPos = start;
          endLLFOutputPos = end;
        }
      }

      /// Get the lifetime of an output of the statement. prepareStatementRead needs to be called first.
      CODI_INLINE Lifetime getLifetime(Lifetime const& stmtId, Identifier const& outputId) {
        CODI_UNUSED(outputId);

        Lifetime r = stmtLivetime[stmtId];
        if (r < invalidId) {
          // Low level function entry
          auto iterPos = std::lower_bound(llfLivetimesId.begin() + curLLFOutputPos,
                                          llfLivetimesId.begin() + endLLFOutputPos, outputId);
          codiAssert(iterPos != llfLivetimesId.end());  // If not found the id references the wrong statement.

          int pos = std::distance(llfLivetimesId.begin(), iterPos);
          r = llfLivetimes[pos];
        }

        return r;
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions for the lifetime analysis of the statement outputs.
      /// @{

      /// Set the lifetime for an output of a statement. The statement should already be finalized.
      CODI_INLINE void setLifetime(Lifetime const& stmtId, Identifier const& outputId, Lifetime const& lifetime) {
        if (stmtLivetime[stmtId] == invalidId) {  // Invalid id indicates single output for statement.
          stmtLivetime[stmtId] = lifetime;
        } else {
          // Low level function.
          codiAssert(stmtLivetime[stmtId] < 0);  // If positive a lifetime was already set.

          auto [start, end] = getLLFRange(stmtLivetime[stmtId]);

          auto iterPos = std::lower_bound(llfLivetimesId.begin() + start, llfLivetimesId.begin() + end, outputId);
          codiAssert(iterPos != llfLivetimesId.end());  // If not found the id references the wrong statement.

          int pos = std::distance(llfLivetimesId.begin(), iterPos);
          llfLivetimes[pos] = lifetime;
        }
      }

      /// Add an output to the current statement.
      CODI_INLINE void addOutputToStatement(Identifier const& id) {
        outputSize += 1;
        llfLivetimesId.push_back(id);
        llfLivetimes.push_back(invalidId);
      }

      /// Finalize the current statement.
      /// Sorts the list of output ids for the statement if it has more than one.
      CODI_INLINE void finalizeStatement() {
        if (0 == outputSize || outputSize > 1) {
          sortCurrentRange();

          // Low level function
          llfLivetimeOffsets.push_back(llfLivetimesId.size());
          stmtLivetime.push_back(generateLLFRange());
        } else {
          // Regular statement
          stmtLivetime.push_back(invalidId);

          // Discard the low level function lifetime.
          llfLivetimesId.pop_back();
          llfLivetimes.pop_back();
        }

        outputSize = 0;
      }

      /// @}

    private:

      /// Get the range for a statement with more than one output.
      CODI_INLINE std::pair<Lifetime, Lifetime> getLLFRange(Lifetime lifetime) {
        Lifetime offset = -lifetime - 2;
        return {llfLivetimeOffsets[offset], llfLivetimeOffsets[offset + 1]};
      }

      /// Create the lookup index for a statement with more than one output.
      CODI_INLINE Lifetime generateLLFRange() {
        return -(Identifier)llfLivetimeOffsets.size();
      }

      /// Sort the output ids of the current statement.
      void sortCurrentRange() {
        // Lifetimes are not set yet, so we just need to sort the identifiers.
        Lifetime start = llfLivetimeOffsets.back();
        Lifetime end = llfLivetimesId.size();

        // Validate assumption
        if (Config::EnableAssert) {
          for (Lifetime i = start; i < end; i += 1) {
            codiAssert(llfLivetimes[i] == invalidId);
          }
        }

        std::sort(llfLivetimesId.begin() + start, llfLivetimesId.end());
        llfLivetimesId.erase(std::unique(llfLivetimesId.begin() + start, llfLivetimesId.end()), llfLivetimesId.end());
        llfLivetimes.resize(llfLivetimesId.size());
      }
  };

  /**
   * @brief Reassigns the identifiers in a tape such that the tape evaluation is optimized.
   *
   * The optimization performs three steps:
   *  1. Analyze the lifetime of the left hand side identifiers for each statement.
   *  2. Based on the lifetime of the statement each statement is given a new identifier.
   *    - If the identifier is used only for a short time it gets a hot identifier. Theses identifiers are from a
   *      special region in the adjoint vector that is accessed quite often.
   *    - Otherwise, the identifier gets a cold identifier. These are taken from the remainder of the available
   *      identifiers.
   *  3. The unused identifiers are removed by shifting the cold identifiers. The adjoint vector usually becomes
   *     smaller.
   *
   *  Since the optimization requires some time, it should only be applied if the tape is evaluated quite often.
   *
   * @tparam T_Tape      Tape tape on which the optimization is applied.
   * @tparam T_Lifetime  The lifetime type. For large tapes this needs to be increased. Signed type required since
   *                     negative values are used for statements with multiple outputs.
   */
  template<typename T_Tape, typename T_Lifetime = int>
  struct IdentifierCacheOptimizerHotCold {
    public:
      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_TAPE);  ///< See IdentifierCacheOptimizerHotCold.
      using Lifetime = CODI_DD(T_Lifetime, int);        ///< See IdentifierCacheOptimizerHotCold.

    private:

      using Real = typename Tape::Real;              ///< See FullTapeInterface.
      using Identifier = typename Tape::Identifier;  ///< See FullTapeInterface.
      using EvalHandle = typename Tape::EvalHandle;  ///< See FullTapeInterface.

      /// Lookup for the lifetime of currently active identifiers.
      using IdLivetimesMap = std::multimap<Lifetime, std::pair<Identifier, Identifier>>;

      Identifier invalidId = {};  ///< See IdentifierInformationTapeInterface.
      Identifier passiveId = {};  ///< See IdentifierInformationTapeInterface.
      Tape& tape;                 ///< Tape that is modified.

      Lifetime hotLiveTimeThreshold = 500;  ///< Threshold for hot identifiers.
      Identifier idMapSize = {};            ///< Maximum size of identifiers.

      IdentifierGenerator<Identifier> generatorHot = {};   ///< Generator for hot identifiers.
      IdentifierGenerator<Identifier> generatorCold = {};  ///< Generator for cold identifiers.

      LifetimeManager<Identifier, Lifetime> lifetimes;  ///< Manager of statement lifetimes.

      /// Status entries of the analysis.
      struct Stats {
          size_t totalHot;
          size_t totalCold;
          size_t total;
          size_t unused;
      };
      Stats stats = {};  ///< Status entries of the analysis.

    public:

      /// Constructor.
      CODI_INLINE IdentifierCacheOptimizerHotCold(Tape& tape)
          : invalidId(tape.getInvalidIndex()), passiveId(tape.getPassiveIndex()), tape(tape), lifetimes(invalidId) {
        Identifier startCold = tape.getIndexManager().getLargestCreatedIndex();
        generatorCold.init(startCold, -1);
        generatorHot.init(1, 1);

        idMapSize = startCold + 1;
      }

      /// Set the threshold for hot variables. Smaller values will be hot and larger ones cold.
      CODI_INLINE void setHotLiveTimeThreshold(Lifetime hotLiveTimeThreshold) {
        this->hotLiveTimeThreshold = hotLiveTimeThreshold;
      }

    private:

      /**
       * @brief Handle the translation of the identifiers.
       *
       * Translates the current identifiers to new ones via a translate map. Each identifier on the left hand side
       * is assigned a new identifier based on its lifetime.
       */
      struct HandleTranslate : public ApplyIdentifierModification<Tape, HandleTranslate> {
          using Base = ApplyIdentifierModification<Tape, HandleTranslate>;  ///< Base class abbreviation.

          IdentifierCacheOptimizerHotCold* parent;  ///< Access general information.
          Lifetime curStmtId = 0;                   ///< Keep track of the current statement.

          std::vector<Identifier> translateMap = {};  ///< Map for id translation.

          IdLivetimesMap currentIdLivetimes = {};  /// Lookup for the lifetimes of currently used identifiers.

          /// Constructor.
          CODI_INLINE HandleTranslate(IdentifierCacheOptimizerHotCold* p) : Base(p->tape), parent(p) {
            translateMap.resize(parent->idMapSize, parent->invalidId);

            translateMap[parent->passiveId] = parent->passiveId;  // 0 -> 0
          }

          /// Add a program input. Needs to be called before the tape is evaluated with the this custom iterator.
          CODI_INLINE void addProgrammInput(Identifier& id) {
            applyToOutput(id);
          }

          /// Translates the id. Should never be called with an untranslated id.
          CODI_INLINE void applyToInput(Identifier& id) {
            Identifier& transId = translateMap[id];

            codiAssert(parent->invalidId != transId);
            id = transId;
          }

          /// Translate an identifiers.
          CODI_INLINE void applyToOutput(Identifier& id) {
            // Early out for passive id.
            if (parent->passiveId == id) {
              if (Config::EnableAssert && parent->lifetimes.isLLFStatement(curStmtId)) {
                parent->lifetimes.getLifetime(curStmtId, id);  // Get lifetime for assert.
              }
              return;
            }

            Lifetime livetime = parent->lifetimes.getLifetime(curStmtId, id);
            Identifier& transId = translateMap[id];

            codiAssert(-1 != livetime);  // There should be no undefined lifetimes.

            if (transId != parent->invalidId) {
              // Ok, duplicated output detected.
            } else {
              // Ok, untranslated id with livetime
              bool isHot = livetime < parent->hotLiveTimeThreshold;

              // Not yet translated generate an entry.
              if (isHot) {
                transId = parent->generatorHot.generate();
              } else {
                transId = parent->generatorCold.generate();
              }

              currentIdLivetimes.emplace(curStmtId + livetime, std::make_pair(id, transId));
            }

            id = transId;
          }

          /// Free all the identifiers, that are no longer used after this statement.
          CODI_INLINE void applyPostInputLogic() {
            using Iter = typename IdLivetimesMap::iterator;
            Iter cur = currentIdLivetimes.begin();
            Iter end = currentIdLivetimes.end();

            // Iterate until we find a live time that is after the current statement id.
            for (; cur != end; cur++) {
              if (cur->first > curStmtId) {
                break;
              }

              bool isHot = parent->generatorHot.isHandledByThis(cur->second.second);
              if (isHot) {
                parent->generatorHot.free(cur->second.second);
              } else {
                parent->generatorCold.free(cur->second.second);
              }
              translateMap[cur->second.first] = parent->invalidId;
            }

            // Remove all freed lifetimes.
            currentIdLivetimes.erase(currentIdLivetimes.begin(), cur);

            parent->lifetimes.prepareStatementRead(curStmtId);
          }

          /// Prepare for the next statement.
          CODI_INLINE void applyPostOutputLogic() {
            curStmtId += 1;
          }
      };

      /**
       * @brief Analysis for the lifetimes of each statement. For statements or low level function with multiple
       * outputs lifetimes are tracked independently.
       *
       * The lifetimes of an identifier is determined by the statement number when it is created and by the statement
       * number of the last use.
       */
      struct HandleHotColdAnalysis : public ApplyIdentifierModification<Tape, HandleHotColdAnalysis> {
          using Base = ApplyIdentifierModification<Tape, HandleHotColdAnalysis>;  ///< Abbreviation for the base.

          IdentifierCacheOptimizerHotCold* parent;  ///< Set lifetimes of identifiers.

          std::vector<Lifetime> idLastUseInStmt = {};  ///< Lookup map for last use.
          std::vector<Lifetime> idCreatedInStmt = {};  ///< Lookup map for creation.

          Lifetime curStmtId = 0;  ///< Counter for statement id.

          /// Constructor.
          CODI_INLINE HandleHotColdAnalysis(IdentifierCacheOptimizerHotCold* p) : Base(p->tape), parent(p) {
            idLastUseInStmt.resize(parent->idMapSize, parent->invalidId);
            idCreatedInStmt.resize(parent->idMapSize, parent->invalidId);
          }

          /// Compute the lifetime of the identifier.
          CODI_INLINE void computeLifetime(Identifier& id) {
            Lifetime& createStmtId = idCreatedInStmt[id];
            Lifetime& lastUseStmtId = idLastUseInStmt[id];

            if (parent->invalidId != createStmtId &&  // ID was created
                parent->invalidId != lastUseStmtId    // ID was actually used
            ) {
              Lifetime livetime = lastUseStmtId - createStmtId;

              parent->lifetimes.setLifetime(createStmtId, id, livetime);
            } else if (parent->invalidId == createStmtId && parent->invalidId == lastUseStmtId) {
              // New identifier, it is used for the first time.
            } else if (parent->invalidId == lastUseStmtId && parent->invalidId != createStmtId) {
              // Identifier was not used.
              parent->lifetimes.setLifetime(createStmtId, id, 0);
            } else {
              CODI_EXCEPTION("Identifier '%d' is used but not created, this is an error in the tape.", (int)id);
            }

            lastUseStmtId = parent->invalidId;  // Reset last use.
          }

          /// Add a program input. Needs to be called before the tape is evaluated with the this custom iterator.
          CODI_INLINE void addProgrammInput(Identifier& id) {
            if (id != parent->passiveId) {
              idCreatedInStmt[id] = curStmtId;
            }

            parent->lifetimes.addOutputToStatement(id);
          }

          /// Updates the last use of the identifier with the current statement number.
          CODI_INLINE void applyToInput(Identifier& id) {
            if (id != parent->passiveId) {
              idLastUseInStmt[id] = curStmtId;
            }
          }

          /// Update the creation of the identifier with the current statement number.
          CODI_INLINE void applyToOutput(Identifier& id) {
            if (id != parent->passiveId) {
              computeLifetime(id);
              idCreatedInStmt[id] = curStmtId;
            }

            parent->lifetimes.addOutputToStatement(id);
          }

          /// Finalizes the statement and prepares for the next statement.
          CODI_INLINE void applyPostOutputLogic() {
            parent->lifetimes.finalizeStatement();

            curStmtId += 1;
          }

          /// Update the lifetime of an output such that it lives longer that the last statement. Called after the tape
          /// iteration process.
          CODI_INLINE void setOutputLifetime(Identifier& id) {
            idLastUseInStmt[id] = curStmtId + parent->hotLiveTimeThreshold + 1;
          }

          /// Update the lifetimes of all identifiers that have not been overwritten yet.
          CODI_INLINE void finalize() {
            Identifier curId = 0;
            for (Lifetime& stmtId : idCreatedInStmt) {
              if (stmtId != parent->invalidId) {
                computeLifetime(curId);
                stmtId = parent->invalidId;
              }
              curId += 1;
            }
          }
      };

      /**
       * @brief Shift the cold identifiers down to remove the gap between hot and cold identifiers.
       */
      struct HandleShift : public ApplyIdentifierModification<Tape, HandleShift> {
          using Base = ApplyIdentifierModification<Tape, HandleShift>;  ///< Abbreviation for the base.

          IdentifierCacheOptimizerHotCold* parent;  ///< Access to hot identifier generator.
          Identifier coldShift;                     ///< How much cold identifiers are shifted.

          /// Constructors.
          CODI_INLINE HandleShift(IdentifierCacheOptimizerHotCold* p)
              : Base(p->tape), parent(p), coldShift(parent->stats.unused) {}

          /// Apply the shift to a cold id.
          CODI_INLINE void applyShift(Identifier& id) {
            if (0 != id && !parent->generatorHot.isHandledByThis(id)) {
              id -= coldShift;
            }
          }

          /// Apply shift to an input identifier.
          CODI_INLINE void applyToInput(Identifier& id) {
            applyShift(id);
          }

          /// Apply shift to an output identifier.
          CODI_INLINE void applyToOutput(Identifier& id) {
            applyShift(id);
          }
      };

      /// Compute the statistics for the optimization.
      CODI_INLINE void updateStats() {
        stats.totalHot = generatorHot.getGeneratedSize() + 1;  // +1 for zero index.
        stats.totalCold = generatorCold.getGeneratedSize();
        stats.total = stats.totalHot + stats.totalCold;
        stats.unused = generatorCold.start - stats.total + 1;  // +1 since start is included in range.
      }

    public:

      /// Perform the tape cache optimization. Se the class description for details.
      template<typename FuncIn, typename FuncOut>
      CODI_NO_INLINE void eval(FuncIn&& iterIn, FuncOut&& iterOut) {
        // Do hot cold analysis
        {
          HandleHotColdAnalysis hotCold = {this};

          // Add inputs as one large low level function.
          iterIn([&](Identifier& id) {
            hotCold.addProgrammInput(id);
          });
          hotCold.applyPostOutputLogic();

          // Analyze the tape.
          tape.iterateForward(hotCold);

          iterOut([&](Identifier& id) {
            hotCold.setOutputLifetime(id);
          });

          // Finalize the analysis.
          hotCold.finalize();
        }

        // Translate tape
        {
          HandleTranslate translate = {this};

          // Handle inputs as large low level function.
          translate.applyPostInputLogic();
          iterIn([&](Identifier& id) {
            translate.addProgrammInput(id);  // Register the starting use of the global inputs.
          });
          translate.applyPostOutputLogic();

          tape.iterateForward(translate);

          iterOut([&](Identifier& id) {
            translate.applyToInput(id);  // Just translate do not generate new translations.
          });
        }

        updateStats();

        // Perform identifier shift
        {
          HandleShift shift = {this};
          tape.iterateForward(shift);

          auto doShift = [&](Identifier& id) {
            shift.applyShift(id);
          };
          iterIn(doShift);
          iterOut(doShift);
        }
      }

      /// Get the new largest created index.
      CODI_INLINE size_t getLargestCreatedIndex() {
        return stats.total - 1;
      }

      /// Write statistics to a stream as a list.
      template<typename Stream>
      void writeStatsVerbose(Stream& out) {
        out << "Hot: " << stats.totalHot << std::endl;
        out << "Cold: " << stats.totalCold << std::endl;
        out << "Total: " << stats.total << std::endl;
        out << "Unused: " << stats.unused << std::endl;
      }

      /// Write the header for the statistics to a stream.
      template<typename Stream>
      void writeStatsHeader(Stream& out) {
        out << "Hot; Cold; Total; Unused;";
      }

      /// Write the data for this optimizer into a row.
      template<typename Stream>
      void writeStatsRow(Stream& out) {
        out << "Hot: " << stats.totalHot << std::endl;
        out << "Cold: " << stats.totalCold << std::endl;
        out << "Total: " << stats.total << std::endl;
        out << "Unused: " << stats.unused << std::endl;
      }
  };
}
