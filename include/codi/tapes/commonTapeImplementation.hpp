/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <algorithm>
#include <type_traits>

#include "../config.h"
#include "../misc/byteDataView.hpp"
#include "../misc/eventSystem.hpp"
#include "../misc/fileIo.hpp"
#include "../misc/macros.hpp"
#include "../misc/temporaryMemory.hpp"
#include "data/dataInterface.hpp"
#include "data/position.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/fullTapeInterface.hpp"
#include "misc/externalFunction.hpp"
#include "misc/lowLevelFunctionEntry.hpp"
#include "misc/vectorAccessInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Interface for the definition of tape types.
   *
   * In CoDiPack, each tape has to define its tape types as a separate struct. As a minimum requirement, tape types have
   * to make the definitions showcased in this interface.
   */
  struct TapeTypesInterface {
    public:

      using Real = CODI_ANY;        ///< Primal computation type, e.g. double.
      using Gradient = CODI_ANY;    ///< Gradient computation type, e.g. double or Direction.
      using Identifier = CODI_ANY;  ///< Identifier for the internal management, e.g. int.

      /// Indicates the storage strategy that will be used by all data vectors. See DataInterface and its
      /// implementations.
      template<typename Chunk, typename Nested>
      using Data = DataInterface<Nested>;

      using NestedData = DataInterface<>;  ///< The root vector of the tape implementation on which the
                                           ///< CommonTapeImplementation builds.
  };

  /**
   * @brief Declares all types used in the CommonTapeImplementation.
   *
   * @tparam T_TapeTypes  Must implement TapeTypesInterface.
   */
  template<typename T_TapeTypes>
  struct CommonTapeTypes {
    public:

      using TapeTypes = CODI_DD(T_TapeTypes, TapeTypesInterface);  ///< See CommonTapeTypes.

      template<typename Chunk, typename Nested>
      using Data = typename TapeTypes::template Data<Chunk, Nested>;  ///< See TapeTypesInterface.

      using NestedData = typename TapeTypes::NestedData;  ///< See TapeTypesInterface.

      /// Token and size data chunk.
      using LowLevelFunctionInfoChunk = Chunk2<Config::LowLevelFunctionToken, Config::LowLevelFunctionDataSize>;
      /// Token and size data for low level functions.
      using LowLevelFunctionInfoData = Data<LowLevelFunctionInfoChunk, NestedData>;

      /// Byte data chunk.
      using LowLevelFunctionByteChunk = Chunk1<char>;
      /// Byte data for low level functions.
      using LowLevelFunctionByteData = Data<LowLevelFunctionByteChunk, LowLevelFunctionInfoData>;

      using Position = typename LowLevelFunctionByteData::Position;  ///< Global position of the tape.
  };

  /**
   * @brief Implementation of all common tape functionality.
   *
   * This basic implementation provides the following functionality:
   *  - external function support with the external functions stored in externalFunctionData,
   *  - tape options gathering,
   *  - activity tracking.
   *
   * It also provides functionality that can be implemented with other functions:
   *  - setter and getter methods,
   *  - non-positional evaluation methods,
   *  - registerOutput,
   *  - TapeValues functions,
   *  - reset functionality.
   *
   * @tparam T_ImplTapeTypes must implement TapeTypesInterface.
   * @tparam T_Impl Type of the full tape implementation.
   */
  template<typename T_ImplTapeTypes, typename T_Impl>
  struct CommonTapeImplementation
      : public FullTapeInterface<typename T_ImplTapeTypes::Real, typename T_ImplTapeTypes::Gradient,
                                 typename T_ImplTapeTypes::Identifier,
                                 typename CommonTapeTypes<T_ImplTapeTypes>::Position> {
    public:

      using ImplTapeTypes = CODI_DD(T_ImplTapeTypes, TapeTypesInterface);  ///< See CommonTapeImplementation.
      using Impl = CODI_DD(T_Impl, CommonTapeImplementation);              ///< See CommonTapeImplementation.

      using Real = typename ImplTapeTypes::Real;              ///< See TapeTypesInterface.
      using Gradient = typename ImplTapeTypes::Gradient;      ///< See TapeTypesInterface.
      using Identifier = typename ImplTapeTypes::Identifier;  ///< See TapeTypesInterface.

      /// See CommonTapeTypes.
      using LowLevelFunctionInfoData = typename CommonTapeTypes<ImplTapeTypes>::LowLevelFunctionInfoData;
      /// See CommonTapeTypes.
      using LowLevelFunctionByteData = typename CommonTapeTypes<ImplTapeTypes>::LowLevelFunctionByteData;
      using Position = typename CommonTapeTypes<ImplTapeTypes>::Position;  ///< See TapeTypesInterface.

      using NestedData = LowLevelFunctionByteData;                         ///< Shorthand.
      using NestedPosition = typename LowLevelFunctionByteData::Position;  ///< Shorthand.

    protected:

      bool active;                       ///< Whether or not the tape is in recording mode.
      std::set<TapeParameters> options;  ///< All options.

      LowLevelFunctionInfoData llfInfoData;  ///< Token and size data for low level functions.
      LowLevelFunctionByteData llfByteData;  ///< Byte data for low level functions.

      Real manualPushLhsValue;             ///< For storeManual, remember the value assigned to the lhs.
      Identifier manualPushLhsIdentifier;  ///< For storeManual, remember the identifier assigned to the lhs.
      size_t manualPushGoal;               ///< Store the number of expected pushes after a storeManual call.
      size_t manualPushCounter;            ///< Count the pushes after storeManual, to identify the last push.

      TemporaryMemory allocator;  ///< Allocator for temporary memory.

      /// Lookup table for low level function.
      static std::vector<LowLevelFunctionEntry<Impl, Real, Identifier>>* lowLevelFunctionLookup;

    private:

      /// External function token is always added first.
      static Config::LowLevelFunctionToken constexpr EXTERNAL_FUNCTION_TOKEN = 0;

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      CODI_INLINE void resetInternal(bool resetAdjoints, AdjointsManagement adjointsManagement,
                                     EventHints::Reset kind) {
        EventSystem<Impl>::notifyTapeResetListeners(cast(), this->getZeroPosition(), kind, resetAdjoints);

        if (resetAdjoints) {
          cast().clearAdjoints(adjointsManagement);
        }

        deleteLowLevelFunctionData(cast().getZeroPosition());

        llfByteData.reset();

        // Requires extra reset since the default vector implementation forwards to resetTo
        cast().indexManager.get().reset();
      }

    protected:

      /// Initialize all manual push data, including the counter. Check that a previous manual store is completed.
      CODI_INLINE void initializeManualPushData(Real const& lhsValue, Identifier const& lhsIndex, size_t size) {
        codiAssert(this->manualPushGoal == this->manualPushCounter);
        if (Config::StatementEvents || Config::EnableAssert) {
          this->manualPushLhsValue = lhsValue;
          this->manualPushLhsIdentifier = lhsIndex;
          this->manualPushCounter = 0;
          this->manualPushGoal = size;
        }
      }

      /// Increment the manual push counter. Check against the declared push goal.
      CODI_INLINE void incrementManualPushCounter() {
        codiAssert(this->manualPushCounter < this->manualPushGoal);

        if (Config::StatementEvents || Config::EnableAssert) {
          this->manualPushCounter += 1;
        }
      }

    protected:

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      TapeValues internalGetTapeValues() const;  ///< Create tape values.

      /// @}

    public:

      /// Constructor
      CommonTapeImplementation()
          : active(false),
            options(),
            llfInfoData(Config::SmallChunkSize),
            llfByteData(Config::ByteDataChunkSize),
            manualPushLhsValue(),
            manualPushLhsIdentifier(),
            manualPushGoal(),
            manualPushCounter(),
            allocator() {
        options.insert(TapeParameters::LLFByteDataSize);
        options.insert(TapeParameters::LLFInfoDataSize);

        if (nullptr == lowLevelFunctionLookup) {
          lowLevelFunctionLookup = new std::vector<LowLevelFunctionEntry<Impl, Real, Identifier>>();

          // Add external function token. So EXTERNAL_FUNCTION_TOKEN is always zero.
          Config::LowLevelFunctionToken token =
              registerLowLevelFunction(ExternalFunctionLowLevelEntryMapper<Impl, Real, Identifier>::create());
          if (token != EXTERNAL_FUNCTION_TOKEN) {
            CODI_EXCEPTION("External function token is not zero.");
          }
        }
      }

      /// Do not allow copy construction.
      CommonTapeImplementation(CommonTapeImplementation const&) = delete;

      /// Do not allow copy assignment.
      CommonTapeImplementation& operator=(CommonTapeImplementation const&) = delete;

      /// Do not allow move construction. Relevant use cases should be covered by \ref swap.
      CommonTapeImplementation(CommonTapeImplementation&&) = delete;

      /// Do not allow move assignment. Relevant use cases should be covered by \ref swap.
      CommonTapeImplementation& operator=(CommonTapeImplementation&&) = delete;

      /*******************************************************************************/
      /// @name Functions from GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::setGradient()
      void setGradient(Identifier const& identifier, Gradient const& gradient,
                       AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        cast().gradient(identifier, adjointsManagement) = gradient;
      }

      /// \copydoc codi::GradientAccessTapeInterface::getGradient()
      Gradient const& getGradient(Identifier const& identifier,
                                  AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const {
        return cast().gradient(identifier, adjointsManagement);
      }

      // Gradient functions are not implemented.

      /// @}
      /*******************************************************************************/
      /// @name Functions from ReverseTapeInterface
      /// @{

      /// \copydoc codi::ReverseTapeInterface::evaluate(AdjointsManagement)
      void evaluate(AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        Impl& impl = cast();

        impl.evaluate(impl.getPosition(), impl.getZeroPosition(), adjointsManagement);
      }

      /// \copydoc codi::ReverseTapeInterface::registerOutput()
      template<typename Lhs>
      void registerOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        cast().template store<Lhs, Lhs>(value, static_cast<ExpressionInterface<Real, Lhs> const&>(value));
        EventSystem<Impl>::notifyTapeRegisterOutputListeners(cast(), value.cast().value(),
                                                             value.cast().getIdentifier());
      }

      /// \copydoc codi::ReverseTapeInterface::setActive()
      void setActive() {
        EventSystem<Impl>::notifyTapeStartRecordingListeners(cast());
        active = true;
      }

      /// \copydoc codi::ReverseTapeInterface::setPassive()
      void setPassive() {
        EventSystem<Impl>::notifyTapeStopRecordingListeners(cast());
        active = false;
      }

      /// \copydoc codi::ReverseTapeInterface::isActive()
      bool isActive() const {
        return active;
      }

      /// \copydoc codi::ReverseTapeInterface::printStatistics()
      template<typename Stream = std::ostream>
      void printStatistics(Stream& out = std::cout) const {
        cast().getTapeValues().formatDefault(out);
      }

      /// \copydoc codi::ReverseTapeInterface::printTableHeader()
      template<typename Stream = std::ostream>
      void printTableHeader(Stream& out = std::cout) const {
        cast().getTapeValues().formatHeader(out);
      }

      /// \copydoc codi::ReverseTapeInterface::printTableRow()
      template<typename Stream = std::ostream>
      void printTableRow(Stream& out = std::cout) const {
        cast().getTapeValues().formatRow(out);
      }

      /// \copydoc codi::ReverseTapeInterface::getTapeValues()
      TapeValues getTapeValues() const {
        TapeValues values = cast().internalGetTapeValues();

        values.addSection("Low level function info data entries");
        llfInfoData.addToTapeValues(values);
        values.addSection("Low level function byte data entries");
        llfByteData.addToTapeValues(values);

        return values;
      }

      /// \copydoc codi::ReverseTapeInterface::reset(bool, AdjointsManagement)
      CODI_INLINE void reset(bool resetAdjoints = true,
                             AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        resetInternal(resetAdjoints, adjointsManagement, EventHints::Reset::Full);
      }

      // clearAdjoints and reset(Position) are not implemented.

      /// @}
      /*******************************************************************************/
      /// @name Functions from DataManagementTapeInterface
      /// @{

      /// \copydoc codi::DataManagementTapeInterface::swap()
      void swap(Impl& other) {
        std::swap(active, other.active);

        llfByteData.swap(other.llfByteData);
      }

      /// \copydoc codi::DataManagementTapeInterface::resetHard()
      void resetHard() {
        Impl& impl = cast();

        // First perform a regular reset.
        resetInternal(false, AdjointsManagement::Automatic, EventHints::Reset::Hard);

        // Then perform the hard resets.
        impl.deleteAdjointVector();

        llfByteData.resetHard();
      }

      /// @}

    private:
      static void writeFunction(const ChunkBase* chunk, FileIo& handle) {
        chunk->writeData(handle);
      }

      static void readFunction(ChunkBase* chunk, FileIo& handle) {
        chunk->readData(handle);
      }

      static void deleteFunction(ChunkBase* chunk) {
        chunk->deleteData();
      }

    public:

      /// @{
      /// \copydoc codi::DataManagementTapeInterface::writeToFile()
      void writeToFile(const std::string& filename) {
        FileIo io(filename, true);

        llfByteData.forEachChunk(writeFunction, true, io);
      }

      /// \copydoc codi::DataManagementTapeInterface::readFromFile()
      void readFromFile(const std::string& filename) {
        FileIo io(filename, false);

        llfByteData.forEachChunk(readFunction, true, io);
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteData()
      void deleteData() {
        llfByteData.forEachChunk(deleteFunction, true);
      }

      /// \copydoc codi::DataManagementTapeInterface::getAvailableParameters()
      std::set<TapeParameters> const& getAvailableParameters() const {
        return options;
      }

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      /// <br><br> Implementation: Handles LLFByteDataSize, LLFInfoDataSize
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::LLFByteDataSize:
            return llfByteData.getDataSize();
            break;
          case TapeParameters::LLFInfoDataSize:
            return llfInfoData.getDataSize();
            break;
          case TapeParameters::ExternalFunctionsSize:
            CODI_WARNING(
                "Tape parameter 'ExternalFunctionsSize' no longer supported. Use 'LLFInfoDataSize' and "
                "'LLFByteDataSize' instead.");
            return 0;
            break;
          default:
            CODI_EXCEPTION("Tried to get undefined parameter for tape.");
            return 0;
            break;
        }
      }

      /// \copydoc codi::DataManagementTapeInterface::hasParameter()
      bool hasParameter(TapeParameters parameter) const {
        return options.cend() != options.find(parameter);
      }

      /// \copydoc codi::DataManagementTapeInterface::setParameter()
      /// <br><br> Implementation: Handles LLFByteDataSize, LLFInfoDataSize
      void setParameter(TapeParameters parameter, size_t value) {
        switch (parameter) {
          case TapeParameters::LLFByteDataSize:
            llfByteData.resize(value);
            break;
          case TapeParameters::LLFInfoDataSize:
            llfInfoData.resize(value);
            break;
          case TapeParameters::ExternalFunctionsSize:
            CODI_WARNING(
                "Tape parameter 'ExternalFunctionsSize' is no longer supported. Use 'LLFInfoDataSize' and "
                "'LLFByteDataSize' instead.");
            break;
          default:
            CODI_EXCEPTION("Tried to set undefined parameter for tape.");
            break;
        }
      }

      // createVectorAccess is not implemented.
      // createVectorAccessCustomAdjoints is not implemented.
      // resizeAdjointVector is not implemented.
      // deleteAdjointVector is not implemented.
      // beginUseAdjointVector is not implemented.
      // endUseAdjointVector is not implemented.

      /// @}
      /*******************************************************************************/
      /// @name Functions from LowLevelFunctionTapeInterface
      /// @{

    protected:

      /// @brief Called by the implementing tapes to store a low level function. The size is reserved and allocated.
      /// The data view is populated with the pointer and can be used to write the data.
      CODI_INLINE void internalStoreLowLevelFunction(Config::LowLevelFunctionToken token, size_t size,
                                                     ByteDataView& dataView) {
        codiAssert((size_t)token < lowLevelFunctionLookup->size());
        if (size >= Config::LowLevelFunctionDataSizeMax) {
          CODI_EXCEPTION(
              "Requested size for low level function is to big. Increase "
              "codi::Config::LowLevelFunctionDataSize or perform a dynamic memory allocation.");
        }

        llfInfoData.reserveItems(1);
        llfByteData.reserveItems(size);

        llfInfoData.pushData(token, size);

        char* dataPointer = nullptr;
        llfByteData.getDataPointers(dataPointer);
        dataView.init(dataPointer, 0, size);
        llfByteData.addDataSize(size);
      }

      /// @brief Called by the implementing tapes during a tape evaluation when a low level function statement has been
      /// reached.
      template<LowLevelFunctionEntryCallKind callType, typename... Args>
      CODI_INLINE static void callLowLevelFunction(Impl& impl, bool forward,
                                                   /* data from low level function byte data vector */
                                                   size_t& curLLFByteDataPos, char* dataPtr,
                                                   /* data from low level function info data vector */
                                                   size_t& curLLFTInfoDataPos,
                                                   Config::LowLevelFunctionToken* const tokenPtr,
                                                   Config::LowLevelFunctionDataSize* const dataSizePtr,
                                                   Args&&... args) {
        if (!forward) {
          curLLFTInfoDataPos -= 1;
          curLLFByteDataPos -= dataSizePtr[curLLFTInfoDataPos];
        }

        size_t endPos = curLLFByteDataPos + dataSizePtr[curLLFTInfoDataPos];
        ByteDataView dataView(dataPtr, curLLFByteDataPos, endPos);

        Config::LowLevelFunctionToken id = tokenPtr[curLLFTInfoDataPos];
        LowLevelFunctionEntry<Impl, Real, Identifier> const& func = (*lowLevelFunctionLookup)[id];
        if (func.template has<callType>()) CODI_Likely {
          func.template call<callType>(&impl, dataView, std::forward<Args>(args)...);

          codiAssert(endPos == dataView.getPosition());
        } else if (LowLevelFunctionEntryCallKind::Delete == callType) CODI_Unlikely {
          // No delete registered. Data is skiped by the curLLFByteDataPos update.
        } else CODI_Unlikely {
          CODI_EXCEPTION("Requested call is not supported for low level function with token '%d'.", (int)id);
        }

        if (forward) {
          curLLFByteDataPos += dataSizePtr[curLLFTInfoDataPos];
          curLLFTInfoDataPos += 1;
        }
      }

    public:

      /// @copydoc LowLevelFunctionTapeInterface::getTemporaryMemory()
      CODI_INLINE TemporaryMemory& getTemporaryMemory() {
        return allocator;
      }

      /// @copydoc LowLevelFunctionTapeInterface::registerLowLevelFunction()
      CODI_INLINE Config::LowLevelFunctionToken registerLowLevelFunction(
          LowLevelFunctionEntry<Impl, Real, Identifier> const& entry) {
        codiAssert(lowLevelFunctionLookup->size() < Config::LowLevelFunctionTokenMaxSize);

        Config::LowLevelFunctionToken token =
            static_cast<Config::LowLevelFunctionToken>(lowLevelFunctionLookup->size());
        lowLevelFunctionLookup->push_back(entry);

        return token;
      }

      // pushLowLevelFunction is not implemented.

      /// @}
      /*******************************************************************************/
      /// @name Functions from ExternalFunctionTapeInterface
      /// @{

      /// \copydoc codi::ExternalFunctionTapeInterface::pushExternalFunction()
      void pushExternalFunction(ExternalFunction<Impl> const& extFunc) {
        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive())) {
          ExternalFunctionLowLevelEntryMapper<Impl, Real, Identifier>::store(cast(), EXTERNAL_FUNCTION_TOKEN, extFunc);
        }
      }

      // registerExternalFunctionOutput is not implemented.

      /// @}
      /*******************************************************************************/
      /// @name Functions from ForwardEvaluationTapeInterface
      /// @{

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      void evaluateForward(AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        Impl& impl = cast();

        impl.evaluateForward(impl.getZeroPosition(), impl.getPosition(), adjointsManagement);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from IdentifierInformationTapeInterface
      /// @{

      /// \copydoc codi::IdentifierInformationTapeInterface::getPassiveIndex()
      Identifier getPassiveIndex() const {
        return IndexManagerInterface<Identifier>::InactiveIndex;
      }

      /// \copydoc codi::IdentifierInformationTapeInterface::getInvalidIndex()
      Identifier getInvalidIndex() const {
        return IndexManagerInterface<Identifier>::InvalidIndex;
      }

      /// \copydoc codi::IdentifierInformationTapeInterface::isIdentifierActive()
      bool isIdentifierActive(Identifier const& index) const {
        return index != cast().getPassiveIndex();
      }

      /// \copydoc codi::IdentifierInformationTapeInterface::deactivateValue()
      template<typename Lhs>
      void deactivateValue(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        value = value.getValue();
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PositionalEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PositionalEvaluationTapeInterface::getPosition()
      Position getPosition() const {
        return llfByteData.getPosition();
      }

      /// \copydoc codi::PositionalEvaluationTapeInterface::getZeroPosition()
      Position getZeroPosition() const {
        return llfByteData.getZeroPosition();
      }

      /// @}

    protected:

      /// Delete all external function data up to `pos`.
      void deleteLowLevelFunctionData(Position const& pos) {
        // Clear external function data.
        auto deleteFunc = [this](
                              /* data from low level function byte data vector */
                              size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
                              /* data from low level function info data vector */
                              size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos,
                              Config::LowLevelFunctionToken* const tokenPtr,
                              Config::LowLevelFunctionDataSize* const dataSizePtr) {
          CODI_UNUSED(endLLFByteDataPos);

          while (curLLFInfoDataPos > endLLFInfoDataPos) {
            callLowLevelFunction<LowLevelFunctionEntryCallKind::Delete>(cast(), false, curLLFByteDataPos, dataPtr,
                                                                        curLLFInfoDataPos, tokenPtr, dataSizePtr);
          }
        };

        llfByteData.template evaluateReverse<1>(cast().getPosition(), pos, deleteFunc);
      }

    public:

      /// @{

      /// \copydoc ::codi::PositionalEvaluationTapeInterface::resetTo
      CODI_INLINE void resetTo(Position const& pos, bool resetAdjoints = true,
                               AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        EventSystem<Impl>::notifyTapeResetListeners(cast(), pos, EventHints::Reset::To, resetAdjoints);

        if (resetAdjoints) {
          Impl& impl = cast();
          impl.clearAdjoints(impl.getPosition(), pos, adjointsManagement);
        }

        deleteLowLevelFunctionData(pos);

        llfByteData.resetTo(pos);
      }

      // clearAdjoints and evaluate are not implemented.

      /// @}
      /*******************************************************************************/
      /// @name Functions from PrimalEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PrimalEvaluationTapeInterface::evaluatePrimal()
      void evaluatePrimal() {
        Impl& impl = cast();

        impl.evaluatePrimal(impl.getZeroPosition(), impl.getPosition());
      }

      /// \copydoc codi::PrimalEvaluationTapeInterface::setPrimal()
      void setPrimal(Identifier const& identifier, Real const& primal) {
        cast().primal(identifier) = primal;
      }

      /// \copydoc codi::PrimalEvaluationTapeInterface::getPrimal()
      Real const& getPrimal(Identifier const& identifier) const {
        return cast().primal(identifier);
      }

      /// @}

    protected:

      /*******************************************************************************/
      /// @name Internal helper functions
      /// @{

      /// Initialize the base class
      void init(typename ImplTapeTypes::NestedData* nested) {
        llfInfoData.setNested(nested);
        llfByteData.setNested(&llfInfoData);
      }

      /// @}
  };

  template<typename ImplTapeTypes, typename Impl>
  std::vector<LowLevelFunctionEntry<Impl, typename ImplTapeTypes::Real, typename ImplTapeTypes::Identifier>>*
      CommonTapeImplementation<ImplTapeTypes, Impl>::lowLevelFunctionLookup = nullptr;
}
