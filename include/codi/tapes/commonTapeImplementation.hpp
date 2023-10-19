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

#include <algorithm>
#include <type_traits>

#include "../config.h"
#include "../misc/byteDataStore.hpp"
#include "../misc/eventSystem.hpp"
#include "../misc/fileIo.hpp"
#include "../misc/macros.hpp"
#include "../misc/temporaryMemoryAllocator.hpp"
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

      using OtherFixedDataChunk = Chunk1<char>;                      ///< Byte data chunk.
      using OtherFixedData = Data<OtherFixedDataChunk, NestedData>;  ///< Byte data for fixed data that is always read.

      using OtherDynamicDataChunk = Chunk1<char>;                            ///< Byte data chunk.
      using OtherDynamicData = Data<OtherDynamicDataChunk, OtherFixedData>;  ///< Byte data for data that is dynamic and
                                                                             ///< based on the fixed data.

      using Position = typename OtherDynamicData::Position;  ///< Global position of the tape.
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

      using OtherFixedData = typename CommonTapeTypes<ImplTapeTypes>::OtherFixedData;      ///< See CommonTapeTypes.
      using OtherDynamicData = typename CommonTapeTypes<ImplTapeTypes>::OtherDynamicData;  ///< See CommonTapeTypes.
      using Position = typename CommonTapeTypes<ImplTapeTypes>::Position;                  ///< See TapeTypesInterface.

      using NestedData = OtherDynamicData;                         ///< Shorthand.
      using NestedPosition = typename OtherDynamicData::Position;  ///< Shorthand.

    protected:

      bool active;                       ///< Whether or not the tape is in recording mode.
      std::set<TapeParameters> options;  ///< All options.

      OtherFixedData otherFixedData;      ///< Byte data for fixed data that is always read.
      OtherDynamicData otherDynamicData;  ///< Byte data for data that is dynamic and based on the fixed data.

      Real manualPushLhsValue;             ///< For storeManual, remember the value assigned to the lhs.
      Identifier manualPushLhsIdentifier;  ///< For storeManual, remember the identifier assigned to the lhs.
      size_t manualPushGoal;               ///< Store the number of expected pushes after a storeManual call.
      size_t manualPushCounter;            ///< Count the pushes after storeManual, to identify the last push.

      TemporaryMemoryAllocator allocator;  ///< Allocator for temporary memory.

      /// Lookup table for low level function.
      static std::vector<LowLevelFunctionEntry<Impl, Real, Identifier>>* lowLevelFunctionLookup;

    private:

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

        deleteExternalFunctionUserData(cast().getZeroPosition());

        otherDynamicData.reset();

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
            otherFixedData(Config::ByteDataChunkSize),
            otherDynamicData(Config::ByteDataChunkSize),
            manualPushLhsValue(),
            manualPushLhsIdentifier(),
            manualPushGoal(),
            manualPushCounter(),
            allocator() {
        options.insert(TapeParameters::OtherFixedDataSize);
        options.insert(TapeParameters::OtherDynamicDataSize);

        if (nullptr == lowLevelFunctionLookup) {
          lowLevelFunctionLookup = new std::vector<LowLevelFunctionEntry<Impl, Real, Identifier>>();

          // Add external function token.
          lowLevelFunctionLookup->push_back(ExternalFunctionLowLevelEntryMapper<Impl, Real, Identifier>::create());
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

        values.addSection("Other fixed data entries");
        otherFixedData.addToTapeValues(values);
        values.addSection("Other variable data entries");
        otherDynamicData.addToTapeValues(values);

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

        otherDynamicData.swap(other.otherDynamicData);
      }

      /// \copydoc codi::DataManagementTapeInterface::resetHard()
      void resetHard() {
        Impl& impl = cast();

        // First perform a regular reset.
        resetInternal(false, AdjointsManagement::Automatic, EventHints::Reset::Hard);

        // Then perform the hard resets.
        impl.deleteAdjointVector();

        otherDynamicData.resetHard();
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

        otherDynamicData.forEachChunk(writeFunction, true, io);
      }

      /// \copydoc codi::DataManagementTapeInterface::readFromFile()
      void readFromFile(const std::string& filename) {
        FileIo io(filename, false);

        otherDynamicData.forEachChunk(readFunction, true, io);
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteData()
      void deleteData() {
        otherDynamicData.forEachChunk(deleteFunction, true);
      }

      /// \copydoc codi::DataManagementTapeInterface::getAvailableParameters()
      std::set<TapeParameters> const& getAvailableParameters() const {
        return options;
      }

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      /// <br><br> Implementation: Handles OtherFixedDataSize, OtherDynamicDataSize
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::OtherFixedDataSize:
            return otherFixedData.getDataSize();
            break;
          case TapeParameters::OtherDynamicDataSize:
            return otherDynamicData.getDataSize();
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
      /// <br><br> Implementation: Handles OtherFixedDataSize, OtherDynamicDataSize
      void setParameter(TapeParameters parameter, size_t value) {
        switch (parameter) {
          case TapeParameters::OtherFixedDataSize:
            return otherFixedData.resize(value);
            break;
          case TapeParameters::OtherDynamicDataSize:
            return otherDynamicData.resize(value);
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

      /// @brief Called by the implementing tapes to store a low level function. The sizes are reserved and allocated.
      /// The data stores are populated with the pointers and can be used to write the data.
      CODI_INLINE void internalStoreLowLevelFunction(size_t fixedSize, size_t dynamicSize, ByteDataStore& fixedData,
                                                     ByteDataStore& dynamicData) {
        otherFixedData.reserveItems(fixedSize);
        otherDynamicData.reserveItems(dynamicSize);

        char* fixedPointer = nullptr;
        char* dynamicPointer = nullptr;
        otherFixedData.getDataPointers(fixedPointer);
        otherDynamicData.getDataPointers(dynamicPointer);
        fixedData.init(fixedPointer, 0, ByteDataStore::Direction::Forward);
        dynamicData.init(dynamicPointer, 0, ByteDataStore::Direction::Forward);
        otherFixedData.addDataSize(fixedSize);
        otherDynamicData.addDataSize(dynamicSize);
      }

      /// @brief Called by the implementing tapes during a tape evaluation when a low level function statement has been
      /// reached.
      template<LowLevelFunctionEntryCallType callType, typename... Args>
      CODI_INLINE static void handleLowLevelFunction(Impl& impl, ByteDataStore::Direction direction,
                                                     /* data from other dynamic data vector */
                                                     size_t& curOtherDynamicDataPos, char const* const otherDynamicPtr,
                                                     /* data from other fixed data vector */
                                                     size_t& curOtherFixedDataPos, char const* const otherFixedPtr,
                                                     Args&&... args) {
        ByteDataStore fixedStore(const_cast<char*>(otherFixedPtr), curOtherFixedDataPos, direction);
        ByteDataStore dynamicStore(const_cast<char*>(otherDynamicPtr), curOtherDynamicDataPos, direction);

        Config::LowLevelFunctionToken magicNumber = *fixedStore.read<Config::LowLevelFunctionToken>(1);
        LowLevelFunctionEntry<Impl, Real, Identifier> const& func = (*lowLevelFunctionLookup)[magicNumber];
        if (func.template has<callType>()) CODI_Likely {
          func.template call<callType>(&impl, fixedStore, dynamicStore, std::forward<Args>(args)...);
        } else if (LowLevelFunctionEntryCallType::Delete != callType) {
          CODI_EXCEPTION("Low level function is not set for token '%d'.", (int)magicNumber);
        }

        curOtherDynamicDataPos = dynamicStore.getPosition();
        curOtherFixedDataPos = fixedStore.getPosition();
      }

    public:

      /// @copydoc LowLevelFunctionTapeInterface::handleLowLevelFunction()
      CODI_INLINE TemporaryMemoryAllocator& getTemporaryMemoryAllocator() {
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
          ExternalFunctionLowLevelEntryMapper<Impl, Real, Identifier>::store(cast(), 0, extFunc);
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
        return otherDynamicData.getPosition();
      }

      /// \copydoc codi::PositionalEvaluationTapeInterface::getZeroPosition()
      Position getZeroPosition() const {
        return otherDynamicData.getZeroPosition();
      }

      /// @}

    protected:

      /// Delete all external function data up to `pos`.
      void deleteExternalFunctionUserData(Position const& pos) {
        // Clear external function data.
        auto deleteFunc =
            [this](
                /* data from other dynamic data vector */
                size_t& curOtherDynamicDataPos, size_t const& endOtherDynamicDataPos, char const* const otherDynamicPtr,
                /* data from other fixed data vector */
                size_t& curOtherFixedDataPos, size_t const& endOtherFixedDataPos, char const* const otherFixedPtr) {
              CODI_UNUSED(endOtherDynamicDataPos);

              while (curOtherFixedDataPos > endOtherFixedDataPos) {
                handleLowLevelFunction<LowLevelFunctionEntryCallType::Delete>(cast(), ByteDataStore::Direction::Reverse,
                                                                              curOtherDynamicDataPos, otherDynamicPtr,
                                                                              curOtherFixedDataPos, otherFixedPtr);
              }
            };

        otherDynamicData.template evaluateReverse<1>(cast().getPosition(), pos, deleteFunc);
      }

    public:

      /// @{

      /// \copydoc codi::PositionalEvaluationTapeInterface::resetTo(Position const&, bool, AdjointsManagement)
      CODI_INLINE void resetTo(Position const& pos, bool resetAdjoints = true,
                               AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        EventSystem<Impl>::notifyTapeResetListeners(cast(), pos, EventHints::Reset::To, resetAdjoints);

        if (resetAdjoints) {
          Impl& impl = cast();
          impl.clearAdjoints(impl.getPosition(), pos, adjointsManagement);
        }

        deleteExternalFunctionUserData(pos);

        otherDynamicData.resetTo(pos);
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
        otherFixedData.setNested(nested);
        otherDynamicData.setNested(&otherFixedData);
      }

      /// @}
  };

  template<typename ImplTapeTypes, typename Impl>
  std::vector<LowLevelFunctionEntry<Impl, typename ImplTapeTypes::Real, typename ImplTapeTypes::Identifier>>*
      CommonTapeImplementation<ImplTapeTypes, Impl>::lowLevelFunctionLookup = nullptr;
}
