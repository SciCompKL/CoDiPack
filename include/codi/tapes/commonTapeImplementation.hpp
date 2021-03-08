#pragma once

#include <algorithm>
#include <type_traits>

#include "../aux/fileIo.hpp"
#include "../aux/macros.hpp"
#include "../config.h"
#include "aux/externalFunction.hpp"
#include "aux/vectorAccessInterface.hpp"
#include "data/dataInterface.hpp"
#include "data/position.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/fullTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Interface definition for required type definitions in the tapes types of a tape.
   */
  struct TapeTypesInterface {
    public:

      using Real = CODI_ANY;        ///< Primal computation type e.g. double
      using Gradient = CODI_ANY;    ///< Gradient computation type e.g. double or Direction
      using Identifier = CODI_ANY;  ///< Identifier for the internal management e.g. int

      /// Declaration for the data vector that is used for the internal storage. See DateInterface implementations.
      template<typename Chunk, typename Nested>
      using Data = DataInterface<Nested>;

      using NestedData = DataInterface<>;  ///< The root vector of the tape implementation on which the
                                           ///< CommonTapeImplementation builds.
  };

  /**
   * @brief Declares all types used in the CommonTapeImplementation.
   *
   * @tparam _TapeTypes  Needs to implement to TapeTypesInterface.
   */
  template<typename _TapeTypes>
  struct CommonTapeTypes {
    public:

      using TapeTypes = CODI_DD(_TapeTypes, TapeTypesInterface);  ///< See CommonTapeTypes

      using NestedData = typename TapeTypes::NestedData;  ///< See TapeTypesInterface.
      template<typename Chunk, typename Nested>
      using Data = typename TapeTypes::template Data<Chunk, Nested>;  ///< See TapeTypesInterface.

      using NestedPosition = typename NestedData::Position;  ///< See TapeTypesInterface.
      using ExternalFunctionChunk =
          Chunk2<ExternalFunctionInternalData, NestedPosition>;  ///< See Data entries for external functions.
      using ExternalFunctionData = Data<ExternalFunctionChunk, NestedData>;  ///< Data vector for external functions.
      using Position = typename ExternalFunctionData::Position;              ///< Global position of the tape.
  };

  /**
   * @brief Implementation of all common tape functionality.
   *
   * This basic implementation provides following functionality:
   *  - External function support with the external functions stored in externalFunctionData
   *  - Tape options gathering
   *  - Activity tracking
   *
   * It also provides functionality that can be implemented with other functions:
   *  - setter and getter methods
   *  - no positional evaluation methods.
   *  - registerOutput
   *  - TapeValues functions.
   *  - reset functionality
   *
   * @tparam _ImplTapeTypes needs to implement TapeTypesInterface;
   * @tparam _Impl Type of the full tape implementation.
   */
  template<typename _ImplTapeTypes, typename _Impl>
  struct CommonTapeImplementation
      : public FullTapeInterface<typename _ImplTapeTypes::Real, typename _ImplTapeTypes::Gradient,
                                 typename _ImplTapeTypes::Identifier,
                                 typename CommonTapeTypes<_ImplTapeTypes>::Position> {
    public:

      using ImplTapeTypes = CODI_DD(_ImplTapeTypes, TapeTypesInterface);  ///< See CommonTapeImplementation.
      using Impl = CODI_DD(
          _Impl, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));  ///< See CommonTapeImplementation.

      using Real = typename ImplTapeTypes::Real;              ///< See TapeTypesInterface.
      using Gradient = typename ImplTapeTypes::Gradient;      ///< See TapeTypesInterface.
      using Identifier = typename ImplTapeTypes::Identifier;  ///< See TapeTypesInterface.
      using NestedData = typename ImplTapeTypes::NestedData;  ///< See TapeTypesInterface.
      using NestedPosition = typename NestedData::Position;   ///< See DataInterface.

      using ExternalFunctionData =
          typename CommonTapeTypes<ImplTapeTypes>::ExternalFunctionData;   ///< See CommonTapeTypes.
      using Position = typename CommonTapeTypes<ImplTapeTypes>::Position;  ///< See TapeTypesInterface.

    protected:

      bool active;                       ///< If tape stores statements or not.
      std::set<TapeParameters> options;  ///< All options.

      ExternalFunctionData externalFunctionData;  ///< Data vector for external function data.

    private:

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

    protected:

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      TapeValues internalGetTapeValues() const;  ///< Create tape values.

      /// @}

    public:

      /// Constructor
      CommonTapeImplementation() : active(false), options(), externalFunctionData(Config::SmallChunkSize) {
        options.insert(TapeParameters::ExternalFunctionsSize);
      }

      /*******************************************************************************/
      /// @name Functions from GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::setGradient()
      void setGradient(Identifier const& identifier, Gradient const& gradient) {
        cast().gradient(identifier) = gradient;
      }

      /// \copydoc codi::GradientAccessTapeInterface::getGradient()
      Gradient const& getGradient(Identifier const& identifier) const {
        return cast().gradient(identifier);
      }

      // gradient functions are not implemented

      /// @}
      /*******************************************************************************/
      /// @name Functions from ReverseTapeInterface
      /// @{

      /// \copydoc codi::ReverseTapeInterface::evaluate()
      void evaluate() {
        Impl& impl = cast();

        impl.evaluate(impl.getPosition(), impl.getZeroPosition());
      }

      /// \copydoc codi::ReverseTapeInterface::registerOutput()
      template<typename Lhs>
      void registerOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        cast().template store<Lhs, Lhs>(value, static_cast<ExpressionInterface<Real, Lhs> const&>(value));
      }

      /// \copydoc codi::ReverseTapeInterface::setActive()
      void setActive() {
        active = true;
      }

      /// \copydoc codi::ReverseTapeInterface::setPassive()
      void setPassive() {
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

        values.addSection("External function entries");
        externalFunctionData.addToTapeValues(values);

        return values;
      }

      /// \copydoc codi::ReverseTapeInterface::reset()
      CODI_INLINE void reset(bool resetAdjoints = true) {
        if (resetAdjoints) {
          cast().clearAdjoints();
        }

        deleteExternalFunctionUserData(cast().getZeroPosition());

        externalFunctionData.reset();
      }

      // clearAdjoints and reset(Position) are not implemented

      /// @}
      /*******************************************************************************/
      /// @name Functions from DataManagementTapeInterface
      /// @{

      /// \copydoc codi::DataManagementTapeInterface::swap()
      void swap(Impl& other) {
        std::swap(active, other.active);

        externalFunctionData.swap(other.externalFunctionData);
      }

      /// \copydoc codi::DataManagementTapeInterface::resetHard()
      void resetHard() {
        Impl& impl = cast();

        impl.reset();
        impl.deleteAdjointVector();

        externalFunctionData.resetHard();
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

        externalFunctionData.forEachChunk(writeFunction, true, io);
      }

      /// \copydoc codi::DataManagementTapeInterface::readFromFile()
      void readFromFile(const std::string& filename) {
        FileIo io(filename, false);

        externalFunctionData.forEachChunk(readFunction, true, io);
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteData()
      void deleteData() {
        externalFunctionData.forEachChunk(deleteFunction, true);
      }

      /// \copydoc codi::DataManagementTapeInterface::getAvailableParameters()
      std::set<TapeParameters> const& getAvailableParameters() const {
        return options;
      }

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      /// <br><br> Implementation: Handles ExternalFunctionsSize
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::ExternalFunctionsSize:
            return externalFunctionData.getDataSize();
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
      /// <br><br> Implementation: Handles ExternalFunctionsSize
      void setParameter(TapeParameters parameter, size_t value) {
        switch (parameter) {
          case TapeParameters::ExternalFunctionsSize:
            externalFunctionData.resize(value);
            break;
          default:
            CODI_EXCEPTION("Tried to set undefined parameter for tape.");
            break;
        }
      }

      // deleteAdjointVector is not implemented

      /// @}
      /*******************************************************************************/
      /// @name Functions from ExternalFunctionTapeInterface
      /// @{

      /// \copydoc codi::ExternalFunctionTapeInterface::pushExternalFunction()
      void pushExternalFunction(ExternalFunction<Impl> const& extFunc) {
        CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive()) {
          externalFunctionData.reserveItems(1);
          externalFunctionData.pushData(extFunc,
                                        externalFunctionData.getPosition().inner);  // TODO: Add getInner zum Interface?
        }
      }

      // registerExternalFunctionOutput is not implemented

      /// @}
      /*******************************************************************************/
      /// @name Functions from ForwardEvaluationTapeInterface
      /// @{

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      void evaluateForward() {
        Impl& impl = cast();

        impl.evaluateForward(impl.getZeroPosition(), impl.getPosition());
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from IdentifierInformationTapeInterface
      /// @{

      /// \copydoc codi::IdentifierInformationTapeInterface::getPassiveIndex()
      Identifier getPassiveIndex() const {
        return IndexManagerInterface<Identifier>::UnusedIndex;
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
        return externalFunctionData.getPosition();
      }

      /// \copydoc codi::PositionalEvaluationTapeInterface::getZeroPosition()
      Position getZeroPosition() const {
        return externalFunctionData.getZeroPosition();
      }

      /// @}

    protected:

      /// Delete all external function data up to `pos`.
      void deleteExternalFunctionUserData(Position const& pos) {
        // clear external function data
        auto deleteFunc = [this](ExternalFunctionInternalData* extFunc, NestedPosition const* endInnerPos) {
          CODI_UNUSED(endInnerPos);

          /* we just need to call the delete function */
          ((ExternalFunction<Impl>*)extFunc)->deleteData(&cast());
        };

        externalFunctionData.forEachReverse(cast().getPosition(), pos, deleteFunc);
      }

    public:

      /// @{

      /// \copydoc codi::PositionalEvaluationTapeInterface::resetTo()
      CODI_INLINE void resetTo(Position const& pos, bool resetAdjoints = true) {
        if (resetAdjoints) {
          Impl& impl = cast();
          impl.clearAdjoints(impl.getPosition(), pos);
        }

        deleteExternalFunctionUserData(pos);

        externalFunctionData.resetTo(pos);
      }

      // clearAdjoints and evaluate are not implemented

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
      void init(NestedData* nested) {
        externalFunctionData.setNested(nested);
      }

      /// Evaluate all external functions from start to end and call `func` for the regions in between.
      template<typename Function, typename... Args>
      CODI_INLINE void internalEvaluateExtFuncPrimal(const Position& start, const Position& end, Function func,
                                                     VectorAccessInterface<Real, Identifier>* vectorAccess,
                                                     Args&&... args) {
        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&](ExternalFunctionInternalData* extFunc, const NestedPosition* endInnerPos) {
          func(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          ((ExternalFunction<Impl>*)extFunc)->evaluatePrimal(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;
        };
        externalFunctionData.forEachForward(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        func(curInnerPos, end.inner, std::forward<Args>(args)...);
      }

      /// Evaluate all external functions from start to end and call `func` for the regions in between.
      template<typename Function, typename... Args>
      CODI_INLINE void internalEvaluateExtFunc(const Position& start, const Position& end, Function func,
                                               VectorAccessInterface<Real, Identifier>* vectorAccess, Args&&... args) {
        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&](ExternalFunctionInternalData* extFunc, const NestedPosition* endInnerPos) {
          func(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          ((ExternalFunction<Impl>*)extFunc)->evaluateReverse(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;
        };
        externalFunctionData.forEachReverse(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        func(curInnerPos, end.inner, std::forward<Args>(args)...);
      }

      /// Evaluate all external functions from start to end and call `func` for the regions in between.
      template<typename Function, typename... Args>
      CODI_INLINE void internalEvaluateExtFuncForward(const Position& start, const Position& end, Function func,
                                                      VectorAccessInterface<Real, Identifier>* vectorAccess,
                                                      Args&&... args) {
        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&](ExternalFunctionInternalData* extFunc, const NestedPosition* endInnerPos) {
          func(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          ((ExternalFunction<Impl>*)extFunc)->evaluateForward(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;
        };
        externalFunctionData.forEachForward(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        func(curInnerPos, end.inner, std::forward<Args>(args)...);
      }

      /// @}
  };
}
