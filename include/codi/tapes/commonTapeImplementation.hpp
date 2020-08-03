#pragma once

#include <algorithm>
#include <type_traits>

#include "../aux/fileIo.hpp"
#include "../aux/macros.hpp"
#include "../config.h"
#include "aux/vectorAccessInterface.hpp"
#include "data/dataInterface.hpp"
#include "data/position.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/fullTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct TapeTypesInterface {
    public:

      using Real = CODI_ANY;
      using Gradient = CODI_ANY;
      using Identifier = CODI_ANY;

      template<typename Chunk, typename Nested>
      using Data = DataInterface<Nested>;

      using NestedData = DataInterface<>;
  };

  template<typename _TapeTypes>
  struct BaseTapeTypes {
    public:

      using TapeTypes = CODI_DECLARE_DEFAULT(_TapeTypes, TapeTypesInterface);

      using NestedData = typename TapeTypes::NestedData;
      template<typename Chunk, typename Nested>
      using Data = typename TapeTypes::template Data<Chunk, Nested>;


      using NestedPosition = typename NestedData::Position;
      using ExternalFunctionChunk = Chunk2<ExternalFunction, NestedPosition>;
      using ExternalFunctionData = Data<ExternalFunctionChunk, NestedData>;
      using Position = typename ExternalFunctionData::Position;
  };

  template<typename _TapeTypes, typename _Impl>
  struct CommonTapeImplementation :
      public FullTapeInterface<
          typename _TapeTypes::Real,
          typename _TapeTypes::Gradient,
          typename _TapeTypes::Identifier,
          typename BaseTapeTypes<_TapeTypes>::Position>
  {
    public:

      using TapeTypes = CODI_DECLARE_DEFAULT(_TapeTypes, TapeTypesInterface);
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(FullTapeInterface<double, double, int, EmptyPosition>));

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using Identifier = typename TapeTypes::Identifier;
      using NestedData = typename TapeTypes::NestedData;
      using NestedPosition = typename NestedData::Position;

      using ExternalFunctionData = typename BaseTapeTypes<TapeTypes>::ExternalFunctionData;
      using Position = typename ExternalFunctionData::Position;

    protected:

      bool active;
      std::set<ConfigurationOption> options;

      ExternalFunctionData externalFunctionData;

    private:

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

    protected:

      /*******************************************************************************
       * Section: Functions expected in the child class
       *
       */

      TapeValues internalGetTapeValues() const;

    public:

      CommonTapeImplementation() :
        active(false),
        options(),
        externalFunctionData(Config::SmallChunkSize)
      {
        options.insert(ConfigurationOption::ExternalFunctionsSize);
      }


      /*******************************************************************************
       * Section: Functions from GradientAccessTapeInterface
       *
       */

      void setGradient(Identifier const& identifier, Gradient const& gradient) {
        cast().gradient(identifier) = gradient;
      }

      Gradient const& getGradient(Identifier const& identifier) const {
        return cast().gradient(identifier);
      }

      // gradient functions are not implemented

      /*******************************************************************************
       * Section: Functions from ReverseTapeInterface
       *
       */

      void evaluate() {
        Impl& impl = cast();

        impl.evaluate(impl.getPosition(), impl.getZeroPosition());
      }

      template<typename Lhs> void registerOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        cast().template store<Lhs, Lhs>(value, static_cast<ExpressionInterface<Real, Lhs> const&>(value));
      }

      void setActive() {active = true; }
      void setPassive() {active = false; }
      bool isActive() const { return active; }

      template<typename Stream = std::ostream> void printStatistics(Stream& out = std::cout) const {
        cast().getTapeValues().formatDefault(out);
      }

      template<typename Stream = std::ostream> void printTableHeader(Stream& out = std::cout) const {
        cast().getTapeValues().formatHeader(out);
      }

      template<typename Stream = std::ostream> void printTableRow(Stream& out = std::cout) const {
        cast().getTapeValues().formatRow(out);
      }

      TapeValues getTapeValues() const {
        TapeValues values = cast().internalGetTapeValues();

        values.addSection("External function entries");
        externalFunctionData.addToTapeValues(values);

        return values;
      }

      CODI_INLINE void reset(bool resetAdjoints = true) {
        if(resetAdjoints) {
          cast().clearAdjoints();
        }

        deleteExternalFunctionUserData(cast().getZeroPosition());

        externalFunctionData.reset();
      }

      // clearAdjoints and reset are not implemented

      /*******************************************************************************
       * Section: Function from DataManagementTapeInterface
       *
       */

      void swap(Impl& other) {
        std::swap(active, other.active);

        externalFunctionData.swap(other.externalFunctionData);
      }

      void resetHard() {
        Impl& impl = cast();

        impl.reset();
        impl.deleteAdjointVector();

        externalFunctionData.resetHard();
      }

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

      void writeToFile(const std::string& filename) {
        FileIo io(filename, true);

        externalFunctionData.forEachChunk(writeFunction, true, io);
      }

      void readFromFile(const std::string& filename) {
        FileIo io(filename, false);

        externalFunctionData.forEachChunk(readFunction, true, io);
      }

      void deleteData() {
        externalFunctionData.forEachChunk(deleteFunction, true);
      }

      std::set<ConfigurationOption> const& getAvailableOptions() const {
        return options;
      }

      size_t getOption(ConfigurationOption option) const {
        switch (option) {
          case ConfigurationOption::ExternalFunctionsSize:
            return externalFunctionData.getDataSize();
            break;
          default:
            CODI_EXCEPTION("Tried to get undefined option for tape.");
            return 0;
            break;
        }
      }

      bool hasOption(ConfigurationOption option) const {
        return options.cend() != options.find(option);
      }

      void setOption(ConfigurationOption option, size_t value) {
        switch (option) {
          case ConfigurationOption::ExternalFunctionsSize:
            externalFunctionData.resize(value);
            break;
          default:
            CODI_EXCEPTION("Tried to set undefined option for tape.");
            break;
        }
      }

      // deleteAdjointVector is not implemented


      /*******************************************************************************
       * Section: Function from ExternalFunctionTapeInterface
       *
       */

      void pushExternalFunction(ExternalFunction const& extFunc) {
        CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive()) {
          externalFunctionData.reserveItems(1);
          externalFunctionData.pushData(extFunc, externalFunctionData.getPosition().inner); // TODO: Add getInner zum Interface?
        }
      }

      // registerExternalFunctionOutput is not implemented

      /*******************************************************************************
       * Section: Function from ForwardEvaluationTapeInterface
       *
       */

      void evaluateForward() {
        Impl& impl = cast();

        impl.evaluateForward(impl.getZeroPosition(), impl.getPosition());
      }


      /*******************************************************************************
       * Section: Function from IdentifierInformationTapeInterface
       *
       */

      Identifier getPassiveIndex() const {
        return IndexManagerInterface<Identifier>::UnusedIndex;
      }

      Identifier getInvalidIndex() const {
        return IndexManagerInterface<Identifier>::InvalidIndex;
      }

      bool isIdentifierActive(Identifier const& index) const {
        return index != cast().getPassiveIndex();
      }

      template<typename Lhs>
      void deactivateValue(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        value = value.getValue();
      }

      /*******************************************************************************
       * Section: Function from PositionalEvaluationTapeInterface
       *
       */

      Position getPosition() const {
        return externalFunctionData.getPosition();
      }

      Position getZeroPosition() const {
        return externalFunctionData.getZeroPosition();
      }

    protected:

      void deleteExternalFunctionUserData(Position const& pos) {
        // clear external function data
        auto deleteFunc = [this] (ExternalFunction* extFunc, NestedPosition const* endInnerPos) {
          CODI_UNUSED(endInnerPos);

          /* we just need to call the delete function */
          extFunc->deleteData(this);
        };

        externalFunctionData.forEachReverse(cast().getPosition(), pos, deleteFunc);
      }

    public:

      CODI_INLINE void resetTo(Position const& pos) {

        Impl& impl = cast();
        impl.clearAdjoints(impl.getPosition(), pos);

        deleteExternalFunctionUserData(pos);

        externalFunctionData.resetTo(pos);
      }

      // clearAdjoints and evaluate are not implemented

      /*******************************************************************************
       * Section: Function from PrimalEvaluationTapeInterface
       *
       */

      void evaluatePrimal() {
        Impl& impl = cast();

        impl.evaluatePrimal(impl.getZeroPosition(), impl.getPosition());
      }

      void setPrimal(Identifier const& identifier, Real const& primal) {
        cast().primal(identifier) = primal;
      }

      Real const& getPrimal(Identifier const& identifier) const {
        return cast().primal(identifier);
      }

    protected:

      void init(NestedData* nested) {
        externalFunctionData.setNested(nested);
      }

      /*******************************************************************************
       * Section: Helper function for accessing and evaluating the external function vector
       *
       */

      template<typename Function, typename ... Args>
      CODI_INLINE void internalEvaluateExtFuncPrimal(const Position& start, const Position &end,
                                 Function func,
                                 VectorAccessInterface<Real, Identifier>* vectorAccess,
                                 Args&&... args){

        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&] (ExternalFunction* extFunc, const NestedPosition* endInnerPos) {

          func(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          extFunc->evaluatePrimal(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;

        };
        externalFunctionData.forEachForward(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        func(curInnerPos, end.inner, std::forward<Args>(args)...);
      }

      template<typename Function, typename ... Args>
      CODI_INLINE void internalEvaluateExtFunc(const Position& start, const Position &end,
                                 Function func,
                                 VectorAccessInterface<Real, Identifier>* vectorAccess,
                                 Args&&... args){

        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&] (ExternalFunction* extFunc, const NestedPosition* endInnerPos) {

          func(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          extFunc->evaluateReverse(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;

        };
        externalFunctionData.forEachReverse(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        func(curInnerPos, end.inner, std::forward<Args>(args)...);
      }

      template<typename Function, typename ... Args>
      CODI_INLINE void internalEvaluateExtFuncForward(const Position& start, const Position &end,
                                 Function func,
                                 VectorAccessInterface<Real, Identifier>* vectorAccess,
                                 Args&&... args){

        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&] (ExternalFunction* extFunc, const NestedPosition* endInnerPos) {

          func(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          extFunc->evaluateForward(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;

        };
        externalFunctionData.forEachForward(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        func(curInnerPos, end.inner, std::forward<Args>(args)...);
      }
  };
}

