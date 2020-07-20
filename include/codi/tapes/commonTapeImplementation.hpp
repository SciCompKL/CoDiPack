#pragma once

#include <algorithm>
#include <type_traits>

#include "../aux/fileIo.hpp"
#include "../aux/macros.h"
#include "../config.h"
#include "aux/vectorAccessInterface.hpp"
#include "data/dataInterface.hpp"
#include "data/chunkVector.hpp"
#include "data/position.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/fullTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct TapeTypesInterface {
    public:

      using Real = ANY;
      using Gradient = ANY;
      using Identifier = ANY;

      using NestedVector = DataInterface<>;
  };

  template<typename _NestedVector>
  struct CommonTapeTypes {
    public:

      using NestedVector = DECLARE_DEFAULT(_NestedVector, TEMPLATE(DataInterface<>));

      using NestedPosition = typename NestedVector::Position;
      using ExternalFunctionChunk = Chunk2<ExternalFunction, NestedPosition>;
      using ExternalFunctionVector = ChunkVector<ExternalFunctionChunk, NestedVector>;
      using Position = typename ExternalFunctionVector::Position;
  };

  template<typename _ImplTapeTypes, typename _Impl>
  struct CommonTapeImplementation :
      public FullTapeInterface<
          typename _ImplTapeTypes::Real,
          typename _ImplTapeTypes::Gradient,
          typename _ImplTapeTypes::Identifier,
          typename CommonTapeTypes<typename _ImplTapeTypes::NestedVector>::Position>
  {
    public:

      using ImplTapeTypes = DECLARE_DEFAULT(_ImplTapeTypes, TapeTypesInterface);
      using Impl = DECLARE_DEFAULT(_Impl, TEMPLATE(FullTapeInterface<double, double, int, EmptyPosition>));

      using Real = typename ImplTapeTypes::Real;
      using Gradient = typename ImplTapeTypes::Gradient;
      using Identifier = typename ImplTapeTypes::Identifier;
      using NestedVector = typename ImplTapeTypes::NestedVector;
      using NestedPosition = typename NestedVector::Position;

      using ExternalFunctionVector = typename CommonTapeTypes<NestedVector>::ExternalFunctionVector;
      using Position = typename ExternalFunctionVector::Position;

    protected:

      bool active;
      std::set<TapeParameters> options;

      ExternalFunctionVector externalFunctionVector;

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
        externalFunctionVector(Config::SmallChunkSize)
      {
        options.insert(TapeParameters::ExternalFunctionsSize);
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
        externalFunctionVector.addToTapeValues(values);

        return values;
      }

      CODI_INLINE void reset(bool resetAdjoints = true) {
        if(resetAdjoints) {
          cast().clearAdjoints();
        }

        deleteExternalFunctionData(cast().getZeroPosition());

        externalFunctionVector.reset();
      }

      // clearAdjoints and reset(Position) are not implemented

      /*******************************************************************************
       * Section: Function from DataManagementTapeInterface
       *
       */

      void swap(Impl& other) {
        std::swap(active, other.active);

        externalFunctionVector.swap(other.externalFunctionVector);
      }

      void resetHard() {
        Impl& impl = cast();

        impl.reset();
        impl.deleteAdjointVector();

        externalFunctionVector.resetHard();
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

        externalFunctionVector.forEachChunk(writeFunction, true, io);
      }

      void readFromFile(const std::string& filename) {
        FileIo io(filename, false);

        externalFunctionVector.forEachChunk(readFunction, true, io);
      }

      void deleteData() {
        externalFunctionVector.forEachChunk(deleteFunction, true);
      }

      std::set<TapeParameters> const& getAvailableOptions() const {
        return options;
      }

      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::ExternalFunctionsSize:
            return externalFunctionVector.getDataSize();
            break;
          default:
            CODI_EXCEPTION("Tried to get undefined parameter for tape.");
            break;
        }
      }

      bool hasParameter(TapeParameters parameter) const {
        return options.cend() != options.find(parameter);
      }

      void setParameter(TapeParameters parameter, size_t value) {
        switch (parameter) {
          case TapeParameters::ExternalFunctionsSize:
            externalFunctionVector.resize(value);
            break;
          default:
            CODI_EXCEPTION("Tried to set undefined parameter for tape.");
            break;
        }
      }

      // deleteAdjointVector is not implemented


      /*******************************************************************************
       * Section: Function from ExternalFunctionTapeInterface
       *
       */

      void pushExternalFunction(ExternalFunction const& extFunc) {
        if(cast().isActive()) {
          externalFunctionVector.reserveItems(1);
          externalFunctionVector.pushData(extFunc, externalFunctionVector.getPosition().inner); // TODO: Add getInner zum Interface?
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
        return externalFunctionVector.getPosition();
      }

      Position getZeroPosition() const {
        return externalFunctionVector.getZeroPosition();
      }

    protected:

      void deleteExternalFunctionData(Position const& pos) {
        // clear external function data
        auto deleteFunc = [this] (ExternalFunction* extFunc, NestedPosition const* endInnerPos) {
          CODI_UNUSED(endInnerPos);

          /* we just need to call the delete function */
          extFunc->deleteData(this);
        };

        externalFunctionVector.forEachReverse(cast().getPosition(), pos, deleteFunc);
      }

    public:

      CODI_INLINE void resetTo(Position const& pos) {

        Impl& impl = cast();
        impl.clearAdjoints(impl.getPosition(), pos);

        deleteExternalFunctionData(pos);

        externalFunctionVector.resetTo(pos);
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

      void init(NestedVector* nested) {
        externalFunctionVector.setNested(nested);
      }

      /*******************************************************************************
       * Section: Helper function for accessing and evaluating the external function vector
       *
       */

      template<typename Function, typename Obj, typename ... Args>
      CODI_INLINE void internalEvaluateExtFuncPrimal(const Position& start, const Position &end,
                                 const Function& func, Obj& obj,
                                 VectorAccessInterface<Real, Identifier>* vectorAccess,
                                 Args&&... args){

        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&] (ExternalFunction* extFunc, const NestedPosition* endInnerPos) {

          (obj.*func)(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          extFunc->evaluatePrimal(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;

        };
        externalFunctionVector.forEachForward(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        (obj.*func)(curInnerPos, end.inner, std::forward<Args>(args)...);
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
        externalFunctionVector.forEachReverse(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        func(curInnerPos, end.inner, std::forward<Args>(args)...);
      }

      template<typename Function, typename Obj, typename ... Args>
      CODI_INLINE void internalEvaluateExtFuncForward(const Position& start, const Position &end,
                                 const Function& func, Obj& obj,
                                 VectorAccessInterface<Real, Identifier>* vectorAccess,
                                 Args&&... args){

        NestedPosition curInnerPos = start.inner;
        auto evalFunc = [&] (ExternalFunction* extFunc, const NestedPosition* endInnerPos) {

          (obj.*func)(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

          extFunc->evaluateForward(&cast(), vectorAccess);

          curInnerPos = *endInnerPos;

        };
        externalFunctionVector.forEachForward(start, end, evalFunc);

        // Iterate over the remainder also covers the case if there have been no external functions.
        (obj.*func)(curInnerPos, end.inner, std::forward<Args>(args)...);
      }
  };
}

