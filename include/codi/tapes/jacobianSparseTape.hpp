/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/helpers/forEachLeafLogic.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../expressions/referenceActiveType.hpp"
#include "../misc/macros.hpp"
#include "../misc/mathUtility.hpp"
#include "../misc/memberStore.hpp"
#include "../traits/computationTraits.hpp"
#include "../traits/expressionTraits.hpp"
#include "commonTapeImplementation.hpp"
#include "data/chunk.hpp"
#include "data/chunkedData.hpp"
#include "indices/indexManagerInterface.hpp"
#include "misc/adjointVectorAccess.hpp"
#include "misc/duplicateJacobianRemover.hpp"
#include "misc/localAdjoints.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct SparseIdentifier {
    public:

      unsigned short chunk;
      unsigned short pos;

      CODI_INLINE SparseIdentifier() = default;

      CODI_INLINE SparseIdentifier(unsigned short chunk, unsigned short pos) : chunk(chunk), pos(pos) {}

      CODI_INLINE bool operator!=(SparseIdentifier const& o) const {
        return chunk != o.chunk || pos != o.pos;
      }

      CODI_INLINE bool operator==(SparseIdentifier const& o) const {
        return chunk == o.chunk && pos == o.pos;
      }

      CODI_INLINE bool operator<(SparseIdentifier const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && pos < o.pos);
      }

      CODI_INLINE bool operator<=(SparseIdentifier const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && pos <= o.pos);
      }

      CODI_INLINE bool operator>(SparseIdentifier const& o) const {
        return o < *this;
      }

      CODI_INLINE bool operator>=(SparseIdentifier const& o) const {
        return o <= *this;
      }
  };

  /**
   * @brief Type definitions for the Jacobian tapes.
   *
   * @tparam T_Real          See TapeTypesInterface.
   * @tparam T_Gradient      See TapeTypesInterface.
   * @tparam T_Data          See TapeTypesInterface.
   * @tparam T_Adjoints      Internal implementation of the adjoint variables.
   */
  template<typename T_Real, typename T_Gradient, template<typename, typename> class T_Data>
  struct JacobianSparseTapeTypes : public TapeTypesInterface {
    public:

      using Real = CODI_DD(T_Real, double);                                              ///< See JacobianSparseTapeTypes.
      using Gradient = CODI_DD(T_Gradient, double);                                      ///< See JacobianSparseTapeTypes.
      template<typename Chunk, typename Nested>
      using Data = CODI_DD(CODI_T(T_Data<Chunk, Nested>),
                           CODI_T(DataInterface<Nested>));  ///< See JacobianSparseTapeTypes.

      using Identifier = SparseIdentifier;

      /// Statement chunk is either \<argument size\> (linear management) or \<lhs identifier, argument size\>
      /// (reuse management).
      using StatementChunk = Chunk2<Identifier, Config::ArgumentSize>;
      using StatementData = Data<StatementChunk, EmptyData>;  ///< Statement data vector.

      using JacobianChunk = Chunk2<Real, Identifier>;           ///< Jacobian chunks is \<Jacobian, rhs index\>.
      using JacobianData = Data<JacobianChunk, StatementData>;  ///< Jacobian data vector.

      using NestedData = JacobianData;  ///< See TapeTypesInterface.
  };

  /**
   * @brief Base class for all standard Jacobian tape implementations.
   *
   * This class provides nearly a full implementation of the FullTapeInterface. There are just a few internal methods
   * left which have to be implemented by the final classes. These methods depend significantly on the index management
   * scheme and are performance critical.
   *
   * Tape evaluations are performed in three steps with two wrapper steps beforehand. Each methods calls the next
   * method:
   * - evaluate
   * - internalEvaluate*
   * - internalEvaluate*_Step1_ExtFunc
   * - internalEvaluate*_Step2_DataExtraction
   * - internalEvaluate*_Step3_EvalStatements
   * The placeholder stands for Reverse, Forward, or Primal.
   *
   * @tparam T_TapeTypes has to implement JacobianSparseTapeTypes.
   */
  template<typename T_TapeTypes>
  struct JacobianSparseTape : public CommonTapeImplementation<T_TapeTypes, JacobianSparseTape<T_TapeTypes>> {
    public:

      /// See JacobianSparseTape.
      using TapeTypes = CODI_DD(
          T_TapeTypes, CODI_T(JacobianSparseTapeTypes<double, double, DefaultChunkedData>));
      using Impl = JacobianSparseTape;

      using Base = CommonTapeImplementation<T_TapeTypes, JacobianSparseTape>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                  ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;          ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;      ///< See TapeTypesInterface.

      using StatementData = typename TapeTypes::StatementData;  ///< See JacobianSparseTapeTypes.
      using JacobianData = typename TapeTypes::JacobianData;    ///< See JacobianSparseTapeTypes.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      using NestedPosition = typename JacobianData::Position;  ///< See JacobianSparseTapeTypes.
      using Position = typename Base::Position;                ///< See TapeTypesInterface.

      template<typename Adjoint>
      using VectorAccess =
          AdjointVectorAccess<Real, Identifier, Adjoint>;  ///< Vector access type generated by this tape.

      /// See GradientAccessTapeInterface.
      using typename GradientAccessTapeInterface<Gradient, Identifier>::ResizingPolicy;

      static bool constexpr AllowJacobianOptimization = true;  ///< See InternalStatementRecordingTapeInterface.
      static bool constexpr HasPrimalValues = false;           ///< See PrimalEvaluationTapeInterface.

      static bool constexpr LinearIndexHandling = true;     ///< See IdentifierInformationTapeInterface.
      static bool constexpr RequiresPrimalRestore = false;  ///< See PrimalEvaluationTapeInterface.

    protected:

#if CODI_RemoveDuplicateJacobianArguments
      DuplicateJacobianRemover<Real, Identifier> jacobianSorter;  ///< Encapsulates jacobianData to remove duplicated
                                                                  ///< Jacobians.
#endif

      EmptyData emptyData;
      StatementData statementData;  ///< Data stream for statement specific data.
      JacobianData jacobianData;    ///< Data stream for argument specific data.

      Gradient emptyAdjoint;
      std::map<Identifier, Gradient> adjoints;  ///< Evaluation vector for AD.

      MemberStore<LinearIndexManager<int>, Impl, false> indexManager;  // TODO: Required for commonTapeImplementation. Refactor there.

    private:

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

    protected:

      CODI_INLINE Identifier createJacobianPosition() {
       typename JacobianData::Position pos = jacobianData.getPosition();
       return Identifier(pos.chunk, pos.data);
     }

     static CODI_INLINE Identifier createStmtPosition(typename StatementData::Position const& pos) {
       return Identifier(pos.chunk, pos.data);
     }

      CODI_INLINE Identifier createStmtPosition() {
       return createStmtPosition(statementData.getPosition());
     }

    public:

      /// Constructor
      JacobianSparseTape()
          : Base(),
#if CODI_RemoveDuplicateJacobianArguments
            jacobianSorter(),
#endif
            statementData(Config::ChunkSize),
            jacobianData(std::max(Config::ChunkSize, Config::MaxArgumentSize)),  // Chunk must be large enough to store
                                                                                 // data for all arguments of one
                                                                                 // statement.
            emptyAdjoint(),
            adjoints(),
            indexManager(0)
      {
        statementData.setNested(&emptyData);
        jacobianData.setNested(&statementData);

        Base::init(&jacobianData);

        Base::options.insert(TapeParameters::AdjointSize);
        Base::options.insert(TapeParameters::JacobianSize);
        Base::options.insert(TapeParameters::StatementSize);
      }

      /*******************************************************************************/
      /// @name Functions from GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&, ResizingPolicy)
      CODI_INLINE Gradient& gradient(Identifier const& identifier,
                                     ResizingPolicy resizingPolicy = ResizingPolicy::CheckAndAdapt) {
        return adjoints[identifier];
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&) const
      CODI_INLINE Gradient const& gradient(Identifier const& identifier) const {
        auto entry = adjoints.find(identifier);
        if (entry == adjoints.cend()) {
          return emptyAdjoint;
        } else {
          return entry->second;
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from InternalStatementRecordingTapeInterface
      /// @{

      /// \copydoc codi::InternalStatementRecordingTapeInterface::initIdentifier()
      template<typename Real>
      CODI_INLINE void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);

        identifier = Identifier();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::destroyIdentifier()
      template<typename Real>
      CODI_INLINE void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);
        CODI_UNUSED(identifier);
      }

      /// @}

    protected:

      /// Pushes Jacobians and indices to the tape.
      struct PushJacobianLogic : public JacobianComputationLogic<PushJacobianLogic> {
        public:
          /// General implementation. Checks for invalid and passive values/Jacobians.
          template<typename Node, typename Jacobian, typename DataVector>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Jacobian jacobianExpr, DataVector& dataVector) {
            Real jacobian = ComputationTraits::adjointConversion<Real>(jacobianExpr);

            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, Identifier() != node.getIdentifier())) {
              if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) {
                if (CODI_ENABLE_CHECK(Config::CheckJacobianIsZero, !RealTraits::isTotalZero(jacobian))) {
                  dataVector.pushData(jacobian, node.getIdentifier());
                }
              }
            }
          }

          /// Specialization for ReferenceActiveType nodes. Delays Jacobian push.
          template<typename Type, typename Jacobian, typename DataVector>
          CODI_INLINE void handleJacobianOnActive(ReferenceActiveType<Type> const& node, Jacobian jacobianExpr,
                                                  DataVector& dataVector) {
            CODI_UNUSED(dataVector);

            Real jacobian = ComputationTraits::adjointConversion<Real>(jacobianExpr);

            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) {
              // Do a delayed push for these leaf nodes, accumulate the jacobian in the local member.
              node.jacobian += jacobian;
            }
          }
      };

      /// Pushes all delayed Jacobians.
      struct PushDelayedJacobianLogic : public ForEachLeafLogic<PushDelayedJacobianLogic> {
        public:

          /// Specialization for ReferenceActiveType nodes. Pushes the delayed Jacobian.
          template<typename Type, typename DataVector>
          CODI_INLINE void handleActive(ReferenceActiveType<Type> const& node, DataVector& dataVector) {
            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier())) {
              if (CODI_ENABLE_CHECK(Config::CheckJacobianIsZero, !RealTraits::isTotalZero(node.jacobian))) {
                dataVector.pushData(node.jacobian, node.getIdentifier());

                // Reset the jacobian here such that it is not pushed multiple times and ready for the next store.
                node.jacobian = Real();
              }
            }
          }

          using ForEachLeafLogic<PushDelayedJacobianLogic>::handleActive;
      };

      /// Push Jacobians and delayed Jacobians to the tape.
      template<typename Rhs>
      CODI_INLINE void pushJacobians(ExpressionInterface<Real, Rhs> const& rhs) {
        PushJacobianLogic pushJacobianLogic;
        PushDelayedJacobianLogic pushDelayedJacobianLogic;

#if CODI_RemoveDuplicateJacobianArguments
        auto& insertVector = jacobianSorter;
#else
        auto& insertVector = jacobianData;
#endif

        pushJacobianLogic.eval(rhs.cast(), Real(1.0), insertVector);
        pushDelayedJacobianLogic.eval(rhs.cast(), insertVector);

#if CODI_RemoveDuplicateJacobianArguments
        jacobianSorter.storeData(jacobianData);
#endif
      }

    public:

      /// @{

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store()
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                             ExpressionInterface<Real, Rhs> const& rhs) {
        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive())) {
          size_t constexpr MaxArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;

          codiAssert(MaxArgs < Config::MaxArgumentSize);

          statementData.reserveItems(1);
          typename JacobianData::InternalPosHandle jacobianStart = jacobianData.reserveItems(MaxArgs);
          Identifier jacobianPos = createJacobianPosition();

          pushJacobians(rhs);

          size_t numberOfArguments = jacobianData.getPushedDataCount(jacobianStart);
          if (CODI_ENABLE_CHECK(Config::CheckEmptyStatements, 0 != numberOfArguments)) {
            statementData.pushData(jacobianPos, (Config::ArgumentSize)numberOfArguments);
            lhs.cast().getIdentifier() = createStmtPosition();

            if (Config::StatementEvents) {
              Real* jacobians;
              Identifier* rhsIdentifiers;
              jacobianData.getDataPointers(jacobianStart, jacobians, rhsIdentifiers);

              EventSystem<Impl>::notifyStatementStoreOnTapeListeners(cast(), lhs.cast().getIdentifier(),
                                                                     rhs.cast().getValue(), numberOfArguments,
                                                                     rhsIdentifiers, jacobians);
            }
          } else {
            lhs.cast().getIdentifier() = Identifier();
          }
        } else {
          lhs.cast().getIdentifier() = Identifier();
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Optimization for copy statements.
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                             LhsExpressionInterface<Real, Gradient, Impl, Rhs> const& rhs) {
        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive())) {
          if (!Config::CopyOptimization) {
            store<Lhs, Rhs>(lhs, static_cast<ExpressionInterface<Real, Rhs> const&>(rhs));
            return;
          } else {
            lhs.cast().getIdentifier() = rhs.cast().getIdentifier();
          }
        } else {
          lhs.cast().getIdentifier() = Identifier();
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Specialization for passive assignments.
      template<typename Lhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, Real const& rhs) {
        lhs.cast().getIdentifier() = Identifier();

        lhs.cast().value() = rhs;
      }

      /// @}
      /*******************************************************************************
       * Protected helper function for ReverseTapeInterface
       */

    protected:

      /// Add a new input to the tape.
      template<typename Lhs>
      CODI_INLINE void internalRegisterInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        statementData.reserveItems(1);
        statementData.pushData(createJacobianPosition(), Config::StatementInputTag);

        value.cast().getIdentifier() = createStmtPosition();
      }

    public:

      /*******************************************************************************/
      /// @name Functions from ReverseTapeInterface
      /// @{

      /// \copydoc codi::ReverseTapeInterface::registerInput()
      template<typename Lhs>
      CODI_INLINE void registerInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        internalRegisterInput(value);
        EventSystem<Impl>::notifyTapeRegisterInputListeners(cast(), value.cast().value(), value.cast().getIdentifier());
      }

      /// \copydoc codi::ReverseTapeInterface::clearAdjoints()
      CODI_INLINE void clearAdjoints() {
        adjoints.clear();
      }

      /// @}

    protected:

      /// Adds data from all streams, the size of the adjoint vector and index manager data.
      CODI_INLINE TapeValues internalGetTapeValues() const {
        std::string name = "CoDi Tape Statistics ( JacobianSparseTape )";
        TapeValues values = TapeValues(name);

        size_t nAdjoints = adjoints.size();
        double memoryAdjoints = static_cast<double>(nAdjoints) * static_cast<double>(sizeof(Gradient));

        values.addSection("Adjoint vector");
        values.addUnsignedLongEntry("Number of adjoints", nAdjoints);
        values.addDoubleEntry("Memory allocated", memoryAdjoints, true, true);

        values.addSection("Statement entries");
        statementData.addToTapeValues(values);
        values.addSection("Jacobian entries");
        jacobianData.addToTapeValues(values);

        return values;
      }

      /******************************************************************************
       * Protected helper function for CustomAdjointVectorEvaluationTapeInterface
       */

    protected:

      /// Performs the AD \ref sec_reverseAD "reverse" equation for a statement.
      template<typename Adjoint>
      CODI_INLINE static void incrementAdjoints(std::map<Identifier, Adjoint>& adjointMap, Adjoint const& lhsAdjoint,
                                                JacobianData& jacobianData,
                                                Identifier jacobianPosition,
                                                Config::ArgumentSize const& numberOfArguments) {
        Real* rhsJacobians;
        Identifier* rhsIdentifiers;
        jacobianData.getDataPointersWithChunk(jacobianPosition.chunk, jacobianPosition.pos, rhsJacobians, rhsIdentifiers);

        if (CODI_ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !RealTraits::isTotalZero(lhsAdjoint))) {
          for(Config::ArgumentSize pos = 0; pos < numberOfArguments; pos += 1) {
            adjointMap[rhsIdentifiers[pos]] += rhsJacobians[pos] * lhsAdjoint;
          }
        }
      }

      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateSparseReverse(
          std::map<Identifier, Adjoint>& seeding,
          StatementData& statementData,
          JacobianData& jacobianData,
          Identifier const& endStmtPos) {

        Identifier* rhsPos;
        Config::ArgumentSize* numberOfJacobians;

        typename std::map<Identifier, Adjoint>::reverse_iterator entry = seeding.rbegin();
        while (entry != seeding.rend() && entry->first >= endStmtPos) {

          Identifier stmtPos = entry->first;
          statementData.getDataPointersWithChunk(stmtPos.chunk, stmtPos.pos - 1, rhsPos, numberOfJacobians); // TODO: Fix of by one error. Identifier() needs to initialize not to zero.
          Config::ArgumentSize const argsSize = *numberOfJacobians;

          if (Config::StatementInputTag != argsSize) {
            Adjoint curAdjoint = entry->second;
            seeding.erase(std::next(entry).base()); // std::map::erase does not support reverse iterators.

            incrementAdjoints(seeding, curAdjoint, jacobianData, *rhsPos, argsSize);

            entry = seeding.rbegin();
          } else {
            ++entry; // TODO: Hack for infinite loop when inputs are registered. Only works if inputs are reqistered at
                     // the start or it will make the evaluation quite slow.
          }
        }
      }


      /// Start for reverse evaluation between external function.
      template<typename Adjoint>
      CODI_NO_INLINE static void internalEvaluateReverse_Step2_DataExtraction(NestedPosition const& start,
                                                                              NestedPosition const& end, Impl& tape,
                                                                              std::map<Identifier, Adjoint>& data,
                                                                              StatementData& statementData,
                                                                              JacobianData& jacobianData) {
        internalEvaluateSparseReverse(data, statementData, jacobianData, createStmtPosition(end.inner));
      }

      // Start for forward evaluation between external function.
      template<typename Adjoint>
      CODI_NO_INLINE static void internalEvaluateForward_Step2_DataExtraction(NestedPosition const& start,
                                                                              NestedPosition const& end, Impl& tape,
                                                                              Adjoint* data,
                                                                              JacobianData& jacobianData) {
        // TODO: Implement.
      }

    public:

      /// @name Functions from CustomAdjointVectorEvaluationTapeInterface
      /// @{

      using Base::evaluate;

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename Adjoint>
      CODI_NO_INLINE void evaluate(Position const& start, Position const& end, std::map<Identifier, Adjoint>& data) {
        VectorAccess<Adjoint> adjointWrapper(NULL); // TODO: Implement.

        EventSystem<Impl>::notifyTapeEvaluateListeners(
            cast(), start, end, &adjointWrapper, EventHints::EvaluationKind::Reverse, EventHints::Endpoint::Begin);

        Base::internalEvaluateReverse_Step1_ExtFunc(
            start, end, JacobianSparseTape::template internalEvaluateReverse_Step2_DataExtraction<Adjoint>,
            &adjointWrapper, cast(), data, statementData, jacobianData);

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &adjointWrapper,
                                                       EventHints::EvaluationKind::Reverse, EventHints::Endpoint::End);
      }

      template<typename Adjoint>
      CODI_NO_INLINE void evaluateSparse(
          Position const& start,
          Position const& end,
          std::map<Identifier, Adjoint>& seeding) {

        internalEvaluateReverse_Step2_DataExtraction(start.inner, end.inner, *this, seeding, statementData, jacobianData);
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename Adjoint>
      CODI_NO_INLINE void evaluateForward(Position const& start, Position const& end, Adjoint* data) {
        VectorAccess<Adjoint> adjointWrapper(data);

        EventSystem<Impl>::notifyTapeEvaluateListeners(
            cast(), start, end, &adjointWrapper, EventHints::EvaluationKind::Forward, EventHints::Endpoint::Begin);

        Base::internalEvaluateForward_Step1_ExtFunc(
            start, end, JacobianSparseTape::template internalEvaluateForward_Step2_DataExtraction<Adjoint>,
            &adjointWrapper, cast(), data, jacobianData);

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &adjointWrapper,
                                                       EventHints::EvaluationKind::Forward, EventHints::Endpoint::End);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from DataManagementTapeInterface
      /// @{

      /// \copydoc codi::DataManagementTapeInterface::swap()
      CODI_INLINE void swap(Impl& other) {
        // Index manager does not need to be swapped, it is either static or swapped with the vector data.
        // Vectors are swapped recursively in the base class.

        adjoints.swap(other.adjoints);

        Base::swap(other);
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteAdjointVector()
      void deleteAdjointVector() {
        adjoints.clear();
      }

      /// \copydoc codi::DataManagementTapeInterface::resizeAdjointVector()
      void resizeAdjointVector() {
        // Do nothing.
      }

      /// \copydoc codi::DataManagementTapeInterface::beginUseAdjointVector()
      void beginUseAdjointVector() {
        // Do nothing.
      }

      /// \copydoc codi::DataManagementTapeInterface::endUseAdjointVector()
      void endUseAdjointVector() {
        // Do nothing.
      }

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::AdjointSize:
            return adjoints.size();
            break;
          case TapeParameters::JacobianSize:
            return jacobianData.getDataSize();
            break;
          case TapeParameters::StatementSize:
            return statementData.getDataSize();
          default:
            return Base::getParameter(parameter);
            break;
        }
      }

      /// \copydoc codi::DataManagementTapeInterface::setParameter()
      void setParameter(TapeParameters parameter, size_t value) {
        switch (parameter) {
          case TapeParameters::AdjointSize:
            CODI_EXCEPTION("Tried to set a get only parameter.");
            break;
          case TapeParameters::JacobianSize:
            jacobianData.resize(value);
            break;
          case TapeParameters::LargestIdentifier:
            CODI_EXCEPTION("Tried to set a get only parameter.");
            break;
          case TapeParameters::StatementSize:
            statementData.resize(value);
            break;
          default:
            Base::setParameter(parameter, value);
            break;
        }
      }

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccess()
      VectorAccess<Gradient>* createVectorAccess() {
        return NULL; // TODO: Implement
      }

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccessCustomAdjoints()
      template<typename Adjoint>
      VectorAccess<Adjoint>* createVectorAccessCustomAdjoints(Adjoint* data) {
        return NULL;  // TODO: Implement
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteVectorAccess()
      void deleteVectorAccess(VectorAccessInterface<Real, Identifier>* access) {
        CODI_UNUSED(access);
         // TODO: Implement
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ExternalFunctionTapeInterface
      /// @{

      /// \copydoc codi::ExternalFunctionTapeInterface::registerExternalFunctionOutput()
      template<typename Lhs>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        internalRegisterInput(value);

        return Real();
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ForwardEvaluationTapeInterface
      /// @{

      using Base::evaluateForward;

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      void evaluateForward(Position const& start, Position const& end) {
        cast().evaluateForward(start, end, adjoints);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ManualStatementPushTapeInterface
      /// @{

      /// \copydoc codi::ManualStatementPushTapeInterface::pushJacobianManual()
      void pushJacobianManual(Real const& jacobian, Real const& value, Identifier const& index) {
        CODI_UNUSED(value);

        cast().incrementManualPushCounter();

        jacobianData.pushData(jacobian, index);

        if (Config::StatementEvents) {
          if (this->manualPushCounter == this->manualPushGoal) {
            // emit statement event
            Real* jacobians;
            Identifier* rhsIdentifiers;
            jacobianData.getDataPointers(jacobianData.reserveItems(0), jacobians, rhsIdentifiers);
            jacobians -= this->manualPushGoal;
            rhsIdentifiers -= this->manualPushGoal;

            EventSystem<Impl>::notifyStatementStoreOnTapeListeners(cast(), this->manualPushLhsIdentifier,
                                                                   this->manualPushLhsValue, this->manualPushGoal,
                                                                   rhsIdentifiers, jacobians);
          }
        }
      }

      /// \copydoc codi::ManualStatementPushTapeInterface::storeManual()
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        codiAssert(size < Config::MaxArgumentSize);

        statementData.reserveItems(1);
        jacobianData.reserveItems(size);

        statementData.pushData(createJacobianPosition(), (Config::ArgumentSize)size);
        lhsIndex = createStmtPosition();

        cast().initializeManualPushData(lhsValue, lhsIndex, size);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PositionalEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PositionalEvaluationTapeInterface::evaluate()
      CODI_INLINE void evaluate(Position const& start, Position const& end) {
        evaluate(start, end, adjoints);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PreaccumulationEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateKeepState()
      void evaluateKeepState(Position const& start, Position const& end) {
        evaluate(start, end);
      }

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateForwardKeepState()
      void evaluateForwardKeepState(Position const& start, Position const& end) {
        evaluateForward(start, end);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PrimalEvaluationTapeInterface
      /// @{

      using Base::evaluatePrimal;

      /// Not implemented, raises an exception.
      void evaluatePrimal(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);

        CODI_EXCEPTION("Accessing primal evaluation of an Jacobian tape.");
      }

      /// Not implemented, raises an exception.
      Real& primal(Identifier const& identifier) {
        CODI_UNUSED(identifier);

        CODI_EXCEPTION("Accessing primal vector of an Jacobian tape.");

        static Real temp;
        return temp;
      }

      /// Not implemented, raises an exception.
      Real primal(Identifier const& identifier) const {
        CODI_UNUSED(identifier);

        CODI_EXCEPTION("Accessing primal vector of an Jacobian tape.");

        return Real();
      }

      /// @}
  };
}
