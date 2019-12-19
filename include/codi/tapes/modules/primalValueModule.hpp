/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "../reverseTapeInterface.hpp"
#include "../../configure.h"
#include "../../tapeTypes.hpp"
#include "../../tools/tapeValues.hpp"
#include "../../typeFunctions.hpp"
#include "../primalTapeExpressions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
    * The module defines the structures stmtVector, indexVector, constantValueVector, primals, primalsSize, primalsIncr.
    * The module defines the static structures InputHandle, CopyHandle and PreaccHandles,
    * The module defines the types PrimalChildVector, StatmentVector, StatementChunk, StmtPosition, IndexVector, IndexChunk,
    * IndexPosition, ConstantValueVector, ConstantValueChunk, ConstantValuePosition.
    *
    * It defines the methods store(Expr), store(const), store(User), pushJacobi, printPrimalValueStatistics from the TapeInterface and ReverseTapeInterface.
    *
    * It defines the methods resizePrimals, checkPrimalsSize, evaluateHandle, evaluateConstantValues, getUsedStatementsSize,
    * getUsedDataEntiresSize, getUsedConstantDataSize, setConstantDataSize, swapPrimalValueModule as interface functions for the including class.
    *
    * It defines the static methods inputHandleFunc, copyHandleFunc, preaccHandleFunc as interface functions for the tape.
    *
    * @tparam TapeTypes  All the types for the tape. Including the calculation type and the vector types.
    * @tparam      Tape  The full tape implementation
    */
  template<typename TapeTypes, typename Tape>
  struct PrimalValueModule: public virtual ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, Tape, typename TapeTypes::Position > {

		private:

    // ----------------------------------------------------------------------
    // All definitions of the module
    // ----------------------------------------------------------------------

      CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

      /** @brief Forward of the type for the handle factory */
      typedef typename TapeTypes::HandleFactory HandleFactory;

      /** @brief The vector for the statement data. */
      typedef typename TapeTypes::StatementVector StatementVector;
      /** @brief The data for the statments */
      typedef typename StatementVector::ChunkType StatementChunk;
      /** @brief The position type of the statment data. */
      typedef typename StatementVector::Position StmtPosition;

      /** @brief The vector for the index data. */
      typedef typename TapeTypes::IndexVector IndexVector;
      /** @brief The data for the indices*/
      typedef typename IndexVector::ChunkType IndexChunk;
      /** @brief The position type of the index data. */
      typedef typename IndexVector::Position IndexPosition;

      /** @brief The vector for the passive data. */
      typedef typename TapeTypes::PassiveValueVector PassiveValueVector;
      /** @brief The data for the passives*/
      typedef typename PassiveValueVector::ChunkType PassiveValueChunk;
      /** @brief The position type of the passive data. */
      typedef typename PassiveValueVector::Position PassiveValuePosition;

      /** @brief The vector for the constant data. */
      typedef typename TapeTypes::ConstantValueVector ConstantValueVector;
      /** @brief The data for the constants*/
      typedef typename ConstantValueVector::ChunkType ConstantValueChunk;
      /** @brief The position type of the constant data. */
      typedef typename ConstantValueVector::Position ConstantValuePosition;

      /** @brief The child vector for the primal data vector. */
      typedef typename StatementVector::NestedVectorType PrimalChildVector;


      /**
       * @brief Cast this class to the full.
       *
       * The full type is able to access all functions from the tape and other modules.
       *
       * @return  The full tape implemenation.
       */
      Tape& cast() {
        return *static_cast<Tape*>(this);
      }

  protected:

      /** @brief The data for the each statements. */
      StatementVector stmtVector;

      /** @brief The data for the indices of each statements. */
      IndexVector indexVector;

      /** @brief The data for the passive values of each statements. */
      PassiveValueVector passiveValueVector;

      /** @brief The data for the constant values of each statements. */
      ConstantValueVector constantValueVector;

      /** @brief The static array of handles for the preaccumulation. */
      const static typename TapeTypes::Handle preaccHandles[MaxStatementIntSize];

      /** @brief The vector with the primal value for each statement. */
      Real* primals;

      /** @brief The current size of the primal vector. */
      Index primalsSize;

      /**
       * @brief The increment of the size of the primal vector.
       *
       * The primals are incremented, when the indices are bigger than the current size.
       */
      Index primalsIncr;

    public:

      /**
       * @brief Default constructor
       */
      PrimalValueModule() :
        stmtVector(DefaultChunkSize),
        indexVector(DefaultChunkSize),
        passiveValueVector(DefaultChunkSize),
        constantValueVector(DefaultChunkSize),
        primals(NULL),
        primalsSize(0),
        primalsIncr(DefaultSmallChunkSize)
      {}

    protected:


      /**
       * @brief Adds information about primal value tape.
       *
       * Adds the number of chunks, the total number of entries, the
       * allocated memory and the used memory for the constant vector, the statement
       * vector and the index vector.
       *
       * @param[in,out] values  The information is added to the values
       */
      void addPrimalValueValues(TapeValues& values) const {

        size_t totalPrimal   = primalsSize;
        size_t sizePrimalEntry = sizeof(Real);
        double memoryAllocPrimal = (double)totalPrimal*(double)sizePrimalEntry* BYTE_TO_MB;

        values.addSection("Primal vector");
        values.addData("Total number", totalPrimal);
        values.addData("Memory allocated", memoryAllocPrimal, true, true);

        values.addSection("Statements");
        values.addStreamData(stmtVector);

        values.addSection("Index entries");
        values.addStreamData(indexVector);

        values.addSection("Passive data entries");
        values.addStreamData(constantValueVector);
      }

      /**
       * @brief Initialize the JacobiModule.
       *
       * Called after all members of the tape have been initialized.
       *
       * @param[in,out] childVector  The child vector for the Jacobi vector.
       */
      void initPrimalValueModule(PrimalChildVector* childVector) {
        stmtVector.setNested(childVector);
        indexVector.setNested(&stmtVector);
        passiveValueVector.setNested(&indexVector);
        constantValueVector.setNested(&passiveValueVector);
      }

      /**
       * @brief Definition of the adjoint vector type in the evaluation process.
       *
       * If CODI_EnableVariableAdjointInterfaceInPrimalTapes is set the type is a general one
       * which allows for arbitrary adjoint vectors but will reduce the performance.
       *
       * If CODI_EnableVariableAdjointInterfaceInPrimalTapes is not set, the adjoint vector
       * is directly accessed but then only the GradientValue definition from the tape
       * is allowed in the reverse interpretation.
       *
       * @tparam AdjointData  The actual type of the adjoint data.
       */
  #if CODI_EnableVariableAdjointInterfaceInPrimalTapes
      template<typename AdjointData>
      using AdjVecType = AdjointInterfacePrimalImpl<Real, Index, AdjointData>;
  #else
      template<typename AdjointData>
      using AdjVecType = AdjointData;
  #endif

      /**
       * @brief The general interface for accessing the adjoint vector without knowing the actual data type.
       *
       * @tparam AdjointData  The actual type of the adjoint data.
       */
      template<typename AdjointData>
      using AdjVecInterface = AdjointInterfacePrimalImpl<Real, Index, AdjointData>;

      /**
       * @brief Helper function for wrapping the adjoint data vector if necessary.
       *
       * @return Either interface or adjointData
       *
       * @param[in,out]   interface  The interface for the wrapped adjointData.
       * @param[in,out] adjointData  The adjoint data vector.
       *
       * @tparam AdjointData  The actual type of the adjoint data.
       */
      template<typename AdjointData>
      AdjVecType<AdjointData>* wrapAdjointVector(AdjVecInterface<AdjointData>& interface, AdjointData* adjointData) {
         CODI_UNUSED(interface);   // To avoid unused warnings
         CODI_UNUSED(adjointData); // To avoid unused warnings

  #if CODI_EnableVariableAdjointInterfaceInPrimalTapes
        return &interface;
  #else
        static_assert(std::is_same<AdjointData, GradientValue>::value,
          "Please enable 'CODI_EnableVariableAdjointInterfacePrimalInPrimalTapes' in order"
          " to use custom adjoint vectors in the primal value tapes.");

        return adjointData;
  #endif
      }

    // ----------------------------------------------------------------------
    // Protected function for the communication with the including class
    // ----------------------------------------------------------------------

      /**
       * @brief Swap the data of the primal value module with the data of the other primal tape module.
       *
       * @param[in] other  The object with the other primal value module.
       */
      void swapPrimalValueModule(Tape& other) {
        std::swap(primals, other.primals);
        std::swap(primalsSize, other.primalsSize);
        std::swap(primalsIncr, other.primalsIncr);
      }

      /**
       * @brief Helper function: Sets the primal vector to a new size.
       *
       * @param[in] size The new size for the primal vector.
       */
      void resizePrimals(const Index& size) {
        Index oldSize = primalsSize;
        primalsSize = size;

        primals = (Real*)realloc((void*)primals, sizeof(Real) * (size_t)primalsSize);

        for(Index i = oldSize; i < primalsSize; ++i) {
          primals[i] = Real();
        }
      }

      /**
       * @brief Helper function: Cheks if the primal vector is big enough. Increases the vector size otherwise.
       */
      CODI_INLINE void checkPrimalsSize() {
        if(primalsSize <= cast().indexHandler.getMaximumGlobalIndex()) {
          Index newSize = 1 + (cast().indexHandler.getMaximumGlobalIndex() + 1) / primalsIncr;
          newSize = newSize * primalsIncr;
          resizePrimals(newSize);
        }
      }

      /**
       * @brief Action that pushes all passive values to the constant value vector.
       *
       * The function is used on the expression templates as an action.
       *
       * @param[in]  data  Not used in the action.
       * @param[in] value  The passive value that is pushed on the vector.
       */
      CODI_INLINE void pushPassive(int data, const PassiveReal& value) {
        CODI_UNUSED(data);
        constantValueVector.setDataAndMove(value);
      }

      /**
       * @brief Action that counts all active values, that is they have a non zero index.
       *
       * @param[in,out] count  The actual count of the active variables.
       * @param[in]     value  Not used in the action.
       * @param[in]     index  The index of the active real.
       */
      CODI_INLINE void countActiveValues(int* count, const Real& value, const Index& index) {
        CODI_UNUSED(value);
        if(0 != index) {
          *count += 1;
        }
      }

      /**
       * @brief Action that adds all indices to the index vector.
       *
       * For passive indices one of the temporary indices is used. e.g 1 to 255.
       * The correspoinding primal value is then pushed to the constant value vector.
       * The constant value is then used in the reverse sweep.
       *
       * @param[in,out] passiveVariableCount  The current count of the  inactive variables.
       * @param[in]                    value  The primal value of the active real.
       * @param[in]                    index  The index of the active real.
       */
      CODI_INLINE void pushIndices(int* passiveVariableCount, const Real& value, const Index& index) {
        Index pushIndex = index;
        if(0 == pushIndex) {
          *passiveVariableCount += 1;
          pushIndex = *passiveVariableCount;
          passiveValueVector.setDataAndMove(value);
        }

        indexVector.setDataAndMove(pushIndex);
      }

    public:

      /**
       * @brief Evaluate one handle in the primal sweep.
       *
       * The function sets the passive values for the evaluation of the expression. Afterwards the primal evaluation is
       * called.
       *
       * The positions are increased such that the next primal handle can be evaluated.
       *
       * @param[in]          funcObj  The function object that performs the reverse AD evaluation of an expression.
       * @param[in]          varSize  The number of variables of the expression.
       * @param[in]        constSize  The number constant variables of the expression.
       * @param[in]   passiveActives  The number of inactive values in the statement.
       * @param[in,out]     indexPos  The position in the index array.
       * @param[in]          indices  The index array.
       * @param[in,out]   passivePos  The position in the passive value array.
       * @param[in]         passives  The passive value array.
       * @param[in,out]  constantPos  The position in the constant value array.
       * @param[in]        constants  The constant value array.
       * @param[in,out] primalVector  The global vector with the primal variables.
       *
       * @tparam FuncObj  The function object that performs the actual evaluation of the expression.
       */
      template<typename FuncObj>
      static CODI_INLINE Real evaluatePrimalHandle(FuncObj funcObj,
                                                   size_t varSize,
                                                   size_t constSize,
                                                   const StatementInt& passiveActives,
                                                   size_t& indexPos,
                                                   Index* &indices,
                                                   size_t& passivePos,
                                                   Real* &passives,
                                                   size_t& constantPos,
                                                   PassiveReal* &constants,
                                                   Real* primalVector) {

        // first restore the primal values of the passive indices
        for(StatementInt i = 0; i < passiveActives; ++i) {
          primalVector[i + 1] = passives[passivePos + i];
        }

        Real result = funcObj(codi::addressof(indices[indexPos]), codi::addressof(constants[constantPos]), primalVector);

        indexPos += varSize;
        constantPos += constSize;
        passivePos += passiveActives;

        return result;
      }

      /**
       * @brief Curry the evaluation of a handle such that the sizes do not need to be stored.
       *
       * This is a helper function that can be stored as a function pointer in order to be able to evaluate a primal
       * handle without knowing the constant sizes for the expression. The function calls evluatePrimalHandle
       *
       * @param[in]   passiveActives  The number of inactive values in the statement.
       * @param[in,out]     indexPos  The position in the index array.
       * @param[in]          indices  The index array.
       * @param[in,out]   passivePos  The position in the passive value array.
       * @param[in]         passives  The passive value array.
       * @param[in,out]  constantPos  The position in the constant value array.
       * @param[in]        constants  The constant value array.
       * @param[in,out] primalVector  The global vector with the primal variables.
       *
       * @tparam Expr  The expression for which the function is created.
       */
      template<typename Expr>
      static CODI_INLINE Real curryEvaluatePrimalHandle(const StatementInt& passiveActives,
                                                        size_t& indexPos,
                                                        Index* &indices,
                                                        size_t& passivePos,
                                                        Real* &passives,
                                                        size_t& constantPos,
                                                        PassiveReal* &constants,
                                                        Real* primalVector) {
        return evaluatePrimalHandle(Expr::template getValue<Index, 0, 0>, ExpressionTraits<Expr>::maxActiveVariables, ExpressionTraits<Expr>::maxConstantVariables,
                                    passiveActives, indexPos, indices, passivePos, passives, constantPos, constants, primalVector);
      }

      /**
       * @brief Evaluate one handle in the reverse sweep.
       *
       * The function sets the primal values in the primal value vector for the inactive values.
       * Then it updates the genaral positions and calls the adjoint function of the handle.
       *
       * @param[in]          funcObj  The function object that performs the reverse AD evaluation of an expression.
       * @param[in]          varSize  The number of variables of the expression.
       * @param[in]        constSize  The number constant variables of the expression.
       * @param[in]              adj  The seed from the lhs of the statement.
       * @param[in]   passiveActives  The number of inactive values in the statement.
       * @param[in,out]     indexPos  The position in the index array.
       * @param[in]          indices  The index array.
       * @param[in,out]   passivePos  The position in the passive value array.
       * @param[in]         passives  The passive value array.
       * @param[in,out]  constantPos  The position in the constant value array.
       * @param[in]        constants  The constant value array.
       * @param[in,out] primalVector  The global vector with the primal variables.
       * @param[in,out]     adjoints  The adjoint vector for the reverse AD evaluation.
       *
       * @tparam FuncObj  The function object that performs the actual evaluation of the expression.
       */
      template<typename FuncObj>
      static CODI_INLINE void evaluateHandle(FuncObj funcObj,
                                             size_t varSize,
                                             size_t constSize,
                                             const PRIMAL_SEED_TYPE& adj,
                                             const StatementInt& passiveActives,
                                             size_t& indexPos,
                                             Index* &indices,
                                             size_t& passivePos,
                                             Real* &passives,
                                             size_t& constantPos,
                                             PassiveReal* &constants,
                                             Real* primalVector,
                                             PRIMAL_ADJOINT_TYPE* adjoints) {
        // first restore the primal values of the passive indices
        passivePos -= passiveActives;
        for(StatementInt i = 0; i < passiveActives; ++i) {
          primalVector[i + 1] = passives[passivePos + i];
        }

        // now update the regular pointers
        indexPos -= varSize;
        constantPos -= constSize;
        ENABLE_CHECK(OptZeroAdjoint, !isTotalZero(adj)){
          funcObj(adj, codi::addressof(indices[indexPos]), codi::addressof(constants[constantPos]), primalVector, adjoints);
        }
      }

      /**
       * @brief Curry the evaluation of a handle such that the sizes do not need to be stored.
       *
       * This is a helper function that can be stored as a function pointer in order to be able to evaluate a
       * handle without knowing the constant sizes for the expression. The function calls evluateHandle
       *
       * @param[in]              adj  The seed from the lhs of the statement.
       * @param[in]   passiveActives  The number of inactive values in the statement.
       * @param[in,out]     indexPos  The position in the index array.
       * @param[in]          indices  The index array.
       * @param[in,out]   passivePos  The position in the passive value array.
       * @param[in]         passives  The passive value array.
       * @param[in,out]  constantPos  The position in the constant value array.
       * @param[in]        constants  The constant value array.
       * @param[in,out] primalVector  The global vector with the primal variables.
       * @param[in,out]     adjoints  The adjoint vector for the reverse AD evaluation.
       *
       * @tparam Expr  The expression for which this helper is instantiated.
       */
      template<typename Expr>
      static CODI_INLINE void curryEvaluateHandle(const PRIMAL_SEED_TYPE& adj,
                                                  const StatementInt& passiveActives,
                                                  size_t& indexPos,
                                                  Index* &indices,
                                                  size_t& passivePos,
                                                  Real* &passives,
                                                  size_t& constantPos,
                                                  PassiveReal* &constants,
                                                  Real* primalVector,
                                                  PRIMAL_ADJOINT_TYPE* adjoints) {
        evaluateHandle(Expr::template evalAdjoint<Index, GradientValue, 0, 0>, ExpressionTraits<Expr>::maxActiveVariables, ExpressionTraits<Expr>::maxConstantVariables,
                       adj, passiveActives, indexPos, indices, passivePos, passives, constantPos, constants, primalVector, adjoints);
      }

      /**
       * @brief Evaluate one handle in the forward sweep.
       *
       * The function sets the primal values in the primal value vector for the inactive values.
       * Then it updates the genaral positions and calls the tangent function of the handle.
       *
       * @param[in]          funcObj  The function object that performs the tangent AD evaluation of an expression.
       * @param[in]          varSize  The number of variables of the expression.
       * @param[in]        constSize  The number constant variables of the expression.
       * @param[in]              adj  The seed from the lhs of the statement.
       * @param[in,out]   lhsAdjoint  The tangent value of the lhs of the expression.
       * @param[in]   passiveActives  The number of inactive values in the statement.
       * @param[in,out]     indexPos  The position in the index array.
       * @param[in]          indices  The index array.
       * @param[in,out]   passivePos  The position in the passive value array.
       * @param[in]         passives  The passive value array.
       * @param[in,out]  constantPos  The position in the constant value array.
       * @param[in]        constants  The constant value array.
       * @param[in,out] primalVector  The global vector with the primal variables.
       * @param[in,out]     adjoints  The adjoint vector for the reverse AD evaluation.
       *
       * @tparam FuncObj  The function object that performs the actual evaluation of the expression.
       */
      template<typename FuncObj>
      static CODI_INLINE Real evaluateForwardHandle(FuncObj funcObj,
                                                    size_t varSize,
                                                    size_t constSize,
                                                    const Real& adj,
                                                    GradientValue& lhsAdjoint,
                                                    const StatementInt& passiveActives,
                                                    size_t& indexPos,
                                                    Index* &indices,
                                                    size_t& passivePos,
                                                    Real* &passives,
                                                    size_t& constantPos,
                                                    PassiveReal* &constants,
                                                    Real* primalVector,
                                                    PRIMAL_ADJOINT_TYPE* adjoints) {

        // first restore the primal values of the passive indices
        for(StatementInt i = 0; i < passiveActives; ++i) {
          primalVector[i + 1] = passives[passivePos + i];
        }

        Real result = funcObj(adj, lhsAdjoint, codi::addressof(indices[indexPos]), codi::addressof(constants[constantPos]), primalVector, adjoints);

        indexPos += varSize;
        constantPos += constSize;
        passivePos += passiveActives;

        return result;
      }

      /**
       * @brief Curry the evaluation of a handle such that the sizes do not need to be stored.
       *
       * This is a helper function that can be stored as a function pointer in order to be able to evaluate a
       * handle without knowing the constant sizes for the expression. The function calls evluateForwardHandle
       *
       * @param[in]              adj  The seed from the lhs of the statement.
       * @param[in,out]   lhsAdjoint  The tangent value of the lhs of the expression.
       * @param[in]   passiveActives  The number of inactive values in the statement.
       * @param[in,out]     indexPos  The position in the index array.
       * @param[in]          indices  The index array.
       * @param[in,out]   passivePos  The position in the passive value array.
       * @param[in]         passives  The passive value array.
       * @param[in,out]  constantPos  The position in the constant value array.
       * @param[in]        constants  The constant value array.
       * @param[in,out] primalVector  The global vector with the primal variables.
       * @param[in,out]     adjoints  The adjoint vector for the reverse AD evaluation.
       *
       * @tparam Expr  The expression for which this helper is instantiated.
       */
      template<typename Expr>
      static CODI_INLINE Real curryEvaluateForwardHandle(const Real& adj,
                                                         GradientValue& lhsAdjoint,
                                                         const StatementInt& passiveActives,
                                                         size_t& indexPos,
                                                         Index* &indices,
                                                         size_t& passivePos,
                                                         Real* &passives,
                                                         size_t& constantPos,
                                                         PassiveReal* &constants,
                                                         Real* primalVector,
                                                         PRIMAL_ADJOINT_TYPE* adjoints) {
        return evaluateForwardHandle(Expr::template evalTangent<Index, GradientValue, 0, 0>, ExpressionTraits<Expr>::maxActiveVariables, ExpressionTraits<Expr>::maxConstantVariables,
                                    adj, lhsAdjoint, passiveActives, indexPos, indices, passivePos, passives, constantPos, constants, primalVector, adjoints);
      }

      /**
       * @brief Set the primal value in the primal value vector.
       *
       * @param[in]  index  The index of the primal value which is set.
       * @param[in] primal  The value which is set into the vector.
       */
      void setPrimalValue(const Index& index, const Real& primal) {
        primals[index] = primal;
      }

      /**
       * @brief Get the primal value from the primal value vector.
       *
       * @param[in]  index  The index of the primal value which is get.
       *
       * @return The value which is get from the vector.
       */
      Real getPrimalValue(const Index& index) const {
        return primals[index];
      }

      /**
       * @brief Get the primal value from the primal value vector.
       *
       * @param[in]  index  The index of the primal value which is get.
       *
       * @return The value which is get from the vector.
       */
      Real& primalValue(const Index& index) {
        return primals[index];
      }

    protected:

      /**
       * @brief Push a copy handle to the tape and all data that is rquired by the handle.
       *
       * The index of the lhs is updated.
       *
       * @param[in]     lhsValue  The value of the active real.
       * @param[in,out] lhsIndex  The index of the active real. The index is renewed in this method.
       * @param[in]     rhsIndex  The index of the rhs value.
       */
      CODI_INLINE void pushCopyHandle(const Real& lhsValue, Index& lhsIndex, const Index& rhsIndex) {
        indexVector.reserveItems(1);
        indexVector.setDataAndMove(rhsIndex);

        cast().pushStmtData(lhsIndex, lhsValue, HandleFactory::template createHandle<CopyExpr<Real>, Tape >(), StatementInt(0));
      }

    public:

    // ----------------------------------------------------------------------
    // Public function from the TapeInterface and ReverseTapeInterface
    // ----------------------------------------------------------------------

      /**
       * @brief Store the indices of the statement on the tape.
       *
       * The indices of the rhs expression are stored on the tape.
       * In addition the constant values of the expression are stored in
       * the constant vector and the passive arguments are also stored in the
       * constant vector. For the statement the handle of the rhs expression
       * is written to the tape.
       *
       * The gradient data of the lhs will get a new index.
       * The primal value of the lhs is set to the primal value of the rhs.
       *
       * @param[out]   lhsValue    The primal value of the lhs. This value is set to the value
       *                           of the right hand side.
       * @param[out]   lhsIndex    The gradient data of the lhs. The index will be updated.
       * @param[in]         rhs    The right hand side expression of the assignment.
       *
       * @tparam Rhs The expression on the rhs of the statement.
       */
      template<typename Rhs>
      CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const Rhs& rhs) {

        static_assert(ExpressionTraits<Rhs>::maxActiveVariables < MaxStatementIntSize, "Expression with to many arguments.");

        ENABLE_CHECK(OptTapeActivity, cast().active){

          int activeCount = 0;
          rhs.valueAction(&activeCount, &Tape::countActiveValues);

          if(0 != activeCount) {
            int passiveVariableNumber = ExpressionTraits<Rhs>::maxActiveVariables - activeCount;

            constantValueVector.reserveItems(ExpressionTraits<Rhs>::maxConstantVariables);
            size_t constantSize = constantValueVector.getChunkPosition();
            rhs.constantValueAction(*this, NULL, &Tape::pushPassive);
            codiAssert(ExpressionTraits<Rhs>::maxConstantVariables == constantValueVector.getChunkPosition() - constantSize);

            indexVector.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
            passiveValueVector.reserveItems(passiveVariableNumber);
            size_t indexSize = indexVector.getChunkPosition();
            int passieveVariableCount = 0;
            rhs.valueAction(&passieveVariableCount, &Tape::pushIndices);
            codiAssert(ExpressionTraits<Rhs>::maxActiveVariables == indexVector.getChunkPosition() - indexSize);
            codiAssert(passieveVariableCount == passiveVariableNumber);

            cast().pushStmtData(lhsIndex, rhs.getValue(), HandleFactory::template createHandle<Rhs, Tape >(), passiveVariableNumber);

            CODI_UNUSED(constantSize);  /* needed to avoid unused variable when the assersts are not enabled. */
            CODI_UNUSED(indexSize);  /* needed to avoid unused variable when the assersts are not enabled. */

  #if CODI_AdjointHandle_Primal
              Index* rhsIndices = NULL;
              PassiveReal* constants = NULL;

              auto posIndex = indexVector.getPosition();
              indexVector.getDataAtPosition(posIndex.chunk, indexSize, rhsIndices);

              auto posPassive = constantValueVector.getPosition();
              constantValueVector.getDataAtPosition(posPassive.chunk, constantSize, constants);

              resizeAdjointsToIndexSize();
              handleAdjointOperation(rhs.getValue(), lhsIndex, ExpressionHandleStore<Real*, Real, Index, Rhs>::getHandle(), passiveVariableNumber, constants, rhsIndices, primals, adjoints);
  #endif
          } else {
            cast().indexHandler.freeIndex(lhsIndex);
          }
        } else {
          cast().indexHandler.freeIndex(lhsIndex);
        }

        /* now set the value of the lhs */
        lhsValue = rhs.getValue();
      }

      /**
       * @brief Optimization for a passive value on the rhs. The lhs index is set to zero.
       *
       * No data is stored in this method.
       *
       * The primal value of the lhs is set to the primal value of the rhs.
       *
       * @param[out]   lhsValue  The primal value of the lhs. This value is set to the value
       *                         of the right hand side.
       * @param[out]   lhsIndex  The gradient data of the lhs. The index will be set to zero.
       * @param[in]         rhs  The right hand side expression of the assignment.
       */
      CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
        cast().indexHandler.freeIndex(lhsIndex);

        lhsValue = rhs;
      }

      /**
       * @brief Manual store routine.
       *
       * Use this routine to add a statement if the corresponding jacobi entries will be manually pushed onto the tape.
       *
       * The Jacobi entries must be pushed immediately after calling this routine using pushJacobiManual.
       *
       * See also the documentation in TapeInterfaceReverse::storeManual.
       *
       * @param[out]   lhsValue  The primal value of the lhs.
       * @param[out]   lhsIndex  The gradient data of the lhs.
       * @param[in]        size  The number of Jacobi entries.
       */
      CODI_INLINE void storeManual(const Real& lhsValue, Index& lhsIndex, StatementInt size) {
        ENABLE_CHECK (OptTapeActivity, cast().active){
          passiveValueVector.reserveItems(size);
          indexVector.reserveItems(size);

          cast().pushStmtData(lhsIndex, lhsValue, preaccHandles[size], size);
        }
      }

      /**
       * @brief Not used in this implementation.
       *
       * @param[in]  data  Not used in this implementation.
       * @param[in] value  Not used in this implementation.
       * @param[in] index  Not used in this implementation.
       */
      template<typename Data>
      CODI_INLINE void pushJacobi(Data& data, const Real& value, const Index& index) {
        CODI_UNUSED(data);
        CODI_UNUSED(value);
        CODI_UNUSED(index);

        codiAssert(false || "Should not be called.");
      }

      /**
       * @brief Not used in this implementation.
       *
       * @param[in]   data  Not used in this implementation.
       * @param[in] jacobi  Not used in this implementation.
       * @param[in]  value  Not used in this implementation.
       * @param[in]  index  Not used in this implementation.
       */
      template<typename Data>
      CODI_INLINE void pushJacobi(Data& data, const Real& jacobi, const Real& value, const Index& index) {
        CODI_UNUSED(data);
        CODI_UNUSED(jacobi);
        CODI_UNUSED(value);
        CODI_UNUSED(index);

        codiAssert(false || "Should not be called.");
      }

      /**
       * @brief Manual jacobi push routine.
       *
       * The indices are stored in the index vector. The jacobies
       * are stored in the constant vector.
       *
       * See also the documentation in TapeReverseInterface::pushJacobiManual.
       *
       * @param[in]   jacobi  Stored in the constant vector.
       * @param[in]    value  Not used in this implementation.
       * @param[in]    index  Stored in the index vector.
       */
      CODI_INLINE void pushJacobiManual(const Real& jacobi, const Real& value, const Index& index) {
        CODI_UNUSED(value);

        passiveValueVector.setDataAndMove(jacobi);
        indexVector.setDataAndMove(index);
      }

      /**
       * @brief Return the number of used stement entries
       * @return The number of used statement entries.
       */
      size_t getUsedStatementsSize() const {
        return stmtVector.getDataSize();
      }

      /**
       * @brief Return the number of used index entries
       * @return The number of used index entries.
       */
      size_t getUsedDataEntriesSize() const {
        return indexVector.getDataSize();
      }

      /**
       * @brief Return the number of used constant data entries
       * @return The number of used constant data entries.
       */
      size_t getUsedConstantDataSize() const {
        return constantValueVector.getDataSize();
      }

      /**
       * @brief Set the number of constant data etnres that should be available.
       *
       * @param[in] constantDataSize  The number of entries that should be available.
       */
      void setConstantDataSize(const size_t& constantDataSize) {
        constantValueVector.resize(constantDataSize);
      }
  };

  /**
   * @brief Creates a preaccumulation handle for a specific size
   *
   * @param[in] size  The size of the preaccumulation handle. This is the size of the index vector and the constant vector.
   */
  #define CREATE_PREACC_HANDLE(size) HandleFactory::template createHandle<PreaccExpr<typename TapeTypes::Real, size>, Tape>()

  /**
   * @brief Instantiation of the preaccumulation handles.
   *
   * @tparam TapeTypes  The definition of the tape types.
   */
  template <typename TapeTypes, typename Tape>
  const typename TapeTypes::Handle PrimalValueModule<TapeTypes, Tape>::preaccHandles[MaxStatementIntSize] = {
    CREATE_PREACC_HANDLE(0), CREATE_PREACC_HANDLE(1), CREATE_PREACC_HANDLE(2), CREATE_PREACC_HANDLE(3), CREATE_PREACC_HANDLE(4),
    CREATE_PREACC_HANDLE(5), CREATE_PREACC_HANDLE(6), CREATE_PREACC_HANDLE(7), CREATE_PREACC_HANDLE(8), CREATE_PREACC_HANDLE(9),
    CREATE_PREACC_HANDLE(10), CREATE_PREACC_HANDLE(11), CREATE_PREACC_HANDLE(12), CREATE_PREACC_HANDLE(13), CREATE_PREACC_HANDLE(14),
    CREATE_PREACC_HANDLE(15), CREATE_PREACC_HANDLE(16), CREATE_PREACC_HANDLE(17), CREATE_PREACC_HANDLE(18), CREATE_PREACC_HANDLE(19),
    CREATE_PREACC_HANDLE(20), CREATE_PREACC_HANDLE(21), CREATE_PREACC_HANDLE(22), CREATE_PREACC_HANDLE(23), CREATE_PREACC_HANDLE(24),
    CREATE_PREACC_HANDLE(25), CREATE_PREACC_HANDLE(26), CREATE_PREACC_HANDLE(27), CREATE_PREACC_HANDLE(28), CREATE_PREACC_HANDLE(29),
    CREATE_PREACC_HANDLE(30), CREATE_PREACC_HANDLE(31), CREATE_PREACC_HANDLE(32), CREATE_PREACC_HANDLE(33), CREATE_PREACC_HANDLE(34),
    CREATE_PREACC_HANDLE(35), CREATE_PREACC_HANDLE(36), CREATE_PREACC_HANDLE(37), CREATE_PREACC_HANDLE(38), CREATE_PREACC_HANDLE(39),
    CREATE_PREACC_HANDLE(40), CREATE_PREACC_HANDLE(41), CREATE_PREACC_HANDLE(42), CREATE_PREACC_HANDLE(43), CREATE_PREACC_HANDLE(44),
    CREATE_PREACC_HANDLE(45), CREATE_PREACC_HANDLE(46), CREATE_PREACC_HANDLE(47), CREATE_PREACC_HANDLE(48), CREATE_PREACC_HANDLE(49),
    CREATE_PREACC_HANDLE(50), CREATE_PREACC_HANDLE(51), CREATE_PREACC_HANDLE(52), CREATE_PREACC_HANDLE(53), CREATE_PREACC_HANDLE(54),
    CREATE_PREACC_HANDLE(55), CREATE_PREACC_HANDLE(56), CREATE_PREACC_HANDLE(57), CREATE_PREACC_HANDLE(58), CREATE_PREACC_HANDLE(59),
    CREATE_PREACC_HANDLE(60), CREATE_PREACC_HANDLE(61), CREATE_PREACC_HANDLE(62), CREATE_PREACC_HANDLE(63), CREATE_PREACC_HANDLE(64),
    CREATE_PREACC_HANDLE(65), CREATE_PREACC_HANDLE(66), CREATE_PREACC_HANDLE(67), CREATE_PREACC_HANDLE(68), CREATE_PREACC_HANDLE(69),
    CREATE_PREACC_HANDLE(70), CREATE_PREACC_HANDLE(71), CREATE_PREACC_HANDLE(72), CREATE_PREACC_HANDLE(73), CREATE_PREACC_HANDLE(74),
    CREATE_PREACC_HANDLE(75), CREATE_PREACC_HANDLE(76), CREATE_PREACC_HANDLE(77), CREATE_PREACC_HANDLE(78), CREATE_PREACC_HANDLE(79),
    CREATE_PREACC_HANDLE(80), CREATE_PREACC_HANDLE(81), CREATE_PREACC_HANDLE(82), CREATE_PREACC_HANDLE(83), CREATE_PREACC_HANDLE(84),
    CREATE_PREACC_HANDLE(85), CREATE_PREACC_HANDLE(86), CREATE_PREACC_HANDLE(87), CREATE_PREACC_HANDLE(88), CREATE_PREACC_HANDLE(89),
    CREATE_PREACC_HANDLE(90), CREATE_PREACC_HANDLE(91), CREATE_PREACC_HANDLE(92), CREATE_PREACC_HANDLE(93), CREATE_PREACC_HANDLE(94),
    CREATE_PREACC_HANDLE(95), CREATE_PREACC_HANDLE(96), CREATE_PREACC_HANDLE(97), CREATE_PREACC_HANDLE(98), CREATE_PREACC_HANDLE(99),
    CREATE_PREACC_HANDLE(100), CREATE_PREACC_HANDLE(101), CREATE_PREACC_HANDLE(102), CREATE_PREACC_HANDLE(103), CREATE_PREACC_HANDLE(104),
    CREATE_PREACC_HANDLE(105), CREATE_PREACC_HANDLE(106), CREATE_PREACC_HANDLE(107), CREATE_PREACC_HANDLE(108), CREATE_PREACC_HANDLE(109),
    CREATE_PREACC_HANDLE(110), CREATE_PREACC_HANDLE(111), CREATE_PREACC_HANDLE(112), CREATE_PREACC_HANDLE(113), CREATE_PREACC_HANDLE(114),
    CREATE_PREACC_HANDLE(115), CREATE_PREACC_HANDLE(116), CREATE_PREACC_HANDLE(117), CREATE_PREACC_HANDLE(118), CREATE_PREACC_HANDLE(119),
    CREATE_PREACC_HANDLE(120), CREATE_PREACC_HANDLE(121), CREATE_PREACC_HANDLE(122), CREATE_PREACC_HANDLE(123), CREATE_PREACC_HANDLE(124),
    CREATE_PREACC_HANDLE(125), CREATE_PREACC_HANDLE(126), CREATE_PREACC_HANDLE(127), CREATE_PREACC_HANDLE(128), CREATE_PREACC_HANDLE(129),
    CREATE_PREACC_HANDLE(130), CREATE_PREACC_HANDLE(131), CREATE_PREACC_HANDLE(132), CREATE_PREACC_HANDLE(133), CREATE_PREACC_HANDLE(134),
    CREATE_PREACC_HANDLE(135), CREATE_PREACC_HANDLE(136), CREATE_PREACC_HANDLE(137), CREATE_PREACC_HANDLE(138), CREATE_PREACC_HANDLE(139),
    CREATE_PREACC_HANDLE(140), CREATE_PREACC_HANDLE(141), CREATE_PREACC_HANDLE(142), CREATE_PREACC_HANDLE(143), CREATE_PREACC_HANDLE(144),
    CREATE_PREACC_HANDLE(145), CREATE_PREACC_HANDLE(146), CREATE_PREACC_HANDLE(147), CREATE_PREACC_HANDLE(148), CREATE_PREACC_HANDLE(149),
    CREATE_PREACC_HANDLE(150), CREATE_PREACC_HANDLE(151), CREATE_PREACC_HANDLE(152), CREATE_PREACC_HANDLE(153), CREATE_PREACC_HANDLE(154),
    CREATE_PREACC_HANDLE(155), CREATE_PREACC_HANDLE(156), CREATE_PREACC_HANDLE(157), CREATE_PREACC_HANDLE(158), CREATE_PREACC_HANDLE(159),
    CREATE_PREACC_HANDLE(160), CREATE_PREACC_HANDLE(161), CREATE_PREACC_HANDLE(162), CREATE_PREACC_HANDLE(163), CREATE_PREACC_HANDLE(164),
    CREATE_PREACC_HANDLE(165), CREATE_PREACC_HANDLE(166), CREATE_PREACC_HANDLE(167), CREATE_PREACC_HANDLE(168), CREATE_PREACC_HANDLE(169),
    CREATE_PREACC_HANDLE(170), CREATE_PREACC_HANDLE(171), CREATE_PREACC_HANDLE(172), CREATE_PREACC_HANDLE(173), CREATE_PREACC_HANDLE(174),
    CREATE_PREACC_HANDLE(175), CREATE_PREACC_HANDLE(176), CREATE_PREACC_HANDLE(177), CREATE_PREACC_HANDLE(178), CREATE_PREACC_HANDLE(179),
    CREATE_PREACC_HANDLE(180), CREATE_PREACC_HANDLE(181), CREATE_PREACC_HANDLE(182), CREATE_PREACC_HANDLE(183), CREATE_PREACC_HANDLE(184),
    CREATE_PREACC_HANDLE(185), CREATE_PREACC_HANDLE(186), CREATE_PREACC_HANDLE(187), CREATE_PREACC_HANDLE(188), CREATE_PREACC_HANDLE(189),
    CREATE_PREACC_HANDLE(190), CREATE_PREACC_HANDLE(191), CREATE_PREACC_HANDLE(192), CREATE_PREACC_HANDLE(193), CREATE_PREACC_HANDLE(194),
    CREATE_PREACC_HANDLE(195), CREATE_PREACC_HANDLE(196), CREATE_PREACC_HANDLE(197), CREATE_PREACC_HANDLE(198), CREATE_PREACC_HANDLE(199),
    CREATE_PREACC_HANDLE(200), CREATE_PREACC_HANDLE(201), CREATE_PREACC_HANDLE(202), CREATE_PREACC_HANDLE(203), CREATE_PREACC_HANDLE(204),
    CREATE_PREACC_HANDLE(205), CREATE_PREACC_HANDLE(206), CREATE_PREACC_HANDLE(207), CREATE_PREACC_HANDLE(208), CREATE_PREACC_HANDLE(209),
    CREATE_PREACC_HANDLE(210), CREATE_PREACC_HANDLE(211), CREATE_PREACC_HANDLE(212), CREATE_PREACC_HANDLE(213), CREATE_PREACC_HANDLE(214),
    CREATE_PREACC_HANDLE(215), CREATE_PREACC_HANDLE(216), CREATE_PREACC_HANDLE(217), CREATE_PREACC_HANDLE(218), CREATE_PREACC_HANDLE(219),
    CREATE_PREACC_HANDLE(220), CREATE_PREACC_HANDLE(221), CREATE_PREACC_HANDLE(222), CREATE_PREACC_HANDLE(223), CREATE_PREACC_HANDLE(224),
    CREATE_PREACC_HANDLE(225), CREATE_PREACC_HANDLE(226), CREATE_PREACC_HANDLE(227), CREATE_PREACC_HANDLE(228), CREATE_PREACC_HANDLE(229),
    CREATE_PREACC_HANDLE(230), CREATE_PREACC_HANDLE(231), CREATE_PREACC_HANDLE(232), CREATE_PREACC_HANDLE(233), CREATE_PREACC_HANDLE(234),
    CREATE_PREACC_HANDLE(235), CREATE_PREACC_HANDLE(236), CREATE_PREACC_HANDLE(237), CREATE_PREACC_HANDLE(238), CREATE_PREACC_HANDLE(239),
    CREATE_PREACC_HANDLE(240), CREATE_PREACC_HANDLE(241), CREATE_PREACC_HANDLE(242), CREATE_PREACC_HANDLE(243), CREATE_PREACC_HANDLE(244),
    CREATE_PREACC_HANDLE(245), CREATE_PREACC_HANDLE(246), CREATE_PREACC_HANDLE(247), CREATE_PREACC_HANDLE(248), CREATE_PREACC_HANDLE(249),
    CREATE_PREACC_HANDLE(250), CREATE_PREACC_HANDLE(251), CREATE_PREACC_HANDLE(252), CREATE_PREACC_HANDLE(253), CREATE_PREACC_HANDLE(254)
  };

  #undef CREATE_PREACC_HANDLE
}
