/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2018 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

/*
 * In order to include this file the user has to define the preprocessor macro CHILD_VECTOR_TYPE, CONSTANT_VECTOR_TYPE,
 * INDEX_VECTOR_TYPE, STMT_VECTOR_TYPE and TAPE_NAME.
 *
 * CHILD_VECTOR_TYPE defines the type of the nested vector for the data vector.
 * CONSTANT_VECTOR_TYPE defines the type of the constant value vector.
 * INDEX_VECTOR_TYPE defines the type of the index vector.
 * STMT_TYPE defines the type of the statement vector.
 * TAPE_NAME defines the name of the including tape.
 *
 * All these macros, except TAPE_NAME, are undefined at the end of the file.
 *
 * The module defines the structures stmtVector, indexVector, constantValueVector, primals, primalsSize, primalsIncr.
 * The module defines the static structures InputHandle, CopyHandle and PreaccHandles,
 * The module defines the types PrimalChildVector, PrimalChildPosition, StatmentVector, StatementChunk, StmtPosition, IndexVector, IndexChunk,
 * IndexPosition, ConstantValueVector, ConstantValueChunk, ConstantValuePosition.
 *
 * It defines the methods store(Expr), store(const), store(User), pushJacobi, printPrimalValueStatistics from the TapeInterface and ReverseTapeInterface.
 *
 * It defines the methods resizePrimals, checkPrimalsSize, evaluateHandle, evaluateConstantValues, getUsedStatementsSize,
 * getUsedDataEntiresSize, getUsedConstantDataSize, setConstantDataSize, swapPrimalValueModule as interface functions for the including class.
 *
 * It defines the static methods inputHandleFunc, copyHandleFunc, preaccHandleFunc as interface functions for the tape.
 */

#ifndef CHILD_VECTOR_TYPE
  #error Please define the type of the child vector
#endif
#ifndef PASSIVE_VECTOR_TYPE
  #error Please define the name of the constant value vector type.
#endif
#ifndef CONSTANT_VECTOR_TYPE
  #error Please define the name of the constant value vector type.
#endif
#ifndef INDEX_VECTOR_TYPE
  #error Please define the name of the index vector type.
#endif
#ifndef STMT_VECTOR_TYPE
  #error Please define the name of the statment vector type.
#endif
#ifndef TAPE_NAME
  #error Please define the name of the tape.
#endif

		private:

  // ----------------------------------------------------------------------
  // All definitions of the module
  // ----------------------------------------------------------------------

    /** @brief The child vector for the primal data vector. */
    typedef CHILD_VECTOR_TYPE PrimalChildVector;

    /** @brief The position type of the primal child vector */
    typedef typename PrimalChildVector::Position PrimalChildPosition;

    /** @brief The vector for the statement data. */
    typedef STMT_VECTOR_TYPE StatementVector;
    /** @brief The data for the statments */
    typedef typename StatementVector::ChunkType StatementChunk;
    /** @brief The position type of the statment data. */
    typedef typename StatementVector::Position StmtPosition;
    /** @brief The data for the each statements. */
    StatementVector stmtVector;

    /** @brief The vector for the index data. */
    typedef INDEX_VECTOR_TYPE IndexVector;
    /** @brief The data for the indices*/
    typedef typename IndexVector::ChunkType IndexChunk;
    /** @brief The position type of the index data. */
    typedef typename IndexVector::Position IndexPosition;
    /** @brief The data for the indices of each statements. */
    IndexVector indexVector;

    /** @brief The vector for the passive data. */
    typedef PASSIVE_VECTOR_TYPE PassiveValueVector;
    /** @brief The data for the passives*/
    typedef typename PassiveValueVector::ChunkType PassiveValueChunk;
    /** @brief The position type of the passive data. */
    typedef typename PassiveValueVector::Position PassiveValuePosition;
    /** @brief The data for the passive values of each statements. */
    PassiveValueVector passiveValueVector;

    /** @brief The vector for the constant data. */
    typedef CONSTANT_VECTOR_TYPE ConstantValueVector;
    /** @brief The data for the constants*/
    typedef typename ConstantValueVector::ChunkType ConstantValueChunk;
    /** @brief The position type of the constant data. */
    typedef typename ConstantValueVector::Position ConstantValuePosition;
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

  private:

#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
    template<typename AdjointData>
    using AdjVecType = AdjointInterfacePrimalImpl<Real, Index, AdjointData>;
#else
    template<typename AdjointData>
    using AdjVecType = AdjointData;
#endif

    template<typename AdjointData>
    using AdjVecInterface = AdjointInterfacePrimalImpl<Real, Index, AdjointData>;

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
  // Private function for the communication with the including class
  // ----------------------------------------------------------------------

    /**
     * @brief Swap the data of the primal value module with the data of the other primal tape module.
     *
     * @param[in] other  The object with the other primal value module.
     */
    void swapPrimalValueModule(TAPE_NAME<TapeTypes>& other) {
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
     *
     * @param[in] size The new size for the primal vector.
     */
    CODI_INLINE void checkPrimalsSize() {
      if(primalsSize <= indexHandler.getMaximumGlobalIndex()) {
        Index newSize = 1 + (indexHandler.getMaximumGlobalIndex() + 1) / primalsIncr;
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
     * @param[in,out] count  The current count of the  inactive variables.
     * @param[in]     value  The primal value of the active real.
     * @param[in]     index  The index of the active real.
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

    private:

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

      pushStmtData(lhsIndex, lhsValue, HandleFactory::template createHandle<CopyExpr<Real>, TAPE_NAME<TapeTypes> >(), StatementInt(0));
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

      ENABLE_CHECK(OptTapeActivity, active){

        int activeCount = 0;
        rhs.valueAction(&activeCount, &TAPE_NAME<TapeTypes>::countActiveValues);

        if(0 != activeCount) {
          int passiveVariableNumber = ExpressionTraits<Rhs>::maxActiveVariables - activeCount;

          constantValueVector.reserveItems(ExpressionTraits<Rhs>::maxConstantVariables);
          size_t constantSize = constantValueVector.getChunkPosition();
          rhs.constantValueAction(*this, NULL, &TAPE_NAME<TapeTypes>::pushPassive);
          codiAssert(ExpressionTraits<Rhs>::maxConstantVariables == constantValueVector.getChunkPosition() - constantSize);

          indexVector.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
          passiveValueVector.reserveItems(passiveVariableNumber);
          size_t indexSize = indexVector.getChunkPosition();
          int passieveVariableCount = 0;
          rhs.valueAction(&passieveVariableCount, &TAPE_NAME<TapeTypes>::pushIndices);
          codiAssert(ExpressionTraits<Rhs>::maxActiveVariables == indexVector.getChunkPosition() - indexSize);
          codiAssert(passieveVariableCount == passiveVariableNumber);

          pushStmtData(lhsIndex, rhs.getValue(), HandleFactory::template createHandle<Rhs, TAPE_NAME<TapeTypes> >(), passiveVariableNumber);

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
          indexHandler.freeIndex(lhsIndex);
        }
      } else {
        indexHandler.freeIndex(lhsIndex);
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
      indexHandler.freeIndex(lhsIndex);

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
      ENABLE_CHECK (OptTapeActivity, active){
        passiveValueVector.reserveItems(size);
        indexVector.reserveItems(size);

        pushStmtData(lhsIndex, lhsValue, preaccHandles[size], size);
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
     * @brief Adds information about primal value tape.
     *
     * Adds the number of chunks, the total number of entries, the
     * allocated memory and the used memory for the constant vector, the statement
     * vector and the index vector.
     *
     * @param[in,out] values  The information is added to the values
     */
    void addPrimalValueValues(TapeValues& values) const {

      size_t nChunksIndex  = indexVector.getNumChunks();
      size_t totalIndex    = indexVector.getDataSize();
      size_t sizeIndexEntry = IndexChunk::EntrySize;
      double memoryUsedIndex = (double)totalIndex*(double)sizeIndexEntry* BYTE_TO_MB;
      double memoryAllocIndex= (double)nChunksIndex*(double)indexVector.getChunkSize()*(double)sizeIndexEntry* BYTE_TO_MB;

      size_t nChunksStmt  = stmtVector.getNumChunks();
      size_t totalStmt    = stmtVector.getDataSize();
      size_t sizeStmtEntry = StatementChunk::EntrySize;
      double memoryUsedStmt = (double)totalStmt*(double)sizeStmtEntry* BYTE_TO_MB;
      double memoryAllocStmt= (double)nChunksStmt*(double)stmtVector.getChunkSize()*(double)sizeStmtEntry* BYTE_TO_MB;

      size_t nChunksPassive  = constantValueVector.getNumChunks();
      size_t totalPassive    = constantValueVector.getDataSize();
      size_t sizePassiveEntry = ConstantValueChunk::EntrySize;
      double memoryUsedPassive = (double)totalPassive*(double)sizePassiveEntry* BYTE_TO_MB;
      double memoryAllocPassive= (double)nChunksPassive*(double)constantValueVector.getChunkSize()*(double)sizePassiveEntry* BYTE_TO_MB;

      size_t totalPrimal   = primalsSize;
      size_t sizePrimalEntry = sizeof(Real);
      double memoryAllocPrimal = (double)totalPrimal*(double)sizePrimalEntry* BYTE_TO_MB;

      values.addSection("Primal vector");
      values.addData("Total number", totalPrimal);
      values.addData("Memory allocated", memoryAllocPrimal, true, true);

      values.addSection("Statements");
      values.addData("Total number", totalStmt);
      values.addData("Number of chunks", nChunksStmt);
      values.addData("Memory used", memoryUsedStmt, true, false);
      values.addData("Memory allocated", memoryAllocStmt, false, true);

      values.addSection("Index entries");
      values.addData("Total number", totalIndex);
      values.addData("Number of chunks", nChunksIndex);
      values.addData("Memory used", memoryUsedIndex, true, false);
      values.addData("Memory allocated", memoryAllocIndex, false, true);

      values.addSection("Passive data entries");
      values.addData("Total number", totalPassive);
      values.addData("Number of chunks", nChunksPassive);
      values.addData("Memory used", memoryUsedPassive, true, false);
      values.addData("Memory allocated", memoryAllocPassive, false, true);
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

#undef CHILD_VECTOR_TYPE
#undef VECTOR_TYPE
