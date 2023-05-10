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

#include <vector>

#include "../../../config.h"
#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../misc/exceptions.hpp"
#include "../../../traits/tapeTraits.hpp"
#include "../../data/direction.hpp"
#include "linearSystemFlags.hpp"
#include "linearSystemInterface.hpp"
#include "linearSystemSpecializationDetection.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   *  Solves Ax=b and registers an external function on the tape which solves the specific AD mode equations.
   *
   *  Forward mode:  x_d  = A_v^-1 * (b_d - A_d * x_v) <br>
   *  Reverse mode: <br>
   *  &emsp;  s    = A_v^T^-1 * x_b <br>
   *  &emsp;  A_b += -x_v * s^T <br>
   *  &emsp;  b_b += s <br>
   *  &emsp;  x_b  = 0 <br>
   *
   *  The hints steer the algorithm, see LinearSystemInterface for details.
   *
   *  See \ref Example_21_Special_handling_of_linear_system_solvers for an example use of this class.
   *
   *  @tparam T_LinearSystem  Implementation of LinearSystemInterface.
   */
  template<typename T_LinearSystem, typename = void>
  struct LinearSystemSolverHandler {
    public:

      /// See LinearSystemSolverHandler.
      using LinearSystem = CODI_DD(T_LinearSystem, CODI_T(LinearSystemInterface<LinearSystemInterfaceTypes>));
      /// See LinearSystemInterfaceTypes.
      using Type = CODI_DD(typename LinearSystem::Type, CODI_DEFAULT_LHS_EXPRESSION);

      using Matrix = typename LinearSystem::Matrix;                      ///< See LinearSystemInterfaceTypes.
      using MatrixReal = typename LinearSystem::MatrixReal;              ///< See LinearSystemInterfaceTypes.
      using MatrixIdentifier = typename LinearSystem::MatrixIdentifier;  ///< See LinearSystemInterfaceTypes.
      using Vector = typename LinearSystem::Vector;                      ///< See LinearSystemInterfaceTypes.
      using VectorReal = typename LinearSystem::VectorReal;              ///< See LinearSystemInterfaceTypes.
      using VectorIdentifier = typename LinearSystem::VectorIdentifier;  ///< See LinearSystemInterfaceTypes.

    private:

      /*******************************************************************************/
      // Additional definitions

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.
      using Gradient = typename Type::Gradient;      ///< See LhsExpressionInterface.

      /// See LhsExpressionInterface.
      using Tape = CODI_DD(typename Type::Tape, CODI_DEFAULT_TAPE);

      /// Vector access definition of the tape.
      using VectorAccess = VectorAccessInterface<Real, Identifier>;

      /// Shortcut for the overload detection.
      using Overloads = LinearSystemSpecializationDetection<LinearSystem>;

      /*******************************************************************************/
      // Implementation of functors for iterators.

      /// Helper for the definition of functors based on the adjoint interface and the dimension.
      struct VectorAccessFunctor {
        public:

          size_t dim;                      ///< Current dimension to access.
          VectorAccess* adjointInterface;  ///< Adjoint interface handle.

          /// Constructor.
          VectorAccessFunctor(size_t dim, VectorAccess* adjointInterface)
              : dim(dim), adjointInterface(adjointInterface) {}
      };

      /// Extract from value into value_v and value_id.
      static void extract(Type const& value, Real& value_v, Identifier& value_id) {
        value_v = value.getValue();
        value_id = value.getIdentifier();
      }

      /// Extract the adjoint with value_id from the interface.
      struct ExtractAdjoint : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_b, Identifier const& value_id) {
            value_b = this->adjointInterface->getAdjoint(value_id, this->dim);
            this->adjointInterface->resetAdjoint(value_id, this->dim);
          }
      };

      /// Set value_v with value.
      static void getOutput(Type const& value, Real& value_v) {
        value_v = value.getValue();
      }

      /// Get the adjoint with value_id from the interface and return it in value_b.
      struct GetAdjoint : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_b, Identifier const& value_id) {
            value_b = this->adjointInterface->getAdjoint(value_id, this->dim);
          }
      };

      /// Get the primal with value_id from the interface and return it in value_v.
      struct GetPrimal : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_v, Identifier const& value_id) {
            value_v = this->adjointInterface->getPrimal(value_id);
          }
      };

      /// Get the primal and adjoint with value_id from the interface and return them in value_v and value_b.
      struct GetPrimalAndGetAdjoint : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_v, Real& value_b, Identifier const& value_id) {
            value_v = this->adjointInterface->getPrimal(value_id);
            value_b = this->adjointInterface->getAdjoint(value_id, this->dim);
          }
      };

      /// Get the primal and tangent with value_id from the interface and return them in value_v and value_b.
      using GetPrimalAndGetTangent = GetPrimalAndGetAdjoint;

      /// Get the tangent with value_id from the interface and return it in value_b.
      using GetTangent = GetAdjoint;

      /// Register value as output of the external function. Set value to value_v and update the value_id.
      static Real registerOutput(Type& value, Real& value_v, Identifier& value_id) {
        value = value_v;
        Real oldTemp = Type::getTape().registerExternalFunctionOutput(value);
        value_id = value.getIdentifier();

        return oldTemp;
      }

      /// Same as registerOutput but also update oldValue.
      static void registerOutputWithPrimal(Type& value, Real& value_v, Identifier& value_id, Real& oldValue) {
        oldValue = registerOutput(value, value_v, value_id);
      }

      /// Set value with value_v.
      static void setOutput(Type& value, Real const& value_v) {
        value = value_v;
      }

      /// Set the tangent at value_id to value_d.
      struct SetTangent : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_d, Identifier const& value_id) {
            this->adjointInterface->resetAdjoint(value_id, this->dim);
            this->adjointInterface->updateAdjoint(value_id, this->dim, value_d);
          }
      };

      /// Set the primal at value_id to value_v.
      struct SetPrimal : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_v, Identifier const& value_id) {
            this->adjointInterface->setPrimal(value_id, value_v);
          }
      };

      /// Set the primal and tangent at value_id to value_v and value_d.
      struct SetPrimalAndSetTangent : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_v, Real& value_d, Identifier const& value_id) {
            this->adjointInterface->setPrimal(value_id, value_v);
            this->adjointInterface->resetAdjoint(value_id, this->dim);
            this->adjointInterface->updateAdjoint(value_id, this->dim, value_d);
          }
      };

      /// Set the primal and tangent at value_id to value_v and value_d. Store the old primal value in oldValue.
      struct SetPrimalAndSetTangentAndUpdateOldPrimal : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_v, Real& value_d, Identifier const& value_id, Real& oldValue) {
            oldValue = this->adjointInterface->getPrimal(value_id);
            this->adjointInterface->setPrimal(value_id, value_v);
            this->adjointInterface->resetAdjoint(value_id, this->dim);
            this->adjointInterface->updateAdjoint(value_id, this->dim, value_d);
          }
      };

      /// Set the primal at value_id to value_v. Store the old primal value in oldValue.
      struct SetPrimalAndUpdateOldPrimals : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_v, Identifier const& value_id, Real& oldValue) {
            oldValue = this->adjointInterface->getPrimal(value_id);
            this->adjointInterface->setPrimal(value_id, value_v);
          }
      };

      /// Update the adjoint at value_id with value_b.
      struct UpdateAdjoint : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Real& value_b, Identifier const& value_id) {
            this->adjointInterface->updateAdjoint(value_id, this->dim, value_b);
          }
      };

      /// Update the adjoint from the dyadic product of x_v and b_b. See LinearSystemSolverHandler for a use case.
      struct UpdateAdjointDyadic : public VectorAccessFunctor {
        public:
          using VectorAccessFunctor::VectorAccessFunctor;

          void operator()(Identifier& mat_id, Real const& x_v, Real const& b_b) {
            Real adjoint = -x_v * b_b;
            this->adjointInterface->updateAdjoint(mat_id, this->dim, adjoint);
          }
      };

      /*******************************************************************************/
      // Detection of constant properties

      /// Only required for primal value tapes.
      static bool constexpr IsPrimalValueTape = TapeTraits::IsPrimalValueTape<Tape>::value;
      static bool constexpr StoreOldPrimals = IsPrimalValueTape & !Tape::LinearIndexHandling;

      /*******************************************************************************/
      // External function handle implementations

      /// Data used in the external functions in CoDiPack.
      struct ExtFuncData {
          MatrixReal* A_v;
          MatrixReal* A_v_trans;
          MatrixIdentifier* A_id;

          VectorIdentifier* b_id;

          VectorReal* x_v;
          VectorIdentifier* x_id;

          VectorReal* oldPrimals;

          LinearSystem lsi;
          LinearSystemSolverHints hints;

          ExtFuncData(LinearSystem lsi, LinearSystemSolverHints hints)
              : A_v(NULL),
                A_v_trans(NULL),
                A_id(NULL),
                b_id(NULL),
                x_v(NULL),
                x_id(NULL),
                oldPrimals(NULL),
                lsi(lsi),
                hints(hints) {}

          ~ExtFuncData() {
            if (NULL != A_v) {
              lsi.deleteMatrixReal(A_v);
            }
            if (NULL != A_v_trans) {
              lsi.deleteMatrixReal(A_v_trans);
            }
            if (NULL != A_id) {
              lsi.deleteMatrixIdentifier(A_id);
            }
            if (NULL != b_id) {
              lsi.deleteVectorIdentifier(b_id);
            }
            if (NULL != x_v) {
              lsi.deleteVectorReal(x_v);
            }
            if (NULL != x_id) {
              lsi.deleteVectorIdentifier(x_id);
            }
            if (NULL != oldPrimals) {
              lsi.deleteVectorReal(oldPrimals);
            }
          }
      };

      /** Reverse mode algorithm
       *  Computes:
       *  s = A^T^-1 * x_b
       *  A_b += -x_v * s^T
       *  b_b += s
       *  x_b = 0
       */
      static void solve_b(Tape* tape, void* d, VectorAccess* adjointInterface) {
        CODI_UNUSED(tape);

        if (!Overloads::SupportsReverseMode()) {
          CODI_EXCEPTION("Missing functionality for linear system reverse mode. iterateDyadic(%d), transposeMatrix(%d)",
                         Overloads::IsDyadicImplemented(), Overloads::IsTransposeImplemented());
        }

        ExtFuncData* data = (ExtFuncData*)d;

        if (!data->hints.test(LinearSystemSolverFlags::ReverseEvaluation)) {
          CODI_EXCEPTION(
              "Linear system reverse mode called without hint 'LinearSystemSolverFlags::ReverseEvaluation'.");
        }

        VectorReal* x_b = data->lsi.createVectorReal(data->x_id);
        VectorReal* s = data->lsi.createVectorReal(data->b_id);

        if (NULL != data->oldPrimals) {
          data->lsi.iterateVector(SetPrimal(0, adjointInterface), data->oldPrimals, data->x_id);
        }

        size_t maxDim = adjointInterface->getVectorSize();
        for (size_t curDim = 0; curDim < maxDim; curDim += 1) {
          data->lsi.iterateVector(ExtractAdjoint(curDim, adjointInterface), x_b, data->x_id);

          data->lsi.solveSystem(data->A_v_trans, x_b, s);

          data->lsi.iterateDyadic(UpdateAdjointDyadic(curDim, adjointInterface), data->A_id, data->x_v, s);
          data->lsi.iterateVector(UpdateAdjoint(curDim, adjointInterface), s, data->b_id);
        }

        data->lsi.deleteVectorReal(x_b);
        data->lsi.deleteVectorReal(s);
      }

      /** Forward mode algorithm
       *  Computes:
       *  x_v = A^-1 * b_v (only for primal value tapes and with the hint RecomputePrimalInForwardEvaluation)
       *  x_d = A^-1 * (b_d - A_d * x_v)
       */
      static void solve_d(Tape* tape, void* d, VectorAccess* adjointInterface) {
        CODI_UNUSED(tape);

        ExtFuncData* data = (ExtFuncData*)d;

        if (!Overloads::SupportsForwardMode()) {
          CODI_EXCEPTION("Missing functionality for linear system forward mode. subtractMultiply(%d)",
                         Overloads::IsSubtractMultiplyImplemented());
        }
        if (!data->hints.test(LinearSystemSolverFlags::ForwardEvaluation)) {
          CODI_EXCEPTION(
              "Linear system forward mode called without hint 'LinearSystemSolverFlags::ForwardEvaluation'.");
        }

        bool const updatePrimals =
            IsPrimalValueTape && data->hints.test(LinearSystemSolverFlags::RecomputePrimalInForwardEvaluation);

        MatrixReal* A_d = data->lsi.createMatrixReal(data->A_id);
        VectorReal* b_v = data->lsi.createVectorReal(data->b_id);  // b_v is also used as a temporary.
        VectorReal* b_d = data->lsi.createVectorReal(data->b_id);
        VectorReal* x_d = data->lsi.createVectorReal(data->x_id);

        size_t maxDim = adjointInterface->getVectorSize();
        for (size_t curDim = 0; curDim < maxDim; curDim += 1) {
          if (0 == curDim && updatePrimals) {
            data->lsi.iterateMatrix(GetPrimalAndGetTangent(curDim, adjointInterface), data->A_v, A_d, data->A_id);
            data->lsi.iterateVector(GetPrimalAndGetTangent(curDim, adjointInterface), b_v, b_d, data->b_id);
          } else {
            data->lsi.iterateMatrix(GetTangent(curDim, adjointInterface), A_d, data->A_id);
            data->lsi.iterateVector(GetTangent(curDim, adjointInterface), b_d, data->b_id);
          }

          if (0 == curDim && updatePrimals) {  // Solve primal system only once and transposed setup only once.

            if (NULL != data->A_v_trans) {
              // Only renew A_v_trans if it already exists.
              data->lsi.deleteMatrixReal(data->A_v_trans);
              data->A_v_trans = data->lsi.transposeMatrix(data->A_v);
            }

            data->lsi.solveSystem(data->A_v, b_v, data->x_v);
          }

          data->lsi.subtractMultiply(b_v, b_d, A_d, data->x_v);  // Use of b_v as temporary.

          data->lsi.solveSystem(data->A_v, b_v /* temporary */, x_d);

          if (updatePrimals) {
            if (NULL != data->oldPrimals) {
              data->lsi.iterateVector(SetPrimalAndSetTangentAndUpdateOldPrimal(curDim, adjointInterface), data->x_v,
                                      x_d, data->x_id, data->oldPrimals);
            } else {
              data->lsi.iterateVector(SetPrimalAndSetTangent(curDim, adjointInterface), data->x_v, x_d, data->x_id);
            }
          } else {
            data->lsi.iterateVector(SetTangent(curDim, adjointInterface), x_d, data->x_id);
          }
        }

        data->lsi.deleteMatrixReal(A_d);
        data->lsi.deleteVectorReal(b_v);
        data->lsi.deleteVectorReal(b_d);
        data->lsi.deleteVectorReal(x_d);
      }

      /** Primal algorithm
       *  Computes:
       *  x_v = A^-1 * b_v (ignores the hint RecomputePrimalInForwardEvaluation)
       */
      static void solve_p(Tape* tape, void* d, VectorAccess* adjointInterface) {
        CODI_UNUSED(tape);

        ExtFuncData* data = (ExtFuncData*)d;

        if (!data->hints.test(LinearSystemSolverFlags::PrimalEvaluation)) {
          CODI_EXCEPTION("Linear system primal mode called without hint 'LinearSystemSolverFlags::PrimalEvaluation'.");
        }

        VectorReal* b_v = data->lsi.createVectorReal(data->b_id);

        data->lsi.iterateMatrix(GetPrimal(0, adjointInterface), data->A_v, data->A_id);
        data->lsi.iterateVector(GetPrimal(0, adjointInterface), b_v, data->b_id);

        data->lsi.solveSystem(data->A_v, b_v, data->x_v);

        if (NULL != data->A_v_trans) {
          // Only renew trans if it already exists.

          data->lsi.deleteMatrixReal(data->A_v_trans);
          data->A_v_trans = data->lsi.transposeMatrix(data->A_v);
        }

        if (NULL != data->oldPrimals) {
          data->lsi.iterateVector(SetPrimalAndUpdateOldPrimals(0, adjointInterface), data->x_v, data->x_id,
                                  data->oldPrimals);
        } else {
          data->lsi.iterateVector(SetPrimal(0, adjointInterface), data->x_v, data->x_id);
        }

        data->lsi.deleteVectorReal(b_v);
      }

      static void deleteData(Tape* tape, void* d) {
        CODI_UNUSED(tape);

        ExtFuncData* data = (ExtFuncData*)d;
        delete data;  /// All vectors and matrices are deleted in the destructor.
      }

    public:

      /**
       * Solves the linear system Ax=b and adds an external function to the tape.
       *
       * The hints steer the memory management and computations of the algorithm.
       */
      void solve(LinearSystem lsi, Matrix* A, Vector* b, Vector* x, LinearSystemSolverHints hints) {
        Tape& tape = Type::getTape();

        MatrixReal* A_v = lsi.createMatrixReal(A);
        MatrixIdentifier* A_id = lsi.createMatrixIdentifier(A);
        VectorReal* b_v = lsi.createVectorReal(b);
        VectorIdentifier* b_id = lsi.createVectorIdentifier(b);
        VectorReal* x_v = lsi.createVectorReal(x);
        VectorIdentifier* x_id = lsi.createVectorIdentifier(x);

        lsi.iterateMatrix(extract, A, A_v, A_id);
        lsi.iterateVector(extract, b, b_v, b_id);

        if (hints.test(LinearSystemSolverFlags::ProvidePrimalSolution)) {
          lsi.iterateVector(getOutput, x, x_v);
        }

        if (Overloads::IsSolvePrimalImplemented()) {
          lsi.solveSystemPrimal(A_v, b_v, x_v);
        } else {
          lsi.solveSystem(A_v, b_v, x_v);
        }

        if (tape.isActive()) {
          MatrixReal* A_v_trans = NULL;
          if (hints.test(LinearSystemSolverFlags::ReverseEvaluation)) {
            A_v_trans = lsi.transposeMatrix(A_v);
          }

          VectorReal* oldPrimals = NULL;
          if (StoreOldPrimals && hints.test(LinearSystemSolverFlags::ReverseEvaluation)) {
            oldPrimals = b_v;  // Reuse b_v here for the primal value handling
            lsi.iterateVector(registerOutputWithPrimal, x, x_v, x_id, oldPrimals);
            b_v = NULL;  // Do not delete b_v
          } else {
            lsi.iterateVector(registerOutput, x, x_v, x_id);
          }

          ExtFuncData* data = new ExtFuncData(lsi, hints);
          if (hints.test(LinearSystemSolverFlags::ForwardEvaluation) ||
              hints.test(LinearSystemSolverFlags::PrimalEvaluation)) {
            data->A_v = A_v;
            A_v = NULL;  // Do not delete A_v
          }
          data->A_v_trans = A_v_trans;
          data->A_id = A_id;
          data->b_id = b_id;
          data->x_v = x_v;
          data->x_id = x_id;
          data->oldPrimals = oldPrimals;

          tape.pushExternalFunction(ExternalFunction<Tape>::create(solve_b, data, deleteData, solve_d, solve_p));

          if (b_v != NULL) {
            lsi.deleteVectorReal(b_v);
          }
          if (A_v != NULL) {
            lsi.deleteMatrixReal(A_v);
          }
        } else {
          lsi.iterateVector(setOutput, x, x_v);

          lsi.deleteMatrixReal(A_v);
          lsi.deleteMatrixIdentifier(A_id);
          lsi.deleteVectorReal(b_v);
          lsi.deleteVectorIdentifier(b_id);
          lsi.deleteVectorReal(x_v);
          lsi.deleteVectorIdentifier(x_id);
        }
      }
  };

#ifndef DOXYGEN_DISABLE
  /// Specialization of LinearSystemSolverHandler for non-AD types.
  /// @tparam T_LinearSystem  Implementation of LinearSystemInterface.
  template<typename T_LinearSystem>
  struct LinearSystemSolverHandler<T_LinearSystem, RealTraits::EnableIfPassiveReal<typename T_LinearSystem::Type>> {
    public:

      /// See LinearSystemSolverHandler.
      using LinearSystem = CODI_DD(T_LinearSystem, CODI_T(LinearSystemInterface<LinearSystemInterfaceTypes>));
      using Matrix = typename LinearSystem::Matrix;  ///< See LinearSystemInterfaceTypes.
      using Vector = typename LinearSystem::Vector;  ///< See LinearSystemInterfaceTypes.

    private:

      using Overloads = LinearSystemSpecializationDetection<LinearSystem>;

    public:

      /** Primal algorithm
       *  Computes:
       *  x_v = A^-1 * b_v
       */
      void solve(LinearSystem lsi, Matrix* A, Vector* b, Vector* x, LinearSystemSolverHints hints) {
        CODI_UNUSED(hints);

        if (Overloads::IsSolvePrimalImplemented()) {
          lsi.solveSystemPrimal(A, b, x);
        } else {
          lsi.solveSystem(A, b, x);
        }
      }
  };

  /// Specialization of LinearSystemSolverHandler for forward mode tapes.
  /// @tparam T_LinearSystem  Implementation of LinearSystemInterface.
  template<typename T_LinearSystem>
  struct LinearSystemSolverHandler<T_LinearSystem,
                                   TapeTraits::EnableIfForwardTape<typename T_LinearSystem::Type::Tape>> {
    public:

      using LinearSystem =
          CODI_DD(T_LinearSystem,
                  CODI_T(LinearSystemInterface<LinearSystemInterfaceTypes>));  ///< See LinearSystemSolverHandler.

      using Type = CODI_DD(typename LinearSystem::Type,
                           CODI_DEFAULT_LHS_EXPRESSION);  ///< See LinearSystemInterfaceTypes.

      using Matrix = typename LinearSystem::Matrix;                      ///< See LinearSystemInterfaceTypes.
      using MatrixReal = typename LinearSystem::MatrixReal;              ///< See LinearSystemInterfaceTypes.
      using MatrixIdentifier = typename LinearSystem::MatrixIdentifier;  ///< See LinearSystemInterfaceTypes.
      using Vector = typename LinearSystem::Vector;                      ///< See LinearSystemInterfaceTypes.
      using VectorReal = typename LinearSystem::VectorReal;              ///< See LinearSystemInterfaceTypes.
      using VectorIdentifier = typename LinearSystem::VectorIdentifier;  ///< See LinearSystemInterfaceTypes.

    private:

      /*******************************************************************************/
      // Additional definitions

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.
      using Gradient = typename Type::Gradient;      ///< See LhsExpressionInterface.

      /// Shortcut for the overload detection.
      using Overloads = LinearSystemSpecializationDetection<LinearSystem>;

      /*******************************************************************************/
      // Implementation of functors for iterators.

      /// Helper for the definition of functors based on the dimension.
      struct DimFunctor {
          size_t dim;

          DimFunctor(size_t dim) : dim(dim) {}
      };

      /// Set value_v with value.
      static void getOutput(Type& value, Real& value_v) {
        value_v = value.getValue();
      }

      /// Get the primals into value_v and tangents into value_d from value.
      struct GetPrimalAndGetTangent : public DimFunctor {
        public:
          using DimFunctor::DimFunctor;

          void operator()(Type const& value, Real& value_v, Real& value_d) {
            value_v = value.getValue();
            value_d = GradientTraits::at(value.getGradient(), this->dim);
          }
      };

      /// Get tangents into value_d from value.
      struct GetTangent : public DimFunctor {
        public:
          using DimFunctor::DimFunctor;

          void operator()(Type const& value, Real& value_d) {
            value_d = GradientTraits::at(value.getGradient(), this->dim);
          }
      };

      /// Set the primals from value_v and tangents from value_d into value.
      struct SetPrimalAndSetTangent : public DimFunctor {
        public:
          using DimFunctor::DimFunctor;

          void operator()(Type& value, Real const& value_v, Real const& value_d) {
            value.value() = value_v;
            GradientTraits::at(value.gradient(), this->dim) = value_d;
          }
      };

      /// Set the tangents from value_d into value.
      struct SetTangent : public DimFunctor {
        public:
          using DimFunctor::DimFunctor;

          void operator()(Type& value, Real const& value_d) {
            GradientTraits::at(value.gradient(), this->dim) = value_d;
          }
      };

    public:

      /** Forward mode algorithm
       *  Computes:
       *  x_v = A^-1 * b_v
       *  x_d = A^-1 * (b_d - A_d * x_v)
       */
      void solve(LinearSystem lsi, Matrix* A, Vector* b, Vector* x, LinearSystemSolverHints hints) {
        CODI_UNUSED(hints);

        MatrixReal* A_v = lsi.createMatrixReal(A);
        MatrixReal* A_d = lsi.createMatrixReal(A);
        VectorReal* b_v = lsi.createVectorReal(b);
        VectorReal* b_d = lsi.createVectorReal(b);
        VectorReal* x_v = lsi.createVectorReal(x);
        VectorReal* x_d = lsi.createVectorReal(x);

        size_t maxDim = GradientTraits::dim<Gradient>();

        if (hints.test(LinearSystemSolverFlags::ProvidePrimalSolution)) {
          lsi.iterateVector(getOutput, x, x_v);
        }

        for (size_t curDim = 0; curDim < maxDim; curDim += 1) {
          if (0 == curDim) {
            lsi.iterateMatrix(GetPrimalAndGetTangent(curDim), A, A_v, A_d);
            lsi.iterateVector(GetPrimalAndGetTangent(curDim), b, b_v, b_d);
          } else {
            lsi.iterateMatrix(GetTangent(curDim), A, A_d);
            lsi.iterateVector(GetTangent(curDim), b, b_d);
          }

          if (0 == curDim) {  // Solve primal system only once.
            // Solve Ax = b
            if (Overloads::IsSolvePrimalImplemented()) {
              lsi.solveSystemPrimal(A_v, b_v, x_v);
            } else {
              lsi.solveSystem(A_v, b_v, x_v);
            }
          }

          // temp(x_d) = b_d - A_d * x
          lsi.subtractMultiply(x_d, b_d, A_d, x_v);

          std::swap(b_d, x_d);  // Move temporary to b_d.

          // Solve A x_d = temp
          lsi.solveSystem(A_v, b_d, x_d);

          if (0 == curDim) {
            lsi.iterateVector(SetPrimalAndSetTangent(curDim), x, x_v, x_d);
          } else {
            lsi.iterateVector(SetTangent(curDim), x, x_d);
          }
        }

        lsi.deleteMatrixReal(A_v);
        lsi.deleteMatrixReal(A_d);
        lsi.deleteVectorReal(b_v);
        lsi.deleteVectorReal(b_d);
        lsi.deleteVectorReal(x_v);
        lsi.deleteVectorReal(x_d);
      }
  };
#endif

  /**
   * Solve Ax=b and add an external function to the tape such that the appropriate AD modes are also solved.
   *
   * @param lsi  The implementation of LinearSystemInterface which defines the linear system solution algorithm.
   * @param A  The matrix.
   * @param b  The right hand side.
   * @param x The solution.
   * @param hints  Hints for the AD algorithm. See function description.
   */
  template<typename LSInterface>
  void solveLinearSystem(LSInterface lsi, typename LSInterface::Matrix& A, typename LSInterface::Vector& b,
                         typename LSInterface::Vector& x,
                         LinearSystemSolverHints hints = LinearSystemSolverHints::ALL()) {
    LinearSystemSolverHandler<LSInterface> handler;
    handler.solve(lsi, &A, &b, &x, hints);
  }
}
