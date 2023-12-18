#pragma once

#include <codi/config.h>

#include <codi/misc/macros.hpp>
#include <codi/tools/lowlevelFunctions/eigenWrappers.hpp>
#include <codi/tools/lowlevelFunctions/generationHelperCoDiPack.hpp>
#include <codi/tools/lowlevelFunctions/lowLevelFunctionCreationUtilities.hpp>

/** \copydoc codi::Namespace */
namespace codi {
  /// Low level function generation for matrixMatrixMultiplication.
  template<Eigen::StorageOptions eigenStore, typename Type>
  struct ExtFunc_matrixMatrixMultiplication {
      /// Abbreviation for vector access interface.
      using AdjointVectorAccess = codi::VectorAccessInterface<typename Type::Real, typename Type::Identifier>*;
      /// Abbreviation for tape.
      using Tape = typename Type::Tape;

      /// Id for this function.
      static codi::Config::LowLevelFunctionToken ID;

      /// Function for forward interpretation.
      CODI_INLINE static void forward(Tape* tape, codi::ByteDataView& dataStore, AdjointVectorAccess adjoints) {
        codi::TemporaryMemory& allocator = tape->getTemporaryMemory();
        codiAssert(allocator.isEmpty());  // No memory should be allocated. We would free it at the end.
        using LLFH = codi::LowLevelFunctionCreationUtilities<2>;

        // Traits for arguments.
        using Trait_A = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_B = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_R = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_n = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_k = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_m = typename LLFH::PassiveStoreTrait<int, uint8_t>;

        // Declare variables.
        typename LLFH::ActivityStoreType activityStore = {};
        typename Trait_A::ArgumentStore A_store = {};
        typename Trait_B::ArgumentStore B_store = {};
        typename Trait_R::ArgumentStore R_store = {};
        typename Trait_n::Store n = {};
        typename Trait_k::Store k = {};
        typename Trait_m::Store m = {};

        bool active_A = false;
        bool active_B = false;

        // Restore data.
        LLFH::restoreActivity(&dataStore, activityStore);
        active_A = LLFH::getActivity(activityStore, 0);
        active_B = LLFH::getActivity(activityStore, 1);
        Trait_n::restore(&dataStore, allocator, 1, true, n);
        Trait_k::restore(&dataStore, allocator, 1, true, k);
        Trait_m::restore(&dataStore, allocator, 1, true, m);
        Trait_A::restore(&dataStore, allocator, n * k, LLFH::createRestoreActions(true, false, active_A, active_B),
                         A_store);
        Trait_B::restore(&dataStore, allocator, k * m, LLFH::createRestoreActions(true, false, active_B, active_A),
                         B_store);
        Trait_R::restore(&dataStore, allocator, n * m, LLFH::createRestoreActions(false, true, false, true), R_store);

        if (Tape::HasPrimalValues) {
          // Get primal values for inputs.
          if (active_A && active_B) {
            Trait_A::getPrimalsFromVector(adjoints, n * k, A_store.identifierIn(), A_store.primal());
          }
          if (active_B && active_A) {
            Trait_B::getPrimalsFromVector(adjoints, k * m, B_store.identifierIn(), B_store.primal());
          }
        }

        for (size_t curDim = 0; curDim < adjoints->getVectorSize(); curDim += 1) {
          // Get input gradients.
          if (active_A) {
            Trait_A::getGradients(adjoints, n * k, false, A_store.identifierIn(), A_store.gradientIn(), curDim);
          }
          if (active_B) {
            Trait_B::getGradients(adjoints, k * m, false, B_store.identifierIn(), B_store.gradientIn(), curDim);
          }
          if (Tape::HasPrimalValues && 0 == curDim) {
            if (!Tape::LinearIndexHandling) {
              // Update old primal values.
              Trait_R::getPrimalsFromVector(adjoints, n * m, R_store.identifierOut(), R_store.oldPrimal());
            }

            // Set new primal values.
            Trait_R::setPrimalsIntoVector(adjoints, n * m, R_store.identifierOut(), R_store.primal());
          }

          // Evaluate forward mode.
          callForward(A_store.primal(), active_A, A_store.gradientIn(), B_store.primal(), active_B,
                      B_store.gradientIn(), R_store.primal(), R_store.gradientOut(), n, k, m);

          Trait_R::setGradients(adjoints, n * m, false, R_store.identifierOut(), R_store.gradientOut(), curDim);
        }

        allocator.free();
      }

      /// Forward function for derivative evaluation.
      CODI_INLINE static void callForward(typename codi::ActiveArgumentStoreTraits<Type*>::Real const* A, bool active_A,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Gradient* A_d_in,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Real const* B, bool active_B,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Gradient* B_d_in,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Real* R,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Gradient* R_d_out,
                                          typename codi::PassiveArgumentStoreTraits<int, uint8_t>::Store n,
                                          typename codi::PassiveArgumentStoreTraits<int, uint8_t>::Store k,
                                          typename codi::PassiveArgumentStoreTraits<int, uint8_t>::Store m) {
        codi::CODI_UNUSED(A, active_A, A_d_in, B, active_B, B_d_in, R, R_d_out, n, k, m);
        if (active_A) {
          mapEigen<eigenStore>(R_d_out, n, m) += mapEigen<eigenStore>(A_d_in, n, k) * mapEigen<eigenStore>(B, k, m);
        }
        if (active_B) {
          mapEigen<eigenStore>(R_d_out, n, m) += mapEigen<eigenStore>(A, n, k) * mapEigen<eigenStore>(B_d_in, k, m);
        }
        mapEigen<eigenStore>(R, n, m) = mapEigen<eigenStore>(A, n, k) * mapEigen<eigenStore>(B, k, m);
      }

      /// Function for reverse interpretation.
      CODI_INLINE static void reverse(Tape* tape, codi::ByteDataView& dataStore, AdjointVectorAccess adjoints) {
        codi::TemporaryMemory& allocator = tape->getTemporaryMemory();
        codiAssert(allocator.isEmpty());  // No memory should be allocated. We would free it at the end.
        using LLFH = codi::LowLevelFunctionCreationUtilities<2>;

        // Traits for arguments.
        using Trait_A = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_B = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_R = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_n = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_k = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_m = typename LLFH::PassiveStoreTrait<int, uint8_t>;

        // Declare variables.
        typename LLFH::ActivityStoreType activityStore = {};
        typename Trait_A::ArgumentStore A_store = {};
        typename Trait_B::ArgumentStore B_store = {};
        typename Trait_R::ArgumentStore R_store = {};
        typename Trait_n::Store n = {};
        typename Trait_k::Store k = {};
        typename Trait_m::Store m = {};

        bool active_A = false;
        bool active_B = false;

        // Restore data.
        LLFH::restoreActivity(&dataStore, activityStore);
        active_A = LLFH::getActivity(activityStore, 0);
        active_B = LLFH::getActivity(activityStore, 1);
        Trait_n::restore(&dataStore, allocator, 1, true, n);
        Trait_k::restore(&dataStore, allocator, 1, true, k);
        Trait_m::restore(&dataStore, allocator, 1, true, m);
        Trait_A::restore(&dataStore, allocator, n * k, LLFH::createRestoreActions(true, false, active_A, active_B),
                         A_store);
        Trait_B::restore(&dataStore, allocator, k * m, LLFH::createRestoreActions(true, false, active_B, active_A),
                         B_store);
        Trait_R::restore(&dataStore, allocator, n * m, LLFH::createRestoreActions(false, true, false, true), R_store);

        if (Tape::HasPrimalValues) {
          if (!Tape::LinearIndexHandling) {
            // Restore old primal values from outputs.
            Trait_R::setPrimalsIntoVector(adjoints, n * m, R_store.identifierOut(), R_store.oldPrimal());
          }

          // Get primal values for inputs.
          if (active_A && active_B) {
            Trait_A::getPrimalsFromVector(adjoints, n * k, A_store.identifierIn(), A_store.primal());
          }
          if (active_B && active_A) {
            Trait_B::getPrimalsFromVector(adjoints, k * m, B_store.identifierIn(), B_store.primal());
          }
        }

        for (size_t curDim = 0; curDim < adjoints->getVectorSize(); curDim += 1) {
          // Get output gradients.
          Trait_R::getGradients(adjoints, n * m, true, R_store.identifierOut(), R_store.gradientOut(), curDim);

          // Evaluate reverse mode.
          callReverse(A_store.primal(), active_A, A_store.gradientIn(), B_store.primal(), active_B,
                      B_store.gradientIn(), R_store.primal(), R_store.gradientOut(), n, k, m);

          if (active_A) {
            Trait_A::setGradients(adjoints, n * k, true, A_store.identifierIn(), A_store.gradientIn(), curDim);
          }
          if (active_B) {
            Trait_B::setGradients(adjoints, k * m, true, B_store.identifierIn(), B_store.gradientIn(), curDim);
          }
        }

        allocator.free();
      }

      /// Reverse function for derivative evaluation.
      CODI_INLINE static void callReverse(typename codi::ActiveArgumentStoreTraits<Type*>::Real const* A, bool active_A,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Gradient* A_b_in,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Real const* B, bool active_B,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Gradient* B_b_in,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Real* R,
                                          typename codi::ActiveArgumentStoreTraits<Type*>::Gradient* R_b_out,
                                          typename codi::PassiveArgumentStoreTraits<int, uint8_t>::Store n,
                                          typename codi::PassiveArgumentStoreTraits<int, uint8_t>::Store k,
                                          typename codi::PassiveArgumentStoreTraits<int, uint8_t>::Store m) {
        codi::CODI_UNUSED(A, active_A, A_b_in, B, active_B, B_b_in, R, R_b_out, n, k, m);
        if (active_A) {
          mapEigen<eigenStore>(A_b_in, n, k) =
              mapEigen<eigenStore>(R_b_out, n, m) * mapEigen<eigenStore>(B, k, m).transpose();
        }
        if (active_B) {
          mapEigen<eigenStore>(B_b_in, k, m) =
              mapEigen<eigenStore>(A, n, k).transpose() * mapEigen<eigenStore>(R_b_out, n, m);
        }
      }

      /// Function for deletion of contents.
      CODI_INLINE static void del(Tape* tape, codi::ByteDataView& dataStore) {
        codi::TemporaryMemory& allocator = tape->getTemporaryMemory();
        codiAssert(allocator.isEmpty());  // No memory should be allocated. We would free it at the end.
        using LLFH = codi::LowLevelFunctionCreationUtilities<2>;

        // Traits for arguments.
        using Trait_A = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_B = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_R = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_n = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_k = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_m = typename LLFH::PassiveStoreTrait<int, uint8_t>;

        // Declare variables.
        typename LLFH::ActivityStoreType activityStore = {};
        typename Trait_A::ArgumentStore A_store = {};
        typename Trait_B::ArgumentStore B_store = {};
        typename Trait_R::ArgumentStore R_store = {};
        typename Trait_n::Store n = {};
        typename Trait_k::Store k = {};
        typename Trait_m::Store m = {};

        bool active_A = false;
        bool active_B = false;

        // Restore data.
        LLFH::restoreActivity(&dataStore, activityStore);
        active_A = LLFH::getActivity(activityStore, 0);
        active_B = LLFH::getActivity(activityStore, 1);
        Trait_n::restore(&dataStore, allocator, 1, true, n);
        Trait_k::restore(&dataStore, allocator, 1, true, k);
        Trait_m::restore(&dataStore, allocator, 1, true, m);
        Trait_A::restore(&dataStore, allocator, n * k, LLFH::createRestoreActions(true, false, active_A, active_B),
                         A_store);
        Trait_B::restore(&dataStore, allocator, k * m, LLFH::createRestoreActions(true, false, active_B, active_A),
                         B_store);
        Trait_R::restore(&dataStore, allocator, n * m, LLFH::createRestoreActions(false, true, false, true), R_store);

        allocator.free();
      }

      /// Store on tape.
      CODI_INLINE static void evalAndStore(Type const* A, Type const* B, Type* R, int n, int k, int m) {
        Tape& tape = Type::getTape();
        codi::TemporaryMemory& allocator = tape.getTemporaryMemory();
        using LLFH = codi::LowLevelFunctionCreationUtilities<2>;

        // Traits for arguments.
        using Trait_A = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_B = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_R = typename LLFH::ActiveStoreTrait<Type*>;
        using Trait_n = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_k = typename LLFH::PassiveStoreTrait<int, uint8_t>;
        using Trait_m = typename LLFH::PassiveStoreTrait<int, uint8_t>;

        // Declare variables.
        typename LLFH::ActivityStoreType activityStore = {};
        typename Trait_A::ArgumentStore A_store = {};
        typename Trait_B::ArgumentStore B_store = {};
        typename Trait_R::ArgumentStore R_store = {};

        // Detect activity.
        bool active_A = Trait_A::isActive(A, n * k);
        bool active_B = Trait_B::isActive(B, k * m);
        bool active = active_A | active_B;

        if (active) {
          // Store function.
          registerOnTape();

          // Count data size.
          size_t dataSize = LLFH::countActivitySize();
          dataSize += Trait_n::countSize(n, 1, true);
          dataSize += Trait_k::countSize(k, 1, true);
          dataSize += Trait_m::countSize(m, 1, true);
          dataSize += Trait_A::countSize(A, n * k, LLFH::createStoreActions(active, true, false, active_A, active_B));
          dataSize += Trait_B::countSize(B, k * m, LLFH::createStoreActions(active, true, false, active_B, active_A));
          dataSize += Trait_R::countSize(R, n * m, LLFH::createStoreActions(active, false, true, false, true));

          // Reserve data.
          codi::ByteDataView dataStore = {};
          tape.pushLowLevelFunction(ID, dataSize, dataStore);

          // Store data.
          LLFH::setActivity(activityStore, 0, active_A);
          LLFH::setActivity(activityStore, 1, active_B);
          LLFH::storeActivity(&dataStore, activityStore);
          Trait_n::store(&dataStore, allocator, n, 1, true);
          Trait_k::store(&dataStore, allocator, k, 1, true);
          Trait_m::store(&dataStore, allocator, m, 1, true);
          Trait_A::store(&dataStore, allocator, A, n * k,
                         LLFH::createStoreActions(active, true, false, active_A, active_B), A_store);
          Trait_B::store(&dataStore, allocator, B, k * m,
                         LLFH::createStoreActions(active, true, false, active_B, active_A), B_store);
          Trait_R::store(&dataStore, allocator, R, n * m, LLFH::createStoreActions(active, false, true, false, true),
                         R_store);
        } else {
          // Prepare passive evaluation.
          Trait_A::store(nullptr, allocator, A, n * k,
                         LLFH::createStoreActions(active, true, false, active_A, active_B), A_store);
          Trait_B::store(nullptr, allocator, B, k * m,
                         LLFH::createStoreActions(active, true, false, active_B, active_A), B_store);
          Trait_R::store(nullptr, allocator, R, n * m, LLFH::createStoreActions(active, false, true, false, true),
                         R_store);
        }

        callPrimal(active, A_store.primal(), active_A, A_store.identifierIn(), B_store.primal(), active_B,
                   B_store.identifierIn(), R_store.primal(), R_store.identifierOut(), n, k, m);

        Trait_R::setExternalFunctionOutput(active, R, n * m, R_store.identifierOut(), R_store.primal(),
                                           R_store.oldPrimal());

        allocator.free();
      }

      /// Primal computation function.
      CODI_INLINE static void callPrimal(bool active, typename codi::ActiveArgumentStoreTraits<Type*>::Real const* A,
                                         bool active_A,
                                         typename codi::ActiveArgumentStoreTraits<Type*>::Identifier const* A_i_in,
                                         typename codi::ActiveArgumentStoreTraits<Type*>::Real const* B, bool active_B,
                                         typename codi::ActiveArgumentStoreTraits<Type*>::Identifier const* B_i_in,
                                         typename codi::ActiveArgumentStoreTraits<Type*>::Real* R,
                                         typename codi::ActiveArgumentStoreTraits<Type*>::Identifier* R_i_out, int n,
                                         int k, int m) {
        codi::CODI_UNUSED(active, A, active_A, A_i_in, B, active_B, B_i_in, R, R_i_out, n, k, m);
        mapEigen<eigenStore>(R, n, m) = mapEigen<eigenStore>(A, n, k) * mapEigen<eigenStore>(B, k, m);

        // User defined activity update.
        if (active) {
          mapEigen<eigenStore>(R_i_out, n, m).setZero();
          if (active_A) {
            mapEigen<eigenStore>(R_i_out, n, m).colwise() += mapEigen<eigenStore>(A_i_in, n, k).rowwise().any();
          }
          if (active_B) {
            mapEigen<eigenStore>(R_i_out, n, m).rowwise() += mapEigen<eigenStore>(B_i_in, k, m).colwise().any();
          }
        }
      }

      /// Register function on tape.
      CODI_INLINE static void registerOnTape() {
        if (codi::Config::LowLevelFunctionTokenInvalid == ID) {
          using Entry = codi::LowLevelFunctionEntry<Tape, typename Type::Real, typename Type::Identifier>;
          ID = Type::getTape().registerLowLevelFunction(Entry(reverse, forward, nullptr, del));
        }
      }
  };

  template<Eigen::StorageOptions eigenStore, typename Type>
  codi::Config::LowLevelFunctionToken ExtFunc_matrixMatrixMultiplication<eigenStore, Type>::ID =
      codi::Config::LowLevelFunctionTokenInvalid;

  /**
   *  Low level function for \f$R = A * B\f$ with
   *   - \f$ R \in \R^{n \times m} \f$
   *   - \f$ A \in \R^{n \times k} \f$
   *   - \f$ B \in \R^{k \times m} \f$
   *
   * @tparam eigenStore One of Eigen::StorageOptions.
   */
  template<Eigen::StorageOptions eigenStore, typename Type>
  void matrixMatrixMultiplication(Type const* A, Type const* B, Type* R, int n, int k, int m) {
    ExtFunc_matrixMatrixMultiplication<eigenStore, Type>::evalAndStore(A, B, R, n, k, m);
  }

  /**
   *  Low level function for \f$R = A * B\f$ with
   *   - \f$ R \in \R^{n \times m} \f$
   *   - \f$ A \in \R^{n \times k} \f$
   *   - \f$ B \in \R^{k \times m} \f$
   */
  template<typename Type>
  void matrixMatrixMultiplicationRowMajor(Type const* A, Type const* B, Type* R, int n, int k, int m) {
    matrixMatrixMultiplication<Eigen::StorageOptions::RowMajor>(A, B, R, n, k, m);
  }
  /**
   *  Low level function for \f$R = A * B\f$ with
   *   - \f$ R \in \R^{n \times m} \f$
   *   - \f$ A \in \R^{n \times k} \f$
   *   - \f$ B \in \R^{k \times m} \f$
   */
  template<typename Type>
  void matrixMatrixMultiplicationColMajor(Type const* A, Type const* B, Type* R, int n, int k, int m) {
    matrixMatrixMultiplication<Eigen::StorageOptions::ColMajor>(A, B, R, n, k, m);
  }
}
