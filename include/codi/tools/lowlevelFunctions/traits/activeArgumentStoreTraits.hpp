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

#include "../../../config.h"
#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../misc/byteDataView.hpp"
#include "../../../misc/macros.hpp"
#include "../../../misc/temporaryMemory.hpp"
#include "../../../tapes/data/position.hpp"
#include "../../../tapes/interfaces/fullTapeInterface.hpp"
#include "../storeAndRestoreActions.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Declares all variables that may be needed to store/restore an active argument which has a pointer type.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.
   */
  template<typename T_Real, typename T_Identifier, typename T_Gradient>
  struct ActiveArgumentPointerStore {
      using Real = CODI_DD(T_Real, double);           ///< See ActiveArgumentPointerStore.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See ActiveArgumentPointerStore.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See ActiveArgumentPointerStore.

      Real* value_v;              ///< Primal value vector.
      Identifier* value_i_in;     ///< Identifier vector of an input value.
      Identifier* value_i_out;    ///< Identifier vector of an output value.
      Gradient* value_deriv_in;   ///< Gradient vector of an input value.
      Gradient* value_deriv_out;  ///< Gradient vector of an output value.

      Real* oldPrimals;  ///< Old primal values in primal value tape setting.

      int passiveValuesCount;  ///< Number of passive values.

      // clang-format off
      Real* value() { return value_v; }                    ///< Get the primal values.
      Identifier* identifierIn() { return  value_i_in; }   ///< Get the input identifiers.
      Identifier* identifierOut() { return value_i_out; }  ///< Get the output identifiers.
      Gradient* gradientIn() { return value_deriv_in; }    ///< Get the input gradients.
      Gradient* gradientOut() { return value_deriv_out; }  ///< Get the output gradients.

      Real* oldPrimal() { return oldPrimals; }  ///< Get old primal values.
      // clang-format on
  };

  /**
   * @brief Declares all variables that may be needed to store/restore an active argument which has a value type.
   *
   * @tparam T_PointerStore  The pointer store implementation.
   */
  template<typename T_PointerStore>
  struct ActiveArgumentValueStore {
      /// See ActiveArgumentValueStore.
      using PointerStore = CODI_DD(T_PointerStore, CODI_T(ActiveArgumentPointerStore<double, int, double>));

      using Real = typename PointerStore::Real;              ///< See ActiveArgumentPointerStore.
      using Identifier = typename PointerStore::Identifier;  ///< See ActiveArgumentPointerStore.
      using Gradient = typename PointerStore::Gradient;      ///< See ActiveArgumentPointerStore.

      PointerStore base;  ///< Declaration of base.

      // clang-format off
      Real& value() { return *base.value(); }                        ///< Get the primal values.
      Identifier& identifierIn() { return  *base.identifierIn(); }   ///< Get the output identifiers.
      Identifier& identifierOut() { return *base.identifierOut(); }  ///< Get the input gradients.
      Gradient& gradientIn() { return *base.gradientIn(); }          ///< Get the output gradients.
      Gradient& gradientOut() { return *base.gradientOut(); }        ///< Get old primal values.

      Real& oldPrimal() { return *base.oldPrimal(); }  ///< See ActiveArgumentPointerStore.
      // clang-format on
  };

  /**
   * @brief Traits for storing active arguments in byte streams.
   *
   * There is one data stream available for storing the data of a type. The size and layout of the data can be
   * arbitrary. The only limitation is the maximum size of the data which is defined by
   * #config::LowLevelFunctionDataSizeMax. This limit holds for the whole external function. If larger data is required
   * it should be allocated with new.
   *
   * \section layout Example data layout for matrix matrix multiplications.
   *
   * The required data for a matrix is the size, the values of the entries, and the identifiers from CoDiPack. The data
   * could then be stored like:
   * \code
   * stream data: | n | m | n * m values | n * m identifiers |
   * \endcode
   * The sizes \c n and \c m are the necessary information for the dynamic data. Therefore, they are stored first.
   * Afterwards the dynamic sized data is stored.
   *
   * \section calls Call process
   *
   * \subsection calls_record Recording
   *
   * The first function, that is usually called, is the #isActive function. It should return true if one CoDiPack
   * entry in the type is active. Afterwards, the #countSize function is called to determine the required byte
   * size for the data stream. This count needs to be exact since the data is preallocated and can not be shortened
   * afterwards. The call to #store initiates the storing of the type in the stream and therefore on the tape. The
   * stored data needs to be smaller or equal to the size returned in #countSize. Finally, the
   * #setExternalFunctionOutput method is called on output arguments. Here, the vectors created in #store need to be
   * updated if necessary. Usually, this are vectors for the identifiers of the output value and old primal values.
   *
   * \subsection calls_eval Tape evaluation
   * During a reverse, forward, primal, etc. evaluation of a tape, the #restore method is called first. It should read
   * all data written in #store.
   *
   * Afterwards the functions #getPrimalsFromVector, #setPrimalsIntoVector, #getGradients and #setGradients are called
   * to populate the vectors created in the restore functions.
   *
   * \subsection store_actions Store actions
   *  - StoreAction::PrimalCreateOnTape : Create a primal vector in the data streams.
   *  - StoreAction::PrimalExtract : Extract the primal form the value. If PrimalCreateOnTape is not requested create
   *                                 a temporary vector for it.
   *  - StoreAction::InputIdentifierCreateAndStore : Create a vector for the input identifiers and store the current
   *                                                 ones from the value.
   *  - StoreAction::OutputIdentifierCreate : Create a vector for the output identifiers. They are populated
   *                                          during the call to #setExternalFunctionOutput after the low level function
   *                                          hast been evaluated.
   *
   * \subsection restore_actions Restore actions
   *  - RestoreAction::PrimalCreate : Create a vector for the primal output values of this argument.
   *  - RestoreAction::PrimalRestore : Restore the primal input values from the data stream.
   *  - RestoreAction::InputIdentifierRestore : Read a vector for the input identifiers from the streams.
   *  - RestoreAction::OutputIdentifierRestore : Read a vector for the output identifiers from the streams.
   *  - RestoreAction::InputGradientCreate : Create a vector for the input gradients. They are then populated/used
   *                                         during a call to #getGradients/#setGradients.
   *  - RestoreAction::OutputGradientCreate : Create a vector for the output gradients. They are then populated/used
   *                                          during a call to #getGradients/#setGradients.
   *
   *
   * @tparam T_T A type that contains CoDiPack values.
   */
  template<typename T_T, typename = void>
  struct ActiveArgumentStoreTraits {
      using T = CODI_DD(T_T, CODI_DEFAULT_LHS_EXPRESSION);  ///< See ActiveArgumentStoreTraits.

      using Real = double;     ///< The type with no CoDiPack values.
      using Identifier = int;  ///< The type for holding the identifiers.
      using Gradient = Real;   ///< The type that can represent the gradient values.

      /// Data for holding all necessary values.
      using ArgumentStore = ActiveArgumentPointerStore<Real, Identifier, Gradient>;

      /**
       * @brief Counts the binary size for the  data stream.
       *
       * This value needs to define the maximum size required to store all data for the type. It should be
       * exact since the allocated memory can not be reduced afterwards.
       *
       * \c actions describe what needs to be done for this argument. \c size is a hint for the implementation, e.g. for
       * pointers the vector size.
       */
      CODI_INLINE static void countSize(size_t& dataSize, T const& value, size_t size, StoreActions const& actions);

      /**
       * @brief Restore the data for this type.
       *
       * \c actions describe what needs to be done for this argument. \c size is a hint for the implementation, e.g. for
       * pointers the vector size. \c data can be used to store data and pointers from the streams.
       */
      CODI_INLINE static void restore(ByteDataView* store, TemporaryMemory& allocator, size_t size,
                                      RestoreActions const& actions, ArgumentStore& data);

      /**
       * @brief Store all data for this type.
       *
       * The amount of data can not be greater than the amount reported by #countSize().
       *
       * \c actions describe what needs to be done for this argument. \c size is a hint for the implementation, e.g. for
       * pointers the vector size. \c data can be used to store data and pointers from the streams.
       */
      CODI_INLINE static void store(ByteDataView* dataStore, TemporaryMemory& allocator, T const& value, size_t size,
                                    StoreActions const& actions, ArgumentStore& data);

      /// Should return true when one element in the type is active.
      CODI_INLINE static bool isActive(T const& value, size_t size);

      /**
       *  Called after the primal evaluation. All active values in \c value need to be registered as outputs of an
       *  external function. \c value needs to be populated with the primal values from \c primal. If \c primal is
       *  \c null then the function can assume that \c value already contains the current values. The identifiers need
       *  to be stored in \c identifiers.
       */
      CODI_INLINE static void setExternalFunctionOutput(bool tapeActive, T& value, size_t size, Identifier& identifier,
                                                        Real& primal, Real& oldPrimals);

      /// Get the primal values from \c data and store them in \c primal.
      CODI_INLINE static void getPrimalsFromVector(VectorAccessInterface<Real, Identifier>* data, size_t size,
                                                   Identifier& identifier, Real& primal);

      /// Set the primal values from \c primal into \c data.
      CODI_INLINE static void setPrimalsIntoVector(VectorAccessInterface<Real, Identifier>* data, size_t size,
                                                   Identifier& identifier, Real& primal);

      /// Get the gradients from \c data and store them in \c gradient.
      CODI_INLINE static void getGradients(VectorAccessInterface<Real, Identifier>* data, size_t size, bool reset,
                                           Identifier& identifier, Gradient& gradient, size_t dim);

      /// Set the gradients from \c gradient into \c data.
      CODI_INLINE static void setGradients(VectorAccessInterface<Real, Identifier>* data, size_t size, bool update,
                                           Identifier& identifier, Gradient& gradient, size_t dim);
  };

#ifndef DOXYGEN_DISABLE

  /// Specialization of ActiveArgumentStoreTraits for value types.
  template<typename T_T>
  struct ActiveArgumentStoreTraits<T_T, typename std::enable_if<!std::is_pointer<T_T>::value>::type> {
      using T = CODI_DD(T_T, CODI_DEFAULT_LHS_EXPRESSION);  ///< See ActiveArgumentStoreTraits.

      using PointerTraits = ActiveArgumentStoreTraits<T*>;  ///< Traits for the pointer version of the type.

      using Real = typename PointerTraits::Real;              ///< See ActiveArgumentStoreTraits.
      using Identifier = typename PointerTraits::Identifier;  ///< See ActiveArgumentStoreTraits.
      using Gradient = typename PointerTraits::Gradient;      ///< See ActiveArgumentStoreTraits.

      using ArgumentStore =
          ActiveArgumentValueStore<typename PointerTraits::ArgumentStore>;  ///< See ActiveArgumentStoreTraits.

      /// @copydoc ActiveArgumentValueStore::countSize()
      CODI_INLINE static void countSize(size_t& dataSize, T const& value, size_t size, StoreActions const& actions) {
        CODI_UNUSED(size);

        PointerTraits::countSize(dataSize, &value, 1, actions);
      }

      /// @copydoc ActiveArgumentValueStore::restore()
      CODI_INLINE static void restore(ByteDataView* store, TemporaryMemory& allocator, size_t size,
                                      RestoreActions const& actions, ArgumentStore& data) {
        CODI_UNUSED_ARG(size);

        PointerTraits::restore(store, allocator, 1, actions, data.base);
      }

      /// @copydoc ActiveArgumentValueStore::store()
      CODI_INLINE static void store(ByteDataView* dataStore, TemporaryMemory& allocator, T const& value, size_t size,
                                    StoreActions const& actions, ArgumentStore& data) {
        CODI_UNUSED_ARG(size);

        PointerTraits::store(dataStore, allocator, value, 1, actions, data.base);
      }

      /// @copydoc ActiveArgumentValueStore::isActive()
      CODI_INLINE static bool isActive(T const& value, size_t size) {
        CODI_UNUSED(size);

        return PointerTraits::isActive(&value, 1);
      }

      /// @copydoc ActiveArgumentValueStore::setExternalFunctionOutput()
      CODI_INLINE static void setExternalFunctionOutput(bool tapeActive, T& value, size_t size, Identifier& identifier,
                                                        Real& primal, Real& oldPrimals) {
        CODI_UNUSED_ARG(size);

        PointerTraits::setExternalFunctionOutput(tapeActive, &value, 1, &identifier, &primal, &oldPrimals);
      }

      /// @copydoc ActiveArgumentValueStore::getPrimalsFromVector()
      CODI_INLINE static void getPrimalsFromVector(VectorAccessInterface<Real, Identifier>* data, size_t size,
                                                   Identifier& identifier, Real& primal) {
        CODI_UNUSED_ARG(size);

        PointerTraits::getPrimalsFromVector(data, 1, &identifier, &primal);
      }

      /// @copydoc ActiveArgumentValueStore::setPrimalsIntoVector()
      CODI_INLINE static void setPrimalsIntoVector(VectorAccessInterface<Real, Identifier>* data, size_t size,
                                                   Identifier& identifier, Real& primal) {
        CODI_UNUSED_ARG(size);

        PointerTraits::setPrimalsIntoVector(data, 1, &identifier, &primal);
      }

      /// @copydoc ActiveArgumentValueStore::getGradients()
      CODI_INLINE static void getGradients(VectorAccessInterface<Real, Identifier>* data, size_t size, bool reset,
                                           Identifier& identifier, Gradient& gradient, size_t dim) {
        CODI_UNUSED_ARG(size);

        PointerTraits::getGradients(data, 1, reset & identifier, &gradient, dim);
      }

      /// @copydoc ActiveArgumentValueStore::setGradients()
      CODI_INLINE static void setGradients(VectorAccessInterface<Real, Identifier>* data, size_t size, bool update,
                                           Identifier& identifier, Gradient& gradient, size_t dim) {
        CODI_UNUSED_ARG(size);

        PointerTraits::setGradients(data, 1, update, &identifier, &gradient, dim);
      }
  };

  /// Specialization of ActiveArgumentStoreTraits for pointers to CoDiPack values.
  template<typename T_T>
  struct ActiveArgumentStoreTraits<T_T*, ExpressionTraits::EnableIfLhsExpression<T_T>> {
      using T = CODI_DD(T_T, CODI_DEFAULT_LHS_EXPRESSION);  ///< See ActiveArgumentStoreTraits.

      using Tape = typename T::Tape;  ///< Tape of active type.

      using Real = typename T::Real;              ///< See ActiveArgumentStoreTraits.
      using Identifier = typename T::Identifier;  ///< See ActiveArgumentStoreTraits.
      using Gradient = Real;  ///< We use the vector access interface so gradients are stored in the same form as the
                              ///< primal.

      using ArgumentStore = ActiveArgumentPointerStore<Real, Identifier, Gradient>;  ///< See ActiveArgumentStoreTraits.

      /// @copydoc ActiveArgumentValueStore::countSize()
      CODI_INLINE static void countSize(size_t& dataSize, T const* value, size_t size, StoreActions const& actions) {
        CODI_UNUSED(dataSize, value);
        if (actions.test(StoreAction::InputIdentifierCreateAndStore)) {
          dataSize += size * sizeof(typename T::Identifier);  // var_i_in
        }
        if (actions.test(StoreAction::PrimalCreateOnTape)) {
          if (Tape::HasPrimalValues) {
            // Primal value tapes only store the passive primal values.
            int passiveIdentifiers = countPassive(value, size);
            dataSize += sizeof(int);  // Number of passive identifiers
            dataSize += passiveIdentifiers * sizeof(typename T::Real);
          } else {
            // Jacobian tape stores full primal values.
            dataSize += size * sizeof(typename T::Real);  // var_v_in
          }
        }
        if (actions.test(StoreAction::OutputIdentifierCreate)) {
          dataSize += size * sizeof(typename T::Identifier);  // var_i_out

          if (Tape::HasPrimalValues && Tape::HasPrimalValues && !Tape::LinearIndexHandling) {
            dataSize += size * sizeof(Real);  // Primal value tapes need to store the old values.
          }
        }
      }

      /// @copydoc ActiveArgumentValueStore::restore()
      CODI_INLINE static void restore(ByteDataView* store, TemporaryMemory& allocator, size_t size,
                                      RestoreActions const& actions, ArgumentStore& data) {
        Real* passiveValues = nullptr;

        if (Tape::HasPrimalValues && actions.test(RestoreAction::PrimalRestore)) {
          data.passiveValuesCount = store->template read<int>();
        }

        if (actions.test(RestoreAction::PrimalRestore)) {
          restoreValue(store, allocator, size, data, passiveValues);
        }
        if (actions.test(RestoreAction::InputIdentifierRestore)) {
          data.value_i_in = store->template read<Identifier>(size);
        }
        if (actions.test(RestoreAction::OutputIdentifierRestore)) {
          data.value_i_out = store->template read<Identifier>(size);
          if (Tape::HasPrimalValues && !Tape::LinearIndexHandling) {
            data.oldPrimals = store->template read<Real>(size);
          }
        }

        if (actions.test(RestoreAction::PrimalRestore)) {
          restorePassiveValues(size, data, passiveValues);
        }

        if (actions.test(RestoreAction::PrimalCreate)) {
          data.value_v = allocator.template alloc<Real>(size);
        }

        if (actions.test(RestoreAction::InputGradientCreate)) {
          data.value_deriv_in = allocator.template alloc<Gradient>(size);
        }
        if (actions.test(RestoreAction::OutputGradientCreate)) {
          data.value_deriv_out = allocator.template alloc<Gradient>(size);
        }
      }

      /// Restores the primal values from the data.
      CODI_INLINE static void restoreValue(ByteDataView* store, TemporaryMemory& allocator, size_t size,
                                           ArgumentStore& data, Real*& passiveValues) {
        if (Tape::HasPrimalValues) {
          // Primal value tapes restore from the tape.
          passiveValues = store->template read<Real>(data.passiveValuesCount);
          if ((size_t)data.passiveValuesCount == size) {
            data.value_v = passiveValues;  // Use the full vector.
          } else {
            // Put the passive primal values at the positions where the identifiers are passive.
            data.value_v = allocator.template alloc<Real>(size);

            // Wait with the restore until the ids are also read. See restorePassiveValues().
          }
        } else {
          // Jacobian tapes read the full vector.
          data.value_v = store->template read<Real>(size);
        }
      }

      /// Copies the passive values into the primal value vector.
      CODI_INLINE static void restorePassiveValues(size_t size, ArgumentStore& data, Real* passiveValues) {
        typename T::Tape& tape = T::getTape();

        if (Tape::HasPrimalValues && (size_t)data.passiveValuesCount != size) {
          // Only restore if we did not use the full vector.
          int passivePos = 0;
          for (size_t i = 0; i < size; i += 1) {
            if (!tape.isIdentifierActive(data.value_i_in[i])) {
              data.value_v[i] = passiveValues[passivePos];
              passivePos += 1;
            }
          }
        }
      }

      /// @copydoc ActiveArgumentValueStore::store()
      CODI_INLINE static void store(ByteDataView* dataStore, TemporaryMemory& allocator, T const* value, size_t size,
                                    StoreActions const& actions, ArgumentStore& data) {
        Real* passiveValues = nullptr;
        if (actions.test(StoreAction::PrimalCreateOnTape)) {
          if (Tape::HasPrimalValues) {
            int passiveIdentifiers = countPassive(value, size);
            dataStore->write(passiveIdentifiers);
            passiveValues = dataStore->template reserve<Real>(passiveIdentifiers);
            data.value_v = allocator.template alloc<Real>(size);
          } else {
            // Jacobian tape stores full primal values.
            data.value_v = dataStore->template reserve<Real>(size);
          }

        } else {
          data.value_v = allocator.template alloc<Real>(size);
        }

        if (actions.test(StoreAction::PrimalExtract)) {
          typename T::Tape& tape = T::getTape();
          int passiveValuesPos = 0;

          for (size_t i = 0; i < size; i += 1) {
            (data.value_v)[i] = value[i].getValue();

            if (actions.test(StoreAction::PrimalCreateOnTape) && Tape::HasPrimalValues &&
                !tape.isIdentifierActive(value[i].getIdentifier())) {
              passiveValues[passiveValuesPos] = value[i].getValue();
              passiveValuesPos += 1;
            }
          }
        }

        if (actions.test(StoreAction::InputIdentifierCreateAndStore)) {
          data.value_i_in = dataStore->template reserve<Identifier>(size);
          for (size_t i = 0; i < size; i += 1) {
            (data.value_i_in)[i] = value[i].getIdentifier();
          }
        }

        if (actions.test(StoreAction::OutputIdentifierCreate)) {
          data.value_i_out = dataStore->template reserve<Identifier>(size);
          for (size_t i = 0; i < size; i += 1) {
            (data.value_i_out)[i] = -1;
          }

          if (Tape::HasPrimalValues && !Tape::LinearIndexHandling) {
            data.oldPrimals = dataStore->template reserve<Real>(size);
          }
        }
      }

      /// @copydoc ActiveArgumentValueStore::isActive()
      CODI_INLINE static bool isActive(T const* value, size_t size) {
        typename T::Tape& tape = T::getTape();
        bool active = true;
        for (size_t i = 0; i < size && active; i += 1) {
          active &= tape.isIdentifierActive(value[i].getIdentifier());
        }
        return active;
      }

      /// @copydoc ActiveArgumentValueStore::setExternalFunctionOutput()
      CODI_INLINE static void setExternalFunctionOutput(bool tapeActive, T* value, size_t size, Identifier* identifier,
                                                        Real* primal, Real* oldPrimals) {
        typename T::Tape& tape = T::getTape();
        for (size_t i = 0; i < size; i += 1) {
          if (nullptr != primal) {
            value[i].setValue(primal[i]);
          }

          if (tapeActive && 0 != identifier[i]) {
            Real oldValue = tape.registerExternalFunctionOutput(value[i]);
            identifier[i] = value[i].getIdentifier();

            if (Tape::HasPrimalValues && !Tape::LinearIndexHandling) {
              oldPrimals[i] = oldValue;
            }
          }
        }
      }

      /// @copydoc ActiveArgumentValueStore::getPrimalsFromVector()
      CODI_INLINE static void getPrimalsFromVector(VectorAccessInterface<Real, Identifier>* data, size_t size,
                                                   Identifier* identifier, Real* primal) {
        typename T::Tape& tape = T::getTape();

        for (size_t i = 0; i < size; i += 1) {
          if (tape.isIdentifierActive(identifier[i])) {
            primal[i] = data->getPrimal(identifier[i]);
          }
        }
      }

      /// @copydoc ActiveArgumentValueStore::setPrimalsIntoVector()
      CODI_INLINE static void setPrimalsIntoVector(VectorAccessInterface<Real, Identifier>* data, size_t size,
                                                   Identifier* identifier, Real* primal) {
        typename T::Tape& tape = T::getTape();

        for (size_t i = 0; i < size; i += 1) {
          if (tape.isIdentifierActive(identifier[i])) {
            data->setPrimal(identifier[i], primal[i]);
          }
        }
      }

      /// @copydoc ActiveArgumentValueStore::getGradients()
      CODI_INLINE static void getGradients(VectorAccessInterface<Real, Identifier>* data, size_t size, bool reset,
                                           Identifier* identifier, Gradient* gradient, size_t dim) {
        for (size_t i = 0; i < size; i += 1) {
          gradient[i] = data->getAdjoint(identifier[i], dim);
          if (reset) {
            data->resetAdjoint(identifier[i], dim);
          }
        }
      }

      /// @copydoc ActiveArgumentValueStore::setGradients()
      CODI_INLINE static void setGradients(VectorAccessInterface<Real, Identifier>* data, size_t size, bool update,
                                           Identifier* identifier, Real* gradient, size_t dim) {
        for (size_t i = 0; i < size; i += 1) {
          if (!update) {
            data->resetAdjoint(identifier[i], dim);
          }
          data->updateAdjoint(identifier[i], dim, gradient[i]);
        }
      }

    private:
      CODI_INLINE static int countPassive(T const* value, size_t size) {
        typename T::Tape& tape = T::getTape();

        int count = 0;
        for (size_t i = 0; i < size; i += 1) {
          if (!tape.isIdentifierActive(value[i].getIdentifier())) {
            count += 1;
          }
        }

        return count;
      }
  };
#endif
}
