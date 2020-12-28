/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include "dummyValue.hpp"
#include "vectorStorage.hpp"
#include "../../typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {


  /**
   * @brief Basic interface definition what the algorithms in the EvaluationHelper, Algorithms, etc. classes use to
   * access the data.
   *
   * @tparam T  The data type of the internal storage.
   */
  template<typename T>
  struct JacobianInterface {

      virtual ~JacobianInterface(){}

      /**
       * @brief Constant access to the specified element of the Jacobian.
       *
       * @param[in] i  Output value of the function. Range: [0, m)
       * @param[in] j  Input value of the function. Range: [0,n)
       *
       * @return Value/reference of the specified location.
       */
      virtual T operator()(const size_t i, const size_t j) const = 0;

      /**
       * @brief Access to the specified element of the hessian.
       *
       * \copydetails operator()(const size_t i, const size_t j) const
       */
      virtual T& operator()(const size_t i, const size_t j) = 0;
  };


  /**
   * @brief A dummy hessian that can be accessed via the call operator.
   */
  struct DummyJacobian {

      /**
       * @brief Returns a dummy value.
       *
       * @param[in] i  Not used
       * @param[in] j  Not used
       * @return A dummy value
       */
      CODI_INLINE DummyValue operator()(const size_t i, const size_t j) const {
        CODI_UNUSED(i);
        CODI_UNUSED(j);

        return DummyValue();
      }
  };

  /**
   * @brief Default Jacobian implementation for algorithms in CoDiPack.
   *
   * The data layout of the Jacobian is described in the \ref FunctionDefinition "Algorithms" documentation.
   * All mathematical symbol names are described there.
   *
   * Structures that implement the same functions can be used a the same places where this Jacobian implementation
   * is used in CoDiPack.
   *
   * @tparam Vec  Either std::vector or std::array. For other types VectorStorage needs to be specialized.
   */
  template <typename Vec>
  struct Jacobian {

    private:
      VectorStorage<Vec> values; /**< Storage for the hessian matrix */

      size_t m; /**< Number of function outputs */
      size_t n; /**< Number of function inputs */

    public:

      /**
       * @brief Inner element of the vector type.
       */
      using T = typename VectorStorage<Vec>::Element;

      /**
       * @brief Create a Jacobian with the given output and input size.
       *
       * @param[in] m  Number of function outputs
       * @param[in] n  Number of function inputs
       */
      explicit Jacobian(size_t m, size_t n) : values(n * m), m(m), n(n) {}

      /**
       * @brief Get the number of function outputs.
       * @return Number of function outputs.
       */
      CODI_INLINE size_t getM() const {
        return m;
      }

      /**
       * @brief Get the number of function inputs.
       * @return Number of function inputs.
       */
      CODI_INLINE size_t getN() const {
        return n;
      }

      /**
       * @brief Constant access to the specified element of the Jacobian.
       *
       * @param[in] i  Output value of the function. Range: [0, m) (Slowest running index)
       * @param[in] j  Input value of the function. Range: [0,n) (Fastest running index)
       *
       * @return Value/reference of the specified location.
       */
      CODI_INLINE T operator()(const size_t i, const size_t j) const {
        return values.data()[computeIndex(i,j)];
      }

      /**
       * @brief Access to the specified element of the hessian.
       *
       * \copydetails operator()(const size_t i, const size_t j) const
       */
      CODI_INLINE T& operator()(const size_t i, const size_t j) {
        return values.data()[computeIndex(i,j)];
      }

      /**
       * @brief Change the size of the Jacobian to the new input and output values.
       *
       * @param[in] m  Number of function outputs
       * @param[in] n  Number of function inputs
       */
      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m*n);
      }

      /**
       * @brief Change the size of the Jacobian to the new input and output values.
       *
       * Keeps  the size of the internal value vector if new shape is smaller than the old one.
       *
       * @param[in] m  Number of function outputs
       * @param[in] n  Number of function inputs
       */
      CODI_INLINE void reshape(size_t const m, size_t const n) {
        if(values.size() < m * n) {
          this->resize(m, n);
        } else {
          this->m = m;
          this->n = n;
        }
      }

      /**
       * @brief The number of entries of the Jacobian.
       *
       * @return The size of the Jacobian.
       */
      CODI_INLINE size_t size() const {
        return m * n;
      }

    protected:

      /**
       * @brief Computes the offset into the data array.
       *
       * @param[in] i  n. Slowest running index
       * @param[in] j  1. Fastest running index
       * @return
       */
      CODI_INLINE size_t computeIndex(const size_t i, const size_t j) const {
        return i * n + j;
      }
  };


  /**
   * @brief Helper structure for the delayed access to the data element.
   *
   * @tparam Impl  The implementation for the data access logic. This type needs to implement setLogic<T>(i,j,v)
   */
  template <typename Impl>
  struct DelayAccessor {

    private:
      size_t i; /**< Output value position. Range: [0, m) */
      size_t j; /**< Input value position. Range: [0, n) */

      Impl& data; /**< Reference to the access logic */

    public:

      /**
       * @brief Default constructor
       * @param[in]    i  Output value position. Range: [0, m)
       * @param[in]    j  Input value position. Range: [0, n)
       * @param[in] data  Reference to the access logic
       */
      DelayAccessor(size_t const i, size_t const j, Impl& data) : i(i), j(j), data(data) {}

      /**
       * @brief Forwards the set operation to the logic implementation.
       *
       * @param[in] v  Value which is set into the logic.
       *
       * @return This
       *
       * @tparam T  The type of the right hand side argument.
       */
      template<typename T>
      DelayAccessor& operator =(T const& v) {
        data.setLogic(i,j, v);

        return *this;
      }

      /**
       * @brief Conversion to the original element
       */
      operator typename Impl::T() const {
        return const_cast<const Impl&>(data).operator()(i,j);
      }
  };

  /**
   * @brief Adds the count for non zero elements in each row.
   *
   * The logic for the class is quite simple. Every time a value is set the implementation checks if the
   * new value is not zero. If so the count for the row is increased by one. The user has to reset the count
   * of non zero entries when the Jacobian is filled again.
   *
   * @tparam Vec Forwarded to Jacobian
   */
  template <typename Vec>
  struct JacobianCountNonZerosRow : public Jacobian<Vec> {

      /**
       * @brief Number of non zero entries per row.
       */
      std::vector<int> nonZerosRow;

    public:
      /** \copydoc Jacobian::T */
      using T = typename Jacobian<Vec>::T;

      /** @brief Type of the delayed accessor object */
      using DelayAcc = DelayAccessor<JacobianCountNonZerosRow<Vec>>;

      /** \copydoc Jacobian::Jacobian */
      explicit JacobianCountNonZerosRow(size_t m, size_t n) : Jacobian<Vec>(m, n), nonZerosRow(m) {}


      /** \copydoc Jacobian::operator()(const size_t i, const size_t j) const */
      CODI_INLINE T operator()(const size_t i, const size_t j) const {
        return Jacobian<Vec>::operator ()(i, j);
      }

      /**
       * @brief Returns an object for the delayed access. This object calls then the access logic which updates the
       * number of non zero elements.
       *
       * \copydetails Jacobian::operator()(const size_t i, const size_t j)
       */
      CODI_INLINE DelayAcc operator()(const size_t i, const size_t j) {
        return DelayAcc(i,j, *this);
      }

      /** \copydoc Jacobian::resize */
      CODI_INLINE void resize(size_t const m, size_t const n) {

        Jacobian<Vec>::resize(m, n);
        nonZerosRow.resize(m);
      }

      /** \copydoc Jacobian::reshape */
      CODI_INLINE void reshape(size_t const m, size_t const n) {
        Jacobian<Vec>::reshape(m, n);
        if(nonZerosRow.size() < m) {
          nonZerosRow.resize(m);
        }
      }
      /**
       * @brief Reference to the number of non zero entries for the specified row.
       *
       * @param[in] i  Output value position. Range: [0, m)
       *
       * @return Reference to the non zero count for the specified row.
       */
      CODI_INLINE int& nonZeroRow(const size_t i) {
        return nonZerosRow[i];
      }

      /**
       * @brief Checks if the element is non zero and updates the counter.
       *
       * @param[in] i  Default
       * @param[in] j  Default
       * @param[in] v  Default
       */
      CODI_INLINE void setLogic(const size_t i, const size_t j, T const& v) {
        if(T() != v) {
          nonZerosRow[i] += 1;
        }
        Jacobian<Vec>::operator ()(i,j) = v;
      }
  };

  /**
   * @brief A wrapper around a Jacobian that converts values when they are set into the Jacobian matrix.
   *
   * Converts the right hand side value with TypeTraits::getBaseValue.
   *
   * @tparam Nested  A type that implements the Jacobian interface.
   */
  template <typename Nested>
  struct JacobianConvertWrapper {

    private:
      Nested& nested;

    public:
      /** \copydoc Jacobian::T */
      using T = typename Nested::T;

      /** @brief Type of the delayed accessor object */
      using DelayAcc = DelayAccessor<JacobianConvertWrapper<Nested>>;

      /**
       * @brief Constructor
       * @param[in] nested Reference to the wrapped Jacobian
       */
      explicit JacobianConvertWrapper(Nested& nested) : nested(nested) {}

      /** \copydoc Jacobian::operator()(const size_t i, const size_t j) const */
      CODI_INLINE T operator()(const size_t i, const size_t j) const {
        return nested(i, j);
      }

      /**
       * @brief Returns an object for the delayed access. This object calls then the access logic which converts the
       * value on the right hand side.
       *
       * \copydetails Jacobian::operator()(const size_t i, const size_t j)
       */
      CODI_INLINE DelayAcc operator()(const size_t i, const size_t j) {
        return DelayAcc(i,j, *this);
      }

      /**
       * @brief Converts the right hand side value with TypeTraits::getBaseValue
       *
       * @param[in] i  Default
       * @param[in] j  Default
       * @param[in] v  Default
       *
       * @tparam SetT  The type of the right hand side value.
       */
      template<typename SetT>
      CODI_INLINE void setLogic(const size_t i, const size_t j, SetT const& v) {
        nested(i, j) = TypeTraits<SetT>::getBaseValue(v);
      }
  };
}
