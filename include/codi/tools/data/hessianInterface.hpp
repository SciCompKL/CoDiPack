#pragma once

#include <vector>

#include "../../aux/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * See \ref sec_namingConventions for the naming conventions.
   *
   * In the access function:
   *  - i: Output value of the function. Range: [0, m)
   *  - j: Input value of the function. First derivative direction. Range: [0,n)
   *  - k: Input value of the function. Second derivative direction. Range: [0,n)
   *
   * @tparam _T  Data type.
   */
  /**
   * @brief General interface for Hessian access in CoDiPack.
   *
   * This interface needs to implemented for data that is given to helper methods which store or read data from
   * a Hessian.
   *
   * See \ref sec_namingConventions for the mathematical nomenclature of the arguments and components.
   *
   * @tparam _T  The data type in the Hessian.
   */
  template<typename _T>
  struct HessianInterface {
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double);  ///< See HessianInterface

      virtual ~HessianInterface() {}  ///< Destructor

      size_t getM() const;  ///< Get size of rows (Output variables)
      size_t getN() const;  ///< Get size of columns  (Input variables)

      /**
       * @brief Value access.
       *
       * @param[in] i  Output value of the function. Range: [0, m)
       * @param[in] j  Input value of the function. First derivative direction. Range: [0,n)
       * @param[in] k  Input value of the function. Second derivative direction. Range: [0,n)
       */
      T operator()(size_t const i, size_t const j, size_t const k) const;

      /**
       * @brief Reference access.
       *
       * \copydetails operator()(size_t const i, size_t const j, size_t const k) const
       */
      T& operator()(size_t const i, size_t const j, size_t const k);

      void resize(size_t const m, size_t const n);  ///< Resize to the new dimensions
      size_t size() const;                          ///< Get total size of the Hessian.
  };
}
