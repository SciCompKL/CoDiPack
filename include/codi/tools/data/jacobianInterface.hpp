#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/aux/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief General interface for Jacobian access in CoDiPack.
   *
   * This interface needs to implemented for data that is given to helper methods which store or read data from
   * a Jacobian.
   *
   * See \ref sec_namingConventions for the mathematical nomenclature of the arguments and components.
   *
   * @tparam _T  The data type in the Jacobian.
   */
  template<typename _T>
  struct JacobianInterface {
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double); ///< See JacobianInterface

      size_t getM() const; ///< Get size of rows (Output variables)
      size_t getN() const; ///< Get size of columns  (Input variables)

      /// Value access, i in [0, ..., m), j in [0, ..., n)
      T operator()(size_t const i, size_t const j) const;

      /// Reference access, i in [0, ..., m), j in [0, ..., n)
      T& operator()(size_t const i, size_t const j);

      void resize(size_t const m, size_t const n); ///< Resize the Jacobian.
      size_t size() const; ///< Get total size of the Jacobian.
  };

  /**
   * Format the Jacobian and write it to the stream.
   *
   * The format is (Matlab):
   * [v11, v12, ..., v1N;
   *  v21, v22, ..., v2N;
   *  ...
   *  vM1, vM2, ..., vMN]
   */
  template<typename Stream, typename Jac, typename = enable_if_base_of<Jac, JacobianInterface<typename Jac::T>>>
  Stream& operator <<(Stream& out, CODI_DD(Jac, CODI_T(JacobianInterface<double>)) const& jacobian) {

    out << "[";
    for(size_t i = 0; i < jacobian.getM(); ++i) {
      if( i != 0) {
        out << " "; // Padding for the '['
      }

      for(size_t j = 0; j < jacobian.getN(); ++j) {
        if(j != 0) {
          out << ", ";
        }
        out << jacobian(i, j);
      }

      if(i + 1 < jacobian.getM()) {
        out << ";\n";
      } else {
        out << "]";
      }
    }

    return out;
  }
}
