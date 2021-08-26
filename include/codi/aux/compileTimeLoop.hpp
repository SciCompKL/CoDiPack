#pragma once

#include <utility>

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Compile time loop evaluation.
   *
   * pos is counted backwards until zero excluding zero.
   *
   * Called range is: (0,pos]
   *
   * @tparam T_pos  Starting value for the loop. Counted downwards.
   */
  template<size_t T_pos>
  struct CompileTimeLoop {
    public:

      static size_t constexpr pos = T_pos;  ///< See CompileTimeLoop.

      /// Func is evaluated with args as func(pos, args...)
      template<typename Func, typename... Args>
      static CODI_INLINE void eval(Func&& func, Args&&... args) {
        func(std::integral_constant<size_t, pos>{}, std::forward<Args>(args)...);

        CompileTimeLoop<pos - 1>::eval(std::forward<Func>(func), std::forward<Args>(args)...);
      }
  };

  /// Termination of loop evaluation.
  template<>
  struct CompileTimeLoop<0> {
    public:

      static size_t constexpr pos = 0;  ///< See CompileTimeLoop.

      /// Nothing is evaluated.
      template<typename... Args>
      static CODI_INLINE void eval(Args&&... args) {
        CODI_UNUSED(args...);
      }
  };
}
