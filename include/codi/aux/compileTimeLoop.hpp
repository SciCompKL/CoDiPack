#pragma once

#include <utility>

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<size_t _pos>
  struct CompileTimeLoop {
    public:

      static size_t constexpr pos = _pos;

      template<typename Func, typename ... Args>
      static void eval(Func&& func, Args&& ... args) {
        func(std::integral_constant<size_t, pos>{}, std::forward<Args>(args)...);

        CompileTimeLoop<pos - 1>::eval(std::forward<Func>(func), std::forward<Args>(args)...);
      }
  };

  template<>
  struct CompileTimeLoop<0> {
    public:

      static size_t constexpr pos = 0;

      template<typename ... Args>
      static void eval(Args&& ... args) {
        CODI_UNUSED(args...);
      }
  };
}
