#pragma once

/** \copydoc codi::Namespace */
namespace codi {

  template <typename Dummy>
  struct StaticDummy {
    public:
      static Dummy dummy; /**< Static value for the dummy */
  };

  template<typename Dummy>
  Dummy StaticDummy<Dummy>::dummy;
}
