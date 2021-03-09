#pragma once

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Dummy>
  struct StaticDummy {
      static Dummy dummy;
  };

  template<typename Dummy>
  Dummy StaticDummy<Dummy>::dummy;
}
