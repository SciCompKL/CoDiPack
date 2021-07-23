#pragma once

/** \copydoc codi::Namespace */
namespace codi {

  /// Static dummy objects for e.g. default reference arguments.
  template<typename Dummy>
  struct StaticDummy {
      static Dummy dummy;  ///< Dummy object.
  };

#ifndef DOXYGEN_DISABLE
  template<typename Dummy>
  Dummy StaticDummy<Dummy>::dummy;
#endif
}
