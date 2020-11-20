#pragma once

#include <new>
#include <utility>

#include "../aux/macros.hpp"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Defines a member, that can either be static or local to the struct.
   *
   * Initialization of the static member is done on a first touch basis.
   *
   * @tparam _Type  The type of the member. Can be anything.
   * @tparam _Parent  The structure where the member is located.
   * @tparam _storeStatic  Define context of the variable.
   */
  template<typename _Type, typename _Parent, bool _storeStatic = false>
  struct MemberStore {
    public:

      ///< See MemberStore
      using Type = _Type; // default declaration breaks auto completion
      using Parent = CODI_DECLARE_DEFAULT(_Parent, CODI_ANY); ///< See MemberStore

      static bool constexpr storeStatic = _storeStatic;  ///< See MemberStore

    private:

      Type member;

    public:

      /// Arguments are forwarded to the constructor of the member
      template<typename ... Args>
      MemberStore(Args&& ... args) : member(std::forward<Args>(args)...) {}

      /// Get a reference to the actual member.
      Type& get() {
        return member;
      }

      /// Get a reference to the actual member.
      Type const& get() const {
        return member;
      }
  };

  /// \copydoc codi::MemberStore
  template<typename _Type, typename _Parent>
  struct MemberStore<_Type, _Parent, true> {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_ANY); ///< See MemberStore
      using Parent = CODI_DECLARE_DEFAULT(_Parent, CODI_ANY);  ///< See MemberStore

      static bool constexpr storeStatic = true;  ///< See MemberStore

    private:

      static char member[sizeof(Type)];
      static bool isInitialized;

    public:

      /// \copydoc codi::MemberStore::MemberStore
      template<typename ... Args>
      MemberStore(Args&& ... args) {
        if(!isInitialized) {
          isInitialized = true;
          new(member) Type(std::forward<Args>(args)...);
        }
      }

      /// \copydoc codi::MemberStore::get()
      CODI_INLINE Type& get() {
        return *((Type*)MemberStore::member);
      }

      /// \copydoc codi::MemberStore::get() const
      CODI_INLINE Type const& get() const {
        return *((Type const*)MemberStore::member);
      }
  };

#ifndef DOXYGEN_DISABLE

  template<typename Type, typename Parent>
  char MemberStore<Type, Parent, true>::member[sizeof(Type)] = {};

  template<typename Type, typename Parent>
  bool MemberStore<Type, Parent, true>::isInitialized = false;
#endif
}
