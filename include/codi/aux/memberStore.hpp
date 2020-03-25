#pragma once

#include <new>
#include <utility>

#include "../aux/macros.h"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename _Parent, bool _storeStatic = false>
  struct MemberStore {

      using Type = _Type; // default declaration breaks auto completion
      using Parent = DECLARE_DEFAULT(_Parent, ANY);

      static bool constexpr storeStatic = _storeStatic;

      Type member;

      template<typename ... Args>
      MemberStore(Args&& ... args) : member(std::forward<Args>(args)...) {}

      Type& get() {
        return member;
      }

      Type const& get() const {
        return member;
      }
  };

  template<typename _Type, typename _Parent>
  struct MemberStore<_Type, _Parent, true> {

      using Type = DECLARE_DEFAULT(_Type, ANY);
      using Parent = DECLARE_DEFAULT(_Parent, ANY);

      static bool constexpr storeStatic = true;

      static char member[sizeof(Type)];
      static bool isInitialized;

      template<typename ... Args>
      MemberStore(Args&& ... args) {
        if(!isInitialized) {
          isInitialized = true;
          new(member) Type(std::forward<Args>(args)...);
        }
      }

      CODI_INLINE Type& get() {
        return *((Type*)MemberStore::member);
      }

      CODI_INLINE Type const& get() const {
        return *((Type const*)MemberStore::member);
      }
  };

  template<typename Type, typename Parent>
  char MemberStore<Type, Parent, true>::member[sizeof(Type)] = {};
  template<typename Type, typename Parent>
  bool MemberStore<Type, Parent, true>::isInitialized = false;
}
