/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <new>
#include <utility>

#include "../config.h"
#include "../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Defines a member that can either be static or local to the struct.
   *
   * Initialization of the static member is done on a first touch basis.
   *
   * @tparam T_Type  The type of the member. Can be anything.
   * @tparam T_Parent  The structure where the member is located.
   * @tparam T_storeStatic  Define context of the variable.
   */
  template<typename T_Type, typename T_Parent, bool T_storeStatic = false>
  struct MemberStore {
    public:

      ///< See MemberStore
      using Type = T_Type;                         // Default declaration breaks auto completion.
      using Parent = CODI_DD(T_Parent, CODI_ANY);  ///< See MemberStore.

      static bool constexpr storeStatic = T_storeStatic;  ///< See MemberStore.

    private:

      Type member;

    public:

      /// Arguments are forwarded to the constructor of the member.
      template<typename... Args>
      MemberStore(Args&&... args) : member(std::forward<Args>(args)...) {}

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
  template<typename T_Type, typename T_Parent>
  struct MemberStore<T_Type, T_Parent, true> {
    public:

      using Type = T_Type;                         ///< See MemberStore.
      using Parent = CODI_DD(T_Parent, CODI_ANY);  ///< See MemberStore.

      static bool constexpr storeStatic = true;  ///< See MemberStore.

    private:

      static char member[sizeof(Type)];
      static bool isInitialized;

    public:

      /// \copydoc codi::MemberStore::MemberStore
      template<typename... Args>
      MemberStore(Args&&... args) {
        if (!isInitialized) {
          isInitialized = true;
          new (member) Type(std::forward<Args>(args)...);
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
