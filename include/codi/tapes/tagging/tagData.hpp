/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../misc/enumBitset.hpp"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Properties for values.
  enum class TagFlags {
    DoNotChange,  ///< DoNotChange: Value can be assigned, but it should not change.
    DoNotUse,     ///< DoNotUse: Value should not be used. That is, it should not be used on the right hand side of an
                  ///< assignment. (Is removed after the value has been overwritten.)
    MaxElement    ///< Maximum number of elements.
  };

  /// Data for a tag.
  template<typename T_Tag>
  struct TagData {
    public:

      using Tag = CODI_DD(T_Tag, int);  ///< See TagData.

      Tag tag;                          ///< Current tag of the value.
      EnumBitset<TagFlags> properties;  ///< Current properties of the value.

      /// Constructor.
      constexpr TagData() : tag(), properties() {}

      /// Constructor.
      TagData(Tag tag) : tag(tag), properties() {}

      /// Operator for satisfying other software.
      TagData& operator+=(TagData const& o) {
        if (*this != o) {
          CODI_EXCEPTION("Operation on different tag objects.");
        }

        return *this;
      }

#if CODI_ImplicitTagConversion
      /// Allow conversion to pure tag.
      operator Tag() const {
        return tag;
      }
#endif
  };

  /// Equal comparison of two tag data objects.
  template<typename Tag>
  bool operator==(TagData<Tag> const& a, TagData<Tag> const& b) {
    return a.tag == b.tag && a.properties == b.properties;
  }

  /// Not equal comparison of two tag data objects.
  template<typename Tag>
  bool operator!=(TagData<Tag> const& a, TagData<Tag> const& b) {
    return a.tag != b.tag || a.properties != b.properties;
  }
}
