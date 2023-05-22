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

#include <bitset>
#include <type_traits>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /* Original idea from http://stackoverflow.com/questions/17350214/using-enum-class-with-stdbitset */

  /**
   * @brief A bitset with enum items as flags.
   *
   * The bitset implementation allows to use an enum class as flags without causing compiler warnings or having to cast
   * the enum elements to integer types.
   *
   * The enum needs to have a 'MaxElement' item and all elements need to be positive. An example definition is:
   * \snippet examples/enumBitset.cpp Enum definition
   * where the overloaded operator | is required to create a bitset out of two flags.
   * An example use case is:
   * \snippet examples/enumBitset.cpp Enum use
   *
   * @tparam T_Enum An enumeration class. It has to have the element 'MaxElement' and only positive values.
   */
  template<typename T_Enum>
  struct EnumBitset {
    public:
      using Enum = T_Enum;  ///< See EnumBitset.

    private:
      using UnderlyingEnumType = typename std::underlying_type<Enum>::type;
      using Bitset = std::bitset<static_cast<UnderlyingEnumType>(Enum::MaxElement)>;

      static UnderlyingEnumType constexpr ALL_VALUE = (1 << static_cast<UnderlyingEnumType>(Enum::MaxElement)) - 1;

      Bitset bitset;

      CODI_INLINE UnderlyingEnumType get_value(Enum v) const {
        return static_cast<UnderlyingEnumType>(v);
      }

      CODI_INLINE EnumBitset(UnderlyingEnumType value) : bitset(value) {}

    public:

      /// Constructor all false.
      CODI_INLINE EnumBitset() : bitset() {}

      /// Constructor which sets one bit to true.
      CODI_INLINE EnumBitset(Enum pos) : bitset() {
        set(pos);
      }

      /// Test if the bit for the enum is set.
      CODI_INLINE bool test(Enum pos) const {
        return bitset.test(get_value(pos));
      }

      /// Reset the bit for the enum to false.
      CODI_INLINE EnumBitset& reset(Enum pos) {
        bitset.reset(get_value(pos));
        return *this;
      }

      /// Flip the bit for the enum.
      CODI_INLINE EnumBitset& flip(Enum pos) {
        bitset.flip(get_value(pos));
        return *this;
      }

      /// Flip all bits.
      CODI_INLINE EnumBitset& flip() {
        bitset.flip();
        bitset &= Bitset(ALL_VALUE);
        return *this;
      }

      /// Set the bit for the enum to true.
      CODI_INLINE EnumBitset& set(Enum pos) {
        bitset.set(get_value(pos));
        return *this;
      }

      /// Reset all bits to false.
      CODI_INLINE EnumBitset& reset() {
        bitset.reset();
        return *this;
      }

      /// Or operation of two bitsets.
      CODI_INLINE EnumBitset& operator|=(EnumBitset const& o) {
        bitset |= o.bitset;
        return *this;
      }

      /// Or operation of the bitset and an enum.
      CODI_INLINE EnumBitset& operator|=(Enum const& pos) {
        return set(pos);
      }

      /// And operation of two bitsets.
      CODI_INLINE EnumBitset& operator&=(EnumBitset const& o) {
        bitset &= o.bitset;
        return *this;
      }

      /// And operation of the bitsets and and enum.
      CODI_INLINE EnumBitset& operator&=(Enum const& pos) {
        return *this &= EnumBitset(pos);
      }

      /// Get the underlying bitset.
      CODI_INLINE Bitset getData() const {
        return bitset;
      }

      /// Constructor for a bitset with all values flagged as true.
      CODI_INLINE static constexpr EnumBitset ALL() {
        return EnumBitset(ALL_VALUE);
      }
  };

  /// Or operation of two bitsets.
  template<typename Enum>
  CODI_INLINE EnumBitset<Enum> operator|(EnumBitset<Enum> const& a, EnumBitset<Enum> const& b) {
    EnumBitset<Enum> r = a;
    r |= b;

    return r;
  }

  /// Or operation of the bitset and an enum.
  template<typename Enum>
  CODI_INLINE EnumBitset<Enum> operator|(EnumBitset<Enum> const& a, Enum b) {
    EnumBitset<Enum> r = a;
    r |= b;

    return r;
  }

  /// Or operation of the bitset and an enum.
  template<typename Enum>
  CODI_INLINE EnumBitset<Enum> operator|(Enum a, EnumBitset<Enum> const& b) {
    return b | a;
  }

  /// And operation of two bitsets.
  template<typename Enum>
  CODI_INLINE EnumBitset<Enum> operator&(EnumBitset<Enum> const& a, EnumBitset<Enum> const& b) {
    EnumBitset<Enum> r = a;
    r &= b;

    return r;
  }

  /// And operation of the bitset and an enum.
  template<typename Enum>
  CODI_INLINE EnumBitset<Enum> operator&(EnumBitset<Enum> const& a, Enum b) {
    EnumBitset<Enum> r = a;
    r &= b;

    return r;
  }

  /// Or operation of the bitset and an enum.
  template<typename Enum>
  CODI_INLINE EnumBitset<Enum> operator&(Enum a, EnumBitset<Enum> const& b) {
    return b & a;
  }

  /// Stream output.
  template<typename Out, typename Enum>
  CODI_INLINE Out& operator<<(Out& out, EnumBitset<Enum> const& b) {
    return out << b.getData();
  }
}
