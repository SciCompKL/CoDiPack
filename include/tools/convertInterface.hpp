/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include <type_traits>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename I, typename I_Impl, typename T, bool isBase>
  struct InterfaceInstImpl;

  template<typename I, typename I_Impl, typename T>
  struct InterfaceInstImpl<I, I_Impl, T, true> {
      I* interface;

      InterfaceInstImpl(T* t) : interface(t) {}

      I* get() {
        return interface;
      }
  };

  template<typename I, typename I_Impl, typename T>
  struct InterfaceInstImpl<I, I_Impl, T, false> {
      I_Impl impl;

      InterfaceInstImpl(T* t) : impl(t) {}

      I* get() {
        return &impl;
      }
  };

  template<typename I, typename I_Impl, typename T>
  struct InterfaceInst {

      typedef InterfaceInstImpl<I, I_Impl, T, std::is_convertible<T*, I*>::value> Impl;

      Impl impl;

      InterfaceInst(T* t) : impl(t) {}

      I* getInterface() {
        return impl.get();
      }
  };
}
