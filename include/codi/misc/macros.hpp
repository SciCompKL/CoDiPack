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

#include "../config.h"

/** @file */

/** \copydoc codi::Namespace */
namespace codi {

  /// Disable unused warnings for an arbitrary number of arguments.
  template<typename... Args>
  void CODI_UNUSED(Args const&...) {}

/// Used in a constexpr context, where using CODI_UNUSED spoils the constexpr.
#define CODI_UNUSED_ARG(arg) /* arg */

/// Check the condition only if the option is true, otherwise the result is always true.
/// Used like if(CODI_ENABLE_CHECK(option, condition)) {...}, option == false means that the if body is always executed.
#define CODI_ENABLE_CHECK(option, condition) (!(option) || (condition))

/// Conversion macro.
#define CODI_TO_STRING2(expression) #expression

/// Conversion macro.
#define CODI_TO_STRING(expression) CODI_TO_STRING2(expression)

/// Check for CPP 14 standard.
#define CODI_IS_CPP14 (201402L <= __cplusplus)

/// Check for CPP 17 standard.
#define CODI_IS_CPP17 (201703L <= __cplusplus)

/*******************************************************************************/
/** @name Default template type declarations
 *  @anchor TemplateDeclarationHelpers
 *
 * These templates are used to employ the design guideline for the default definitions of template
 * arguments. See \ref TemplateDeclaration for Details.
 *
 * @{
 */

/**
 * CODI_IDE can be defined to use the default declaration of type names. This enables auto completion in the IDEs.
 *
 * Every using declaration in all CoDiPack classes should declare its variables as:
 *  using TYPE = CODI_DECLARE_DEFAULT(T_TYPE, Default);
 */
#if CODI_IDE
  #define CODI_DECLARE_DEFAULT(Type, Default) Default
#else
  #define CODI_DECLARE_DEFAULT(Type, Default) Type
#endif

/// Abbreviation for CODI_DECLARE_DEFAULT
// Does not work in QT if CODI_DECLARE_DEFAULT is used
#if CODI_IDE
  #define CODI_DD(Type, Default) Default
#else
  #define CODI_DD(Type, Default) Type
#endif

/// Used in default declarations of expression templates.
#define CODI_ANY int

#ifndef DOXYGEN_DISABLE
  /// Placeholer for the implementation of an interface.
  struct ImplProxy {};
#endif
/// Used in interface declarations to indicate the type of the implementing class.
#define CODI_IMPLEMENTATION ImplProxy

/// Expand template types in preprocessor macros.
#define CODI_TEMPLATE(...) __VA_ARGS__

/// Abbreviation for CODI_TEMPLATE
#define CODI_T(...) CODI_TEMPLATE(__VA_ARGS__)

/// Used in interface declarations for types that have to be defined in the specializations.
#define CODI_UNDEFINED void

/// Used in interface declarations for variables that have to be defined in the specializations.
#define CODI_UNDEFINED_VALUE false

#if CODI_IDE
  #define CODI_STATIC_ASSERT(cond, message) /* Do not check in IDE mode */
#else
  #define CODI_STATIC_ASSERT(cond, message) static_assert(cond, message)
#endif

#if CODI_IDE
  /// Proxy definition for an ActiveType in the real traits.
  struct ActiveTypeProxy {
      using Real = double;
      using Identifier = int;
  };

  /// Declaration of the default full tape interface.
  #define CODI_DEFAULT_TAPE FullTapeInterface<double, double, int, EmptyPosition>

  /// Declaration of a default parallel tape interface.
  #define CODI_DEFAULT_PARALLEL_TAPE CODI_UNION<CODI_DEFAULT_TAPE, EditingTapeInterface<EmptyPosition>>

  /// Declaration of the default lhs expression interface.
  #define CODI_DEFAULT_LHS_EXPRESSION LhsExpressionInterface<double, double, CODI_DEFAULT_TAPE, CODI_ANY>
#endif

#ifndef DOXYGEN_DISABLE
  /// Creates a union of interface definitions.
  template<typename First, typename... Tail>
  struct CODI_UNION : public First, public CODI_UNION<Tail...> {};
#endif

  /// Creates a union of interface definitions.
  template<typename First>
  struct CODI_UNION<First> : public First {};

/// @}

/// Wrap a function in a function object. Used for performance optimizations.
#define CODI_WRAP_FUNCTION(NAME, FUNC)        \
  struct NAME {                               \
    public:                                   \
      /** Empty */                            \
      template<typename... Args>              \
      void operator()(Args&&... args) const { \
        FUNC(std::forward<Args>(args)...);    \
      }                                       \
  }

/// Wrap a function in a function object. Used for performance optimizations.
#define CODI_WRAP_FUNCTION_TEMPLATE(NAME, FUNC)   \
  template<typename... TT>                        \
  struct NAME {                                   \
    public:                                       \
      /** Empty */                                \
      template<typename... Args>                  \
      void operator()(Args&&... args) const {     \
        FUNC<TT...>(std::forward<Args>(args)...); \
      }                                           \
  }
}
