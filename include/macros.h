/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 *
 * Originally based on Adept 1.0 (http://www.met.rdg.ac.uk/clouds/adept/)
 * released under GPL 3.0 (Copyright (C) 2012-2013 Robin Hogan and the University of Reading).
 */

#pragma once

/**
 * @brief Create a store by value or store by reference member based on the setting of the expression template.
 *
 * The constant member storeAsReference is looked of from the type named by Name and the
 * either a reference or a value is stored.
 *
 * @param   Name  The name of the type for which a member declaration is crated.
 */
#define CODI_CREATE_STORE_TYPE(Name) \
  typename std::conditional<Name::storeAsReference, const Name &, const Name>::type


/**
 * @brief Combine two preprocessor variables into one.
 *
 * This helper routine is needed in order two expand the arguments.
 *
 * @param A First parameter that is combined
 * @param B Second parameter that is combined
 */
#define COMBINE2(A,B) A ## B

/**
 * @brief Combine two preprocessor variables into one.
 *
 * @param A First parameter that is combined
 * @param B Second parameter that is combined
 */
#define COMBINE(A,B) COMBINE2(A,B)
