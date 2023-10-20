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

#include "../../config.h"
#include "../../misc/enumBitset.hpp"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Action flags for activities during a restore. See \ref restore_actions for a detailed description.
  enum class RestoreAction {
    PrimalCreate,
    PrimalRestore,
    InputIdentifierRestore,
    OutputIdentifierRestore,
    InputGradientCreate,
    OutputGradientCreate,
    MaxElement
  };

  /// Action set for a variable during a store. See \ref restore_actions for a detailed description.
  using RestoreActions = EnumBitset<RestoreAction>;

#define ENUM RestoreAction
#include "../../misc/enumOperations.tpp"

  /// Action flags for activities during a store. See \ref store_actions for a detailed description.
  enum class StoreAction {
    PrimalCreateOnTape,
    PrimalExtract,
    InputIdentifierCreateAndStore,
    OutputIdentifierCreate,
    MaxElement
  };

  /// Action set for a variable during a store. See \ref store_actions for a detailed description.
  using StoreActions = EnumBitset<StoreAction>;

#define ENUM StoreAction
#include "../../misc/enumOperations.tpp"
}
