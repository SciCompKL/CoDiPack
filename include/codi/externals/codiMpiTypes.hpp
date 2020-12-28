/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include <type_traits>

#include "codiForwardMediPackTypes.hpp"
#include "codiMediPackTypes.hpp"

#include "../tapes/tapeTraits.hpp"

#include <medi/ampi/typeDefault.hpp>
#include <medi/ampi/types/indexTypeHelper.hpp>

template<typename CoDiType,
         typename ToolImpl =
            typename std::conditional<
              codi::isForwardTape<typename CoDiType::TapeType>::value,
              CoDiPackForwardTool<CoDiType>,
              CoDiPackTool<CoDiType>
           >::type
         >
struct CoDiMpiTypes {

    using MPIType = medi::MpiTypeDefault<ToolImpl>;
    using Tool = ToolImpl;

  private:

    MPI_Datatype codiMpiType;
    MPI_Datatype modifiedMpiType;
    MPI_Datatype primalMpiType;
    MPI_Datatype adjointMpiType;

    ToolImpl adTool;

  public:

    MPIType* MPI_TYPE;
    medi::AMPI_Datatype MPI_INT_TYPE;

    CoDiMpiTypes() :
      codiMpiType(createByteType(sizeof(CoDiType))),
      modifiedMpiType(codiMpiType),
      primalMpiType(createByteType(sizeof(typename CoDiType::Real))),
      adjointMpiType(primalMpiType),
      adTool(primalMpiType, adjointMpiType),
      MPI_TYPE(nullptr),
      MPI_INT_TYPE()
    {
      MPI_TYPE = new MPIType(&adTool, codiMpiType, modifiedMpiType);
      MPI_INT_TYPE = ToolImpl::OpHelper::createIntType(MPI_TYPE);
    }

    ~CoDiMpiTypes() {
      ToolImpl::OpHelper::freeIntType(MPI_INT_TYPE);

      if(nullptr != MPI_TYPE) {
        delete MPI_TYPE;
        MPI_TYPE = nullptr;
      }

      MPI_Type_free(&codiMpiType);
      MPI_Type_free(&primalMpiType);
    }

  private:

    static MPI_Datatype createByteType(size_t size) {
      MPI_Datatype type;
      MPI_Type_contiguous(size, MPI_BYTE, &type);
      MPI_Type_commit(&type);

      return type;
    }
};
