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

#include <medi/medi.hpp>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../traits/tapeTraits.hpp"
#include "codiForwardMeDiPackTool.hpp"
#include "codiReverseMeDiPackTool.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief MPI datatype implementation for CoDipack types in the type wrapper of MeDiPack.
   *
   * See \ref Example_13_MPI_communication for an example.
   *
   * Use the member MPI_TYPE as the type for the communication in MeDiPack wrapped MPI routines or
   * MPI_INT_TYPE for pairs of CoDiPack and an int.
   *
   * @tparam T_Type  CoDiPack active type. Must implement LhsExpressionInterface.
   * @tparam T_Tool  Actual tool implementation of the MeDiPack interface.
   */
  template<typename T_Type,
           typename T_Tool = typename std::conditional<codi::TapeTraits::IsForwardTape<typename T_Type::Tape>::value,
                                                       CoDiPackForwardTool<T_Type>, CoDiPackReverseTool<T_Type> >::type>
  struct CoDiMpiTypes {
    public:

      /// See CoDiMpiTypes.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using Tool = CODI_DD(T_Tool, medi::ADToolInterface);  ///< See CoDiMpiTypes.

      using MPIType = medi::MpiTypeDefault<Tool>;  ///< MeDiPack default implementation.

    private:

      MPI_Datatype codiMpiType;
      MPI_Datatype modifiedMpiType;
      MPI_Datatype primalMpiType;
      MPI_Datatype adjointMpiType;

      Tool adTool;

    public:

      MPIType* MPI_TYPE;                 ///< MPI_Datatype for the specified CoDiPack type.
      medi::AMPI_Datatype MPI_INT_TYPE;  ///< MPI_Datatype for the specified CoDiPack type and an int.

      /// Constructor
      CoDiMpiTypes()
          : codiMpiType(createByteType(sizeof(Type))),
            modifiedMpiType(codiMpiType),
            primalMpiType(createByteType(sizeof(typename Type::Real))),
            adjointMpiType(primalMpiType),
            adTool(primalMpiType, adjointMpiType),
            MPI_TYPE(nullptr),
            MPI_INT_TYPE() {
        MPI_TYPE = new MPIType(&adTool, codiMpiType, modifiedMpiType);
        MPI_INT_TYPE = Tool::OpHelper::createIntType(MPI_TYPE);
      }

      /// Destructor
      ~CoDiMpiTypes() {
        Tool::OpHelper::freeIntType(MPI_INT_TYPE);

        if (nullptr != MPI_TYPE) {
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
}
