#pragma once

#include <medi/medi.hpp>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/tapeTraits.hpp"
#include "codiForwardMediPackTool.hpp"
#include "codiReverseMediPackTool.hpp"

/** \copydoc codi::Namespace */
namespace codi {


  /**
   * @brief Mpi datatype implementation for CoDipack types in the type wrapper of MeDiPack.
   *
   * See \ref Example_13_MPI_communication for an example.
   *
   * Use the member MPI_TYPE as the type for the communication in MeDiPack wrapped MPI routines or
   * MPI_INT_TYPE for pairs of CoDiPack and an int.
   *
   * @tparam _Type  CoDiPack active type. Needs to implement LhsExpressionInterface.
   * @tparam _Tool  Actual tool implementation of the MeDiPack interface.
   */
  template<typename _Type,
           typename _Tool =
              typename std::conditional<
                codi::TapeTraits::IsForwardTape<typename _Type::Tape>::value,
                CoDiPackForwardTool<_Type>,
                CoDiPackReverseTool<_Type>
             >::type
           >
  struct CoDiMpiTypes {
    public:

      /// See CoDiMpiTypes
      using Type = CODI_DD(_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Tool = CODI_DD(_Tool, medi::ADToolInterface); ///< See CoDiMpiTypes

      using MPIType = medi::MpiTypeDefault<Tool>; ///< MeDiPack default implementation

    private:

      MPI_Datatype codiMpiType;
      MPI_Datatype modifiedMpiType;
      MPI_Datatype primalMpiType;
      MPI_Datatype adjointMpiType;

      Tool adTool;

    public:

      MPIType* MPI_TYPE; ///< MPI_Datatype for the specified CoDiPack type.
      medi::AMPI_Datatype MPI_INT_TYPE;  ///< MPI_Datatype for the specified CoDiPack type and an int.

      /// Constructor
      CoDiMpiTypes() :
        codiMpiType(createByteType(sizeof(Type))),
        modifiedMpiType(codiMpiType),
        primalMpiType(createByteType(sizeof(typename Type::Real))),
        adjointMpiType(primalMpiType),
        adTool(primalMpiType, adjointMpiType),
        MPI_TYPE(nullptr),
        MPI_INT_TYPE()
      {
        MPI_TYPE = new MPIType(&adTool, codiMpiType, modifiedMpiType);
        MPI_INT_TYPE = Tool::OpHelper::createIntType(MPI_TYPE);
      }

      /// Destructor
      ~CoDiMpiTypes() {
        Tool::OpHelper::freeIntType(MPI_INT_TYPE);

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
}
