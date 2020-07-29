#pragma once

#include <medi/ampi/types/indexTypeHelper.hpp>
#include <medi/ampi/typeDefault.hpp>

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../traits/tapeTraits.hpp"
#include "codiForwardMediPackTool.hpp"
#include "codiReverseMediPackTool.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type,
           typename _Tool =
              typename std::conditional<
                codi::isForwardTape<typename _Type::Tape>::value,
                CoDiPackForwardTool<_Type>,
                CoDiPackReverseTool<_Type>
             >::type
           >
  struct CoDiMpiTypes {
    public:

      using Type = DECLARE_DEFAULT(_Type, TEMPLATE(LhsExpressionInterface<double, double, ANY, ANY>));
      using Tool = DECLARE_DEFAULT(_Tool, medi::ADToolInterface);

      using MPIType = medi::MpiTypeDefault<Tool>;

    private:

      MPI_Datatype codiMpiType;
      MPI_Datatype modifiedMpiType;
      MPI_Datatype primalMpiType;
      MPI_Datatype adjointMpiType;

      Tool adTool;

    public:

      MPIType* MPI_TYPE;
      medi::AMPI_Datatype MPI_INT_TYPE;

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
