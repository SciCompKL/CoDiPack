function(setVar name doc type)
  set(OPTION_NAME CODI_${name})

  # Find default
  string(REGEX MATCH "#define ${OPTION_NAME} ([a-zA-Z0-9_]+)" _ ${configFile})
  set(OPTION_DEFAULT ${CMAKE_MATCH_1})

  if(DEFINED OPTION_DEFAULT)
    if(${type} STREQUAL "BOOL")
      # Convert default to ON and OFF.
      if(OPTION_DEFAULT)
        set(OPTION_DEFAULT "ON")
      else()
        set(OPTION_DEFAULT "OFF")
      endif()
    endif()

    # Set option
    set(${OPTION_NAME} ${OPTION_DEFAULT} CACHE ${type} ${doc})
    mark_as_advanced(${OPTION_NAME})

    # Add as compile option if not default
    if(NOT ${${OPTION_NAME}} STREQUAL ${OPTION_DEFAULT})
      if(${type} STREQUAL "BOOL")
        set(OPTION_VALUE $<BOOL:${${OPTION_NAME}}>)
      else()
        set(OPTION_VALUE ${${OPTION_NAME}})
      endif()
      target_compile_definitions(${CODIPACK_NAME} INTERFACE ${OPTION_NAME}=${OPTION_VALUE})
    endif()
  else()
    message(SEND_ERROR "Default for ${OPTION_NAME} not found.")
  endif()
endfunction()


include(${CMAKE_CURRENT_LIST_DIR}/codipack-include.cmake)

add_library(${CODIPACK_NAME} INTERFACE IMPORTED)

file(READ "${CODIPACK_INCLUDE_DIR}/codi/config.h" configFile)
setVar(ADWorkflowEvents "Enable AD workflow events, also known as Tape* events. Enabled by default." BOOL)
setVar(AvoidedInlines "Do not inline functions like evaluate()." BOOL)
setVar(AnnotateBranchLikelihood "Annotate branches with likely or unlikely, e.g., for if and else." STRING)
setVar(CheckEmptyStatements "Tapes push statements only if at least one Jacobian was pushed." BOOL)
setVar(CheckExpressionArguments "Check for invalid arguments to expressions like division by zero." BOOL)
setVar(CheckJacobianIsZero "Ignore Jacobians that are zero in Jacobian based tapes." BOOL)
setVar(CheckTapeActivity "Makes it possible to ignore certain code parts. If turned of everything will be recorded." BOOL)
setVar(CheckZeroIndex "Ignore active types that are not dependent on any input value in Jacobian tapes." BOOL)
setVar(ChunkSize "Default size of chunks (ChunkBase) used in ChunkedData in reverse tape implementations." STRING)
setVar(CopyOptimization "Do not store copy statements like a = b\; if the identity handler allows it." BOOL)
setVar(EnableAssert "Enables asserts in CoDiPack for consistency checking." BOOL)
setVar(EnableEigen "Enable Eigen specific implementations." BOOL)
setVar(EnableMPI "Add MPI and MeDiPack specific headers." BOOL)
setVar(EnableOpDiLib "Add OpDiLib specific headers. Requires codi::Config::EnableOpenMP == true." BOOL)
setVar(EnableOpenMP "Add OpenMP specific headers." BOOL)
setVar(ForcedInlines "Force inlining instead of using the heuristics from the compiler." BOOL)
setVar(IgnoreIntelNoInlineWarning "Disables warnings of the sort:  warning #2196: routine is both \"inline\" and \"noinline\"." BOOL)
setVar(IgnoreInvalidJacobians "Ignore invalid Jacobians like NaN or Inf." BOOL)
setVar(ImplicitConversion "This will give a warning every time an implicit conversion is instantiated. This" BOOL)
setVar(ImplicitConversionWarning "Warn about implicit conversions in the code." BOOL)
setVar(IndexEvents "Enable index management events. Disabled by default." BOOL)
setVar(OverflowCheck "Check in the index manager if an overflow occurred." BOOL)
setVar(PreaccEvents "Enable preaccumulation events. Disabled by default." BOOL)
setVar(RemoveDuplicateJacobianArguments "Extra pass in Jacobian tapes that combines arguments with the same identifier." BOOL)
setVar(ReversalZeroesAdjoints "With a linear index management, control if adjoints are set to zero during reversal." BOOL)
setVar(SkipZeroAdjointEvaluation "Do not perform a reverse evaluation of a statement if the seeding adjoint is zero." BOOL)
setVar(SmallChunkSize "Default smaller size of chunks (ChunkBase) used in ChunkedData in reverse tape implementations." STRING)
setVar(SortIndicesOnReset "Reuse index tapes will sort their indices on a reset." BOOL)
setVar(StatementEvents "Enable statement events. Disabled by default." BOOL)
setVar(VariableAdjointInterfaceInPrimalTapes "Allow custom adjoint vector in primal values tapes." BOOL)

set_target_properties(${CODIPACK_NAME} PROPERTIES
  INTERFACE_COMPILE_FEATURES ${CODIPACK_CXX_VERSION}
  INTERFACE_INCLUDE_DIRECTORIES "${CODIPACK_INCLUDE_DIR}"
)
