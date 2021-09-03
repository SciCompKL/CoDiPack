Migration from 1.9 to 2.0 {#Migration_1_9_2_0}
=======

### General usage

In general we tried to keep the main user interface intact. For a simple _Hello World_ case only the access to the tape
has changed. 

CoDiPack 1.9:
```{.cpp}
  using Real = codi::RealReverse;
  using Tape = typename Real::TapeType;
  Tape& tape = Real::getGlobalTape();
```

CoDiPack 2.0:
```{.cpp}
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;  // TapeType -> Tape
  Tape& tape = Real::getTape();      // getGlobalTape() -> getTape()
```

Overall, the rewrite of CoDiPack was done to improve the maintainability of the existing functionality and to support
future features. For all internal classes, it was checked if the naming of the class properly hints at the use case the
class wants to solve. The same has been done for function names and type declarations. Therefore, nearly everything in
the internal structure of CoDiPack has changed and a complete list would be hard to assemble. The most notable changes
are:
 - Type declarations inside of classes have dropped the `Type` suffix, especially true for type definitions of template
   arguments in classes.
 - `getGradientData`, `GradientData` have been renamed to `getIdentifier`, `Identifier`. This change reflects the
   general meaning of what the gradient data always has been, the identifier of the active types for the reverse tapes.
   This naming makes now less sence for the forward mode, but in general it should improve the understandability.

Incomplete list of major user interface changes.
| CoDiPack 1.9 | CoDiPack 2.0 |
|:--------|:--------|
| ActiveReal | ActiveType |
| getGradientData | getIdentifier |
| TapeType | Tape |
| GradientData | Identifier |
| GradientValue | Gradient |
| getGlobalTape | getTape |

### Mpi communication

The MeDiPack implementation moved now to the final destination in `codi/tools/mpi/codiMpiTypes.hpp`. The structure
`CoDiMpiTypes` moved into the `codi` namespace. \ref Example_13_MPI_communication provides a demonstration of the full
use case.

The correct definition is now:
```{.cpp}
#include "codi/tools/mpi/codiMpiTypes.hpp"
using MpiTypes = codi::CoDiMpiTypes<codi::RealReverse>;
```

### External functions

The interface of the codi::ExternalFunctionHelper is nearly the same, only the passive evaluation function got renamed.
The manual push of external functions changed a lot, please see \ref Example_11_External_function_user_data for an
overview of the new interface.

| CoDiPack 1.9 | CoDiPack 2.0 |
|:--------|:--------|
| callPassiveFunc | callPrimalFuncWithADType |

### Configuration options

Configuration options are renamed to have better indicate their meaning. In general, the `Opt`, `Enable`, `Disable` prefixes
are dropped. See below a full list of the changes.

| CoDiPack 1.9 | CoDiPack 2.0 |
|:--------|:--------|
| CODI_UseForcedInlines | CODI_ForcedInlines |
| CODI_UseAvoidedInlines | CODI_AvoidedInlines |
| CODI_OptIgnoreInvalidJacobians | CODI_IgnoreInvalidJacobians |
| CODI_OptJacobiIsZero | CODI_CheckJacobianIsZero |
| CODI_OptCheckZeroIndex | CODI_CheckZeroIndex |
| CODI_OptCheckEmptyStatements | CODI_CheckEmptyStatements |
| CODI_OptTapeActivity | CODI_CheckTapeActivity |
| CODI_ZeroAdjointReverse | CODI_ReversalZeroesAdjoints |
| CODI_OptZeroAdjoint | CODI_SkipZeroAdjointEvaluation |
| CODI_DisableAssignOptimization | CODI_CopyOptimization |
| CODI_EnableImplicitConversion | CODI_ImplicitConversion |
| CODI_DisableImplicitConversionWarning | CODI_ImplicitConversionWarning=false |
| CODI_DisableSortIndicesOnReset | CODI_SortIndicesOnReset |
| CODI_EnableVariableAdjointInterfaceInPrimalTapes | CODI_VariableAdjointInterfaceInPrimalTapes |
| CODI_EnableCombineJacobianArguments | CODI_RemoveDuplicateJacobianArguments |
| CODI_DisableCalcGradientSpecialization | removed |
| CODI_AdjointHandle_Jacobi | future release |
| CODI_AdjointHandle_Jacobi_Reverse | future release |
| CODI_AdjointHandle_Primal | future release |
| CODI_AdjointHandle_Tangent | future release |
| CODI_IndexHandle | future release |
| DisableIntelNoInlineWarning | CODI_IgnoreIntelNoInlineWarning |


