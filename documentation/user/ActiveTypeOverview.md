Active type overview {#ActiveTypeOverview}
=======

This is a short selection overview about which CoDiPack type you should use. Details can be found at \ref ActiveTypeList.
A plus in the table means that this is supported or a beneficial type, a minus that it is not supported or not beneficial.

| Type | C-like memory operations | Vector mode | Large inactive parts | Hessian computation |
|:-----|:-----:|:-----:|:-----:|:-----:|
| codi::RealReverse | + | - | + | - |
| codi::RealReverseIndex | - | + | + | - |
| codi::RealReversePrimal | + | - | - | + |
| codi::RealReversePrimalIndex | - | + | - | + |

__Legend:__
 - __C-like memory operations:__ Operations like `memcpy` can be used on the CoDiPack type. Assuming that the correct
                                 size is used (e.g. `memcpy(y, x, n * sizeof(codi::RealReverse))`.
 - __Vector mode:__ Additional memory consumpition for the vector mode is quite small. Vector modes can either be run
                    with the *Vec [variant](@ref ActiveTypeList) of the type (@ref Tutorial_04_Vector_mode_AD) or with
                    the codi::CustomAdjointVectorHelper (@ref Example_02_Custom_adjoint_vector_evaluation).
 - __Large inactive parts:__ CoDiPack performs an [online activity analysis](@\ref ActivityAnalysis). If variables do not depend on the input
                             variables - registered with calls to `registerInput` (independent variables) - then the
                             sensitivites are not propagated for these _passive_ variables. Jacobian taping approaches
                             pefrom usually better in such a case.
 - __Hessian computation:__ For a full Hessian computation the seeding for the second order sensitivites need to be
                            changed. Primal value tapes can simply reevaluate these seedings which is usually faster
                            than rerunning the program. See TODO for an automaic Hessian computation.
