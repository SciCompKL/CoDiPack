Tutorial 4 - Vector mode AD {#Tutorial_04_Vector_mode_AD}
=======

**Goal:** Vector mode with the forward and reverse mode AD CoDiPack types.

**Prerequisite:** \ref Tutorial_01_Forward_mode_AD, \ref Tutorial_02_Reverse_mode_AD

**Function:** \ref func_simpleNto2
\snippet tutorials/Tutorial_04_Vector_mode_AD.cpp Function

**Full code:**
\snippet tutorials/Tutorial_04_Vector_mode_AD.cpp Tutorial 4 - Vector mode AD

The vector mode of AD is introduced by using matrices for the adjoint values and tangent values \f$\bar y\f$,
\f$\bar x\f$, \f$\dot x\f$ and \f$\dot y\f$. If the vector mode dimension is denoted by \f$d\f$, the matrices have the
dimension \f$\bar Y \in \R^{m \times d}\f$, \f$\bar X \in \R^{n \times d}\f$, \f$\dot X \in \R^{n \times d}\f$ and
\f$\dot Y \in \R^{m \times d}\f$. The basic idea is that derivatives for multiple directions are computed with one tape
evaluation or forward mode run.

With CoDiPack, the vector mode is enabled by exchanging the type used for the gradient computation. A fixed size vector
is implemented in the [Direction](@ref codi::Direction) class. The `codi::Real*Vec` types use this class to provide the
vector mode. Only the vector dimension needs to be provided by the user, e.g. `codi::RealReverseVec<5>`.

This is nearly the only change one has to make when switching to the vector mode in CoDiPack. In order to access the
elements of the vector mode, the Direction class offers array access operators. Gradients can also be initialized with
initializer lists that are also supported by the Direction class.

The code example above highlights the necessary changes for the vector mode with the comments *Step 1* to *Step 3*. In step 1,
the CoDiPack type is defined as a vector type. Step 2 shows the seeding for the vector directions. In step 3, the
Jacobian matrix is extracted.

#### Performance notes ####

CoDiPack offers only a fixed size vector mode that has to be defined at compile time. This decision against providing
a vector mode that can be defined at run time was made deliberately. In the forward mode, all operations have to check
if the vectors have the correct size and need to allocate new memory accordingly. This reduces the performance of a
forward vector mode. In addition, if the vector size is known a priori, the compiler can optimize the forward and
reverse vector mode even further by using e.g. SIMD instructions.

In the reverse mode, it is possible to introduce a run time decision on the vector size by using the custom evaluation
methods described in \ref Example_08_Vector_helper_interface_access and \ref Example_02_Custom_adjoint_vector_evaluation






