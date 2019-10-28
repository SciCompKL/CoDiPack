Tutorials {#tutorialPage}
=======

The CoDiPack tutorials are structured in three parts. The first part are the beginner (B) tutorials where we use some
helper tools in CoDiPack to compute the derivative, Jacobians and Hessians of some functions.

The second section of the tutorials are the regular CoDiPack tutorials. Here, all details of the usage of CoDiPack are
explained step step and the user has the full control about how the derivatives are computed.

The third section are the advanced (A) tutorials which introduce some additional helpers, which can be used to reduce
the required memory and the runtime of the tape evaluations.

The beginner CoDiPack tutorials:
  - \subpage TutorialB1 "Tutorial B1": Evaluation helper introduction: Compute Jacbioans and Hessians of functions
   - \subpage TutorialB1_1 "Tutorial B1.1": Function objects for the the evaluation helper: How more generalized function objects can be defined.
   - \subpage TutorialB1_2 "Tutorial B1.2": Run time optimization for the evaluation helper: How to evaluated a function object at several different points.
  - \subpage TutorialB2 "Tutorial B2": Tape helper introduction
    - Simple tape management and algorithms for computation of gradients, Jacobians and Hessians

The CoDiPack tutorials:
  - \subpage Tutorial1 "Tutorial 1": Explains the basic usage for the RealForward type.
  - \subpage Tutorial2 "Tutorial 2": Explains the basic usage for the RealReverse type.
  - \subpage Tutorial3 "Tutorial 3": Compute the full Jacobi matrix of a function with the RealForward type.
  - \subpage Tutorial4 "Tutorial 4": Compute the full Jacobi matrix of a function with the RealReverse type.
  - \subpage Tutorial5 "Tutorial 5": ReferenceActiveReal demonstration
  - \subpage Tutorial6 "Tutorial 6": Vector mode demonstration
  - \subpage Tutorial7 "Tutorial 7": Higher order derivatives
   - \subpage Tutorial7_1 "Tutorial 7.1": Higher order derivatives with the template interface
   - \subpage Tutorial7_2 "Tutorial 7.2": Higher order derivatives without a helper interface
  - \subpage Tutorial8 "Tutorial 8": Several different gradient computations / reset of reverse tapes
  - \subpage Tutorial9 "Tutorial 9": Gradient evaluations for different program configurations

The advanced CoDiPack tutorials:
  - \subpage TutorialA1 "Tutorial A1": External functions
    - Can be used for handling libraries or to optimize the derivative computation.
  - \subpage TutorialA2 "Tutorial A2": Preaccumulation
    - Reduce tape memory in a non intrusive way.
  - \subpage TutorialA3 "Tutorial A3": Manual storing of statements
    - Apply very local and small tape optimizations
    - Handle small functions where the derivative is known
  - \subpage TutorialA4 "Tutorial A4": Advanced vector mode usage
    - \subpage TutorialA4_1 "Tutorial A4.1": OpenMP reverse mode evaluation
  - \subpage TutorialA5 "Tutorial A5": Custom Jacobian and Hessian storage
