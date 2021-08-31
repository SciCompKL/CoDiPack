Active type definitions {#ActiveTypeList}
=======

CoDiPack defines several different tapes and options for these tapes, which results in a large range of possible tape
configurations. Each tape is bound to an active type which has to be used in the application for the differentiation.
The most used types are:
 - codi::RealForward : AD forward mode implementation, see \ref sec_forwardAD.
 - codi::RealReverse : AD reverse mode implementation, see \ref sec_reverseAD.
 
These two types are configured such that they work mostly out of the box and the user does not need to worry that much
about modifying his application for the type.

The other codi::RealReverse* type definitions have some common postfixes. Each of these postfixes changes some aspects
of the tape implementation. If they are not used, the default is used instead.

The most common reverse AD types are:
 - codi::RealReverse : Jacobian taping with a linear index management scheme. Data is allocted automatically.
 - codi::RealReverseIndex : Jacobian taping with a reuse index management scheme. Data is allocted automatically.
 - codi::RealReversePrimal : Primal value taping with a linear index management scheme. Data is allocted automatically.
 - codi::RealReversePrimalIndex : Primal value taping with a reuse index management scheme. Data is allocted automatically.

The current postfixes are:
 - Index: Uses a resue index management scheme with copy optimizations. See \ref IndexManagers for details.
   - Absence: Uses a linear index management scheme with copy optimizations.
 - Primal: Uses a primal value taping approach. See \ref TapingStrategy for details.
   - Absence: Uses a Jacobian taping approach.
 - Unchecked: Memory allocation has to be done beforehand. CoDiPack does not check in theses tapes if memory for the
              next statement is available.
   - Absence: Memory is allocated automatically in a chunked fashion.
 - Vec: Fixed vector mode evaluations.
   - Absence: Scalar reverse and forward evaluations.
 - Gen: Generalized template definitions that allow to modify other aspects of the type
   - Real: The primal computation type, default: double.
   - Gradient: The computation type for the gradients, default: Real.
   - IndexManager: The manger for the identifiers. See \ref IndexManagers.
   - Index: The identifier type for the linear index managers.
   - StatementEvaluator: How statements are stored for primal value types. See \ref StatementEvaluators.
 
Index managers {#IndexManagers}
-------

Index managers handle the distribution of the identifiers used by CoDiPack. Each identifier is coupled to a variable and
depending on the manager, different strategies are used for the creation, deletion and reuse of identifiers.
There are currently three index managers available in CoDiPack. 

codi::LinearIndexManager distributes the identifiers in an increasing fashion. Each identifier is used only once during
the tape recording. A reset on the tape will also reset the identifiers and all values need to be registered again. This
index manager is compatible with C-like memory operations (e.g. memcpy). Tapes which use this index manager do not
need to store statements which perform a copy operation. It is the default index manager in CoDiPack.

codi::ReuseIndexManager allows identifiers to be reused. The lifetime of an identifier is determined by the lifetime of
the associated variable (value). If the value is overwritten or the variable is freed, then the lifetime of the
identifier ends and it can be reused. In case of an overwrite a new identifier is generated. This index manager is
currently not used directly in CoDiPack. The application cannot use C-like memory operations. Copy operations need to
be stored by the tape.

codi::MultiUseIndexManager is an extension of the codi::ReuseIndexManager. It allows that copy operations do not need to
be stored by the tape. The application still cannot use C-like memory operations. It is the default index manager for
types that use the 'Index' postfix.

Statement evaluators {#StatementEvaluators}
-------

The statement evaluators define how the expressions are stored in primal value tapes. Depending on the implementation
the support for different forward and reverse evaluations might vary. The number of jumps the CPU has to do might
also be affected by the implementation.

codi::ReverseStatementEvaluator stores the function pointer for the reverse evaluation of the expression in the tape.
Therefore it supports only reverse tape evaluations.

codi::DirectStatementEvaluator creates static handles in the binaries and stores the pointers to theses handles in the
tape. This requires one additional address lookup by the CPU. It supports all evaluation modes of the tape.

codi::InnerStatementEvaluator uses the same strategy as the codi::DirectStatementEvaluator for the handle creation, but
it shifts the boundary between the tape evaluation and the statement evaluation towards the statement evaluation. This
allows the compiler to optimize also for the general setup of the statement evaluation (e.g. copying passive values).
