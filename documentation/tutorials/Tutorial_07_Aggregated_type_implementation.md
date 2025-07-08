Tutorial 7 - Aggregated type implementation {#Tutorial_07_Aggregated_type_implementation}
============

**Goal:** Learn how to add your own aggregated types to CoDiPack.

**Prerequisite:** \ref Tutorial_02_Reverse_mode_AD

**Function:** Norm of two scaled vector.
\snippet tutorials/Tutorial_07_Aggregated_type_implementation.cpp Function

**Full code:**
\snippet tutorials/Tutorial_07_Aggregated_type_implementation.cpp Tutorial 7 - Aggregated type implementation

#### Introduction ####

An aggregated type, in the sense of CoDiPack, is a structure that can be expressed by a set of `double` values. E.g.
`std::complex<double>` can be represented by two double values. If aggregated types are added to the CoDiPack expression
tree, then memory and runtime of differentiated programs can be reduced.


#### Defining the basic types and the active vector type ####

In this tutorial we want to add specialized expressions for the Eigen `Vector4` type. Therefore, we decalare first the 
basic types as:
\snippet tutorials/Tutorial_07_Aggregated_type_implementation.cpp Base declaration

It is now possible to declare the vector type that integrates into the CoDiPack expressions. The simplest way is to
extend from #codi::AggregatedActiveType. This base class implements all the necessary interfaces and functions for the
compatibility with the CoDiPack expressions. In order to do this, #codi::AggregatedActiveType needs some information
about the used type, which is in our case `Vector4`. This is done by the specialization of 
#codi::RealTraits::AggregatedTypeTraits. Here, the array constructor, the array access methods, and the number of
elements for the type are defined. If the type can be described as an array of floating point values, then the default
implementation #codi::RealTraits::ArrayAggregatedTypeTraitsBase can be used. In addition we need to specify a transposed
operation for the type. All together this is done with:
\snippet tutorials/Tutorial_07_Aggregated_type_implementation.cpp Vector4 specializations
\snippet tutorials/Tutorial_07_Aggregated_type_implementation.cpp ActiveVector4 definition
Currently we only define a constructor and the array access for the `ActiveVector4`.

#### Definition of operations ####
For our function example, we need two operations the scalar multiplication and the vector addition. Both can be defined
by implementing the codi::BinaryOperation interface. This interface defines the primal operation and the gradient
computations. For the real valued case, the operation implementations can always assume that the Jacobian has the same
type as the primal result. For the vector operations, this is no longer the case. Therfore, we also need to declare the
Jacobian types in the operation implementations:
\snippet tutorials/Tutorial_07_Aggregated_type_implementation.cpp Operation implementation

#### Definition of member operations ####
If the base type contains member operations, then these need to be specified to. These operations could be specified in
the implementation of ActiveVector4, but then they would not be availalbe on expressions. For example, it would
not be possible to write `(v1 + v2).norm()`, since `(v1 + v2)` has the type #codi::ComputeExpression. It is therefore
necessary to inject the member operations into all expressions. This can be done by a specialization of 
#codi::ExpressionMemberOperations for `Vector4`:
\snippet tutorials/Tutorial_07_Aggregated_type_implementation.cpp Member operation implementation

#### Advantages and further hints ####
With the finished implementation we can now call the function defined above. We will do this once with an unspecialized
type (`using Vector4WithActiveType = Eigen::Matrix<ActiveNumber, 4, 1>`) and once with the ActiveVector4. In the test we
measure the tape memory that is recorded during the function evaluation. The result is:
~~~
Running example with 'Vector4WithActiveType' vector type. No specialization are used for Eigen vector.
-------------------------------------
CoDi Tape Statistics ( JacobianLinearTape )
-------------------------------------
  Total memory used      :     648.00 B
-------------------------------------
Statement entries
-------------------------------------
  Total number           :         20
-------------------------------------
Jacobian entries
-------------------------------------
  Total number           :         39
-------------------------------------

Running example with 'ActiveVector4' vector type. The specializations are used for Eigen vector.
-------------------------------------
CoDi Tape Statistics ( JacobianLinearTape )
-------------------------------------
  Total memory used      :     129.00 B
-------------------------------------
-------------------------------------
Statement entries
-------------------------------------
  Total number           :          1
-------------------------------------
Jacobian entries
-------------------------------------
  Total number           :         10
-------------------------------------
~~~
With the specialization only 129 bytes are used which is a reduction of about 80%. Since we integrated the ActiveVector4
into the CoDiPack expressions, the tape no longer records all intermediate steps in the function evaluation. 
The statement in the function is seen by CoDiPack as one big expression which reduces the number of statment and
Jacobian entries.

Currently, all uses of Vector4 in an application would need to be exchanged with ActiveVector4. A declaration like
`Eigen::Matrix<ActiveNumber, 4, 1>` would not automatically use the specialization. An automatic use can be achieved by
a direct specialization of `Eigen::Matrix` for a CoDiPack type. E.g:
~~~{.cpp}
template<>
struct Eigen::Matrix<ActiveNumber, 4, 1, Eigen::AutoAlign, 4, 1> :
    public codi::AggregatedActiveType<Vector4, ActiveNumber, ActiveVector4> {
    
  ... // Same as ActiveVector4.
};
~~~
With this specialization both runs of the test function have the same reduced memory consumption.

#### Notes on implementation ####

The current implementation of aggregated types is tailored to the implementation of `std::complex`. The current
expression framework of CoDiPack has problems in handling operations where the operators have not the commutativity
property. It is also problematic to handle operations where the AD reverse mode can not be expressed as a Jacobian times
a bar value. These cases can be handled, but their implementation is quite cumbersome. We will extend the expression
framework in the near future such that these cases can be handled in an improved way.
