Expressions {#Expressions}
=======

The expressions in CoDiPack help to perform operations on the statement level instead of the operator level.
In an operator overloading framework, the general design of the overloads is like
\f[
  \text{Real} \circ \text{Real} \rightarrow \text{Real},
\f]
where \f$\text{Real}\f$ is the computation type and \f$\circ\f$ the operator.
Expressions change this concept such that the return type of the operator is no longer the result itself but a structure
that contains information about the operator. An example with two different types \f$A\f$ and \f$B\f$ and an expression
context is
\f[
  \text{Expr}_A \circ \text{Expr}_B \rightarrow \text{Expr}_{A \circ B}.
\f]
The final return value is now an expression which contains information about the operation, which is indicated by the
suffix \f$A \circ B\f$.

This allows us to move from the differentiation of single elemental operations to the differentiation of full statements with right hand sides that are composed of multiple elemental operations.

On the equation
\f[
  w = ((a + b) * (c - d))^2
\f]
the pure operator overloading approach would create three intermediate statements
\f{align*}{
  t_1 = & a + b \\
  t_2 = & c - d \\
  t_3 = & t_1 * t_2 \\
  w = & t_3^2 \eqdot
\f}
The expression template approach creates one large expression that is represented by the structure
\f[
  \text{SQUARE} < \text{MULT} < \text{ADD} < \text{Real}, \text{Real}>, \text{SUB}<\text{Real}, \text{Real}> > > \eqdot
\f]
It can be used to evaluate the result \f$w\f$ and access additional information about the operations. In CoDiPack
the functions for each expression are kept rather small. The only directly callable method is `getValue` which evaluates
the expression. For all other operations the traversal logic classes have to be used.

The current most used expression implementation are:
 - [ActiveType](@ref codi::ActiveType): LhsExpression implementation for a concrete value.
 - [BinaryExpression](@ref codi::BinaryExpression): Expression with two arguments.
 - [ConstantExpression](@ref codi::ConstantExpression): Represents a constant value in the expression tree like 4.0 or
                                                        other non CoDiPack types.
 - [UnaryExpression](@ref codi::UnaryExpression): Expression with one argument.
 - [LhsExpressionInterface](@ref codi::LhsExpressionInterface): General interface for expressions that are lvalues.

Expression traversal (custom operations on expressions) {#customExpressionLogic}
-------

Custom operations on expressions in CoDiPack are all implemented with the [TraversalLogic](@ref codi::TraversalLogic)
or [CompileTimeTraversalLogic](@ref codi::CompileTimeTraversalLogic) classes. The first one enables custom logic in a
runtime context. Variables and results can be stored in the implementation. The CompileTimeTraversalLogic allows for the
computation of results in a compile time context. Both classes contain functions for the visit of links, nodes and
leaves in the expression graph. The terms are explained in the following picture:
```
  ┌─┐   ┌─┐   ┌─┐   ┌─┐
  │a│   │b│   │c│   │d│  # Nodes without arguments are leaf nodes, they are primary objects like floating point
  └┬┘   └┬┘   └┬┘   └┬┘  # constants or ActiveType values.
   │     │     │     │
   └──┬──┘     └──┬──┘   # Links go from a parent to a child and describe an argument relation. The child is used as an
      │           │      # argument by the parent.
     ┌┴┐         ┌┴┐
     │+│         │-│     # Nodes with arguments are regular nodes, they are intermediate objects.
     └┬┘         └┬┘
      └─────┬─────┘
            │
           ┌┴┐
           │*│
           └┬┘
            │
           ┌┴┐
           │²│            # The root node describes the full expression. It is used to initialize the traversal logic.
           └─┘
```
Nodes in the traversal logic are all intermediate operations and leaves are the values on which the operations
are evalauted. Links indicate the relations between nodes. Child nodes serve as arguments to their respective parent node.

The graph is traversed with a depth-first search (DFS). On each node the first link is visited, followed by the
child of the link. Recursively, the child visits all links and so on. The above example would have the following walk:
```
Node ²
Link ², *
Node *
Link *, +
Node +
Link +, a
Leaf a
Link +, b
Leaf b
Link *, -
...
```
Information can be provided or stored in the implementation or propagated via the arguments of the function calls. The
default implementation forwards all arguments from the initial call.

Specializations for some common use cases are available:
 - [ForEachLeafLogic](@ref codi::ForEachLeafLogic): Default traversal of the tree with no logic in the nodes and links.
         Calls the function [handleActive](@ref codi::ForEachLeafLogic::handleActive()) for all left hand side
         expression objects. [handleConstant](@ref codi::ForEachLeafLogic::handleConstant()) is called for all constant
         expression objects.
 - [JacobianComputationLogic](@ref codi::JacobianComputationLogic): Evaluates the [reverse AD logic](@ref sec_forwardAD)
          for each link. The first user argument is expected to be the seed value for the output of the expression. For each
          left hand side expression object the method
          [handleJacobianOnActive](@ref codi::JacobianComputationLogic::handleJacobianOnActive()) is called. In this case, the first
          user argument contains the derivative of the origin node with respect to the current leaf times the
          initial seeding.

The [CompileTimeTraversalLogic](@ref codi::CompileTimeTraversalLogic) works in the same way as the
TraversalLogic and results mostly in the same order of calls. The only difference is that binary expressions
perform a reduction of the results from the two children.
