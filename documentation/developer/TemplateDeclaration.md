Template declarations {#TemplateDeclaration}
=======

Template declarations for structures and classes follow a special layout in CoDiPack. There are three rules any
developer needs to follow:
 - Template arguments begin with a `T_`.
 - Template arguments are redeclared inside the class with using without the `T_` prefix.
 - A default declaration of the template type needs to be done.
 
The first two rules are simple and an example is:
```{.cpp}
template<typename T_T>
struct Test {
  public:
    using T = T_T;
};
```

The third rule is in place to allow for auto completion in IDEs on template arguments. In your IDE you can use the
preprocessor flag `-DCODI_IDE=1` to enable this feature. In order for the IDEs to be able to know about the template
type the default type declarations are required. CoDiPack defines preprocessor helpers to declare these in the file
[macros.hpp](@ref TemplateDeclarationHelpers). All declarations have to use `CODI_DECLARE_DEFAULT` or short `CODI_DD`.
The first argument is the template argument and the second argument is the default type declaration. Since the
preprocessor engine does not know about templates, template argument lists are parsed as arguments to the preprocessor
macro. It is therefore necessary to wrap these declarations in `CODI_TEMPLATE` or short `CODI_T` macros.
While using defining the defaults for `CODI_IDE` the clangd compatibility should be kept in mind. This may make it
necessary to not declare a default template argument in special cases.

Example declarations are
```{.cpp}
  using Tape = CODI_DECLARE_DEFAULT(T_Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, EmptyPosition>));
  using Real = CODI_DD(T_Real, double);
  using Operation = CODI_DD(CODI_T(T_Operation<Real>), CODI_T(BinaryOperation<Real>));
  using Chunk = CODI_DD(T_Chunk, CODI_T(Chunk1<CODI_ANY>));
```



