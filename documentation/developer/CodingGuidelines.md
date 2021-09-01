Coding guidelines {#CodingGuidelines}
===========================

Indentation
------
 - Use only spaces.
 - Tab is 2 spaces long.
 - Declarations spread over multiple lines need to be indented once more than the adjacent block.

Classes/Structures
------
 - Always use `public/private/protected`.
 - `public/private/protected` is indented.
 - Members, functions, etc. are indented again (once from `public/private/proteced`).
 - Template arguments:
  * Start with 'T_'.
  * Template arguments are made available by `using` declarations without the 'T_'.
  * There, a default type needs to be declared using `DECLARE_DEFAULT` (see \ref TemplateDeclaration for details).

Argument/member/variable declarations:
------
 - const declaration is on the right hand side (e.g. `T const`).
 - Pointer/reference symbol is part of the type (e.g. `T& a`).

Namespaces
------
 - Everything is indented once.

for/if/while/switch
------
 - Space between keyword and bracket (e.g. `while (true)`).
 - `case` statements are indented once.
 - Code inside `case` statements is indented again (once from the `case` statement itself).
 - Curly brackets are mandatory.

Functions
------
 - Declaration order: `static CODI_INLINE <ret> <name>`
 - Single-line declaration:
  * Curly brackets are on the same line.
 - Multiline declaration:
  * Arguments are indented twice.
  * Closing bracket and opening curly bracket are on a new line and not indented (not enforced yet).

Files
------
 - Includes are sorted by name.
 - First external includes, then internal includes.
 - Includes are relative to the file location.
 - At most 120 characters in a line.

Documentation
------
 - The documentation should be ASCII readable.
 - Documentation of all classes, members, functions and types (no parameters).
 - Template parameters are documented for classes.
 - Interface declarations need to have a detailed documentation, especially for the class itself. This should contain
   * the abstraction concept,
   * usage of the interface and
   * examples.

   Other functions in the interface should be documented according to their exposal to the user.
 - Implementation of interfaces can copy the documentation from the interface. Only important details need to be added.
   * Short example (preferred):
```
/// \copydoc GradientAccessTapeInterface::getGradient
```
   * Long example:
```
/// \copydoc GradientAccessTapeInterface::getGradient <br><br>
/// Implementation: Details
```
   * Long example:
```
/** \copydoc GradientAccessTapeInterface::getGradient <br><br>
 * Implementation: Details
 */
```

Example
----
```
// External includes first, sorted by name.
#include <iostream>
#include <sstream>

// Internal includes second, sorted by name.
#include "../tapes/interfaces/reverseTapeInterface.hpp"
#include "../traits/realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  // One level of indentation corresponds to one tab.
  // One tabs consists of 2 spaces.
  // Everything in a namespace is indented once.

  struct Test {
    public: // Always use access specifier, access specifiers are indented.

      int v; // Members are indented once again.
  };

  template<typename T_T> // Template parameters are declared with a 'T_' prefix
  struct Test2 {
    public:
      using T = DECLARE_DEFAULT(T_T, int); // Template arguments are made available by `using` declarations without
                                           // the 'T_'. A default type must be declared.                               ^
                                           //                                                                          |
                                           //                                           at most 120 characters in a line

  };

  struct Test3 {
    public:
      int constexpr m1;  // const declarations are on the right hand side.
      int const* m2;      // const declarations are on the right hand side.

      CODI_INLINE int const* func1(int const& offset) {  // const declarations are on the right hand side.
        return &m2[offset];
      }

      CODI_INLINE int func2(int& a) {  // Pointer/reference symbol is part of the type.
        return ++a;
      }
  };

  CODI_INLINE void func(int& i) {
    if (i > 0) {            // Space between if and bracket.
      while (i != 0) {      // Space between while and bracket.
        i -= 1;
      }
    } else {
      switch (i) {          // Space between switch and bracket.
        case -1:            // case statements are indented once.
          i = 10;           // case bodies are indented once more.
          break;
        case -2:
          i = 100;
          break;
        default:
          i = 0;
          break;
      }
    }

    for (int j = 0; j < 10; j += 1) { // Space between for and bracket.
      i += j;
    }                                 // Always use curly brackets.
  }

  CODI_INLINE int func(
      int a,  // Arguments of multiline function declarations are indented twice.
      int b,
      int c
  ) { // Closing bracket is on a new line.
    return a + b + c;
  }
}
```
