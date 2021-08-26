Coding guidelines {#CodingGuidelines}
===========================

Identation
------
 - Only spaces
 - Tab is 2 spaces long
 - Declarations spread over multiple lines need to be indented once more than the adjacent block

Classess/Structures
------
 - Always use `public/private/protected`.
 - `public/private/protected` is indented.
 - Members, functions, etc. are indented again (once from `public/private/proteced`).
 - Template arguments:
  * Start with '_'.
  * Template arguments are made available with `using` declarations without the '_'.
  * There, a default type needs to be declared using `DECLARE_DEFAULT` (See \ref TemplateDeclaration for details).

Argument/member/variable declarations:
------
 - const declaration is on the right hand side (e.g. `T const`).
 - Pointer/reference symbol is part of the type (e.g. `T& a`).

Namespaces
------
 - Everything is indented once.

for/if/while/switch
------
 - Space between bracket (e.g. `while (true)`).
 - `case` statements are indented once.
 - Code inside `case` statements is indented again (once from the `case` statement itself).
 - Curly brackets are mandatory.

Functions
------
 - Declaration order: `static CODI_INLINE <ret> <name>`
 - Single line declaration:
  * Curly brackets are on the same line.
 - Multi line declaration:
  * Arguments are indented twice.
  * Closing bracket and opening curly bracket are on a new line not indented. (Not enforced yet)

Files
------
 - Includes are sorted by name.
 - First external includes.
 - Second internal includes.
 - Includes are relative to the file location.
 - 120 maximum characters in a line.

Example
----
```
// External includes first, sorted by name
#include <iostream>
#include <sstream>

// Internal includes second, sorted by name
#include "../tapes/interfaces/reverseTapeInterface.hpp"
#include "../traits/realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  // Indentation in namespace is one tab
  // Tabs are only spaces and 2 spaces long

  struct Test {
    public: // Always use access specifier, access specifier are indented

      int v; // Member variables are indented once again.
  };

  template<typename T_T> // Template parameters are declared with a 'T_' prefix
  struct Test2 {
    public:
      using T = DECLARE_DEFAULT(T_T, int); // Template arguments are made available with `using` declarations without 
                                           // the 'T_'. A default type needs to be declared                            ^
                                           //                                                                          |
                                           //                                            120 maximum character in a line

  };

  struct Test3 {
    public:
      int constexpre m1;  // const declarations are on the right hand side
      int const* m2;      // const declarations are on the right hand side

      CODI_INLINE int const* func1(int const& offset) {  // const declarations are on the right hand side
        return &m2[offset];
      }

      CODI_INLINE int func2(int& a) {  // Pointer/reference symbol is part of the type
        return ++a;
      }
  };

  CODI_INLINE void func(int& i) {
    if (i > 0) {            // Space between bracket in if
      while (i != 0) {      // Space between bracket in while
        i -= 1;
      }
    } else {
      switch (i) {          // Space between bracket in switch
        case -1:            // case statements are indented once
          i = 10;           // case bodies are indented once more
          break;
        case -2:
          i = 100;
          break;
        default:
          i = 0;
          break;
      }
    }

    for (int j = 0; j < 10; j += 1) { // Space between bracket in while
      i += j;
    }                                 // Always use curly brackets
  }

  CODI_INLINE int func(
      int a,  // Arguments of a multi line function declrations are indented twice
      int b,
      int c,
  ) { // Closing bracket is on a new line
    return a + b + c;
  }
}
```
