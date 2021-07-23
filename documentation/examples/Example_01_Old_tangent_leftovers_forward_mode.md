Example 1 - Old tangent leftovers in forward mode AD {#Example_01_Old_tangent_leftovers_forward_mode}
=======

**Goal:** Learn how side effects can effect the tangent computation.

**Prerequisite:** \ref Tutorial_01_Forward_mode_AD

**Function:**
\snippet examples/Example_01_Old_tangent_leftovers_forward_mode.cpp Function

**Full code:**
\snippet examples/Example_01_Old_tangent_leftovers_forward_mode.cpp Example 1 - Old tangent leftovers forward mode

The computational path of the function is changed via the parameter `updateGlobal`. If this parameter is `true` then
`x` also enters the computation of the global variable. If the parameter is `false` then the global value is seen as a
constant with respect to AD.

The three different calls demonstrate the error in the tangent computation. The first and third call are correct, the
second one is wrong. Here, we fix the issue by directly resetting the tangent of the `global` variable. An other option
would be to call `func(x, true)` again, with a tangent value of zero for `x`.

#### Notes ####

This problem is not specific to CoDiPack. Nearly all tapeless forward AD tools suffer from this problem. One possible
final fix is to implement a Gradient type, that tags tangents with an epoch. The user would then need to mange the
global valid epoch and the Gradient type can ignore tangents from older epochs.
