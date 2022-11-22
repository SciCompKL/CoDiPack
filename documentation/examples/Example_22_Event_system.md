Example 22 - Event system {#Example_22_Event_system}
=======

**Goal:** Use CoDiPack's event system to gain insight into the AD workflow.

**Prequesties:** \ref Tutorial_02_Reverse_mode_AD

**Function:**
\snippet examples/Example_22_Event_system.cpp Function

CoDiPack features an event system that can be used, for example,
- to gain insight into the AD workflow,
- to debug AD, from the level of tapes to the level of individual statements,
- to insert custom code into the AD workflow,
- to monitor the performance of AD constructs.

For each event, custom callback functions can be registered. CoDiPack invokes these callbacks as the corresponding
events occur.

We begin with the definition of callbacks for the AD workflow.

\snippet examples/Example_22_Event_system.cpp AD Workflow callback definitions

In `onTapeRegisterInput`, we make use of the possibility to associate custom data with the callback.

The callbacks in this example are written for [codi::RealReverse](@ref codi::RealReverse) and its corresponding tape.
As the event system is tape specific, signatures have to be adapted if other tapes are used. You can refer to the event
test system for general, templated callbacks.

We complement the AD workflow callbacks by callbacks for individual statements.

\snippet examples/Example_22_Event_system.cpp Statement callback definitions

AD workflow events are enabled by default. Statement events have to be enabled by compiling with the flag
`-DCODI_StatementEvents`, which enables the switch  [Config::StatementEvents](@ref codi::Config::StatementEvents). There
are also events for preaccumulation and index management, both with corresponding flags and switches.

At the beginning of the usual AD workflow, we have to register the callbacks.

\snippet examples/Example_22_Event_system.cpp Callback registration

Note that `onTapeRegisterInput` is registered twice, once with custom data and once without.
The required callback signatures can be looked up in [codi::EventSystemBase](@ref codi::EventSystemBase) and
[codi::EventSystem](@ref codi::EventSystem).
They should also be displayed by IDE code completion for the corresponding `register*` call.
The following output is produced by the example.

~~~~{.txt}
TapeStartRecording
TapeRegisterInput value 1 identifier 1
TapeRegisterInput value 1 identifier 1
  custom data 42
TapeRegisterInput value 2 identifier 2
TapeRegisterInput value 2 identifier 2
  custom data 42
TapeRegisterInput value 3 identifier 3
TapeRegisterInput value 3 identifier 3
  custom data 42
TapeRegisterInput value 4 identifier 4
TapeRegisterInput value 4 identifier 4
  custom data 42
TapeRegisterInput value 5 identifier 5
TapeRegisterInput value 5 identifier 5
  custom data 42
StatementStoreOnTape lhsIdentifier 6 newValue 1 numActiveVariables 1
  1 1;
StatementStoreOnTape lhsIdentifier 7 newValue 1 numActiveVariables 1
  1 1;
StatementStoreOnTape lhsIdentifier 8 newValue 3 numActiveVariables 2
  6 1; 2 1;
StatementStoreOnTape lhsIdentifier 9 newValue 2 numActiveVariables 2
  7 2; 2 1;
StatementStoreOnTape lhsIdentifier 10 newValue 6 numActiveVariables 2
  8 1; 3 1;
StatementStoreOnTape lhsIdentifier 11 newValue 6 numActiveVariables 2
  9 3; 3 2;
StatementStoreOnTape lhsIdentifier 12 newValue 10 numActiveVariables 2
  10 1; 4 1;
StatementStoreOnTape lhsIdentifier 13 newValue 24 numActiveVariables 2
  11 4; 4 6;
StatementStoreOnTape lhsIdentifier 14 newValue 15 numActiveVariables 2
  12 1; 5 1;
StatementStoreOnTape lhsIdentifier 15 newValue 120 numActiveVariables 2
  13 5; 5 24;
StatementStoreOnTape lhsIdentifier 16 newValue 15 numActiveVariables 1
  14 1;
TapeRegisterOutput value 15 identifier 16
StatementStoreOnTape lhsIdentifier 17 newValue 120 numActiveVariables 1
  15 1;
TapeRegisterOutput value 120 identifier 17
TapeStopRecording
TapeEvaluate begin from [0, 0, [0, 20, [0, 17, 17]]] to [0, 0, [0, 0, [0, 0, 0]]]
StatementEvaluate lhsIdentifier 17 numAdjoints 1
  2
StatementEvaluate lhsIdentifier 16 numAdjoints 1
  1
StatementEvaluate lhsIdentifier 15 numAdjoints 1
  2
StatementEvaluate lhsIdentifier 14 numAdjoints 1
  1
StatementEvaluate lhsIdentifier 13 numAdjoints 1
  10
StatementEvaluate lhsIdentifier 12 numAdjoints 1
  1
StatementEvaluate lhsIdentifier 11 numAdjoints 1
  40
StatementEvaluate lhsIdentifier 10 numAdjoints 1
  1
StatementEvaluate lhsIdentifier 9 numAdjoints 1
  120
StatementEvaluate lhsIdentifier 8 numAdjoints 1
  1
StatementEvaluate lhsIdentifier 7 numAdjoints 1
  240
StatementEvaluate lhsIdentifier 6 numAdjoints 1
  1
TapeEvaluate end from [0, 0, [0, 20, [0, 17, 17]]] to [0, 0, [0, 0, [0, 0, 0]]]
f(1 .. 5) = (15, 120)
df/dx (1 .. 5) [1 2]^T = (241 121 81 61 49)
TapeReset position [0, 0, [0, 0, [0, 0, 0]]] clear adjoints 1
~~~~

Note two calls of onTapeRegisterInput per input. There is output for the ten statements in `func` as well as two copy
statements in the course of `registerOutput` calls. All of them are complemented by `StatementEvaluate` events. The
statement events are surrounded by AD workflow events.

Here is the full code for the example.

\snippet examples/Example_22_Event_system.cpp Example 22 - Event system
