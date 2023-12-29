Tutorials & Examples {#TutorialsAndExamples}
=======

Tutorials provide a detailed explanation of the most common CoDiPack features. They try to explain nearly everything so
that also beginners should have no problem to apply AD to a code.

Examples provide in most cases only the code for using the necessary features. Some examples have some helpful remarks
or pointers to other features.

|Tutorial | |
|:--------|:--------|
| \subpage Tutorial_01_Forward_mode_AD "" | Explains the usage of the codi::RealForward type. |
| \subpage Tutorial_02_Reverse_mode_AD "" | Explains the usage of the codi::RealReverse type. |
| \subpage Tutorial_03_Full_jacobian_computation "" | Full Jacobian computation in CoDiPack with the forward and reverse mode. |
| \subpage Tutorial_04_Vector_mode_AD "" | Vector mode examples with CoDiPack. |
| \subpage Tutorial_05_Repeated_tape_recordings "" | How different tapes are recorded in in CoDiPack. |
| \subpage Tutorial_06_Higher_order_types_helper_acces "" | Describes how higher order AD types can be constructed and used. |


| Example | |
|:--------|:--------|
| \subpage Example_01_Old_tangent_leftovers_forward_mode "" | Shows possible errors if the computational path is changed. |
| \subpage Example_02_Custom_adjoint_vector_evaluation "" | How custom types can be used in an reverse evaluation, on tapes that are already recorded. |
| \subpage Example_03_Positional_tape_evaluations "" | Demonstrates how to evaluate only parts of a tape. |
| \subpage Example_04_Higher_order_types_helper_access_compile_time "" | Example of the higher order AD types accessed with compile time constructs. |
| \subpage Example_05_Higher_order_types_direct_access "" | Example of higher order AD types accessed with the basic CoDiPack data functions. |
| \subpage Example_06_Forward_tape_evaluation "" | Demonstrates how a tape can be evaluated in a forward AD mode. |
| \subpage Example_07_Primal_tape_evaluation "" | Demonstrates how primal value tapes can be reevaluated for a different point without recording a new tape. |
| \subpage Example_08_Vector_helper_interface_access "" | Using a generalized interface for the custom vector access. |
| \subpage Example_09_OpenMP_reverse_evaluation "" | Shows how OpenMP can be used to evaluate the same tape concurrently with multiple threads. |
| \subpage Example_10_External_function_helper "" | Ease of access structure for adding custom function to the tape. |
| \subpage Example_11_External_function_user_data "" | How user data can be added to external functions. |
| \subpage Example_12_Manual_statement_creation "" | Describes how custom derivatives for small statements can be added to the tape. |
| \subpage Example_13_MPI_communication "" | Demonstrates how MPI constructs can be handled with CoDiPack types. |
| \subpage Example_14_ReferenceActiveType "" | Shows how the codi::ReferenceActiveType class is used. |
| \subpage Example_15_Preaccumulation_of_code_parts "" | Provides an example of memory reduction through preaccumulation |
| \subpage Example_16_TapeHelper "" | Demonstrates a simpler interface for the CoDiPack types. |
| \subpage Example_17_EvaluationHelper "" | Automatic computation of Hessians and Jacobians for function objects. (EvaluationHelper) |
| \subpage Example_18_EvaluationHelper_function_objects "" | Function object kinds for the EvaluationHelper (Example 17). |
| \subpage Example_19_EvaluationHelper_handle_creation "" | Handle creation for the EvaluationHelper (Example 17). |
| \subpage Example_20_Aggregated_active_type_handling "" | Generalized data extraction of aggregated active types in external functions. |
| \subpage Example_21_Special_handling_of_linear_system_solvers "" | Add special handling linear system solvers. |
| \subpage Example_22_Event_system "" | Use CoDiPack's event system to gain insight into the AD workflow. |
| \subpage Example_23_OpenMP_Parallel_Codes "" | Use CoDiPack together with OpDiLib for the differentiation of OpenMP parallel codes. |
| \subpage Example_24_Enzyme_external_function_helper "" | Adding Enzyme-differentiated functions to the CoDiPack tapes. |

The graph shows how the tutorials and examples are connected. Usually it is better to understand first the prerequisites
of a tutorial/example before reading the actual example.

\dotfile user/TutorialsGraph.dot
