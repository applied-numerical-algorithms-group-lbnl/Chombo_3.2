# Chombo_3.2/releasedExamples/Proto
* Proto is our library that we use for performance portability.
* Chombo can use proto to run solvers on accellerators (GPUs).
* This set of examples is designed to show how to use Chombo with proto.

# Chombo's AMRMultiGrid with Proto
* This is the easiest way to  make your Chombo code use Proto for performance portability.
* Chombo's AMRMultiGrid handles the Martin-Cartwright algorithm dance.
* Proto handles device data and executes stencil operations.

# Rules of engagement
* Proto knows nothing about Chombo.
* Interactions with  proto are  mediated through a plain old data  interface.
* These examples rely heavily upon namespaces.
* To use this stuff, you have to make your code compile with NAMESPACE=TRUE.

# Procedure for porting your Chombo code to Proto.
* This directory holds a bunch of very small examples.
* This will hold proto-based solvers for all Chombo (non-EB) elliptic operators (EB users should consider Chombo4).
* This will also have all the tools for the conversion.
* Each proto operator example will also have a counterpart that uses the standard Chombo operator.
* Choose which solver you would like to switch to using proto operators.
* Compare its new operator example to the standard Chombo operator example.
* Use the tools provided to switch to the new solver for your own code.

# Directories
* common holds code used by several examples
* test_harness holds the machinery to do and compare multiple runs of this stuff.
* data_transfer tests moving data between proto and Chombo data types.
* amr_helmoltz uses AMRMultiGrid and Proto_Helmholtz_Op to solve the constant-coefficient Helmholtz equation.
* _old_helmoltz does the same thing with AMRPoissonOp so the two can be easily compared.
* amr_conductivity uses AMRMultiGrid and Proto_Conductivity_Op to solve the variable-coefficient Helmholtz equation.
* _old_conductivity  does the same thing with VCAMRPoissonOp so the two can be easily compared.
* amr_viscous_tensor uses AMRMultiGrid and Proto_Viscous_Tensor_Op to solve the variable-coefficient Viscous Tensor equation.
* _old_viscous_tensor  does the same thing with ViscousTensorOp so the two can be easily compared.
* amr_resisitivity uses AMRMultiGrid and Proto_Resistivity_Op to solve the variable-coefficient magnetic resistivity equation.
* _old_resistivity  does the same thing with ResistivityOp so the two can be easily compared.
* _proto_sg_poisson is a Chombo port of a proto example to demonstrate Proto syntax.

