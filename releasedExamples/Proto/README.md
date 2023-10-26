# Chombo_3.2/releasedExamples/Proto
* A set of examples to show how to use Chombo with proto for performance portability.


# Chombo's AMRMultiGrid with Proto
* This is the easiest way to  make your Chombo code use Proto for performance portability.
* Chombo's AMRMultiGrid handles the Martin-Cartwright algorithm dance.
* Proto handles device data and executes stencil operations.

# Rules of engagement
* Proto knows nothing about Chombo.
* The interaction proto is mediated through a plain old data  interface.
* Make your code compile with NAMESPACE=TRUE.

# Procedure for porting your Chombo code to Proto.
* This directory holds a bunch of very small examples.
* This will hold solvers for all Chombo (non-EB) elliptic operators (EB users should consider Chombo4).
* This will also have all the tools for the conversion.
* Each proto operator example will also have a counterpart that uses the standard Chombo operator.
* Choose which solver you would like to switch to using proto operators.
* Compare its new operator example to the standard Chombo operator example.
* Use the tools provided to switch to the new solver for your own code.

# Directories
* Poisson is a Chombo port of a proto example (single grid Poisson solver using multigrid).
* common holds code used by several examples
* data_transfer tests moving data between proto and Chombo data types.
* amr_helmoltz uses AMRMultiGrid and proto to solve the constant-coefficient Helmholtz equation.
* _old_helmoltz does the same thing with AMRPoissonOp so the two can be easily compared.


# We plan on also porting in 2023:
* Variable coefficient Helmholtz (VCAMRPoissonOp). 
* The viscous tensor equation (ViscousTensorOp). 
* The magnetic resistivity  equation (ResitivityOp). 

