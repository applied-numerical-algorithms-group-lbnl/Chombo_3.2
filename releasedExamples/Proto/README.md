# Chombo3.2/releasedExamples/Proto
* A set of examples to show how to use Chombo with proto for performance portability.

# Expectations for developers:
* These codes are written with NAMESPACE=TRUE.
* The interaction proto is mediated through a plain old data (POD) interface.

# AMRMultiGrid with Proto
* One of the easiest way to effectively migrate to  proto.
* AMRMultiGrid handles the Martin-Cartwright algorithm dance.
* Proto handles  device data and executes stencil operations.


# Directories
* Poisson is a Chombo port of a proto example (single grid Poisson solver using multigrid).
* common holds code used by several examples
* data_transfer tests moving data between proto and Chombo data types.
* amr_helmoltz uses AMRMultiGrid and proto to solve the constant-coefficient Helmholtz equation.
* _old_helmoltz does the same thing with AMRPoissonOp so the two can be easily compared


# We plan on also porting in 2023:
* Variable coefficient Helmholtz (VCAMRPoissonOp), 
* The viscous tensor equation (ViscousTensorOp), 
* The magnetic resistivity  equation (ResitivityOp), 

