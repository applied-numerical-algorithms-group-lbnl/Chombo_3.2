# Chombo_3.2/releasedExamples/Proto/amr_viscous_tensor
* amr_viscous_tensor uses AMRMultiGrid and Proto_FastVTO to solve the variable-coefficient viscous tensor equation.
* There are some obvious optimizations that Proto_FastVTO will have that Proto_Viscous_Tensor_Op does not.
* Yes, Proto_FastVTO is a dumb name and it will go away Pretty Soon (tm).
* The plan is for Proto_FastVTO to become Proto_Viscous_Tensor_Op and for the old operator to go away.


# Files
* amr_fast_vto.cpp : main file that runs solver.
* _inputs/case_0.inputs : inputs for single level with domain decomposition solves 
* _inputs/case_1.inputs : simple AMR case

