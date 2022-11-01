#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"

#include "CONSTANTS.H"
#include "memusage.H"
#include "Proto.H"
#include "ProtoTestCommon.H"



/**
End pirated bits.   Back to ChomboLand.
**/
  
shared_ptr<pr_amr_poisson>  >
getAMRPoissonSolver(const shared_ptr<pr_amrgrid>& pr_hierarchy)
{
  const auto& hierarchy = *pr_hierarchy;
  unsigned int domsize = hierarchy[0].problemDomain.sizes()[0];
  Real dx  = 1./Real(domsize);
  std::array<double, DIM> dxarray;
  dxarray.fill(dx);
  typedef BoxOp_Laplace<double> POISSON_OP;
  
  shared_ptr<AMRSolver_FASMultigrid<POISSON_OP, double> >  retval =
    new AMRSolver_FASMultigrid<POISSON_OP, double>(hierarchy, dxarray);
  
  return retval;
}
/**
 **/
int runFASSolver()
{
  //Proto has to drive the grid generation
  //get the grid hierarchy
  shared_ptr<pr_amr_grid>   pr_hierarchy =  getProtoAMRGrid();
  //get the chombo equivalent of the grid hierarchy
  shared_ptr<ch_amr_grid>   ch_hierarchy =  getChomboAMRGrid(pr_hierarchy);
  
  //get the poisson solver
  shared_ptr<pr_amr_solver> amrpoisPtr = getAMRPoissonSolver(pr_hierarchy);
  
}

/**
 **/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = runFASSolver();
    status += solverStatus;
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
