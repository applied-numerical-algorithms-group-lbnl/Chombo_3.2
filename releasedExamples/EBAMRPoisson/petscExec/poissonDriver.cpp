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
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "PoissonUtilities.H"

#include "EBFABView.H"
#include "EBDebugDump.H"

#include "EBAMRPoissonOp.H"
#include "EBLevelDataOps.H"
#include "BaseBCValue.H"
#include "BaseDomainBC.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "NeumannPoissonEBBC.H"

#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLevelDataOps.H"
#include "EBSimpleSolver.H"
#include "BiCGStabSolver.H"
#include "EBEllipticLoadBalance.H"
#include "EBLoadBalance.H"
#include "LoadBalance.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"
#include "EBPoissonPetscSolver.H"
#include "memusage.H"
#include "memtrack.H"


/********/
void solve(const PoissonParameters&  a_params)
{
  pout() << "ignoring params numLevels and just doing single level solve"  << endl;
  Vector<Box> boxes;
  Vector<int> procs;
  Real alpha, beta;
  int blockfactor, maxboxsize;
  ParmParse pp;
  pp.get("block_factor" , blockfactor);
  pp.get("max_grid_size", maxboxsize);
  pp.get("alpha"        , alpha);
  pp.get("beta"         , beta);
  domainSplit(a_params.coarsestDomain, boxes, maxboxsize, blockfactor);
  LoadBalance(procs, boxes);
  
  DisjointBoxLayout grids(boxes, procs, a_params.coarsestDomain);
  EBISLayout        ebisl;
  int numGhost = 4;
  const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(ebisl, grids, a_params.coarsestDomain, numGhost);

  Vector<DisjointBoxLayout> vgrids(1, grids);
  Vector<EBISLayout>        vebisl(1, ebisl);
  EBCellFactory cellfact(ebisl);

  LevelData<EBCellFAB> phi(grids, 1, a_params.ghostPhi, cellfact);
  LevelData<EBCellFAB> rhs(grids, 1, a_params.ghostRHS, cellfact);

  EBLevelDataOps::setVal(phi, 0.0);
  EBLevelDataOps::setVal(rhs, 1.0);



  //create the solver
  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > opfactory;
  //not  important in this context
  int numprecond = 0; int  relaxtype = 0;
  Real time = 0; 
  getEBAMRPFactory(opfactory, vgrids, vebisl, a_params, numprecond, relaxtype, time, alpha, beta);

  EBAMRPoissonOp*  rawptr = dynamic_cast<EBAMRPoissonOp*>(opfactory->AMRnewOp(a_params.coarsestDomain));
  if(rawptr == NULL)
    { 
      MayDay::Error("dynamic cast failed");
    }
  RefCountedPtr<EBAMRPoissonOp> opptr(rawptr);

  pout() << "defining solver with homogeneous bcs " << endl;
  EBPoissonPetscSolver solver(opptr);

  pout() << "solving " << endl;

  solver.solve(phi, rhs);

#ifdef CH_USE_HDF5
  bool fileOut;
  pp.get("do_file_output", fileOut);
  if (fileOut)
    {
      pout() << "outputting solution to file" << endl;
      string filephi = string("phi.hdf5");
      string filerhs = string("rhs.hdf5");
      writeEBLevelname(&phi, filephi.c_str());
      writeEBLevelname(&rhs, filerhs.c_str());
    }
#endif
  //checking the answer
  //now let us compute the residual and output its norm
  LevelData<EBCellFAB> resid;
  opptr->create(resid, rhs);
  opptr->residual(resid, phi, rhs, false);
  
  int comp = 0;
  Real volume;
  for(int inorm = 0; inorm < 3; inorm++)
    {
      Real norm = EBArith::norm(volume, resid, grids, ebisl,  comp, inorm);
      if(inorm == 0)
        {
          pout() << "final maximum \t norm of the residual = " << norm << endl;
        }
      else
        {
          pout() << "final " << inorm << "th \t norm of the residual = " << norm << endl;
        }
    }
}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {
    CH_TIME("main_stuff");
    PetscInt ierr;
    ierr = PetscInitialize(&argc, &argv, "./petscrc",PETSC_NULL); CHKERRQ(ierr); 
    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters params;

    //read params from file
    getPoissonParameters(params);

    //define geometry from given params
    definePoissonGeometry(params);

    //solve the stinking problem and output everything
    solve(params);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();


    ierr = PetscFinalize(); CHKERRQ(ierr); 
  }    // End scoping trick

  dumpmemoryatexit();
#ifdef CH_MPI
  pout() << "dumping timers" << endl;
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

}

