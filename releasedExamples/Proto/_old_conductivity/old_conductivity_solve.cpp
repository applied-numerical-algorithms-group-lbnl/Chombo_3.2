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
#include "ParmParse.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "BRMeshRefine.H"
#include "FABView.H"
#include "DebugDump.H"
#include "VCLocalFunctions.H"
#include "MultilevelLinearOp.H"
#include "BiCGStabSolver.H"
#include "AMRIO.H"



namespace Chombo
{

/******/
  int
  poissonSolve(Vector<       LevelData< FArrayBox >* >& a_phi,
               const Vector< LevelData< FArrayBox >* >& a_rhs,
               const Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(nlevels);
               const Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(nlevels);
               const Vector< DisjointBoxLayout >      & a_grids,
               const VCPoissonParameters              & a_params)
  {
    int nlevels = a_params.numLevels;

    // set up solver
    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > opFactory
      = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
      (defineOperatorFactory(a_grids, vectDomains, aCoef, bCoef, a_params));

    int lBase = 0;

    BiCGStabSolver<Vector<LevelData<FArrayBox>* > > solver;

    bool homogeneousBC = false;
    solver.define(&mlOp, homogeneousBC);
    solver.m_verbosity = a_params.verbosity;
    solver.m_normType = 0;
    solver.m_eps = tolerance;
    solver.m_imax = max_iter;

    solver.solve(a_phi, a_rhs);

    int exitStatus = solver.m_exitStatus;
    // note that for AMRMultiGrid, success = 1.
    exitStatus -= 1;
    return exitStatus;

  }
  int 
  defineData(Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_phi,
             Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_rhs,
             Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_aCoef,
             Vector< RefCountedPtr< LevelData<   FluxBox > > >& a_bCoef;
             const Vector< DisjointBoxLayout >                & a_grids,
             const VCPoissonParameters                        & a_params)
  {
    a_phi  .resize(a_grids.size());
    a_rhs  .resize(a_grids.size());
    a_aCoef.resize(a_grids.size());
    a_bCoef.resize(a_grids.size());
    
    for(int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      a_phi  [ilev] = RefCountedPtr< LevelData< FArrayBox > > >( new   LevelData< FArrayBox > > (grids[ilev], 1, IntVect::Unit));
      a_rhs  [ilev] = RefCountedPtr< LevelData< FArrayBox > > >( new   LevelData< FArrayBox > > (grids[ilev], 1, IntVect::Zero));
      a_aCoef[ilev] = RefCountedPtr< LevelData< FArrayBox > > >( new   LevelData< FArrayBox > > (grids[ilev], 1, IntVect::Zero));
      a_bCoef[ilev] = RefCountedPtr< LevelData<   FluxBox > > >( new   LevelData<   FluxBox > > (grids[ilev], 1, IntVect::Zero));
    }  ///end loop voer levels
  }    //end function defineData

  /***************/
  void outputData(const Vector< LevelData< FArrayBox >* >&   a_phi,
                  const Vector< LevelData< FArrayBox >* >&   a_rhs,
                  const Vector< DisjointBoxLayout       >&   a_grids,
                  const VCPoissonParameters              &   a_params)
  {
#ifdef CH_USE_HDF5

    string         phiFileName(          "phi.hdf5");
    Vector<string> phiVarNames(1, string("phi"));
    string         rhsFileName(          "rhs.hdf5");
    Vector<string> rhsVarNames(1, string("rhs"));

    Real fakeTime = 1.0;
    Real fakeDt = 1.0;
    WriteAMRHierarchyHDF5(phiFileName,  a_grids,
                          a_phi, phiVarNames,
                          a_params.coarsestDomain.domainBox(),
                          a_params.coarsestDx,
                          fakeDt, fakeTime,
                          a_params.refRatio,
                          a_params.numLevels);
    WriteAMRHierarchyHDF5(rhsFileName,  a_grids,
                          a_phi, rhsVarNames,
                          a_params.coarsestDomain.domainBox(),
                          a_params.coarsestDx,
                          fakeDt, fakeTime,
                          a_params.refRatio,
                          a_params.numLevels);

#endif

  }


  int local_test_harness()
  {
    VCPoissonParameters param;
    Vector<DisjointBoxLayout> grids;

    //read params from file
    getPoissonParameters(param);
    int nlevels = param.numLevels;
    Vector<RefCountedPtr<LevelData< FArrayBox > > > phi;
    Vector<RefCountedPtr<LevelData< FArrayBox > > > rhs;
    Vector<RefCountedPtr<LevelData< FArrayBox > > > aCoef;
    Vector<RefCountedPtr<LevelData<   FluxBox > > > bCoef;

    setGrids(grids,  param);

    defineData(phi, rhs, aCoef, bCoef,  nlevels, grids, param);
    poissonSolve(phi, rhs, grids,  param);
    outputData(  phi, rhs, grids, param);

    // clear memory
    for (int level = 0; level<phi.size(); level++)
    {
      if (phi[level] != NULL)
      {
        delete phi[level];
        phi[level] = NULL;
      }
      if (rhs[level] != NULL)
      {
        delete rhs[level];
        rhs[level] = NULL;
      }
    }

    return 0;
  }
} //end namespace Chombo



int main(int argc, char* argv[])
{
  int status = 0;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    status = Chombo::local_test_harness();
  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
