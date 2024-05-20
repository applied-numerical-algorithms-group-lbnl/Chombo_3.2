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
#include "BoxIterator.H"
#include "memusage.H"
#include "BCFunc.H"
#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CH_HDF5.H"



///
/**
  This simple example is meant to show how to use MPI_Comm_Split to 
  have more than one Chombo solving running simultaneously.
  The following code is a bunch of stuff adapted from AMRPoisson/execCell example
  hardwired to be single level.
 */
void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  a_values[0]=0.;
}

void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{

  if (!a_domain.domainBox().contains(a_state.box()))
  {
    Box valid = a_valid;
    for (int idir =0; idir <CH_SPACEDIM; idir++)
    {
      Box ghostBoxLo = adjCellBox(valid, idir, Side::Lo, 1);
      Box ghostBoxHi = adjCellBox(valid, idir, Side::Hi, 1);
      if (!a_domain.domainBox().contains(ghostBoxLo))
      {
        DiriBC(a_state,
               valid,
               a_dx,
               true,
               ParseValue,
               idir,
               Side::Lo,
               1);
      }

      if (!a_domain.domainBox().contains(ghostBoxHi))
      {
        DiriBC(a_state,
               valid,
               a_dx,
               true,
               ParseValue,
               idir,
               Side::Hi,
               1);
      }
    }
  }
}

void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
            Vector<ProblemDomain>& a_amrDomains,
            Vector<int>& a_refRatios,
            Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setRHS");

  for (int lev=0; lev<=a_finestLevel; lev++)
  {
    LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
    const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();

    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& thisRhs = levelRhs[levelDit];
      thisRhs.setVal(1.0);
    } // end loop over grids on this level
  } // end loop over levels
}



void
setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver,
            LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
            const Vector<DisjointBoxLayout>& a_amrGrids,
            const Vector<ProblemDomain>& a_amrDomains,
            const Vector<int>& a_refRatios,
            const Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setupSolver");

  ParmParse ppSolver("setupSolver");

  int numLevels = a_finestLevel+1;

  AMRPoissonOpFactory opFactory;

  // solving poisson problem here
  Real alpha =0.0;
  Real beta = 1.0;

  opFactory.define(a_amrDomains[0],
                   a_amrGrids,
                   a_refRatios,
                   a_amrDx[0],
                   &ParseBC, alpha, beta);

  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

  a_amrSolver->define(a_amrDomains[0], castFact,
                      &a_bottomSolver, numLevels);

  // multigrid solver parameters
  int numSmooth, numMG, maxIter, verbosity;
  Real eps, hang;
  ppSolver.get("num_smooth"    , numSmooth);
  ppSolver.get("num_mg"        , numMG    );
  ppSolver.get("max_iterations", maxIter  );
  ppSolver.get("tolerance"     , eps      );
  ppSolver.get("hang"          , hang     );
  ppSolver.get("verbosity"     , verbosity);

  Real normThresh = 1.0e-30;
  a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                   numMG, maxIter, eps, hang, normThresh);
  a_amrSolver->m_verbosity = verbosity;
}
  
int
runSolver(const DisjointBoxLayout  & a_grid,
          const ProblemDomain      & a_domain,
          const Real               & a_dx,
          int       a_icolor) //included in file names if positive
{
  CH_TIME("runSolver");

  ParmParse ppMain("runSolver");

  pout() << "runSolver: making vector structures" << endl;
  // single level only solve
  Vector<DisjointBoxLayout> amrGrids(1, a_grid);
  Vector<ProblemDomain>   amrDomains(1, a_domain);
  Vector<Real> amrDx(1, a_dx);
  Vector<int> refRatios(1, 2);
  int finestLevel = 0;

  pout() << "runSolver: defining the solver " << endl;
  // initialize solver
  AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;
  amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
  BiCGStabSolver<LevelData<FArrayBox> > bottomSol;
  bottomSol.m_verbosity  = 0;
  setupSolver(amrSolver, bottomSol, amrGrids, amrDomains,
              refRatios, amrDx, finestLevel);

  pout() << "runSolver: allocate solution and RHS" << endl;
  int numLevels = amrGrids.size();
  Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);

  for (int lev=0; lev<=finestLevel; lev++)
  {
    const DisjointBoxLayout& levelGrids = amrGrids[lev];
    phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
    rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
  }

  pout() << "runSolver: setting the rhs" << endl;
  setRHS(rhs, amrDomains, refRatios, amrDx, finestLevel );

  pout() << "runSolver: running amrmultigrid::solve" << endl;
  bool zeroInitialGuess = true;
  bool print            = true;
  bool forceHomogeneous = false;
  amrSolver->solve(phi, rhs, finestLevel, 0,
                   zeroInitialGuess, forceHomogeneous, print);

  
  pout() << "runSolver: clean up " << endl;
  for (int lev=0; lev<phi.size(); lev++)
  {
    delete phi[lev];
    delete rhs[lev];
  }
  delete amrSolver;

  return 0;
}

int
runColoredSolvers()
{
  ParmParse pp("runColoredSolvers");
  int nx, maxBoxSize, procsPerColor;
  pp.get("nx"        , nx         );
  pp.get("maxBoxSize", maxBoxSize );
  pp.get("procsPerColor" , procsPerColor  );

  pout() << "runColoredSolvers:nx         =" << nx         << endl;
  pout() << "runColoredSolvers:maxBoxSize =" << maxBoxSize << endl;
  pout() << "runColoredSolvers:procsPerColor  =" << procsPerColor  << endl;
  IntVect ivlo =        IntVect::Zero;
  IntVect ivhi = (nx-1)*IntVect::Unit;
  Box domain(ivlo, ivhi);
  Vector<Box> boxes;
  Vector<int> procs;
  domainSplit(domain, boxes,  maxBoxSize);
  double dx = 1./nx;
  LoadBalance(procs, boxes);

  ///
  // Get the rank and size in the original communicator
  int world_rank, world_size;
  MPI_Comm_rank(Chombo_MPI::comm, &world_rank);
  MPI_Comm_size(Chombo_MPI::comm, &world_size);

  ///actual timed test
  if(1)
  {
    CH_TIME("runColoredSolvers: actual test");
    // Determine color based on original rank
    int icolor = world_rank / procsPerColor;

    MPI_Comm color_comm;
    MPI_Comm_split(Chombo_MPI::comm, icolor, world_rank, &color_comm);
  
    int color_rank, color_size;
    MPI_Comm_rank(color_comm, &color_rank);
    MPI_Comm_size(color_comm, &color_size);
  
    pout() << "runColoredSolvers: world rank = " << world_rank << endl;
    pout() << "runColoredSolvers: world size = " << world_size << endl;
    pout() << "runColoredSolvers: color rank = " << color_rank << endl;
    pout() << "runColoredSolvers: color size = " << color_size << endl;
    pout() << "runColoredSolvers: proc color = " << icolor     << endl;
  

    DisjointBoxLayout dblWorld(boxes, procs, Chombo_MPI::comm);
    DisjointBoxLayout dblColor(boxes, procs, color_comm);
    int nColor = world_size/procsPerColor;
    pout() << "runColoredSolvers: running world solver Ncolor times" << endl;
    {
      CH_TIME("running world solver Ncolor times");
      for(int isolve = 0; isolve < nColor; isolve++)
      {
        pout() << "isolve = " << isolve << endl;
        // negative color indicates world communicator (for I/O)
        runSolver(dblWorld, domain, dx,  -4586);  
      }
    }
    pout() << "runColoredSolvers: running colored solvers" << endl;
    {
      CH_TIME("colored solves");
      runSolver(dblColor, domain, dx,   icolor);
    }
  }

  return 0;
}


/// init mpi, init parmparse.  call test.  dump timers. finalize mpi.
int main(int argc, char* argv[])
{
  ///init mpi
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 4586;
  if (argc < 2)
  {
    cerr << " usage " << argv[0] << " <input_file_name> " << endl;
    return status;
  }

  ///init parmparse using the ancient formula
  char* in_file = argv[1];
  ParmParse  pp(argc-2,argv+2,NULL,in_file);

  /// call test
  status = runColoredSolvers();
  if(status != 0)
  {
    cerr << "runColoredSolvers returned error code " << status << endl;
    return status;
  }
  
#ifdef CH_MPI
  dumpmemoryatexit();
  ///dump timers
  CH_TIMER_REPORT();
  ///finalize mpi
  MPI_Finalize();
#endif

  return(status);
}
