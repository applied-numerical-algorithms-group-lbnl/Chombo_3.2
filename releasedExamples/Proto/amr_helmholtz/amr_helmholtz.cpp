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

#include "PrChUtilities.H"  //lives in releasedExamples/Proto/common
#include "Proto_Helmholtz_Op.H"

namespace Chombo
{
  class local_test_harness
  {
  public:
    ///no data --- just static functions and typedefs 
    local_test_harness()
    {;}

    ~local_test_harness()
    {;}

    typedef        LevelData<FArrayBox>      ch_ldf;
    typedef        ProblemDomain             ch_dom;
    typedef Proto::ProblemDomain             pr_dom;
    typedef Proto::LevelBoxData<double, 1>   pr_lbd;
    
    static void
    setRHS(Vector<ch_ldf* >   a_rhs,
           Vector<ch_dom  >&  a_amrDomains,
           Vector<int     >&  a_refRatios,
           Vector<Real    >&  a_amrDx,
           int a_finestLevel)
    {
      CH_TIME("setRHS");

      for (int lev=0; lev<=a_finestLevel; lev++)
      {
        ch_ldf& levelRhs = *(a_rhs[lev]);
        const ch_dbl& levelGrids = levelRhs.getBoxes();

        // rhs is cell-centered...
        RealVect ccOffset = 0.5*a_amrDx[lev]*RealVect::Unit;

        DataIterator levelDit = levelGrids.dataIterator();
        for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          FArrayBox& thisRhs = levelRhs[levelDit];
          thisRhs.setVal(1.);
        }
      }
    }


  
    static void
    setupSolver(AMRMultiGrid<ch_ldf > *a_amrSolver,
                LinearSolver<ch_ldf >& a_bottomSolver,
                const Vector<ch_dbl>& a_amrGrids,
                const Vector<ch_dom>& a_amrDomains,
                const Vector<int>& a_refRatios,
                const Vector<Real>& a_amrDx,
                int a_finestLevel)
    {
      CH_TIME("setupSolver");

      ParmParse ppSolver("solver");


      AMRPoissonOpFactory opFactory;

      Real alpha =4586.;
      Real beta = 4586.;
      ppSolver.get("alpha", alpha);
      ppSolver.get("beta" , beta) ;

//  int numLevels = a_finestLevel+1;
//  opFactory.define(a_amrDomains[0],
//                   a_amrGrids,
//                   a_refRatios,
//                   a_amrDx[0],
//                   &ParseBC, alpha, beta);
//
//  AMRLevelOpFactory<ch_ldf >& castFact = (AMRLevelOpFactory<ch_ldf >& ) opFactory;
//
//  a_amrSolver->define(a_amrDomains[0], castFact,
//                      &a_bottomSolver, numLevels);

      // multigrid solver parameters
      int numSmooth, numMG, maxIter;
      Real eps, hang;
      ppSolver.get("num_smooth", numSmooth);
      ppSolver.get("num_mg",     numMG);
      ppSolver.get("max_iterations", maxIter);
      ppSolver.get("tolerance", eps);
      ppSolver.get("hang",      hang);

      Real normThresh = 1.0e-30;
      a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                       numMG, maxIter, eps, hang, normThresh);
      a_amrSolver->m_verbosity = 3;

      // optional parameters
      ppSolver.query("num_pre", a_amrSolver->m_pre);
      ppSolver.query("num_post", a_amrSolver->m_post);
      ppSolver.query("num_bottom", a_amrSolver->m_bottom);
    }

    static int runSolver()
    {
      CH_TIME("runSolver");

      int status = 0;
      ParmParse ppMain("main");

      // set up grids&
      Vector<ch_dbl> amrGrids;
      Vector<ch_dom> amrDomains;
      Vector<int> refRatios;
      Vector<Real> amrDx;
      int finestLevel;

      PrChUtilities<1>::setupLLCornerGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);
      for(int ilev = 0; ilev < amrGrids.size(); ilev++)
      {
        pout() << "ilev = " << ilev << ", grids = " << endl;
        amrGrids[ilev].print();
      }
      AMRMultiGrid<ch_ldf > *amrSolver;
      amrSolver = new AMRMultiGrid<ch_ldf >();
      BiCGStabSolver<ch_ldf > bottomSolver;
      bottomSolver.m_verbosity = 0;
      setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
                  refRatios, amrDx, finestLevel);

      int numLevels = amrGrids.size();
      Vector<ch_ldf* > phi_ch(numLevels, NULL);
      Vector<ch_ldf* > rhs(numLevels, NULL);

      for (int lev=0; lev<=finestLevel; lev++)
      {
        const ch_dbl& levelGrids = amrGrids[lev];
        phi_ch[lev] = new ch_ldf(levelGrids, 1, IntVect::Unit);
        rhs[lev] = new ch_ldf(levelGrids, 1, IntVect::Zero);
      }

      setRHS(rhs, amrDomains, refRatios, amrDx, finestLevel );

      amrSolver->solve(phi_ch, rhs, finestLevel, 0, true); //bool zeroInitialGuess = true;

#ifdef CH_USE_HDF5

      string fname("phi.hdf5");
      Vector<string> varNames(1, string("phi"));

      Real bogusVal = 4586.0;

      WriteAMRHierarchyHDF5(fname,
                            amrGrids,
                            phi,
                            varNames,
                            amrDomains[0].domainBox(),
                            amrDx[0],
                            bogusVal,
                            bogusVal,
                            refRatios,
                            numLevels);

#endif // end if HDF5

      // clean up
      for (int lev=0; lev<phi.size(); lev++)
      {
        delete phi[lev];
        delete rhs[lev];
      }

      delete amrSolver;

      return status;
    }
  };
}
/*****/
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
      std::cerr<< " usage " << argv[0] << " <input_file_name> " << std::endl;
      exit(0);
    }
    char* in_file = argv[1];
    Chombo::ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = Chombo::local_test_harness::runSolver();
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
