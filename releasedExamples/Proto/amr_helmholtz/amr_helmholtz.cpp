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
    typedef        AMRLevelOpFactory<pr_lbd> ch_op_fact_pr;   
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
    solveForPhi(Vector<ch_ldf* >   a_rhs_ch,
                Vector<ch_ldf* >   a_phi_ch,
                Vector<ch_dom  >&  a_amrDomains,
                Vector<int     >&  a_refRatios,
                Vector<Real    >&  a_amrDx,
                int a_finestLevel)
    {
      CH_TIME("solveForPhi");
      ///the solver declaration has to change 
      shared_ptr<AMRMultiGrid<pr_lbd > > amr_solver_ptr(new AMRMultiGrid<   pr_lbd > ());
      shared_ptr<LinearSolver<pr_lbd>  > bott_solve_ptr(new BiCGStabSolver< pr_lbd > ());
      bottomSolver.m_verbosity = 0;
      ParmParse ppSolver("solver");
      Real alpha = 4586.;
      Real beta  = 4586.;
      ppSolver.get("alpha", alpha);
      ppSolver.get("beta" , beta) ;

      shared_ptr<ch_op_fact_pr> solver_factory_ptr =
        PrChUtilities<1>::getProtoHelmholtzOpFactory(a_amrDomains[0], a_refRatios, a_amrDx[0], domainBC, alpha, beta);

      PrChUtilities<1>::setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
                                    refRatios, amrDx, finestLevel, solver_factory_ptr);

      int numLevels = amrGrids.size();

      PrChUtilities<1>::copyToDevice(rhs_pr, a_rhs_ch);
      
      amrSolver->solve(phi_pr, rhs_pr, finestLevel, 0, true); //bool zeroInitialGuess = true;
      
      PrChUtilities<1>::copyToHost  (a_phi_ch, phi_pr);
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

      Vector<ch_ldf* > phi_ch(numLevels, NULL);
      Vector<ch_ldf* > rhs_ch(numLevels, NULL);

      for (int lev=0; lev<=finestLevel; lev++)
      {
        const ch_dbl& levelGrids = amrGrids[lev];
        phi_ch[lev] = new ch_ldf(levelGrids, 1, IntVect::Unit);
        rhs_ch[lev] = new ch_ldf(levelGrids, 1, IntVect::Zero);
      }

      setRHS(rhs_ch, amrDomains, refRatios, amrDx, finestLevel );

      solveForPhi(phi_ch, rhs_ch, amrDomains, refRatios, amrDx, finestLevel );

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
