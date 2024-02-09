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
#include "DebuggingTools.H"

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
    typedef        DisjointBoxLayout         ch_dbl;
    typedef Proto::DisjointBoxLayout         pr_dbl;
    typedef Proto::Point                     pr_pt;

    ///
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

    ///
    static void
    solveForPhi(Vector<ch_ldf* >   a_phi_ch,
                Vector<ch_ldf* >   a_rhs_ch,
                Vector<ch_dbl  >&  a_amr_grids,
                Vector<ch_dom  >&  a_amr_domains,
                Vector<int     >&  a_ref_ratios,
                Vector<Real    >&  a_amrDx,
                int a_finestLevel)
    {
      CH_TIME("solveForPhi");
      ///the solver declaration has to change 
      shared_ptr<AMRMultiGrid<pr_lbd > > amr_solver_ptr(new AMRMultiGrid<   pr_lbd > ());
      shared_ptr<LinearSolver<pr_lbd>  > bott_solve_ptr(new BiCGStabSolver< pr_lbd > ());

      ParmParse ppSolver("solver");
      Real alpha = 4586.;
      Real beta  = 4586.;
      ppSolver.get("alpha", alpha);
      ppSolver.get("beta" , beta) ;
      string domain_bc;
      ppSolver.get("domain_bc", domain_bc);
      shared_ptr<ch_op_fact_pr> solver_factory_ptr =
        HelmholtzUtilities::getProtoHelmholtzOpFactory(a_amr_domains[0],
                                                       a_ref_ratios,
                                                       a_amr_grids,
                                                       a_amrDx[0],
                                                       domain_bc, alpha, beta);

      PrChUtilities<1>::setupSolver(amr_solver_ptr, bott_solve_ptr, a_amr_grids, a_amr_domains,
                                    a_ref_ratios, a_amrDx, solver_factory_ptr);

      ///define the proto versions of the data
      Vector<pr_lbd*>  phi_pr(a_amr_grids.size(), NULL);
      Vector<pr_lbd*>  rhs_pr(a_amr_grids.size(), NULL);
      for(int ilev = 0; ilev < a_amr_grids.size(); ilev++)
      {
        shared_ptr<pr_dbl> layout_ptr =
          PrChUtilities<1>::getProtoLayout(a_amr_grids[ilev]);
        pr_pt ghost_phi = ProtoCh::getPoint(a_phi_ch[ilev]->ghostVect());
        pr_pt ghost_rhs = ProtoCh::getPoint(a_rhs_ch[ilev]->ghostVect());
        phi_pr[ilev] = new pr_lbd(*layout_ptr, ghost_phi);
        rhs_pr[ilev] = new pr_lbd(*layout_ptr, ghost_rhs);
      }

      //get the rhs onto the device
      for(int ilev = 0; ilev < a_amr_grids.size(); ilev++)
      {
        PrChUtilities<1>::copyToDevice(*rhs_pr[ilev], *a_rhs_ch[ilev]);
      }

      //solve for phi on the device
      bool zeroInitialGuess = true;
      amr_solver_ptr->solve(phi_pr, rhs_pr, a_finestLevel, 0, zeroInitialGuess); 

      //get phi back onto the host
      for(int ilev = 0; ilev < a_amr_grids.size(); ilev++)
      {
        PrChUtilities<1>::copyToHost(*a_phi_ch[ilev], *phi_pr[ilev]);
      }

      //clean up device data
      for(int ilev = 0; ilev < a_amr_grids.size(); ilev++)
      {
        delete(phi_pr[ilev]);
        delete(rhs_pr[ilev]);
      }
    }

  

    ///
    static int runSolver()
    {
      CH_TIME("runSolver");


      // set up grids&
      Vector<ch_dbl> amrGrids;
      Vector<ch_dom> amrDomains;
      Vector<int> refRatios;
      Vector<Real> amrDx;

      int finestLevel;
      PrChUtilities<1>::setupLLCornerGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);
      int numLevels = amrGrids.size();
      if(numLevels != finestLevel+1)
      {
        Chombo::MayDay::Error("runSolver: logic error");
      }
      
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        pout() << "ilev = " << ilev << ", grids = " << endl;
        amrGrids[ilev].print();
      }
      
      Vector<ch_ldf* > phi_ch(numLevels, NULL);
      Vector<ch_ldf* > rhs_ch(numLevels, NULL);
      for (int lev=0; lev< numLevels; lev++)
      {
        const ch_dbl& levelGrids = amrGrids[lev];
        phi_ch[lev] = new ch_ldf(levelGrids, 1, 2*IntVect::Unit);
        rhs_ch[lev] = new ch_ldf(levelGrids, 1,   IntVect::Zero);
      }

      setRHS(rhs_ch, amrDomains, refRatios, amrDx, finestLevel );
      solveForPhi(phi_ch, rhs_ch, amrGrids, amrDomains, refRatios, amrDx, finestLevel );

#ifdef CH_USE_HDF5
      //write the answer to file
      string fname("phi.hdf5");
      Vector<string> varNames(1, string("phi"));

      Real bogusVal = 4586.0;

      WriteAMRHierarchyHDF5(fname,
                            amrGrids,
                            phi_ch,
                            varNames,
                            amrDomains[0].domainBox(),
                            amrDx[0],
                            bogusVal,
                            bogusVal,
                            refRatios,
                            numLevels);

#endif // end if HDF5

      // clean up
      for (int lev=0; lev< phi_ch.size(); lev++)
      {
        delete phi_ch[lev];
        delete rhs_ch[lev];
      }


      return 0;
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
    Chombo::pout() << "Running amr_helmholtz for DIM= " <<  Chombo::SpaceDim << endl;
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
