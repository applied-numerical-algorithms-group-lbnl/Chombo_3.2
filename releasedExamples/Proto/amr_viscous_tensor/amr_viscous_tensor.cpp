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
#include "Proto_Conductivity_Op.H"
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

    typedef        LevelData<FArrayBox>          ch_ldf_cell;
    typedef        LevelData<FluxBox>            ch_ldf_flux;
    typedef        ProblemDomain                 ch_dom;
    typedef Proto::ProblemDomain                 pr_dom;
    typedef Proto::LevelBoxData<double, 1>       pr_lbd_sca;
    typedef Proto::LevelBoxData<double, DIM>     pr_lbd_vec;
    typedef        AMRLevelOpFactory<pr_lbd_vec> ch_op_fact_pr;
    typedef        DisjointBoxLayout             ch_dbl;
    typedef Proto::DisjointBoxLayout             pr_dbl;
    typedef Proto::Point                         pr_pt;

    ///
    static void
    setRHS(Vector<ch_ldf_cell* >   a_rhs,
           Vector<ch_dom  >&  a_amrDomains,
           Vector<int     >&  a_refRatios,
           Vector<Real    >&  a_amrDx,
           int a_finestLevel)
    {
      CH_TIME("setRHS");

      for (int lev=0; lev<=a_finestLevel; lev++)
      {
        ch_ldf_cell& levelRhs = *(a_rhs[lev]);
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
    }//end function setRHS

    ///
    static void
    defineAndSetCoefficients(RefCountedPtr< ch_ldf_cell> & a_aco,
                             RefCountedPtr< ch_ldf_flux> & a_eta,
                             RefCountedPtr< ch_ldf_flux> & a_lam,
                             ch_dbl                      & a_grids,
                             double a_aco_val, double a_eta_val, double a_lam_val)
    {
      CH_TIME("defineAndSetCoefficients");

      a_aco = RefCountedPtr<ch_ldf_cell>(new ch_ldf_cell(a_grids, 1, IntVect::Zero));
      a_eta = RefCountedPtr<ch_ldf_flux>(new ch_ldf_flux(a_grids, 1, IntVect::Zero));
      a_lam = RefCountedPtr<ch_ldf_flux>(new ch_ldf_flux(a_grids, 1, IntVect::Zero));
      DataIterator dit = a_grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        (*a_aco)[dit[ibox]].setVal(a_aco_val);
        (*a_eta)[dit[ibox]].setVal(a_eta_val);
        (*a_lam)[dit[ibox]].setVal(a_lam_val);
      }
    } //end function defineAndSetCoefficients
    
    ///
    static void
    solveForPhi(Vector<ch_ldf_cell* >   a_phi_ch,
                Vector<ch_ldf_cell* >   a_rhs_ch,
                Vector<ch_dbl  >&  a_amr_grids,
                Vector<ch_dom  >&  a_amr_domains,
                Vector<int     >&  a_ref_ratios,
                Vector<Real    >&  a_amrDx)
    {
      CH_TIME("solveForPhi");
      int numLevels = a_amr_grids.size();
      
      Chombo::ParmParse pp("viscous_op");
      double aco_val = 1; double eta_val = 1; double lam_val = 1.;
      pp.get("acoef_value" , aco_val);
      pp.get("eta_value"   , eta_val);
      pp.get("lambda_value", lam_val);
      ///the solver declaration has to change because amrmultigrid is templated on data type
      shared_ptr<AMRMultiGrid<pr_lbd_vec > > amr_solver_ptr(new AMRMultiGrid<   pr_lbd_vec > ());
      shared_ptr<LinearSolver<pr_lbd_vec>  > bott_solve_ptr(new BiCGStabSolver< pr_lbd_vec > ());
      Vector<RefCountedPtr< ch_ldf_cell > > aco(numLevels);
      Vector<RefCountedPtr< ch_ldf_flux > > eta(numLevels);
      Vector<RefCountedPtr< ch_ldf_flux > > lam(numLevels);

      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        
        defineAndSetCoefficients(aco[ilev], eta[ilev], lam[ilev],
                                 a_amr_grids[ilev], aco_val, eta_val, lam_val);
      }
      
      ParmParse ppSolver("solver");
      Real alpha = 4586.;
      Real beta  = 4586.;
      ppSolver.get("alpha", alpha);
      ppSolver.get("beta" , beta) ;
      string domain_bc;
      ppSolver.get("domain_bc", domain_bc);
      ch_dom coarsest_dom = a_amr_domains[0];
      shared_ptr<ch_op_fact_pr> solver_factory_ptr =
        ViscousTensorUtilities::getProtoViscousTensorOpFactory(coarsest_dom, //a_amr_domains[0],
                                                               a_ref_ratios,
                                                               a_amr_grids,
                                                               a_amrDx[0],
                                                               aco, eta, lam,  
                                                               domain_bc, alpha, beta);

      PrChUtilities<DIM>::setupSolver(amr_solver_ptr, bott_solve_ptr, a_amr_grids,
                                      a_amr_domains, a_ref_ratios, a_amrDx,
                                      solver_factory_ptr);

      ///define the proto versions of the data
      Vector<pr_lbd_vec*>   phi_pr(numLevels, NULL);
      Vector<pr_lbd_vec*>   rhs_pr(numLevels, NULL);
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        shared_ptr<pr_dbl> layout_ptr =
          PrChUtilities<1>::getProtoLayout(a_amr_grids[ilev]);
        pr_pt ghost_phi = ProtoCh::getPoint(a_phi_ch[ilev]->ghostVect());
        pr_pt ghost_rhs = ProtoCh::getPoint(a_rhs_ch[ilev]->ghostVect());
        phi_pr[ilev] = new pr_lbd_vec(*layout_ptr, ghost_phi);
        rhs_pr[ilev] = new pr_lbd_vec(*layout_ptr, ghost_rhs);
      }

      //get the rhs onto the device
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        PrChUtilities<DIM>::copyToDevice(*rhs_pr[ilev], *a_rhs_ch[ilev]);
      }

      //solve for phi on the device
      bool zeroInitialGuess = true;  int lbase = 0; int lmax = numLevels-1; bool print = true;
      amr_solver_ptr->solve(phi_pr, rhs_pr, lmax, lbase, zeroInitialGuess, print); 

      //get phi back onto the host
      for(int ilev = 0; ilev < a_amr_grids.size(); ilev++)
      {
        PrChUtilities<DIM>::copyToHost(*a_phi_ch[ilev], *phi_pr[ilev]);
      }

      //clean up device data
      for(int ilev = 0; ilev < a_amr_grids.size(); ilev++)
      {
        delete(phi_pr[ilev]);
        delete(rhs_pr[ilev]);
      }
    }//end function solveForPhi
  

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
      
      Vector<ch_ldf_cell* > phi_ch(numLevels, NULL);
      Vector<ch_ldf_cell* > rhs_ch(numLevels, NULL);

      for (int lev=0; lev< numLevels; lev++)
      {
        const ch_dbl& levelGrids = amrGrids[lev];
        phi_ch[lev] = new ch_ldf_cell(levelGrids, DIM, IntVect::Unit);
        rhs_ch[lev] = new ch_ldf_cell(levelGrids, DIM, IntVect::Zero);
      }

      setRHS(rhs_ch, amrDomains, refRatios, amrDx, finestLevel );
      solveForPhi(phi_ch, rhs_ch, amrGrids, amrDomains, refRatios, amrDx);

#ifdef CH_USE_HDF5

      string fname("phi.hdf5");
      Vector<string> varNames(DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        varNames[idir] = string("phi_comp_") + std::to_string(idir);
      }
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
    } // end function runSolver
  };  // end class local_test_harness
}     // end namespace Chombo
/*****/
int main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  if (argc < 2)
  {
    std::cerr<< " usage " << argv[0] << " <input_file_name> " << std::endl;
    exit(0);
  }
  char* in_file = argv[1];
  Chombo::ParmParse  pp(argc-2,argv+2,NULL,in_file);
  int status = Chombo::local_test_harness::runSolver();
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}