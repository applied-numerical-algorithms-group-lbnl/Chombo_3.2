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
/****************/
PROTO_KERNEL_START 

unsigned int setPhiQuadF(int                     a_pt[DIM],
                         Proto::Var<double, 1>   a_phi,
                         double                  a_aco,
                         double                  a_bco,
                         double                  a_cco,
                         double                  a_dx,
                         int                     a_direction)
{
  double xval = a_dx*(double(a_pt[a_direction]) + 0.5);
  double phival = a_aco*xval*xval + a_bco*xval + a_cco;
  a_phi(0) = phival;
  return 0;

}
PROTO_KERNEL_END(setPhiQuadF, setPhiQuad) 

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

    typedef        LevelData<FArrayBox>      ch_ldf_cell;
    typedef        LevelData<FluxBox>        ch_ldf_flux;
    typedef        ProblemDomain             ch_dom;
    typedef Proto::ProblemDomain             pr_dom;
    typedef Proto::LevelBoxData<double, 1>   pr_lbd;
    typedef        AMRLevelOpFactory<pr_lbd> ch_op_fact_pr;
    typedef        DisjointBoxLayout         ch_dbl;
    typedef Proto::DisjointBoxLayout         pr_dbl;
    typedef Proto::Point                     pr_pt;
    typedef Proto::Box                       pr_box;
    typedef Proto::BoxData<double, 1>        pr_box_data;

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
    defineAndSetCoefficients(RefCountedPtr< ch_ldf_cell> & a_acoef,
                             RefCountedPtr< ch_ldf_flux> & a_bcoef,
                             ch_dbl     & a_grids,
                             double a_aval, double a_bval)
    {
      CH_TIME("defineAndSetCoefficients");

      a_acoef = RefCountedPtr<ch_ldf_cell>(new ch_ldf_cell(a_grids, 1, IntVect::Zero));
      a_bcoef = RefCountedPtr<ch_ldf_flux>(new ch_ldf_flux(a_grids, 1, IntVect::Zero));
      DataIterator dit = a_grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        (*a_acoef)[dit[ibox]].setVal(a_aval);
        (*a_bcoef)[dit[ibox]].setVal(a_bval);
      }
    } //end function defineAndSetCoefficients

    ///
    static int
    testFluxRegister(Vector<ch_dbl  >&  a_amr_grids,
                     Vector<ch_dom  >&  a_amr_domains,
                     Vector<int     >&  a_ref_ratios,
                     Vector<Real    >&  a_amrDx)
    {
      typedef PrCh_Tools::FluxRegister<1> prch_flux_reg;
      typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data;
      typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >   pr_xfac_data;
      typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >   pr_yfac_data;
      typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >   pr_zfac_data;
      ///define the proto versions of the data
      int numLevels = a_amr_domains.size();
      vector<shared_ptr< pr_dbl> >    grids_pr(numLevels);
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        grids_pr[ilev] = PrChUtilities<1>::getProtoLayout(a_amr_grids[ilev]);
      }
      
      /**
         Outline:
         1. Define a bunch of fluxes at each level. 
         2. Set them to constant != 0, 
         3. Do the reflux dance and check to see if the difference is zero
      **/
      double flux_val = 4.586;
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        if(ilev < (numLevels-1))
        {
          pr_dbl& grids_coar = (*grids_pr[ilev  ]);
          pr_dbl& grids_fine = (*grids_pr[ilev+1]);
          int ref_ratio = a_ref_ratios[ilev];
          double dx_coar = a_amrDx[ilev];
          shared_ptr<pr_cell_data>   solun_incr(new pr_cell_data(grids_coar, pr_pt::Zeros()));
          shared_ptr<pr_xfac_data>   xflux_coar(new pr_xfac_data(grids_coar, pr_pt::Zeros()));
          shared_ptr<pr_xfac_data>   xflux_fine(new pr_xfac_data(grids_fine, pr_pt::Zeros()));
          shared_ptr<pr_yfac_data>   yflux_coar(new pr_yfac_data(grids_coar, pr_pt::Zeros()));
          shared_ptr<pr_yfac_data>   yflux_fine(new pr_yfac_data(grids_fine, pr_pt::Zeros()));
#if DIM==3      
          shared_ptr<pr_zfac_data>   zflux_coar(new pr_zfac_data(grids_coar, pr_pt::Zeros()));
          shared_ptr<pr_zfac_data>   zflux_fine(new pr_zfac_data(grids_fine, pr_pt::Zeros()));
#else
          shared_ptr<pr_zfac_data>   zflux_coar;
          shared_ptr<pr_zfac_data>   zflux_fine;
#endif
          solun_incr->setVal(0.);
          xflux_coar->setVal(flux_val);
          xflux_fine->setVal(flux_val);
          yflux_coar->setVal(flux_val);
          yflux_fine->setVal(flux_val);
#if DIM==3         
          zflux_coar->setVal(flux_val);
          zflux_fine->setVal(flux_val);
#endif
          prch_flux_reg flux_register(grids_coar, grids_fine, ref_ratio);
          bool verbose = true;
          flux_register.setCoarFlux(  xflux_coar, yflux_coar, zflux_coar, verbose);
          flux_register.setFineFlux(  xflux_fine, yflux_fine, zflux_fine, verbose);
          flux_register.reflux(*solun_incr, dx_coar, verbose);
          double remainderNorm = solun_incr->absMax();
          double tol = 1.0e-10;
          if(remainderNorm > tol)
          {
            pout() << "testFluxRegister: ||remainder|| = " << remainderNorm << ", tol = "<< tol << endl;
            return -4586;
          }
        } //end if ilev < numLevels-1
      }//   end loop over levels
      return 0;
    }//end function testFluxRegister

    ///
    /**
       Sets initial guess to phi to a polynomial.
       Does not make much sense outside debugging.
     **/
    static void
    setInitialPhi(Vector<pr_lbd*>   a_phi_pr,   Vector<double> a_amrDx)
    {
      ///phi = a x^2 + b x + c
      ///this is a debugging thing so I am just going to set them
      double acoef = 0., bcoef = 1., ccoef = 4.586; int direction = 0;
      for(int ilev = 0; ilev < a_phi_pr.size(); ilev++)
      {
        auto dit   = a_phi_pr[ilev]->begin();
        auto grids = a_phi_pr[ilev]->layout();
        double dx = a_amrDx[ilev];
        for(int ibox = 0; ibox < dit.localSize(); ibox++)
        {
          pr_box_data& phi = (*a_phi_pr[ilev])[dit[ibox]];
          pr_box valid     =             grids[dit[ibox]];
          Proto::forallInPlace_i(setPhiQuad, valid, phi,
                                 acoef, bcoef, ccoef, dx, direction);
        }
      }
      return;
    }
    ///
    static int
    solveForPhi(Vector<ch_ldf_cell* >   a_phi_ch,
                Vector<ch_ldf_cell* >   a_rhs_ch,
                Vector<ch_dbl  >&  a_amr_grids,
                Vector<ch_dom  >&  a_amr_domains,
                Vector<int     >&  a_ref_ratios,
                Vector<Real    >&  a_amrDx)
    {
      CH_TIME("solveForPhi");
      int numLevels = a_amr_grids.size();

      ///opening act slipped in here because it was convenient
      int ifluxReg = testFluxRegister(a_amr_grids, a_amr_domains, a_ref_ratios, a_amrDx);
      if(ifluxReg != 0)
      {
        pout() << "solveForPhi: problem in testFluxRegister" << endl;
        return ifluxReg;
      }
      else
      {
        pout() << "solveForPhi: testFluxRegister passed" << endl;
      }
      ///the solver declaration has to change because amrmultigrid is templated on data type
      shared_ptr<AMRMultiGrid<pr_lbd > > amr_solver_ptr(new AMRMultiGrid<   pr_lbd > ());
      shared_ptr<LinearSolver<pr_lbd>  > bott_solve_ptr(new BiCGStabSolver< pr_lbd > ());
      Vector<RefCountedPtr< ch_ldf_cell > > acoef(numLevels);
      Vector<RefCountedPtr< ch_ldf_flux > > bcoef(numLevels);

      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        double aval = 1; double bval = 1;
        defineAndSetCoefficients(acoef[ilev], bcoef[ilev],
                                 a_amr_grids[ilev], aval, bval);
      }
      
      ParmParse ppSolver("solver");
      Real alpha = 4586.;
      Real beta  = 4586.;
      ppSolver.get("alpha", alpha);
      ppSolver.get("beta" , beta) ;
      string domain_bc;
      ppSolver.get("domain_bc", domain_bc);
      shared_ptr<ch_op_fact_pr> solver_factory_ptr =
        HelmholtzUtilities::getProtoConductivityFactory(a_amr_domains[0],
                                                        a_ref_ratios,
                                                        a_amr_grids,
                                                        a_amrDx[0],
                                                        acoef, bcoef, 
                                                        domain_bc, alpha, beta);

      PrChUtilities<1>::setupSolver(amr_solver_ptr, bott_solve_ptr, a_amr_grids,
                                    a_amr_domains, a_ref_ratios, a_amrDx,
                                    solver_factory_ptr);

      ///define the proto versions of the data
      Vector<pr_lbd*>   phi_pr(numLevels, NULL);
      Vector<pr_lbd*>   rhs_pr(numLevels, NULL);
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        shared_ptr<pr_dbl> layout_ptr =
          PrChUtilities<1>::getProtoLayout(a_amr_grids[ilev]);
        pr_pt ghost_phi = ProtoCh::getPoint(a_phi_ch[ilev]->ghostVect());
        pr_pt ghost_rhs = ProtoCh::getPoint(a_rhs_ch[ilev]->ghostVect());
        phi_pr[ilev] = new pr_lbd(*layout_ptr, ghost_phi);
        rhs_pr[ilev] = new pr_lbd(*layout_ptr, ghost_rhs);
      }

      //get the rhs onto the device
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        PrChUtilities<1>::copyToDevice(*rhs_pr[ilev], *a_rhs_ch[ilev]);
      }

      //solve for phi on the device
#if 0
      ///standard code where initial guess for phi is zero
      bool zeroInitialGuess = true;  int lbase = 0; int lmax = numLevels-1; bool print = true;
      amr_solver_ptr->solve(phi_pr, rhs_pr, lmax, lbase, zeroInitialGuess, print); 
#else
      ///code where I set phi to something I know and so I can look at initial residuals
      setInitialPhi(phi_pr, a_amrDx);
      bool zeroInitialGuess = false;  int lbase = 0; int lmax = numLevels-1; bool print = true;
      amr_solver_ptr->solve(phi_pr, rhs_pr, lmax, lbase, zeroInitialGuess, print);
#endif      

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
      return 0;
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
      IntVect ivghost =  4.*IntVect::Unit;
      for (int lev=0; lev< numLevels; lev++)
      {
        const ch_dbl& levelGrids = amrGrids[lev];
        phi_ch[lev] = new ch_ldf_cell(levelGrids, 1, ivghost);
        rhs_ch[lev] = new ch_ldf_cell(levelGrids, 1, ivghost);
      }

      setRHS(rhs_ch, amrDomains, refRatios, amrDx, finestLevel );
      int iret = solveForPhi(phi_ch, rhs_ch, amrGrids, amrDomains, refRatios, amrDx);
      if(iret != 0)
      {
        pout() << "runSolver: solveForPhi returned "  << iret << endl;
        return iret;
      }

#ifdef CH_USE_HDF5
#if 0      
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

#endif // end if
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
  Chombo::pout() << "Running amr_conductivity for DIM= " <<  Chombo::SpaceDim << endl;
  int status = Chombo::local_test_harness::runSolver();
  
#ifdef CH_MPI
  Chombo::dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
