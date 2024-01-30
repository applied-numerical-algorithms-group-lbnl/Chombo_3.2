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
#include "PrChUtilities.H"
#include "RelaxSolver.H"
#include "Vector.H"
#include "DebuggingTools.H"
#include "AMRIO.H"



namespace Chombo
{

  ///  just static functions and typedefs
  class viscous_utils
  {
  public:
    //no data to initialize
    viscous_utils()
    {}
    //no pointers to delete
    ~viscous_utils()
    {}
    
    typedef LevelData<FArrayBox> ch_ldf_cell;
    typedef LevelData<FluxBox>   ch_ldf_flux;
    typedef DisjointBoxLayout    ch_dbl; 
    typedef AMRLevelOpFactory<ch_ldf_cell >      ch_amr_levelop;
    typedef LinearSolver<LevelData<FArrayBox> >  ch_base_solver;
    ///
    static void
    poissonSolve(Vector<RefCountedPtr<ch_ldf_cell>       > & a_phi,
                 const Vector<RefCountedPtr<ch_ldf_cell> > & a_rhs,
                 RefCountedPtr< ch_amr_levelop>            & a_opFactory,
                 const Chombo::Vector< ch_dbl >            & a_grids,
                 const VCPoissonParameters                 & a_params)
    {
      string which_solver;
      ParmParse ppSolver("solver");
      // set up multigrid solver parameters
      int numSmooth, numMG, maxIter;
      Real eps, hang;

    
      ppSolver.get("num_smooth", numSmooth);
      ppSolver.get("num_mg",     numMG);
      ppSolver.get("max_iterations", maxIter);
      ppSolver.get("tolerance", eps);
      ppSolver.get("hang",      hang);
      Real normThresh = 1.0e-10;
      string bottom_flag;
      ppSolver.get("bottom_flag", bottom_flag);



      shared_ptr<ch_base_solver> bottom_solver;
      if(bottom_flag == string("relax"))
      {
        int num_bottom;
        ppSolver.get("num_bottom", num_bottom);
        RelaxSolver<ch_ldf_cell >* relax_raw_ptr = new RelaxSolver<ch_ldf_cell >();
        relax_raw_ptr->m_imax = num_bottom;
        ch_base_solver* base_raw_ptr =  static_cast<ch_base_solver*>(relax_raw_ptr);
        bottom_solver = shared_ptr<ch_base_solver>(base_raw_ptr);
      }
      else if(bottom_flag == string("bicgstab"))
      {
        BiCGStabSolver<ch_ldf_cell >* bicgstab_raw_ptr = new BiCGStabSolver<ch_ldf_cell >();
        ch_base_solver * base_raw_ptr = static_cast<ch_base_solver*>(bicgstab_raw_ptr);
        bottom_solver = shared_ptr<ch_base_solver>(base_raw_ptr);
      }
      else
      {
        MayDay::Error("poissonSolve: bogus bottom_flag");
      }
    
      AMRMultiGrid<  ch_ldf_cell > amr_mg_solver;
      amr_mg_solver.define(a_params.m_coarsestDomain,
                           (*a_opFactory),
                           &(*bottom_solver),
                           a_params.m_numLevels);
    
      amr_mg_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                                        numMG, maxIter, eps, hang, normThresh);

      amr_mg_solver.m_verbosity = 3;
      Vector<ch_ldf_cell*> raw_phi = PrChUtilities<1>::getRawVectorCh(a_phi);
      Vector<ch_ldf_cell*> raw_rhs = PrChUtilities<1>::getRawVectorCh(a_rhs);
      
      amr_mg_solver.solve(raw_phi, raw_rhs, a_params.m_maxLevel, 0, true);

    }//end function poissonSolve
  
    ///call define but do not fill data
    static void
    defineData(Vector< RefCountedPtr<ch_ldf_cell> >& a_phi,
               Vector< RefCountedPtr<ch_ldf_cell> >& a_rhs,
               Vector< RefCountedPtr<ch_ldf_cell> >& a_aco,
               Vector< RefCountedPtr<ch_ldf_flux> >& a_eta,
               Vector< RefCountedPtr<ch_ldf_flux> >& a_lam,
               const Vector< ch_dbl >              & a_grids)
    {
      a_phi.resize(a_grids.size());
      a_rhs.resize(a_grids.size());
      a_aco.resize(a_grids.size());
      a_eta.resize(a_grids.size());
      a_lam.resize(a_grids.size());
    
      for(int ilev = 0; ilev < a_grids.size(); ilev++)
      {
        a_phi[ilev] = RefCountedPtr<ch_ldf_cell>( new ch_ldf_cell(a_grids[ilev],DIM, IntVect::Unit));
        a_rhs[ilev] = RefCountedPtr<ch_ldf_cell>( new ch_ldf_cell(a_grids[ilev],DIM, IntVect::Zero));
        a_aco[ilev] = RefCountedPtr<ch_ldf_cell>( new ch_ldf_cell(a_grids[ilev], 1 , IntVect::Zero));
        a_eta[ilev] = RefCountedPtr<ch_ldf_flux>( new ch_ldf_flux(a_grids[ilev], 1 , IntVect::Zero) );
        a_lam[ilev] = RefCountedPtr<ch_ldf_flux>( new ch_ldf_flux(a_grids[ilev], 1 , IntVect::Zero));
      }  ///end loop over levels
    }    //end function defineData

    /***************/
    static void
    outputData(const Vector< RefCountedPtr<ch_ldf_cell > >&   a_phi,
               const Vector< RefCountedPtr<ch_ldf_cell > >&   a_rhs,
               const Vector< ch_dbl >                     &   a_grids,
               const VCPoissonParameters                  &   a_params)
    {
#ifdef CH_USE_HDF5
      /// a good UI is a joy forever
      Vector<ch_ldf_cell*> raw_phi = PrChUtilities<1>::getRawVectorCh(a_phi);
      Vector<ch_ldf_cell*> raw_rhs = PrChUtilities<1>::getRawVectorCh(a_rhs);
      string         phiFileName("phi.hdf5");
      string         rhsFileName("rhs.hdf5");
    
      Vector<string> phiVarNames(DIM);
      Vector<string> rhsVarNames(DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        phiVarNames[idir] = string("vel_comp_") + std::to_string(idir);
        rhsVarNames[idir] = string("rhs_comp_") + std::to_string(idir);
      }
      
      Real fakeTime = 1.0;
      Real fakeDt = 1.0;
      WriteAMRHierarchyHDF5(phiFileName,  a_grids,
                            raw_phi, phiVarNames,
                            a_params.m_coarsestDomain.domainBox(),
                            a_params.m_coarsestDx,
                            fakeDt, fakeTime,
                            a_params.m_refRatio,
                            a_params.m_numLevels);
      WriteAMRHierarchyHDF5(rhsFileName,  a_grids,
                            raw_rhs, rhsVarNames,
                            a_params.m_coarsestDomain.domainBox(),
                            a_params.m_coarsestDx,
                            fakeDt, fakeTime,
                            a_params.m_refRatio,
                            a_params.m_numLevels);

#endif
    } //end function ouputData

    ///
    static void
    setACoef(Vector<RefCountedPtr< ch_ldf_cell > >  & a_aCoefVec,
             const VCPoissonParameters              & a_params)
    {

      Chombo::ParmParse pp("viscous_op");
      double aco_val;
      pp.get("acoef_value", aco_val);
      for(int ilev = 0; ilev < a_params.m_numLevels; ilev++)
      {
        DataIterator dit = a_aCoefVec[ilev]->dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& aCoef = (*a_aCoefVec[ilev])[dit];
          aCoef.setVal(aco_val);
        }   // end loop over grids
      }     // end loop over levels
    }       //end function setACoef

    /// 
    static void
    setEtaAndLambda(Vector< RefCountedPtr< ch_ldf_flux > > & a_eta,
                    Vector< RefCountedPtr< ch_ldf_flux > > & a_lam,
                    const VCPoissonParameters              & a_params)
    {
      Real etaVal, lamVal;
      Chombo::ParmParse pp("viscous_op");
      pp.get("eta_value"   , etaVal);
      pp.get("lambda_value", lamVal);
      for(int ilev = 0; ilev < a_params.m_numLevels; ilev++)
      {
        DataIterator dit = a_eta[ilev]->dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& thisEta = (*a_eta[ilev])[dit];
          FluxBox& thisLam = (*a_lam[ilev])[dit];
          for (int face_dir=0; face_dir < SpaceDim; face_dir++)
          {
            FArrayBox& etafab = thisEta[face_dir];
            FArrayBox& lamfab = thisLam[face_dir];
            etafab.setVal(etaVal);
            lamfab.setVal(lamVal);
          }   // end loop over face direction
        }     // end loop over boxes
      }       // end loop over levels
    }//.end function setBCoef


    /**
       crude blob code
       if(rad^2 < a_rad^2) rhs= a_value;
       else rhs = 0;
     **/
    static double
    rhs_blob(const IntVect& a_loc, double a_dx, double a_cent_loc, double a_rad_sq, double a_value)
    {
      double dist_sq = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        double x_loc       = a_dx*(double(a_loc[idir]) + 0.5);
        double x_minus_xno = x_loc - a_cent_loc;
        dist_sq += (x_minus_xno)*(x_minus_xno);
      }
      double retval = 4586.;
      if(dist_sq < a_rad_sq)
      {
        retval = a_value;
      }
      else
      {
        retval = 0.;
      }
      return retval;     
    }
             
    /********/
    static void
    setRHS(Vector< RefCountedPtr<ch_ldf_cell > >   a_rhs,
           Vector<Real>                            a_amrDx)
    {
      
      double cent_loc = 0.5;
      double rad_sq   = 0.001;
      double rhs_val  = 0.2;
      for(int ilev = 0; ilev < a_rhs.size(); ilev++)
      {
        auto& rhs_lev = *(a_rhs[ilev]);
        double dx = a_amrDx[ilev];
        for (DataIterator dit = rhs_lev.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& rhs_fab = rhs_lev[dit()];
          Box        rhs_box = rhs_fab.box();
          for(BoxIterator bit(rhs_box); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            double rhsval = rhs_blob(iv, dx, cent_loc, rad_sq, rhs_val);
            for(int ivar = 0; ivar < rhs_fab.nComp(); ivar++)
            {
              rhs_fab(iv, ivar) = rhsval;
            } //end loop over variables
          }   //end loop over cells
        }     //end loop over boxes    
      }       //end loop over levels
      return;
    }///end function setRHS
    static  int
    local_test_harness()
    {
      
      Vector<RefCountedPtr<ch_ldf_cell> > phi;
      Vector<RefCountedPtr<ch_ldf_cell> > rhs;
      Vector<RefCountedPtr<ch_ldf_cell> > aCoef;
      Vector<RefCountedPtr<ch_ldf_flux> > eta;
      Vector<RefCountedPtr<ch_ldf_flux> > lambda;

      Vector<DisjointBoxLayout> amrGrids;
      Vector<ProblemDomain> amrDomains;
      Vector<int> refRatios;
      Vector<Real> amrDx;
      int finestLevel;

      RealVect tag_location;
      ParmParse pp("test_harness");
      Vector<Real> tag_location_vec(DIM);
      pp.getarr("grid_tag_location", tag_location_vec, 0, DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        tag_location[idir] = tag_location_vec[idir];
      }
      PrChUtilities<1>::setupGridsWithTagLocation(amrGrids, amrDomains, refRatios, amrDx, finestLevel, tag_location);
      for(int ilev = 0; ilev < amrGrids.size(); ilev++)
      {
        pout() << "ilev = " << ilev << ", grids = " << endl;
        amrGrids[ilev].print();
      }

      VCPoissonParameters param = 
        getVCParameters(finestLevel, amrDx[0], amrDomains[0], refRatios);
    

      defineData(phi, rhs, aCoef, eta, lambda,  amrGrids);
      setRHS(rhs, amrDx);
      setACoef(aCoef, param);
      setEtaAndLambda(eta, lambda, param);

      RefCountedPtr<AMRLevelOpFactory<ch_ldf_cell >  >  
        op_factory = VCLocalFunctions::getViscousTensorFactory(aCoef, eta, lambda, amrGrids, param);

      poissonSolve(phi, rhs, op_factory, amrGrids, param);
      outputData(  phi, rhs,             amrGrids, param);

      return 0;
    }
  }; //end class viscous_utils
} //end namespace Chombo



int main(int argc, char* argv[])
{
  /**
     We are not using petsc so we can avoid the petsc/mpi init/finalize dance.
  **/
#ifdef CH_MPI   
  MPI_Init(&argc, &argv);
#endif
  if (argc < 2)
  {
    cerr << " usage " << argv[0] << " <input_file_name> " << endl;
    exit(0);
  }

  char* inFile = argv[1];
  Chombo::ParmParse pp(argc-2,argv+2,NULL,inFile);

  Chombo::pout() << "Running old_viscous_tensor for DIM= " <<  Chombo::SpaceDim << endl;
  int status = Chombo::viscous_utils::local_test_harness();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
