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
       This is from AMRPoisson/execViscousTensor example (the Gaussian blob case).
     **/
    static double
    rhs_blob(const IntVect& a_cell, double a_dx)
    {
      //old code:
      //       int numGaussians = 3;
      //        Vector<RealVect> center(numGaussians,RealVect::Zero);
      //        Vector<Real> scale(numGaussians, 1.0);
      //        Vector<Real> strength(numGaussians, 1.0);
      //doing this with POD for port.
#define BLOB_NUM_GAUSS 3
      double    cente[BLOB_NUM_GAUSS][DIM];
      double    scale[BLOB_NUM_GAUSS];
      double    stren[BLOB_NUM_GAUSS];
      ///Gaussian 0
      stren[0]    = 1.0;
      scale[0]    = 1.0e-2;
      cente[0][0] = 0.25;
      cente[0][1] = 0.25;
#if DIM==3      
      cente[0][2] = 0.25;
#endif      
      ///Gaussian 1
      stren[1] = 3.0;
      scale[1] = 1.0e-2;
      cente[1][0] = 0.50;
      cente[1][1] = 0.75;
#if DIM==3      
      cente[1][2] = 0.75;
#endif
      ///Gaussian 2
      stren[2] = 2.0;
      scale[2] = 1.0e-2;
      cente[2][0] = 0.75; 
      cente[2][1] = 0.50; 
      cente[2][2] = 0.50; 

      double x_loc[DIM];
      for(int idir = 0; idir < DIM; idir++)
      {
        x_loc[idir] = a_dx*(double(a_cell[idir]) + 0.5);
      }

      double retval = 4586;
      retval = 0; // this is additive so 0 is important
      for(int igauss = 0; igauss < BLOB_NUM_GAUSS; igauss++)
      {
        double dist[DIM];
        double rad_sq = 0;
        for(int idir = 0; idir < DIM; idir++)
        {
          dist[idir] = x_loc[idir] - cente[igauss][idir];
          rad_sq += dist[idir]*dist[idir];
        }

        double val = stren[igauss]*exp(-rad_sq/scale[igauss]);
        retval += val;
      }   // end loop over Gaussians
      return retval;     
#undef BLOB_NUM_GAUSS 
    }
             
    /********/
    static void
    setRHS(Vector< RefCountedPtr<ch_ldf_cell > >   a_rhs,
           Vector<Real>                            a_amrDx)
    {
      
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
            double rhsval = rhs_blob(iv, dx);
            for(int ivar = 0; ivar < rhs_fab.nComp(); ivar++)
            {
              rhs_fab(iv, ivar) = rhsval;
            } //end loop over variables
          }   //end loop over cells
        }     //end loop over boxes    
      }       //end loop over levels
      return;
    }///end function setRHS

    ///
    /**
       From releasedExamples/AMRPoisson/execViscousTensor
    **/
    static void
    setupGrids(Vector<DisjointBoxLayout>& a_amrGrids,
               Vector<ProblemDomain>& a_amrDomains,
               Vector<int>& a_refRatios,
               Vector<Real>& a_amrDx,
               int& a_finestLevel)
    {
      a_finestLevel = 0;
      ParmParse ppGrids("setupGrids");

      // get grid generation parameters
      int maxLevel, maxBoxSize, blockFactor;
      Real fillRatio;

      ppGrids.get("max_level", maxLevel);

      ppGrids.get("max_box_size",maxBoxSize);

      ppGrids.get("block_factor", blockFactor);

      ppGrids.get("fillRatio", fillRatio);

      // note that there only need to be numLevels-1 refinement ratios
      a_refRatios.resize(maxLevel);
      ppGrids.getarr("ref_ratio", a_refRatios, 0, maxLevel);

      Vector<int>  is_periodic_int;
      bool is_periodic[SpaceDim];
      ppGrids.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
      for (int dir=0; dir<SpaceDim; dir++)
      {
        is_periodic[dir] = (is_periodic_int[dir] == 1);
      }

      IntVect numCells;
      Vector<int> incells(SpaceDim);
      ppGrids.getarr("num_cells", incells, 0, SpaceDim);
      numCells = IntVect(D_DECL6(incells[0],incells[1],incells[2],
                                 incells[3],incells[4],incells[5]) );

      RealVect domainSize = RealVect::Unit;
      if (ppGrids.contains("domain_size"))
      {
        Vector<Real> insize(SpaceDim);
        ppGrids.getarr("domain_size", insize, 0, SpaceDim);
        domainSize = RealVect(D_DECL6(insize[0],insize[1],insize[2],
                                      insize[3],insize[4],insize[5]) );
      }

      // resize dataholders
      int maxNumLevels = maxLevel +1;
      a_amrGrids.resize(maxNumLevels);
      a_amrDomains.resize(maxNumLevels);
      a_amrDx.resize(maxNumLevels,-1);
      a_finestLevel = 0;

      // assumes dx=dy=dz
      a_amrDx[0] = domainSize[0]/numCells[0];

      IntVect domLo = IntVect::Zero;
      IntVect domHi  = numCells - IntVect::Unit;

      ProblemDomain baseDomain(domLo, domHi, is_periodic);
      a_amrDomains[0] = baseDomain;

      // set up refined domains, etc
      for (int lev=1; lev<= maxLevel; lev++)
      {
        a_amrDomains[lev] = a_amrDomains[lev-1];
        a_amrDomains[lev].refine(a_refRatios[lev-1]);
        a_amrDx[lev] = a_amrDx[lev-1]/a_refRatios[lev-1];
      }

      Vector<Vector<Box> > vectBoxes(maxLevel+1);

      // local scope. for base-level grid generation
      {
        // generate base level grids

        domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

        Vector<int> procAssign(vectBoxes[0].size(), 0);

        LoadBalance(procAssign, vectBoxes[0]);

        DisjointBoxLayout baseGrids(vectBoxes[0], procAssign, baseDomain);

        a_amrGrids[0] = baseGrids;
      }


      if (maxLevel > 0)
      {
        // tag on grad(rhs)
        int bufferSize = 1;
        BRMeshRefine meshGen(a_amrDomains[0],
                             a_refRatios,
                             fillRatio,
                             blockFactor,
                             bufferSize,
                             maxBoxSize);

        // to be used by MeshRefine...
        Vector<Vector<Box> > oldMeshes(maxLevel+1);
        oldMeshes[0] = vectBoxes[0];
        for (int lev=1; lev<oldMeshes.size(); lev++)
        {
          oldMeshes[lev].push_back(a_amrDomains[lev].domainBox());
        }

        Real refineThresh;
        ppGrids.get("refine_threshold", refineThresh);

        Real threshSqr = refineThresh*refineThresh;

        bool moreLevels = true;
        while (moreLevels)
        {
          // tag based on grad(rhs)
          // first need to allocate RHS
          Vector<RefCountedPtr<ch_ldf_cell> > tempRHS(a_finestLevel+1);
          for (int lev=0; lev<= a_finestLevel; lev++)
          {
            // note that we add a ghost cell to simplify gradients
            tempRHS[lev] = RefCountedPtr<ch_ldf_cell>(new ch_ldf_cell(a_amrGrids[lev],
                                                                      SpaceDim,
                                                                      IntVect::Unit));
          }

          double max_grad_sq = 0;
          setRHS(tempRHS, a_amrDx);
          Vector<IntVectSet> tags(a_finestLevel+1);
          for (int lev=0; lev<a_finestLevel+1; lev++)
          {
            const DisjointBoxLayout& levelGrids = a_amrGrids[lev];
            const LevelData<FArrayBox>& levelRHS = *tempRHS[lev];
            IntVectSet& levelTags = tags[lev];

            // compute mag(gradient)
            DataIterator dit = levelGrids.dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
            {
              const FArrayBox& rhsFab = levelRHS[dit];
              // local storage foer gradient
              FArrayBox gradFab(levelGrids[dit],SpaceDim);
              gradFab.setVal(0.0);
              Real thisGrad;

              BoxIterator bit(levelGrids[dit]);
              for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv=bit();
                //
                double sum_grad_sq = 0; //
                for (int dir=0; dir<SpaceDim; dir++)
                {
                  // use mag(undivided gradient)
                  IntVect hi = iv + BASISV(dir);
                  IntVect lo = iv - BASISV(dir);
                  thisGrad = rhsFab(hi,dir) - rhsFab(lo,dir);
                  gradFab(iv, dir) = (thisGrad*thisGrad);
                  sum_grad_sq += (thisGrad*thisGrad);
                }
                if(sum_grad_sq > max_grad_sq)
                {
                  max_grad_sq = sum_grad_sq;
                }
              } // end loop over cells

              //gradFab now has mag(grad*dx)^2

              // tag where mag(gradient) > tolerance^2
              for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                for (int comp=0; comp<gradFab.nComp(); comp++)
                {
                  double gradval = gradFab(iv, comp);
                  if (gradval > threshSqr)
                  {
                    levelTags |= iv;
                  }
                } // end loop over grad comps
              } // end loop over cells
            } // end loop over grids on this level

          } // end loop over levels
          Chombo::pout() << "max sum(rhs gradient) on this proc = "  << max_grad_sq << endl;


          // call meshRefine.
          for (int lev=1; lev<=a_finestLevel; lev++)
          {
            oldMeshes[lev] = vectBoxes[lev];
          }

          int topLevel = a_finestLevel;
          int newFinestLevel =  meshGen.regrid(vectBoxes,
                                               tags,
                                               0,
                                               topLevel,
                                               oldMeshes);


          // define new grids if necessary and test to see if we're done
          if (newFinestLevel > a_finestLevel)
          {
            a_finestLevel = newFinestLevel;

            // setup new grid hierarchy
            for (int lev=1; lev<=a_finestLevel; lev++)
            {
              Vector<int> procAssign(vectBoxes[lev].size(),0);
              LoadBalance(procAssign, vectBoxes[lev]);
              DisjointBoxLayout levelGrids(vectBoxes[lev],
                                           procAssign,
                                           a_amrDomains[lev]);
              a_amrGrids[lev] = levelGrids;
            }
          }
          else
          {
            moreLevels = false;
          }

          if (a_finestLevel == maxLevel)
          {
            moreLevels = false;
          }

        } // end while (moreLevels)

      }

      // fill in remaining levels with empty DisjointBoxLayouts
      if(maxLevel != (a_finestLevel+1))
      {
        MayDay::Error("I think underrefining will not work here");
      }

    }
    
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

      ParmParse pp("test_harness");
      setupGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);
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
