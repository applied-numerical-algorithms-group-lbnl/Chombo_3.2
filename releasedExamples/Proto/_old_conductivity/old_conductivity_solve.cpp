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
#include "AMRIO.H"



namespace Chombo
{


  
  ///
  void
  poissonSolve(Vector<              RefCountedPtr<   LevelData< FArrayBox > >  > & a_phi,
               const Vector<        RefCountedPtr<   LevelData< FArrayBox > >  > & a_rhs,
               RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >          & a_opFactory,
               const Vector< DisjointBoxLayout >                                 & a_grids,
               const VCPoissonParameters                                         & a_params)
  {
    string which_solver;
    ParmParse ppSolver("solver");
    ppSolver.get("which_solver", which_solver);
    // set up multigrid solver parameters
    int numSmooth, numMG, maxIter;
    Real eps, hang;

    ppSolver.get("num_smooth", numSmooth);
    ppSolver.get("num_mg",     numMG);
    ppSolver.get("max_iterations", maxIter);
    ppSolver.get("tolerance", eps);
    ppSolver.get("hang",      hang);
    Real normThresh = 1.0e-10;

    if(which_solver == string("AMRMultiGrid"))
    {
      BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
      AMRMultiGrid<  LevelData<FArrayBox> > amr_mg_solver;
      amr_mg_solver.define(a_params.m_coarsestDomain,
                           (*a_opFactory),
                           &bottomSolver,
                           a_params.m_numLevels);
    
      amr_mg_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                                        numMG, maxIter, eps, hang, normThresh);

      amr_mg_solver.m_verbosity = 3;
      Vector<LevelData<FArrayBox>*> raw_phi = PrChUtilities<1>::getRawVectorCh(a_phi);
      Vector<LevelData<FArrayBox>*> raw_rhs = PrChUtilities<1>::getRawVectorCh(a_rhs);
      
      amr_mg_solver.solve(raw_phi, raw_rhs, a_params.m_maxLevel, 0, true);
    } //end use amrmultigrid
    else if(which_solver == string("BiCGStab"))
    {
      Vector<RealVect>           vectDx(a_params.m_numLevels);
      Vector<ProblemDomain> vectDomains(a_params.m_numLevels);
      
      RealVect dxLev     =      a_params.m_coarsestDx*RealVect::Unit;
      ProblemDomain domLev(     a_params.m_coarsestDomain);
      for (int ilev = 0; ilev < a_params.m_numLevels; ilev++)
      {
        vectDomains[ilev] = domLev;
        vectDx     [ilev] =  dxLev;
        if(ilev < (a_params.m_numLevels-1))
        {
          domLev.refine  (a_params.m_refRatio[ilev]);
          dxLev   /= Real(a_params.m_refRatio[ilev]);
        }
      }
      
      int lBase = 0;
      MultilevelLinearOp<FArrayBox> mlOp;

      mlOp.m_num_mg_iterations = numMG;
      mlOp.m_num_mg_smooth = numSmooth;
      int preCondSolverDepth = -1;
      mlOp.m_preCondSolverDepth = preCondSolverDepth;


      mlOp.define(a_grids, a_params.m_refRatio, vectDomains,
                  vectDx, a_opFactory, lBase);

      BiCGStabSolver<Vector<LevelData<FArrayBox>* > > solver;

      bool homogeneousBC = false;
      solver.define(&mlOp, homogeneousBC);
      solver.m_verbosity = 5;///bicgstab does not print much if you use lower
      solver.m_normType = 0;
      solver.m_eps = eps;
      solver.m_imax = maxIter;

      Vector<LevelData<FArrayBox>*> raw_phi = PrChUtilities<1>::getRawVectorCh(a_phi);
      Vector<LevelData<FArrayBox>*> raw_rhs = PrChUtilities<1>::getRawVectorCh(a_rhs);
      solver.solve(raw_phi, raw_rhs);
    } //end using bicgstab
    else
    {
      MayDay::Error("poissonSolver: bogus which_solver");
    }//end switch on solver type
  }//end function poissonSolve
  
  ///call define but do not fill data
  void
  defineData(Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_phi,
             Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_rhs,
             Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_aCoef,
             Vector< RefCountedPtr< LevelData<   FluxBox > > >& a_bCoef,
             const Vector< DisjointBoxLayout >                & a_grids)
  {
    a_phi  .resize(a_grids.size());
    a_rhs  .resize(a_grids.size());
    a_aCoef.resize(a_grids.size());
    a_bCoef.resize(a_grids.size());
    
    for(int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      a_phi  [ilev] = RefCountedPtr< LevelData< FArrayBox > >( new   LevelData< FArrayBox > (a_grids[ilev], 1, IntVect::Unit));
      a_rhs  [ilev] = RefCountedPtr< LevelData< FArrayBox > >( new   LevelData< FArrayBox > (a_grids[ilev], 1, IntVect::Zero));
      a_aCoef[ilev] = RefCountedPtr< LevelData< FArrayBox > >( new   LevelData< FArrayBox > (a_grids[ilev], 1, IntVect::Zero));
      a_bCoef[ilev] = RefCountedPtr< LevelData<   FluxBox > >( new   LevelData<   FluxBox > (a_grids[ilev], 1, IntVect::Zero));
    }  ///end loop over levels
  }    //end function defineData

  /***************/
  void
  outputData(const Vector< RefCountedPtr<LevelData< FArrayBox > > >&   a_phi,
             const Vector< RefCountedPtr<LevelData< FArrayBox > > >&   a_rhs,
             const Vector< DisjointBoxLayout       >&   a_grids,
             const VCPoissonParameters              &   a_params)
  {
#ifdef CH_USE_HDF5

    /// a good UI is a joy forever
    Vector<LevelData<FArrayBox>*> raw_phi = PrChUtilities<1>::getRawVectorCh(a_phi);
    Vector<LevelData<FArrayBox>*> raw_rhs = PrChUtilities<1>::getRawVectorCh(a_rhs);
    string         phiFileName(          "phi.hdf5");
    Vector<string> phiVarNames(1, string("phi"));
    string         rhsFileName(          "rhs.hdf5");
    Vector<string> rhsVarNames(1, string("rhs"));

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

  /// sets acoef(x, y, z) = x
  void
  setACoef(Vector<RefCountedPtr< LevelData<FArrayBox> > >  & a_aCoefVec,
           const VCPoissonParameters                       & a_params)
  {

    for(int ilev = 0; ilev < a_params.m_numLevels; ilev++)
    {
      Real dx = a_params.m_coarsestDx;
      DataIterator dit = a_aCoefVec[ilev]->dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& aCoef = (*a_aCoefVec[ilev])[dit];
        for(BoxIterator boxit(aCoef.box()); boxit.ok(); ++boxit)
        {
          auto iv = boxit();
          RealVect pos = VCLocalFunctions::cellLocation(iv, dx);
          aCoef(iv, 0) = pos[0];
        } // end loop over cells
      }   // end loop over grids
      if(ilev < a_params.m_maxLevel)
      {
        dx /= a_params.m_refRatio[ilev];
      }   // end if on a coarse level
    }     // end loop over levels
  }       //end function setACoef

  /// sets bcoef(x, y, z) = x + y + z
  static void
  setBCoef(Vector< RefCountedPtr< LevelData<FluxBox> > >  & a_bCoef,
           const VCPoissonParameters                      & a_params)
  {
    for(int ilev = 0; ilev < a_params.m_numLevels; ilev++)
    {
      Real  dx = a_params.m_coarsestDx;
      DataIterator dit = a_bCoef[ilev]->dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        FluxBox& thisBCoef = (*a_bCoef[ilev])[dit];
        for (int face_dir=0; face_dir < SpaceDim; face_dir++)
        {
          FArrayBox& bcofab = thisBCoef[face_dir];
          const Box& dirBox = bcofab.box();
          RealVect pos;
          for(BoxIterator boxit(dirBox); boxit.ok(); ++boxit)
          {
            IntVect   iv  = boxit();
            RealVect  pos = VCLocalFunctions::faceLocation(iv, dx, face_dir);
            //bcoef = x + y + z in the amrpoissonop example
            Real sum = pos.sum();
            bcofab(iv,  0) = sum;
          } // end loop over cells
        }   // end loop over face direction
      }     // end loop over boxes
      if(ilev < a_params.m_maxLevel)
      {
        dx /= a_params.m_refRatio[ilev];
      }     // end if on a coarse level
    }       // end loop over levels
  }//.end function setBCoef

  /********/
  void setRHSConst(Vector< RefCountedPtr<LevelData<FArrayBox> > >   a_rhs,
                   Real                                             a_value)
  {
    for(int ilev = 0; ilev < a_rhs.size(); ilev++)
    {
      auto& rhs_lev = *(a_rhs[ilev]);
      for (DataIterator dit = rhs_lev.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& thisRHS = rhs_lev[dit()];
        thisRHS.setVal(a_value);
      }
    }
  }///end function setRHSConst

  int local_test_harness()
  {
    Vector<RefCountedPtr<LevelData< FArrayBox > > > phi;
    Vector<RefCountedPtr<LevelData< FArrayBox > > > rhs;
    Vector<RefCountedPtr<LevelData< FArrayBox > > > aCoef;
    Vector<RefCountedPtr<LevelData<   FluxBox > > > bCoef;

    Vector<DisjointBoxLayout> amrGrids;
    Vector<ProblemDomain> amrDomains;
    Vector<int> refRatios;
    Vector<Real> amrDx;
    int finestLevel;

    PrChUtilities<1>::setupLLCornerGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);
    for(int ilev = 0; ilev < amrGrids.size(); ilev++)
    {
      pout() << "ilev = " << ilev << ", grids = " << endl;
      amrGrids[ilev].print();
    }

    VCPoissonParameters param = 
      getVCParameters(finestLevel, amrDx[0], amrDomains[0], refRatios);
    

    defineData(phi, rhs, aCoef, bCoef,   amrGrids);
    setRHSConst(rhs, 1.0);
    setACoef(aCoef, param);
    setBCoef(bCoef, param);

    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> >  >  
      op_factory = VCLocalFunctions::getConductivityFactory(aCoef, bCoef, amrGrids, param);

    poissonSolve(phi, rhs, op_factory, amrGrids, param);
    outputData(  phi, rhs,             amrGrids, param);

    return 0;
  }
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

  int status = Chombo::local_test_harness();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
