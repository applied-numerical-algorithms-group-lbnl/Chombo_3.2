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
#include "DebuggingTools.H"
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
    // set up multigrid solver parameters
    int numSmooth, numMG, maxIter;
    Real eps, hang;

    ppSolver.get("num_smooth", numSmooth);
    ppSolver.get("num_mg",     numMG);
    ppSolver.get("max_iterations", maxIter);
    ppSolver.get("tolerance", eps);
    ppSolver.get("hang",      hang);
    Real normThresh = 1.0e-10;


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

  }//end function poissonSolve
  
  ///call define but do not fill data
  void
  defineData(Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_phi,
             Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_rhs,
             Vector< RefCountedPtr< LevelData< FArrayBox > > >& a_aco,
             Vector< RefCountedPtr< LevelData<   FluxBox > > >& a_eta,
             Vector< RefCountedPtr< LevelData<   FluxBox > > >& a_lam,
             const Vector< DisjointBoxLayout >                & a_grids)
  {
    a_phi.resize(a_grids.size());
    a_rhs.resize(a_grids.size());
    a_aco.resize(a_grids.size());
    a_eta.resize(a_grids.size());
    a_lam.resize(a_grids.size());
    
    for(int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      a_phi[ilev] = RefCountedPtr< LevelData< FArrayBox > >( new   LevelData< FArrayBox > (a_grids[ilev],DIM, IntVect::Unit));
      a_rhs[ilev] = RefCountedPtr< LevelData< FArrayBox > >( new   LevelData< FArrayBox > (a_grids[ilev],DIM, IntVect::Zero));
      a_aco[ilev] = RefCountedPtr< LevelData< FArrayBox > >( new   LevelData< FArrayBox > (a_grids[ilev], 1 , IntVect::Zero));
      a_eta[ilev] = RefCountedPtr< LevelData<   FluxBox > >( new   LevelData<   FluxBox > (a_grids[ilev], 1 , IntVect::Zero));
      a_lam[ilev] = RefCountedPtr< LevelData<   FluxBox > >( new   LevelData<   FluxBox > (a_grids[ilev], 1 , IntVect::Zero));
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

  /// sets acoef(x, y, z) = x
  void
  setACoef(Vector<RefCountedPtr< LevelData<FArrayBox> > >  & a_aCoefVec,
           const VCPoissonParameters                       & a_params)
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

  /// sets bcoef(x, y, z) = x + y + z
  static void
  setEtaAndLambda(Vector< RefCountedPtr< LevelData<FluxBox> > >  & a_eta,
                  Vector< RefCountedPtr< LevelData<FluxBox> > >  & a_lam,
                  const VCPoissonParameters                      & a_params)
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
    Vector<RefCountedPtr<LevelData<   FluxBox > > > eta;
    Vector<RefCountedPtr<LevelData<   FluxBox > > > lambda;

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
    

    defineData(phi, rhs, aCoef, eta, lambda,  amrGrids);
    setRHSConst(rhs, 1.0);
    setACoef(aCoef, param);
    setEtaAndLambda(eta, lambda, param);

    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> >  >  
      op_factory = VCLocalFunctions::getViscousTensorFactory(aCoef, eta, lambda, amrGrids, param);

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