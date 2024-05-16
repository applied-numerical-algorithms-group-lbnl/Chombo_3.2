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
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "FABView.H"
#include "DebugDump.H"
#include "SSAUtilities.H"
#include "SigmaCS.H"
#include "defineSigmaCS.H"
#include "BiCGStabSolver.H"
#include "PicardSolver.H"
#include "ConstitutiveRelation.H"

#include "UsingNamespace.H"

#ifdef CH_Linux
// Should be undefined by default
//#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

using std::cerr;

enum rhsTypes {zeroRHS = 0,
               hydrostaticRHS, 
               numRHSTypes};

enum muTypes {unityMu = 0,
              linearMu,
              GlensLawMu,
              numMuTypes};

// define the test problem
int muType = GlensLawMu;
//int muType = unityMu;
int solverType = multigrid;
//int solverType = BiCGStab;
//int thicknessType = constantThickness;
int thicknessType = sinusoidalH;
int basalType = constantZb;
//int basalType = xInclineZb;
int phiType = constZero;
//int rhsType = zeroRHS;
int rhsType = hydrostaticRHS;


void
setBasalShear(LevelData<FArrayBox>& a_beta,
              const PoissonParameters& a_params,
              const RealVect& a_dx)
{
  RealVect pos;
  int num;
  DataIterator dit = a_beta.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& beta = a_beta[dit];
      ForAllXBNN(Real, beta, beta.box(), 0, beta.nComp());
      {
        num = nR;
        D_TERM(pos[0]=a_dx[0]*(iR+0.5);,
               pos[1]=a_dx[1]*(jR+0.5);,
               pos[2]=a_dx[2]*(kR+0.5));
        //betaR = pos[0];
        // constant-coefficient
        betaR = 1000.0;
        //betaR = 0.0;
      }EndFor;
    } // end loop over grids
}



void getRHS(LevelData<FArrayBox>& a_rhs, 
            const LevelData<SigmaCS>& a_coordSys, 
            const PoissonParameters& a_params)
{
  // ice density from Pattyn(2003) in kg/m^3
  Real rho = 910.0;
  // gravitational acceleration in m/s^2
  Real g = 9.81;
    

  if (rhsType == zeroRHS)
    {
      // initial stab at things -- rhs = 0
      DataIterator dit = a_rhs.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          a_rhs[dit].setVal(0.0);
        }
    } 
  else if (rhsType == hydrostaticRHS)
    {
      const DisjointBoxLayout& grids = a_rhs.getBoxes();
      // need a ghost cell here in order to be able to compute gradients
      LevelData<FArrayBox> surfaceZ(grids, 1, IntVect::Unit);
      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const SigmaCS& thisCS = a_coordSys[dit];
          const FArrayBox& thisZb = thisCS.getBaseHeight();
          const FArrayBox& thisH = thisCS.getH();
          const RealVect& dx = thisCS.dx();

          FArrayBox& thisZsurf = surfaceZ[dit];
          thisZsurf.copy(thisZb, thisZsurf.box());
          thisZsurf.plus(thisH, thisZsurf.box(), 0, 0, 1);
          // now take grad(z_s)
          FArrayBox& thisRhs = a_rhs[dit];
          BoxIterator bit(grids[dit]);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              for (int dir=0; dir<2; dir++)
                {
                  // 2D implementation
                  if (SpaceDim != 2)
                    {
                      MayDay::Error("testPicardSolver getRHS only implemented in 2D");
                    }
                  IntVect offsetVect = BASISV(dir);
                  IntVect hi = iv + offsetVect;
                  IntVect lo = iv - offsetVect;
                  
                  Real grad = (thisZsurf(hi,0) - thisZsurf(lo,0))/dx[dir];
                  
                  thisRhs(iv,dir) = rho*g*thisH(iv,0)*grad;
                } // end loop over directions
            } // end loop over cells in this box
        } // end loop over boxes
    } // end if hydrostatic RHS
  else
    { 
      MayDay::Error("bad RHS type");
    }
}


/******/
int doSolve(Vector<LevelData<FArrayBox>* >& a_phi,
            Vector<LevelData<FArrayBox>* >& a_rhs,
            const Vector< DisjointBoxLayout >&   a_grids,
            Vector<RefCountedPtr<LevelData<SigmaCS> > >& a_coordSys,
            const PoissonParameters&             a_params)
{
  ParmParse pp;

  Real initialGuess = 1.0;     
  pp.query("initial_guess", initialGuess);

  int nlevels = a_params.numLevels;
  a_phi.resize(nlevels);
  a_rhs.resize(nlevels);
  Vector<LevelData<FArrayBox>* > beta(nlevels);

  Vector<ProblemDomain> vectDomains(nlevels);
  Vector<RealVect> vectDx(nlevels);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  IntVect sigmaCSGhost = IntVect::Unit;
          
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      // rhs and phi have 2 components each
      a_rhs[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 2,IntVect::Zero);
      a_phi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 2,IntVect::Unit);
      beta[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero);
      // ncomp for SigmaCS is meaningless, but required by LevelData API
      a_coordSys[ilev] = RefCountedPtr<LevelData<SigmaCS> >(new LevelData<SigmaCS>(a_grids[ilev], 1, sigmaCSGhost));
      vectDomains[ilev] = domLev;
      vectDx[ilev] = dxLev;

      setBasalShear(*beta[ilev], a_params, dxLev);

      RealVect dx = dxLev*RealVect::Unit;
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_phi[ilev])[dit()].setVal(initialGuess);

          Box ghostBox(a_grids[ilev][dit]);
          ghostBox.grow(sigmaCSGhost);
          SigmaCS& localCS = (*a_coordSys[ilev])[dit];
          localCS.define(ghostBox,dx);
        }
        
      RealVect domainSize = a_params.domainSize;
      defineSigmaCS(*a_coordSys[ilev], domainSize, thicknessType, basalType);

      getRHS (*a_rhs[ilev], *a_coordSys[ilev], a_params);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
   }



  // set up solver
  ConstitutiveRelation* constRel = NULL;
  if (muType == unityMu)
    {
      constMuRelation* constMu = new constMuRelation;
      constMu->setConstVal(1.0);
      constRel = static_cast<ConstitutiveRelation*>(constMu);
    }
  else if (muType == GlensLawMu)
    {
      GlensFlowRelation* glensMu = new GlensFlowRelation;
      constRel = static_cast<ConstitutiveRelation*>(glensMu);
    }
  else
    {
      MayDay::Error("invalid Constitutive Relation");
    }

  PicardSolver solver;
  solver.setSolverType(solverType);
  solver.define(a_params.coarsestDomain,
                constRel,
                a_grids,
                a_params.refRatio,
                a_params.coarsestDx,
                &ParseBC,
                nlevels);

  solver.setVerbosity(a_params.verbosity);

  int lBase = 0;
  int maxLevel = nlevels -1;
  int exitStatus = solver.solve(a_phi, a_rhs, beta, a_coordSys, 
                                lBase, maxLevel);
                                

  // clean up local storage
  delete constRel;

  for (int lev=0; lev<beta.size(); lev++)
    {
      delete beta[lev];
      beta[lev] = NULL;
    }
  return exitStatus;

}

/***************/
void outputData(const Vector<LevelData<FArrayBox>* >&   a_phi,
                const Vector<LevelData<FArrayBox>* >&   a_rhs,
                const Vector< DisjointBoxLayout >&       a_grids,
                const Vector<RefCountedPtr<LevelData<SigmaCS> > >& a_coordSys,
                const PoissonParameters&                 a_params)
{
#ifdef CH_USE_HDF5

#if CH_SPACEDIM==2
    string fileName("PicardTestOut.");
#else
    string fileName("PicardTestOut.");
#endif

    int nPhiComp = a_phi[0]->nComp();
    int nRhsComp = a_rhs[0]->nComp();
    int totalComp = nPhiComp + nRhsComp;

    Vector<string> phiNames(totalComp);
    // hardwire to single-component
    CH_assert(totalComp == 4);
    phiNames[0] = "xVel";
    phiNames[1] = "yVel";
    phiNames[2] = "rhs_x";
    phiNames[3] = "rhs_y";


    CH_assert(a_phi.size() == a_rhs.size());
    Vector<LevelData<FArrayBox>* > tempData(a_phi.size(), NULL);
    for (int level=0; level<a_phi.size(); level++)
      {
        tempData[level] = new LevelData<FArrayBox>(a_grids[level], totalComp);
        Interval phiComps(0, nPhiComp-1);
        Interval rhsComps(nPhiComp, totalComp-1);
        a_phi[level]->copyTo(a_phi[level]->interval(),
                             *tempData[level], phiComps);
        a_rhs[level]->copyTo(a_rhs[level]->interval(),
                             *tempData[level], rhsComps);
      }


    // need to use pointers instead of refCountedPtr's
    Vector<const LayoutData<SigmaCS>* > coords(a_coordSys.size(), NULL);
    for (int lev=0; lev<coords.size(); lev++)
      {
        coords[lev] = dynamic_cast<const LayoutData<SigmaCS>* >(&(*a_coordSys[lev]));
      }

    Real fakeTime = 1.0;
    Real fakeDt = 1.0;

    WriteSigmaMappedAMRHierarchyHDF5(fileName, a_grids,
                                tempData, phiNames,
                                coords,
                                a_params.coarsestDomain.domainBox(),
                                fakeDt, fakeTime,
                                a_params.refRatio,
                                a_params.numLevels);

    // clean up temporary storage
    for (int level=0; level<a_phi.size(); level++)
      {
        delete tempData[level];
        tempData[level] = NULL;
      }
#endif

}




int main(int argc, char* argv[])
{
  int status = 0;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#if 0
    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
#endif

    char* inFile;
    if (argc > 2)
      {
        inFile = argv[1];
      }
    else
      {
        pout() << "using picardSolver.inputs" << endl;
        inFile = "picardSolver.inputs";
      }        
    
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters param;
    Vector<DisjointBoxLayout> grids;

    //read params from file
    getPoissonParameters(param);
    int nlevels = param.numLevels;
    Vector<LevelData<FArrayBox>* > phi(nlevels, NULL);
    Vector<LevelData<FArrayBox>* > rhs(nlevels, NULL);
    Vector<RefCountedPtr<LevelData<SigmaCS> > > vectCoordSys(nlevels);

    setGrids(grids,  param);
    

    status = doSolve(phi, rhs, grids,  vectCoordSys, param);

    int dofileout = false;
    pp.get("write_output", dofileout);
    if(dofileout == 1)
      {
        outputData(phi, rhs, grids, vectCoordSys, param);
      }

    // clear memory
    for (int level = 0; level<phi.size(); level++)
      {
        if (phi[level] != NULL)
          {
            delete phi[level];
            phi[level] = NULL;
          }
        if (rhs[level] != NULL)
          {
            delete rhs[level];
            rhs[level] = NULL;
          }
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
