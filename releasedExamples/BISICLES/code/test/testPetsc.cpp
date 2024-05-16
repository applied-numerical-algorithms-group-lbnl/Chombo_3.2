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
#include "computeNorm.H"
#include "ViscousTensorOp.H"
#include "PetscSolver.H"
#include "CONSTANTS.H"

#include "FABView.H"


#include "UsingNamespace.H"

enum etaTypes {unityEta = 0,
               linearEta,linearyEta,linearxEta,
               sinusoidalEta,
               numEtaTypes};

enum lambdaTypes {zeroLambda = 0,
                  twiceEta,
                  numLambdaTypes};


enum aTypes {zeroA = 0,
             unityA,
             sinusoidalA,
             numATypes};

enum phiTypes {constZero = 0,
               constOne,
               constU,
               constV,
               linearPhi,
               sinusoidalPhi,
               numPhiTypes};

bool writePlots = true;

//int aType = zeroA;
int aType = unityA;
//int aType = sinusoidalA;

//int etaType = unityEta;
int etaType = sinusoidalEta;
//int etaType = linearxEta;

//int lambdaType = zeroLambda;
int lambdaType = twiceEta;

//int phiType = constZero;
//int phiType = constOne;
//int phiType = constU;
//int phiType = constV;
int phiType = sinusoidalPhi;

Real diffTol = 1.0e-8;

void noOpBC(FArrayBox& a_state,
            const Box& a_vlaid,
            const ProblemDomain& a_domain,
            Real a_dx, 
            bool a_homogeneous)
{
  // does nothing, at lest for the moment
}

// put this BCHolder class in the global scope
BCHolder testBC(&noOpBC);


void
setEta(LevelData<FluxBox>& a_eta,
       Real a_dx)
{
  DataIterator dit = a_eta.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisEta = a_eta[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& dirFlux = thisEta[dir];
          const Box& dirBox = dirFlux.box();
          // this sets up a vector which is 0 in the dir
          // direct and 0.5 in the other (cell-centered) directions
          RealVect offsets = BASISREALV(dir);
          RealVect pos;
          offsets -= RealVect::Unit;
          offsets *= -0.5;
          int n;
          ForAllXBNN(Real, dirFlux, dirBox, 0, dirFlux.nComp())
            {
              n = nR;
              D_TERM(pos[0] = a_dx*(iR+offsets[0]);,
                     pos[1] = a_dx*(jR+offsets[1]);,
                     pos[2] = a_dx*(kR+offsets[2]));
              if (etaType == unityEta)
                {
                  dirFluxR = 1.0;
                }
              else if (etaType == linearEta)
                {
                  dirFluxR = D_TERM(pos[0], +pos[1], +pos[2]);
                }
              else if (etaType == linearxEta)
                {
                  dirFluxR = pos[0];
                }
              else if (etaType == linearyEta)
                {
                  dirFluxR = pos[1];
                }
              else if (etaType == sinusoidalEta)
                {
                  dirFluxR = 2.0*cos(4.0*Pi*pos[0])*sin(4.0*Pi*pos[1]);
                }
              else
                {
                  MayDay::Error("bad etaType");
                }
            }EndFor
       } // end loop over directions
    }
}

void
setLambda(LevelData<FluxBox>& a_lambda,
          const LevelData<FluxBox>& a_eta,
          Real a_dx)
{
  DataIterator dit = a_lambda.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FluxBox& thisEta = a_eta[dit];
      FluxBox& thisLambda = a_lambda[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& dirLambda = thisLambda[dir];
          const FArrayBox& dirEta = thisEta[dir];

          if (lambdaType == zeroLambda)
            {
              dirLambda.setVal(0.0);
            }
          else if (lambdaType == twiceEta)
            {
              // mirror BISICLES and set lambda = 2*eta
              dirLambda.copy(dirEta);
              dirLambda *= 2.0;
            }
        } // end loop over directions
    }
}



void
setA(LevelData<FArrayBox>& a_Acoef,
     Real a_dx)
{
  const DisjointBoxLayout grids = a_Acoef.getBoxes();
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisA = a_Acoef[dit];

      RealVect offsets = 0.5*RealVect::Unit;
      RealVect pos;

      
      BoxIterator bit(thisA.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          
          IntVect iv = bit();
          D_TERM(pos[0] = a_dx*(iv[0]+offsets[0]);,
                 pos[1] = a_dx*(iv[1]+offsets[1]);,
                 pos[2] = a_dx*(iv[2]+offsets[2]));

          if (aType == zeroA)
            {
                     thisA(iv,0) = 0.0;
            }
          else if (aType == unityA)
            {
              thisA(iv,0) = 1.0;
            }
          else if (aType == sinusoidalA)
            {
              // offset to keep A positive
              thisA(iv,0) = 2.0 + sin(two*Pi*pos[0])*cos(two*Pi*pos[1]);
            }
          else
            {
              MayDay::Error("Bad aType");
            }

        } // end loop over cells
    }
}

void
setPhi(LevelData<FArrayBox>& a_phi,
       Real a_dx)
{
  // this only sets Phi on the valid cells, since
  // we're trying to test to be sure that BC's are being properly 
  // handled
  const DisjointBoxLayout grids = a_phi.getBoxes();
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];
      const Box& gridBox = grids[dit];

      RealVect offsets = 0.5*RealVect::Unit;
      RealVect pos;

      
      BoxIterator bit(gridBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          
          IntVect iv = bit();
          D_TERM(pos[0] = a_dx*(iv[0]+offsets[0]);,
                 pos[1] = a_dx*(iv[1]+offsets[1]);,
                 pos[2] = a_dx*(iv[2]+offsets[2]));

          if (phiType == constZero)
            {
              D_TERM(
                     thisPhi(iv,0) = 0.0;,
                     thisPhi(iv,1) = 0.0;,
                     thisPhi(iv,2) = 0.0; )
            }
          else if (phiType == constOne)
            {
              D_TERM(
                     thisPhi(iv,0) = 1.0;,
                     thisPhi(iv,1) = 1.0;,
                     thisPhi(iv,2) = 0.0; )
            }
          else if (phiType == constU)
            {
              D_TERM(
                     thisPhi(iv,0) = 1.0;,
                     thisPhi(iv,1) = 0.0;,
                     thisPhi(iv,2) = 0.0; )
            }
          else if (phiType == constV)
            {
              D_TERM(
                     thisPhi(iv,0) = 0.0;,
                     thisPhi(iv,1) = 1.0;,
                     thisPhi(iv,2) = 0.0; )
            }
          else if (phiType == sinusoidalPhi)
            {
              D_TERM(
                     thisPhi(iv,0) = sin(two*Pi*pos[0])*cos(4*Pi*pos[1]);,
                     thisPhi(iv,1) = cos(two*Pi*pos[1])*sin(2*Pi*pos[1]);,
                     thisPhi(iv,2) = sin(two*Pi*(pos[0] + pos[1])); )
            }                     
          else
            {
              MayDay::Error("Bad phiType");
            }

        } // end loop over cells
    }
}


ViscousTensorOp*
defineSolvers(
#ifdef CH_USE_PETSC
              PetscSolverViscousTensor<LevelData<FArrayBox> > *& a_petscSolverPtr,
#endif
              const DisjointBoxLayout& a_levelGrids, 
              const DisjointBoxLayout& a_crseGrids,                   
              const DisjointBoxLayout& a_fineGrids,
              const ProblemDomain& a_domain, 
              int a_nRefCrse,
              int a_nRefFine,
              Real a_dx, 
              Real a_dxCrse)
{
  
  // first, allocate storage for coefficients
  RefCountedPtr<LevelData<FArrayBox> > aCoef(new LevelData<FArrayBox>(a_levelGrids, 1, IntVect::Zero));
  RefCountedPtr<LevelData<FluxBox> >  etaCoef(new LevelData<FluxBox>(a_levelGrids, 1, IntVect::Zero));
  RefCountedPtr<LevelData<FluxBox> > lambdaCoef(new LevelData<FluxBox>(a_levelGrids, 1, IntVect::Zero));
  
  // set coefficient values
  setA(*aCoef, a_dx);
  setEta(*etaCoef, a_dx);
  setLambda(*lambdaCoef, *etaCoef, a_dx);
  
  Real alpha = 1.0;
  Real beta = 1.0;

  // define operators
  ViscousTensorOp* VTopPtr = new ViscousTensorOp(a_levelGrids,
                                                 a_fineGrids, a_crseGrids, 
                                                 etaCoef, lambdaCoef, aCoef,
                                                 alpha, beta, 
                                                 a_nRefFine, a_nRefCrse,
                                                 a_domain, a_dx,
                                                 a_dxCrse,
                                                 testBC);

#ifdef CH_USE_PETSC
  //  a_petscSolverPtr = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
  a_petscSolverPtr->define(VTopPtr, true );
  // a_petscSolverPtr->define( alpha, beta, &(*aCoef), &(*etaCoef), &(*lambdaCoef) );
  //a_petscSolverPtr->formMatrix(a_petscSolverPtr->m_mat);
  a_petscSolverPtr->setInitialGuessNonzero();
#endif

  return VTopPtr;
}


int testPetsc()
{
  int status = 0;

  // this code will test the PetscViscousTensorOp vs. 
  // the normal Chombo ViscousTensorOp for various configurations.
  // success means that the two have identical stencils and so produce 
  // identical results from applyOp.

  // first test -- domain-spanning periodic box -- first do single box, then 
  // 2x2, etc.
  {
    int domSize = 16;
    Box domainBox(IntVect::Zero, (domSize-1)*IntVect::Unit);
    ProblemDomain domain(domainBox);
    for (int dir=0; dir<SpaceDim; dir++)
      {
        domain.setPeriodic(dir,true);
      }

    Real dx = 1.0/domSize;

    for (int numBoxes=1; numBoxes<5; numBoxes*=2)
      {
        int maxBoxSize = domSize/numBoxes;
        
        Vector<Box> vectBox;
        domainSplit(domain, vectBox, maxBoxSize);

        Vector<int> procAssign(vectBox.size(), 0);
#ifdef CH_MPI
        LoadBalance(procAssign, vectBox);
#endif
        DisjointBoxLayout levelGrids(vectBox, procAssign, domain);

        // allocate operators
        ViscousTensorOp* VTopPtr = NULL;
#ifdef CH_USE_PETSC
        PetscSolverViscousTensor<LevelData<FArrayBox> >* petscSolverPtr;
#endif

        // no coarse or fine levels for this test
        DisjointBoxLayout bogusDL;
        int nRefCrse = -1;
        int nRefFine = -1;
        Real dxCrse = dx*nRefCrse;
        
        
#ifdef CH_USE_PETSC
        petscSolverPtr = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
#endif
        
        VTopPtr = defineSolvers(
#ifdef CH_USE_PETSC
                                petscSolverPtr,
#endif
                                levelGrids, 
                                bogusDL,
                                bogusDL,
                                domain, 
                                nRefCrse,
                                nRefFine,
                                dx, 
                                dxCrse);
        
        // initialize storage for phi and L(phi)
        LevelData<FArrayBox> phi(levelGrids, SpaceDim, IntVect::Unit);
        LevelData<FArrayBox> phiPetsc(levelGrids, SpaceDim, IntVect::Unit);        
        LevelData<FArrayBox> Lphi(levelGrids, SpaceDim, IntVect::Zero);
        LevelData<FArrayBox> LphiPetsc(levelGrids, SpaceDim, IntVect::Zero);
        LevelData<FArrayBox> diff(levelGrids, SpaceDim, IntVect::Zero);
        
        setPhi( phi, dx );
        
        VTopPtr->applyOp( Lphi, phi, true );
#ifdef CH_USE_PETSC
        // calling solve initializes the matrix, leaves it in a usable state. 
        petscSolverPtr->solve(phiPetsc, Lphi);
        int ierr = petscSolverPtr->applyOp( LphiPetsc, phi );
        if (ierr != 0)
          {
            pout() << "petsc solver applyOp call returned nonzero error code: " 
                   << ierr << endl;
          }
#endif
  
        DataIterator dit = diff.dataIterator();

        for (dit.begin(); dit.ok(); ++dit)
          {
            diff[dit].copy(Lphi[dit]);
            diff[dit].minus(LphiPetsc[dit],0,0,SpaceDim);
          }
        
        // compute norms of diff
        RealVect L1diff, L2diff, maxDiff;
        DisjointBoxLayout* nullPtr = NULL;
        
        for (int dir=0; dir<SpaceDim; dir++)
          {
            Interval comps(dir,dir);
            int normType = 0;
            maxDiff[dir] = computeNorm(diff, nullPtr,
                                       nRefFine, dx, 
                                       comps, normType);
            
            if (maxDiff[dir] > diffTol)
              {
                pout() << "max(diff)[" << dir << "] = " << maxDiff[dir] << endl;
                status++;
              }

            normType = 1;
            L1diff[dir] = computeNorm(diff, nullPtr,
                                      nRefFine, dx, 
                                      comps, normType);

            if (L1diff[dir] > diffTol)
              {
                pout() << "L1(diff)[" << dir << "] = " << L1diff[dir] << endl;
                status++;
              }

            
            normType = 2;
            L2diff[dir] = computeNorm(diff, nullPtr,
                                      nRefFine, dx, 
                                      comps, normType);

            if (L2diff[dir] > diffTol)
              {
                pout() << "L2(diff)[" << dir << "] = " << L2diff[dir] << endl;
                status++;
              }
            
          }
        pout() << "baseLevel test, nGrids = " << numBoxes << ": L1(diff) = " << L1diff << endl;
        pout() << "                         "  << "  L2(diff) = " << L2diff << endl;
        pout() << "                         " << numBoxes << ": max(diff) = " << maxDiff << endl;
        
        // delete storage
        delete VTopPtr;
#ifdef CH_USE_PETSC
        delete petscSolverPtr;
#endif

      } // end loop over numBoxes
  }


  // now test refined-mesh case
  {
    int domSize = 56;
    Box domainBox(IntVect::Zero, (domSize-1)*IntVect::Unit);
    ProblemDomain domain(domainBox);
    for (int dir=0; dir<SpaceDim; dir++)
      {
        domain.setPeriodic(dir,true);
      }

    Real dx = 1.0/domSize;
    
    // increasing amounts of complexity with increasing number of boxes

    Vector<Box> vectBox(8);
    int numBoxes=0; 
    {
      Box gridBox;
      IntVect lo = 8*IntVect::Unit;
      IntVect hi = 15*IntVect::Unit;
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    {
      Box gridBox;
      IntVect lo = 16*IntVect::Unit;
      IntVect hi(D_DECL(39,23,7));
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    {
      Box gridBox;
      IntVect lo(D_DECL(40,8,0));
      IntVect hi(D_DECL(47,15,7));
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    {
      Box gridBox;
      IntVect lo(D_DECL(16,24,0));
      IntVect hi(D_DECL(23,31,7));
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    {
      Box gridBox;
      IntVect lo(D_DECL(32,24,0));
      IntVect hi(D_DECL(39,31,7));
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    {
      Box gridBox;
      IntVect lo(D_DECL(16,32,0));
      IntVect hi(D_DECL(39,39,7));
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    {
      Box gridBox;
      IntVect lo(D_DECL(8,40,0));
      IntVect hi(D_DECL(15,47,7));
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    {
      Box gridBox;
      IntVect lo(D_DECL(40,40,0));
      IntVect hi(D_DECL(47,47,7));
      gridBox.define(lo,hi);
      vectBox[numBoxes++] = (gridBox);
    }
    
    Vector<int> procAssign(vectBox.size(), 0);
#ifdef CH_MPI
    LoadBalance(procAssign, vectBox);
#endif
    DisjointBoxLayout levelGrids(vectBox, procAssign, domain);
    
    // allocate operators
    ViscousTensorOp* VTopPtr = NULL;
#ifdef CH_USE_PETSC
    PetscSolverViscousTensor<LevelData<FArrayBox> >* petscSolverPtr;
#endif
    
    // no coarse or fine levels for this test
    DisjointBoxLayout bogusDL;
    int nRefCrse = 2;
    int nRefFine = -1;
    Real dxCrse = dx*nRefCrse;

    // set up a coarse level in order to compute full multilevel operator
    ProblemDomain crseDomain = domain;
    crseDomain.coarsen(nRefCrse);
    Vector<Box> crseBoxes(1, crseDomain.domainBox());
    Vector<int> crseProcAssign(1,0);
    DisjointBoxLayout crseGrids(crseBoxes, crseProcAssign, crseDomain);

#ifdef CH_USE_PETSC
    petscSolverPtr = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
#endif
    
    VTopPtr = defineSolvers(
#ifdef CH_USE_PETSC
                            petscSolverPtr,
#endif
                            levelGrids, 
                            crseGrids,
                            bogusDL,
                            domain, 
                            nRefCrse,
                            nRefFine,
                            dx, 
                            dxCrse);
    
    // initialize storage for phi and L(phi)
    LevelData<FArrayBox> phi(levelGrids, SpaceDim, IntVect::Unit);
    LevelData<FArrayBox> phiPetsc(levelGrids, SpaceDim, IntVect::Unit);    
    LevelData<FArrayBox> LphiAMR(levelGrids, SpaceDim, IntVect::Zero);
    LevelData<FArrayBox> Lphi(levelGrids, SpaceDim, IntVect::Zero);
    LevelData<FArrayBox> LphiPetsc(levelGrids, SpaceDim, IntVect::Zero);
    LevelData<FArrayBox> diff(levelGrids, SpaceDim, IntVect::Zero);
    LevelData<FArrayBox> diffAMR(levelGrids, SpaceDim, IntVect::Zero);
    
    LevelData<FArrayBox> phiCrse(crseGrids, SpaceDim, IntVect::Unit);
    DataIterator crseDit = crseGrids.dataIterator();
    for (crseDit.begin(); crseDit.ok(); ++crseDit)
      {
        phiCrse[crseDit].setVal(0.0);
      }


    setPhi( phi, dx );

    // compute AMR operator with crse values set to zero
    VTopPtr->AMROperatorNF(LphiAMR, phi, phiCrse, true);

    VTopPtr->homogeneousCFInterp(phi);
    VTopPtr->applyOp( Lphi, phi, true );
#ifdef CH_USE_PETSC
    petscSolverPtr->solve(phiPetsc, Lphi );    
    int ierr = petscSolverPtr->applyOp( LphiPetsc, phi );
    if (ierr != 0)
      {
        pout() << "petsc solver applyOp call returned nonzero error code: " 
               << ierr << endl;
      }
#endif
    
    DataIterator dit = diff.dataIterator();
    
    for (dit.begin(); dit.ok(); ++dit)
      {
        diff[dit].copy(Lphi[dit]);
        diffAMR[dit].copy(Lphi[dit]);
        diff[dit].minus(LphiPetsc[dit],0,0,SpaceDim);
        diffAMR[dit].minus(LphiAMR[dit],0,0,SpaceDim);
      }

    if (writePlots)
      {
        writeLevelname(&Lphi, "chomboOperator.2d.hdf5");
        writeLevelname(&LphiPetsc, "petscOperator.2d.hdf5");
        writeLevelname(&diff, "chomboPetscDiff.2d.hdf5");
      }
    
    // compute norms of diff
    RealVect L1diff, L2diff, maxDiff;
    RealVect L1diffAMR, L2diffAMR, maxDiffAMR;
    DisjointBoxLayout* nullPtr = NULL;
    
    for (int dir=0; dir<SpaceDim; dir++)
      {
        Interval comps(dir,dir);
        int normType = 0;
        maxDiff[dir] = computeNorm(diff, nullPtr,
                                   nRefFine, dx, 
                                   comps, normType);

        maxDiffAMR[dir] = computeNorm(diffAMR, nullPtr,
                                      nRefFine, dx, 
                                      comps, normType);
        
        if (maxDiff[dir] > diffTol)
          {
            pout() << "max(diff)[" << dir << "] = " << maxDiff[dir] << endl;
            status++;
          }

        if (maxDiffAMR[dir] > diffTol)
          {
            pout() << "AMR vs. homogeneous CFInterp -- max(diff)["
                   << dir << "] = " << maxDiff[dir] << endl; 
            status++;
          }
        
        normType = 1;
        L1diff[dir] = computeNorm(diff, nullPtr,
                                  nRefFine, dx, 
                                  comps, normType);
        
        L1diffAMR[dir] = computeNorm(diffAMR, nullPtr,
                                     nRefFine, dx, 
                                     comps, normType);
        
        if (L1diff[dir] > diffTol)
          {
            pout() << "L1(diff)[" << dir << "] = " << L1diff[dir] << endl;
            status++;
          }
        if (L1diffAMR[dir] > diffTol)
          {
            pout() << "AMR vs. homogeneous CFInterp -- L1(diff)[" << dir << "] = " << L1diff[dir] << endl;
            status++;
          }
        
        
        normType = 2;
        L2diff[dir] = computeNorm(diff, nullPtr,
                                  nRefFine, dx, 
                                  comps, normType);

        L2diffAMR[dir] = computeNorm(diffAMR, nullPtr,
                                     nRefFine, dx, 
                                     comps, normType);
        
        if (L2diffAMR[dir] > diffTol)
          {
            pout() << "AMR vs. homogeneous CFInterp -- L2(diff)[" << dir << "] = " << L2diff[dir] << endl;
            status++;
          }
        
        if (L2diff[dir] > diffTol)
          {
            pout() << "L2(diff)[" << dir << "] = " << L2diff[dir] << endl;
            status++;
          }
        
      }

    pout() << "refined-level test, nGrids = " << numBoxes << ": L1(diff) = " << L1diff << endl;
    pout() << "                         "  << "  L2(diff) = " << L2diff << endl;
    pout() << "                         " << numBoxes << ": max(diff) = " << maxDiff << endl;
    
    // delete storage
    delete VTopPtr;
#ifdef CH_USE_PETSC
    delete petscSolverPtr;
#endif
    
  }
  
  if (status == 0)
    {
      pout() << "testPetsc: PASSED all tests!" << endl;
    }
  else
    {
      pout() << "testPetsc: FAILED at least one test!" << endl;
    }
  return status;
}


int main(int argc, char* argv[])
{
  int ierr = 0;
#ifdef CH_USE_PETSC
  // use same .petscrc as in the main code
  ierr = PetscInitialize(&argc, &argv,"../exec2D/.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
#endif

  // Scoping trick
  int retVal;
  {
    retVal = testPetsc();
  }// End scoping trick

  retVal += ierr;

#ifdef CH_USE_PETSC
  ierr = PetscFinalize(); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
#endif

  return retVal;
}
