#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;


#include "parstream.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "ConstitutiveRelation.H"
#include "LevelSigmaCS.H"
#include "FABView.H"
#include "CONSTANTS.H"
#include "IceConstants.H"

#include "defineLevelSigmaCS.H"

#include "computeNorm.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testConstitutive" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
//static bool verbose = false ;
static bool verbose = true ;

static Real tolerance = 2.0e-7;
static Real rateTolerance = 0.1;

#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif

// probtype enum
enum velTypeEnum{zeroVel = 0,
                 sinusoidalxVel,
                 sinusoidalyVel,
                 simpleSinusoidalVel,
                 transverseSinusoidalxVel,
                 transverseSinusoidalyVel,
                 fullSinusoidalxVel,
                 fullSinusoidalyVel,
                 sinusoidalVel,
                 num_probtype};

//int veltype = transverseSinusoidalyVel;
//int veltype = fullSinusoidalyVel;
int veltype = sinusoidalVel;

int thicknesstype = constantThickness;
int basalType = sinusoidalZb;

///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for( int i = 1 ; i < argc ; ++i )
    {
      if( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              //argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              //argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return;
}


// helper function
Real getConvergenceRate(Real fineErr, Real crseErr)
{
  Real rate = log(Abs(crseErr/fineErr))/log(2.0);
  return rate;
}



RealVect pointVelocity(const RealVect& x, const RealVect& domainSize)
{  
  RealVect vel;

  if (veltype == zeroVel)
    {
      vel = IntVect::Zero;
    }
  else if (veltype == sinusoidalxVel)
    {
      vel = RealVect::Zero;
      vel[0] = sin(2*Pi*x[0]/domainSize[0]);
    }
  else if (veltype == sinusoidalyVel)
    {
      vel = RealVect::Zero;
      vel[1] = sin(2*Pi*x[1]/domainSize[0]);
    }
  else if (veltype == simpleSinusoidalVel)
    {
      vel[0] = sin(2*Pi*x[0]/domainSize[0]);
      vel[1] = sin(2*Pi*x[1]/domainSize[0]);
    }
  else if (veltype == transverseSinusoidalxVel)
    {
      vel = RealVect::Zero;
      vel[0] = sin(2*Pi*x[1]/domainSize[1]);
    }
  else if (veltype == transverseSinusoidalyVel)
    {
      vel = RealVect::Zero;
      vel[1] = sin(2*Pi*x[0]/domainSize[0]);
    }
  else if (veltype == fullSinusoidalxVel)
    {
      vel = RealVect::Zero;
      vel[0] = D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                      *sin(2*Pi*x[1]/domainSize[1]),
                      *sin(2*Pi*x[2]/domainSize[2]));
    }
  else if (veltype == fullSinusoidalyVel)
    {
      vel = RealVect::Zero;
      vel[1] = D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                      *sin(2*Pi*x[1]/domainSize[1]),
                      *sin(2*Pi*x[2]/domainSize[2]));
    }
  else if (veltype == sinusoidalVel)
    {
      vel = RealVect::Unit;
      vel *= D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                    *sin(2*Pi*x[1]/domainSize[1]),
                    *sin(2*Pi*x[2]/domainSize[2]));
    }
  else
    {
      MayDay::Error("bad veltype");
    }
  return vel;
}


Real pointEpsSqr(const RealVect& x, const RealVect& domainSize)
{
  Real epsSqr = 1.23456e10; 

  Real dudx;
  Real dudy;
  Real dvdx;
  Real dvdy;

  if (veltype == zeroVel)
    {
      dudx = 0;
      dudy = 0;
      dvdx = 0;
      dvdy = 0;
    }
  else if (veltype == sinusoidalxVel)
    {
      dudx = (2*Pi/domainSize[0])*cos(2*Pi*x[0]/domainSize[0]);
      
      dudy = 0;
      
      dvdx = 0;
      
      dvdy = 0;
    }
  else if (veltype == sinusoidalyVel)
    {
      dudx = 0;
      
      dudy = 0;
      
      dvdx = 0;
      
      dvdy = (2*Pi/domainSize[1])*cos(2*Pi*x[1]/domainSize[0]);
    }
  else if (veltype == simpleSinusoidalVel)
    {
      dudx = (2*Pi/domainSize[0])*cos(2*Pi*x[0]/domainSize[0]);
      
      dudy = 0;
      
      dvdx = 0;
      
      dvdy = (2*Pi/domainSize[1])*cos(2*Pi*x[1]/domainSize[0]);
    }
  else if (veltype == transverseSinusoidalxVel)
    {
      dudx = 0;
      
      dudy = (2*Pi/domainSize[1])*cos(2*Pi*x[1]/domainSize[1]);
      
      dvdx = 0;
      
      dvdy = 0;
    }
  else if (veltype == transverseSinusoidalyVel)
    {
      dudx = 0;
      
      dudy = 0;
      
      dvdx = (2*Pi/domainSize[0])*cos(2*Pi*x[0]/domainSize[0]);
      
      dvdy = 0;
    }
  else if (veltype == fullSinusoidalxVel)
    {
      dudx = (2*Pi/domainSize[0])*D_TERM(cos(2*Pi*x[0]/domainSize[0]),
                                         *sin(2*Pi*x[1]/domainSize[1]),
                                         *sin(2*Pi*x[2]/domainSize[2]));
      
      dudy = (2*Pi/domainSize[1])*D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                                         *cos(2*Pi*x[1]/domainSize[1]),
                                         *sin(2*Pi*x[2]/domainSize[2]));
      
      dvdx = 0;
      
      dvdy = 0;
    }
  else if (veltype == fullSinusoidalyVel)
    {
      dudx = 0.0;
      
      dudy = 0.0;
      
      dvdx = (2*Pi/domainSize[0])*D_TERM(cos(2*Pi*x[0]/domainSize[0]),
                                         *sin(2*Pi*x[1]/domainSize[1]),
                                         *sin(2*Pi*x[2]/domainSize[2]));
      
      dvdy = (2*Pi/domainSize[1])*D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                                         *cos(2*Pi*x[1]/domainSize[1]),
                                         *sin(2*Pi*x[2]/domainSize[2]));
    }
  else if (veltype == sinusoidalVel)
    {
      dudx = (2*Pi/domainSize[0])*D_TERM(cos(2*Pi*x[0]/domainSize[0]),
                                         *sin(2*Pi*x[1]/domainSize[1]),
                                         *sin(2*Pi*x[2]/domainSize[2]));

      if (SpaceDim > 1)
        {      
          dudy = (2*Pi/domainSize[1])*D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                                             *cos(2*Pi*x[1]/domainSize[1]),
                                             *sin(2*Pi*x[2]/domainSize[2]));
          
          dvdx = (2*Pi/domainSize[0])*D_TERM(cos(2*Pi*x[0]/domainSize[0]),
                                             *sin(2*Pi*x[1]/domainSize[1]),
                                             *sin(2*Pi*x[2]/domainSize[2]));
          
          dvdy = (2*Pi/domainSize[1])*D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                                             *cos(2*Pi*x[1]/domainSize[1]),
                                             *sin(2*Pi*x[2]/domainSize[2]));
        }
    }
  else
    {
      MayDay::Error("bad veltype for pointEpsSqr");
    }
      
  if (SpaceDim == 1) 
    {
      epsSqr = dudx*dudx;
    }
  else if (SpaceDim == 2) 
    {
      //epsSqr = (dudx*dudx) +(dvdy*dvdy) +(dudx + dvdy)*(dudx + dvdy) +0.5*(dudy + dvdx)*(dudy + dvdx);
      
      //slc: i think this is the correct strain rate invariant (see e.g MacAyeal 1996)
      epsSqr = dudx*dudx + dvdy*dvdy + dudx*dvdy + 0.25 * (dudy+dvdx)*(dudy+dvdx);
      
    }
  else
    {
      MayDay::Error("exact strain invariant not implemented yet for DIM= SPACEDIM");
    }
 return epsSqr;
}



void initVel(LevelData<FArrayBox>& a_velocity, 
             LevelSigmaCS& a_sigmaCoords,
             const RealVect& a_domainSize)
{
  if (verbose)
    {
      if (veltype == zeroVel)
        {
          pout() << "using zero velocity field" << endl;
        }
      else if (veltype == sinusoidalxVel )
        {
          pout() << "using sinusoidal x-velocity field" << endl;
        }
      else if (veltype == sinusoidalyVel)
        {
          pout() << "using sinusoidal y-velocity field" << endl;
        }
      else if (veltype == simpleSinusoidalVel )
        {
          pout() << "using simple sinusoidal velocity field" << endl;
        }
      else if (veltype == transverseSinusoidalxVel )
        {
          pout() << "using tranverse-sinusoidal x-velocity field" << endl;
        }
      else if (veltype == transverseSinusoidalyVel )
        {
          pout() << "using tranverse-sinusoidal y-velocity field" << endl;
        }
      else if (veltype == fullSinusoidalxVel )
        {
          pout() << "using full-sinusoidal x-velocity field" << endl;
        }
      else if (veltype == fullSinusoidalyVel )
        {
          pout() << "using full-sinusoidal y-velocity field" << endl;
        }
      else if (veltype == sinusoidalVel )
        {
          pout() << "using full-sinusoidal velocity field" << endl;
        }
      else 
        {
          pout() << "undefined velocity type!" << endl;
        }
    } // end if verbose

  const RealVect& dx = a_sigmaCoords.dx();
  DataIterator dit=a_velocity.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& velfab = a_velocity[dit];
      
      BoxIterator ccBit(velfab.box());
      RealVect ccOffset = 0.5*RealVect::Unit;
      for (ccBit.begin(); ccBit.ok(); ++ccBit)
        {
          IntVect iv = ccBit();
          RealVect mappedLoc(iv);
          mappedLoc += ccOffset;
          mappedLoc *= dx;
          RealVect realLoc = a_sigmaCoords.realCoord(mappedLoc, dit());
          RealVect thisVel = pointVelocity(realLoc, a_domainSize);

          // note that velfab has 2 components regardless of
          // dimensionality, since it's the horizontal velocity          
          velfab(iv,0) = thisVel[0];
          if (SpaceDim > 1)
            {
              velfab(iv,1) = thisVel[1];
            }
        }
    } // end loop over boxes
}

void computeEpsSqrError(LevelData<FArrayBox>& a_ccError, 
                        LevelData<FluxBox>& a_fcError, 
                        LevelData<FArrayBox>& a_ccEpsSqr,
                        LevelData<FluxBox>& a_fcEpsSqr,
                        LevelSigmaCS& a_sigmaCoords,
                        const RealVect& a_domainSize)
{
  const RealVect& dx = a_sigmaCoords.dx();
  DataIterator dit=a_ccError.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // first do cell-centered data
      FArrayBox& ccErr = a_ccError[dit];
      FArrayBox& ccEps = a_ccEpsSqr[dit];
      BoxIterator ccBit(ccErr.box());
      RealVect ccOffset = 0.5*RealVect::Unit;
      for (ccBit.begin(); ccBit.ok(); ++ccBit)
        {
          IntVect iv = ccBit();
          RealVect mappedLoc(iv);
          mappedLoc += ccOffset;
          mappedLoc *= dx;
          RealVect realLoc = a_sigmaCoords.realCoord(mappedLoc, dit());

          Real exactEpsSqr = pointEpsSqr(realLoc, a_domainSize);
          //ccErr(iv,0) = ccEps(iv,0) - exactEpsSqr;
          ccErr(iv,0) = exactEpsSqr;
        }
      ccErr.minus(ccEps);
      
      // now do face-centered data
      for (int dir=0; dir<SpaceDim; dir++)
        {
          RealVect offset(0.5*RealVect::Unit);
          offset[dir] = 0.0;

          FArrayBox& fcerr = a_fcError[dit][dir];
          FArrayBox& fcEpsSqr = a_fcEpsSqr[dit][dir];
          BoxIterator fcBit(fcEpsSqr.box());
          for (fcBit.begin(); fcBit.ok(); ++fcBit)
            {
              IntVect iv = fcBit();
              RealVect mappedLoc(iv);
              mappedLoc += offset;
              mappedLoc *= dx;
              RealVect realLoc = a_sigmaCoords.realCoord(mappedLoc, dit());

              Real exactEpsSqr = pointEpsSqr(realLoc, a_domainSize);
              fcerr(iv,0) =  exactEpsSqr;                
            }
          fcerr.minus(fcEpsSqr);
        } // end loop over face directions
    } // end loop over boxes


}


int
testConstitutive();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testConstitutive();

  if( status == 0 )
    pout() << indent << "ConstitutiveRelation test" << " passed." << endl ;
  else
    pout() << indent << "ConstitutiveRelation test" << " failed with return code "
         << status << endl ;



  if( status == 0 )
    pout() << indent << pgmname << " passed." << endl ;
  else
    pout() << indent << pgmname << " failed with return code " << status << endl ;


#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int
testConstitutive() 
{
  int status = 0;

  
  IntVect numCells(32*IntVect::Unit);
  if (CH_SPACEDIM == 3)
    {
      numCells[0] = 8;
    }
  IntVect hiVect = numCells - IntVect::Unit;
  Box domainBox(IntVect::Zero, hiVect);
  ProblemDomain entireDomain(domainBox);
  for (int dir=1; dir<SpaceDim; dir++)
    {
      entireDomain.setPeriodic(dir,true);
    }

  RealVect domainSize(10.0*RealVect::Unit);
  if (CH_SPACEDIM == 3)
    {
      domainSize[0] = 1.0;
    }
  RealVect dx = domainSize;
  dx /= numCells;

  //int maxBoxSize = 8;
  int maxBoxSize = 32;
  Vector<Box> gridBoxes;
  domainSplit(domainBox, gridBoxes, maxBoxSize); 
  
  Vector<int> procAssign(gridBoxes.size(), 0);
  LoadBalance(procAssign, gridBoxes);

  Real Aval = 1.0e-18; // typical value
  Real time_scale = 1.0; // |u| in m/a
  // test zero-velocity, constant temperature
  {
   
    
    // allocate storage
    DisjointBoxLayout grids(gridBoxes, procAssign, entireDomain);

    IntVect ghostVect = IntVect::Unit;    
    LevelSigmaCS sigmaCoords(grids, dx, ghostVect);
    // now loop over grids and define coordSys on each grid
    
    RealVect basalSlope = RealVect::Zero;
    defineLevelSigmaCS(sigmaCoords, domainSize, thicknesstype, basalType,
                       basalSlope);

    LevelData<FArrayBox> A(grids, 1, ghostVect);
    LevelData<FluxBox> faceA(grids, 1, ghostVect);
    int numVelComps = min(SpaceDim, 2);
    LevelData<FArrayBox> horizontalVel(grids, numVelComps, ghostVect);
    // for now, no coarser level
    LevelData<FArrayBox>* crseVelPtr = NULL;
    int nRefCrse = -1;

    DataIterator dit = grids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        A[dit].setVal(Aval);
        faceA[dit].setVal(Aval);
        horizontalVel[dit].setVal(0);
      }

    // cell- and face-centered mu:
    LevelData<FArrayBox> cellMu(grids, 1);
    LevelData<FluxBox> faceMu(grids,1);

    // ConstitutiveRelation object
    GlensFlowRelation constitutiveRelation;
    ArrheniusRateFactor rateFactor(time_scale);
    Real epsSqrZero = 1e-30;
    constitutiveRelation.setParameters(3.0 , epsSqrZero, 0.0);


    Real exponent = -1.0/3.0;
    Real exactMu0 = 0.5*pow(Aval,exponent);
    // zero velocity field -- this is the effect of the
    // small parameter epsSqr_0
     Real exactMu = exactMu0*pow(epsSqrZero,exponent);
        
    IntVect muGhost = IntVect::Zero;
    // cell-centered Mu
    constitutiveRelation.computeMu(cellMu,
                                   horizontalVel, time_scale,
                                   crseVelPtr,
                                   nRefCrse,
                                   A,
                                   sigmaCoords,
				   entireDomain,
                                   muGhost);

    // face-centered Mu
    constitutiveRelation.computeFaceMu(faceMu,
                                       horizontalVel, time_scale,
                                       crseVelPtr,
                                       nRefCrse,
                                       faceA,
                                       sigmaCoords,
				       entireDomain,
                                       muGhost);

    
    // compute errors
    Real maxCellErr = 0.0;
    RealVect maxFaceErr = RealVect::Zero;
    for (dit.begin(); dit.ok(); ++dit)
      {
        cellMu[dit] -= exactMu;
        // relative error
        cellMu[dit] /= exactMu;
        Real thisMax = cellMu[dit].norm(0);
        if (thisMax > maxCellErr) maxCellErr = thisMax;
        for (int dir=0; dir<SpaceDim; dir++)
          {
            faceMu[dit][dir] -= exactMu;
            faceMu[dit][dir] /= exactMu;
            thisMax = faceMu[dit][dir].norm(0);
            if (thisMax > maxFaceErr[dir]) maxFaceErr[dir] = thisMax;
          }
      }
    
    pout() << "constant Temp and zero velocity test: " << endl;
    pout() << "      max relative cell err = " << maxCellErr 
           << ",   max relative face err = ";
    for (int dir=0; dir<SpaceDim; dir++)
      {
        pout() << maxFaceErr[dir] << "   ";
      }
    pout() << endl;
    pout() << endl;                                           

    // test for acceptable error
    if (maxCellErr > tolerance) 
      {
        pout() << "  CELL ERR FAILS!" << endl;
        status += 1;
      }
    for (int dir=0; dir<SpaceDim; dir++)
      if (maxFaceErr[dir] > tolerance)
        {
          pout() << "  FACE ERR ON FACE " << dir << " FAILS!" << endl;
          status += 10;
        }
  }
  
  
  // do mesh refinement testing 

  // storage for epsSqr errors
  Real L1ErrSave;
  Real L2ErrSave;
  Real MaxErrSave;

  RealVect maxFCErrSave(RealVect::Zero);

  int numRefine = 4;
  for (int refinement = 0; refinement <numRefine; refinement++)
    {
      if (verbose)  pout() << "refinement " << refinement << endl;

      // allocate storage
      DisjointBoxLayout grids(gridBoxes, procAssign, entireDomain);

      IntVect ghostVect = IntVect::Unit;      
      LevelSigmaCS sigmaCoords(grids, dx, ghostVect);

      RealVect basalSlope = RealVect::Zero;
      defineLevelSigmaCS(sigmaCoords, domainSize, thicknesstype, basalType,
                         basalSlope);
      
      // set up data -- note that velocity is 2 components
      // in 2- and 3-D, since it's only the horizontal 
      // velocity
      int numVelComp = min(SpaceDim, 2);
      LevelData<FArrayBox> vel(grids,numVelComp, ghostVect);
      LevelData<FArrayBox> theta(grids, 1, ghostVect);
      // no coarse-level vel
      LevelData<FArrayBox>* crseVelPtr = NULL;
      int nRefCrse = -1;

      initVel(vel, sigmaCoords, domainSize);

      // now compute gradients
      Interval comps(0,0);
      Interval derivDirections(0,SpaceDim-1);
      LevelData<FArrayBox> ccEpsSqr(grids, 1);
      LevelData<FluxBox> fcEpsSqr(grids, 1);

      // constitutive relation object
      GlensFlowRelation constitutiveRelation;
      ArrheniusRateFactor rateFactor(time_scale);
      constitutiveRelation.setParameters(3.0 , 1.0e-30, 0.0);
      
      IntVect epsSqrGhost = IntVect::Zero;
      constitutiveRelation.computeStrainRateInvariant(ccEpsSqr,
                                                      vel,
                                                      crseVelPtr,
                                                      nRefCrse,
                                                      sigmaCoords,
                                                      epsSqrGhost);
      
      // face-centered Epsilon-squared
      constitutiveRelation.computeStrainRateInvariantFace(fcEpsSqr,
                                                          vel,
                                                          crseVelPtr,
                                                          nRefCrse,
                                                          sigmaCoords,
                                                          epsSqrGhost);
      
      // now evaluate error
      LevelData<FArrayBox> ccError(grids, 1);
      LevelData<FluxBox> fcError(grids, 1);
      computeEpsSqrError(ccError, fcError, ccEpsSqr, fcEpsSqr, 
                         sigmaCoords, domainSize);

      // now compute norms of errors
      Real L1Err = 0;
      Real L2Err = 0;
      Real MaxErr = 0;
      
      RealVect maxFCErr(RealVect::Zero);
      
      DisjointBoxLayout* fakeFineGrids = NULL;
      int nRefFine = -1;
      Interval ErrorComps(0,0);
      L1Err = computeNorm(ccError, fakeFineGrids,
                          nRefFine, dx[0], ErrorComps,
                          1);


      L2Err = computeNorm(ccError, fakeFineGrids,
                          nRefFine, dx[0], ErrorComps,
                          2);

      
      MaxErr = computeNorm(ccError, fakeFineGrids,
                           nRefFine, dx[0], ErrorComps,
                           0);
          
      // loop over faces and evaluate max(err)
      RealVect localMaxErr = RealVect::Zero;
      DataIterator dit = fcError.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& faceErr = fcError[dit];
          for (int faceDir = 0; faceDir<SpaceDim; faceDir++)
            {
              Real thisNorm = faceErr[faceDir].norm(0,0, 1);
              if (thisNorm > localMaxErr[faceDir]) localMaxErr[faceDir] = thisNorm;

            } // end loop over face directions
        } // end loop over grids

      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          maxFCErr[faceDir] = localMaxErr[faceDir];
        }
    
      // write errors
      if (verbose)
        {
          pout() << "CC Strain invariant Errors: " << endl;
          

          pout() << "L1 = " << L1Err
                 << ",   L2 = " << L2Err
                 << ",   max = " << MaxErr << endl;
          
          Real L1Rate, L2Rate, MaxRate;

          if (refinement > 0)
            {
              L1Rate = getConvergenceRate(L1Err, L1ErrSave);
              L2Rate = getConvergenceRate(L2Err, L2ErrSave);
              MaxRate = getConvergenceRate(MaxErr, MaxErrSave); 
              // compute convergence rates
              pout() << "   rate: L1 = " << L1Rate 
                     << "     L2 = " << L2Rate
                     << "    max = " << MaxRate;

              if (L1Rate < 2.0-rateTolerance) 
                {
                  pout() << endl;
                  pout() << "L1 rate FAILS!" << endl;
                  status += 100;
                }

              if (L2Rate < 2.0-rateTolerance) 
                {
                  pout() << endl;
                  pout() << "L2 rate FAILS!" << endl;
                  status += 100;
                }


              if (MaxRate < 2.0-rateTolerance) 
                {
                  pout() << endl;
                  pout() << "MaxNorm rate FAILS!" << endl;
                  status += 100;
                }            
            }
          pout() << endl;
          
          pout() << "FC Strain invariant Errors: " << endl;
          for (int faceDir=0; faceDir<SpaceDim; faceDir++)
            {
              pout() << "   Face dir = " << faceDir << ":" 
                     << "  max(Err) = " << maxFCErr[faceDir];
            }
          pout() << endl;

          RealVect maxFCrate;
          if (refinement > 0)
            {
              for (int faceDir = 0; faceDir<SpaceDim; faceDir++)
                {
                  maxFCrate[faceDir] = getConvergenceRate(maxFCErr[faceDir],
                                                          maxFCErrSave[faceDir]);
                }
              
              pout() << "           rate: " 
                D_TERM(<< " max(xFaceErr) = " << maxFCrate[0],
                       << "  max(yFaceErr) = " << maxFCrate[1],
                       << "  max(zFaceErr) = " << maxFCrate[2])
                     << endl;

              // test for acceptable convergence
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  if (maxFCrate[dir] < 2.0-rateTolerance) 
                    {
                      pout() << "  " << dir
                             << "-direction maxnorm convergence test FAILS!"
                             << endl;
                      status += 1000;
                    }
                }
            }        
        } // end if verbose
      
      // save errors
      L1ErrSave = L1Err;
      L2ErrSave = L2Err;
      MaxErrSave = MaxErr;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          maxFCErrSave[dir] = maxFCErr[dir];
        }
      

      bool writePlotFiles = true;
      if (writePlotFiles)
        {
          // 2 velocity components, epsSqr, err(epsSqr)
          int nPlotComp = 4;
          LevelData<FArrayBox> plotData(grids, nPlotComp);
          for (dit.begin(); dit.ok(); ++dit)
            {
              //const SigmaCS& thisCS = sigmaCoords[dit];
              //const FArrayBox& DeltaFactors = thisCS.deltaFactors();
              const FArrayBox& thisVel = vel[dit];
              //plotData[dit].copy(DeltaFactors,0,0,SpaceDim-1);
              plotData[dit].copy(thisVel,0,0,thisVel.nComp());
              const FArrayBox& thisCCeps = ccEpsSqr[dit];
              plotData[dit].copy(thisCCeps,0,2,1);
              const FArrayBox& thisErr = ccError[dit];
              plotData[dit].copy(thisErr, 0, 3, 1);
            }
          
          string fileRoot("testConstitutive.");
          // use unigrid version here
          dit.reset();
          WriteSigmaMappedUGHDF5(fileRoot, grids, plotData, sigmaCoords,
                                 domainBox);
        } // end if writing plot files
      
      // set up for next refinement -- refine by factor of 2
      for (int boxNo=0; boxNo<gridBoxes.size(); boxNo++)
        {
          gridBoxes[boxNo].refine(2);
        }
      entireDomain.refine(2);
      domainBox.refine(2);
      dx *= 0.5;

      
    } // end loop over refinements
      
  return status;
      
}
