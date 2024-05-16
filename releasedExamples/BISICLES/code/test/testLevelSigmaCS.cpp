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
#include "LevelSigmaCS.H"
#include "FABView.H"
#include "CONSTANTS.H"
#include "LevelMappedDerivatives.H"
#include "computeNorm.H"
#include "defineLevelSigmaCS.H"
#include "VieliPayneIBC.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testLevelSigmaCS" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
//static bool verbose = false ;
static bool verbose = true ;

#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif

static Real rateTolerance = 0.02;

// probtype enum
enum phitypeEnum{sinusoidalPhi = 0,
                 linearInZ,
                 num_probtype};

int phitype = sinusoidalPhi;
//int phitype = linearInZ;

//int heightType = constantHeight;
int heightType = doubleZb;
//int heightType = circle;

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

Real getConvergenceRate(Real fineErr, Real crseErr)
{
  Real rate = log(Abs(crseErr/fineErr))/log(2.0);
  return rate;
}

Real pointVal(const RealVect& x, const RealVect& domainSize)
{  
  Real val;

  if (phitype == sinusoidalPhi)
    {
      val = D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                   *sin(2*Pi*x[1]/domainSize[1]),
                   *sin(2*Pi*x[2]/domainSize[2]));
    }
  else if (phitype ==  linearInZ)
    {
      val = x[0];
    }
  else
    {
      MayDay::Error("bad phitype");
    }
  return val;
}


RealVect pointGrad(const RealVect& x, const RealVect& domainSize)
{
  RealVect grad; 

  for (int derivDir=0; derivDir<SpaceDim; derivDir++)
    {
      grad[derivDir] = 2*Pi/domainSize[derivDir];
      for (int dir=0; dir<SpaceDim; dir++)
        { 
          if (dir == derivDir) 
            {
              grad[derivDir] *=  cos(2*Pi*x[dir]/domainSize[dir]);
            }
          else 
            {
              grad[derivDir] *= sin(2*Pi*x[dir]/domainSize[dir]);
            }
        } // end loop over directions
    } // end loop over deriv directions
  return grad;
}



void initData(LevelData<FArrayBox>& a_ccData, 
              LevelData<FluxBox>& a_fcData, 
              LevelSigmaCS& a_sigmaCoords,
              const RealVect& a_domainSize)
{
  const RealVect& dx = a_sigmaCoords.dx();
  DataIterator dit=a_ccData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {

     
      // first do cell-centered data
      FArrayBox& ccfab = a_ccData[dit];
      BoxIterator ccBit(ccfab.box());
      RealVect ccOffset = 0.5*RealVect::Unit;
      for (ccBit.begin(); ccBit.ok(); ++ccBit)
        {
          IntVect iv = ccBit();
          RealVect mappedLoc(iv);
          mappedLoc += ccOffset;
          mappedLoc *= dx;
          RealVect realLoc = a_sigmaCoords.realCoord(mappedLoc,dit());

          ccfab(iv,0) = pointVal(realLoc, a_domainSize);
        }

      // now do face-centered data
      for (int dir=0; dir<SpaceDim; dir++)
        {
          RealVect offset(0.5*RealVect::Unit);
          offset[dir] = 0.0;

          FArrayBox& fcfab = a_fcData[dit][dir];
          BoxIterator fcBit(fcfab.box());
          for (fcBit.begin(); fcBit.ok(); ++fcBit)
            {
              IntVect iv = fcBit();
              RealVect mappedLoc(iv);
              mappedLoc += offset;
              mappedLoc *= dx;
              RealVect realLoc = a_sigmaCoords.realCoord(mappedLoc,dit());

              fcfab(iv,0) = pointVal(realLoc, a_domainSize);
            }
        } // end loop over face directions
    } // end loop over boxes
      
}

void computeGradError(LevelData<FArrayBox>& a_ccError, 
                      LevelData<FluxBox>& a_fcError, 
                      LevelData<FluxBox>& a_fcError2, 
                      LevelData<FArrayBox>& a_ccGrad,
                      LevelData<FluxBox>& a_fcGrad,
                      LevelData<FluxBox>& a_fcGrad2,
                      LevelSigmaCS& a_sigmaCoords,
                      const RealVect& a_domainSize)
{
  const RealVect& dx = a_sigmaCoords.dx();
  DataIterator dit=a_ccError.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // first do cell-centered data
      FArrayBox& ccErr = a_ccError[dit];
      FArrayBox& ccGrad = a_ccGrad[dit];
      BoxIterator ccBit(ccErr.box());
      RealVect ccOffset = 0.5*RealVect::Unit;
      for (ccBit.begin(); ccBit.ok(); ++ccBit)
        {
          IntVect iv = ccBit();
          RealVect mappedLoc(iv);
          mappedLoc += ccOffset;
          mappedLoc *= dx;
          RealVect realLoc = a_sigmaCoords.realCoord(mappedLoc, dit());

          RealVect exactGrad = pointGrad(realLoc, a_domainSize);
          for (int gradDir=0; gradDir<SpaceDim; gradDir++)
            {
              ccErr(iv,gradDir) = ccGrad(iv,gradDir) - exactGrad[gradDir];
            }
        }

      // now do face-centered data
      for (int dir=0; dir<SpaceDim; dir++)
        {
          RealVect offset(0.5*RealVect::Unit);
          offset[dir] = 0.0;

          FArrayBox& fcerr = a_fcError[dit][dir];
          FArrayBox& fcerr2 = a_fcError2[dit][dir];
          FArrayBox& fcgrad = a_fcGrad[dit][dir];
          FArrayBox& fcgrad2 = a_fcGrad2[dit][dir];
          BoxIterator fcBit(fcgrad.box());
          for (fcBit.begin(); fcBit.ok(); ++fcBit)
            {
              IntVect iv = fcBit();
              RealVect mappedLoc(iv);
              mappedLoc += offset;
              mappedLoc *= dx;
              RealVect realLoc = a_sigmaCoords.realCoord(mappedLoc,
                                                         dit());

              RealVect exactGrad = pointGrad(realLoc, a_domainSize);
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  fcerr(iv,dir) = fcgrad(iv,dir) - exactGrad[dir];
                  fcerr2(iv,dir) = fcgrad2(iv,dir) - exactGrad[dir];
                }
            }
        } // end loop over face directions
    } // end loop over boxes


}

int
testLevelSigmaCS();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testLevelSigmaCS();

  if( status == 0 )
    pout() << indent << "LevelSigmaCS test" << " passed." << endl ;
  else
    pout() << indent << "LevelSigmaCS test" << " failed with return code "
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
testLevelSigmaCS() 
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

  {
    // quick test of floating mask
    Real thickness = 500;
    RealVect basalSlope = RealVect::Zero;
    basalSlope[0] = -0.005;
    Real originElevation = -250.0;
    Real seaLevel = 0.0;

    VieliPayneIBC VPibc;
    VPibc.setParameters(thickness, basalSlope, 
                        originElevation, seaLevel);
    

    RealVect vpDomainSize = 100000.0*RealVect::Unit;
    RealVect vpDx = vpDomainSize;
    vpDx /= numCells;
      
    // vieli-Payne only periodic in y-direction
    ProblemDomain vpdomain(entireDomain);
    vpdomain.setPeriodic(0,false);
    VPibc.define(vpdomain, dx[0]);
    
    // allocate storage
    DisjointBoxLayout grids(gridBoxes, procAssign, entireDomain);
    LevelSigmaCS coordsys(grids, vpDx);
  
    Real time = 0.0;
    // crse level info -- empty for now
    LevelSigmaCS* crsePtr = NULL;
    int nRefCrse = -1;
    VPibc.initializeIceGeometry(coordsys, vpDx, vpDomainSize, time,
                                crsePtr, nRefCrse);

    // need to make this a pass-fail test at some point...
  }

  // do mesh refinement testing 

  // storage for err(grad(phi))
  RealVect L1GradErrSave(RealVect::Zero);
  RealVect L2GradErrSave(RealVect::Zero);
  RealVect MaxGradErrSave(RealVect::Zero);

  Vector<RealVect> maxFCGradErrSave(SpaceDim, RealVect::Zero);
  Vector<RealVect> maxFCGradErrSave2(SpaceDim, RealVect::Zero);

  int numRefine = 4;
  for (int refinement = 0; refinement <numRefine; refinement++)
    {
      if (verbose)  pout() << "refinement " << refinement << endl;

      // allocate storage
      DisjointBoxLayout grids(gridBoxes, procAssign, entireDomain);

      IntVect ghostVect = IntVect::Unit;
      
      LevelSigmaCS sigmaCoords(grids, dx, ghostVect);
      // now loop over grids and define coordSys on each grid
      
      RealVect basalSlope = RealVect::Zero;
      Real thicknessScale = 1.0;

      defineLevelSigmaCS(sigmaCoords, domainSize, heightType, basalType,
                         basalSlope, thicknessScale);

      // set up data
      IntVect ccGhostVect(ghostVect);
      ccGhostVect += IntVect::Unit;
      LevelData<FArrayBox> ccData(grids, 1, ccGhostVect);
      LevelData<FluxBox> fcData(grids, 1, ghostVect);

      initData(ccData, fcData, sigmaCoords, domainSize);

      // now compute gradients
      Interval comps(0,0);
      Interval derivDirections(0,SpaceDim-1);
      LevelData<FArrayBox> ccGrad(grids, SpaceDim, ghostVect);
      LevelData<FluxBox> fcGrad(grids, SpaceDim);
      LevelData<FluxBox> fcGrad2(grids, SpaceDim);

      computeCCDerivatives(ccGrad,
                           ccData,
                           sigmaCoords,
                           comps,
                           derivDirections,
                           ghostVect);

      // fc derivatives of ccData next
      computeFCDerivatives(fcGrad,
                           ccData,
                           sigmaCoords,
                           comps, derivDirections);
      // fc derivatives of ccData next, using alternate function
      computeFCDerivatives(fcGrad2,
                           ccData,
                           ccGrad,
                           sigmaCoords,
                           comps, derivDirections);                         
      
      
      // now evaluate error
      LevelData<FArrayBox> ccError(grids, SpaceDim);
      LevelData<FluxBox> fcError(grids, SpaceDim);
      LevelData<FluxBox> fcError2(grids, SpaceDim);
      computeGradError(ccError, fcError, fcError2,
                       ccGrad, fcGrad, fcGrad2,
                       sigmaCoords, domainSize);

      // now compute norms of errors
      RealVect L1GradErr(RealVect::Zero);
      RealVect L2GradErr(RealVect::Zero);
      RealVect MaxGradErr(RealVect::Zero);
      
      Vector<RealVect> maxFCGradErr(SpaceDim, RealVect::Zero);
      Vector<RealVect> maxFCGradErr2(SpaceDim, RealVect::Zero);
      
      DisjointBoxLayout* fakeFineGrids = NULL;
      int nRefFine = -1;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          Interval ErrorComps(dir,dir);
          L1GradErr[dir] = computeNorm(ccError, fakeFineGrids,
                                       nRefFine, dx[0], ErrorComps,
                                       1);


          L2GradErr[dir] = computeNorm(ccError, fakeFineGrids,
                                       nRefFine, dx[0], ErrorComps,
                                       2);


          MaxGradErr[dir] = computeNorm(ccError, fakeFineGrids,
                                        nRefFine, dx[0], ErrorComps,
                                        0);

          // loop over faces and evaluate max(err)
          RealVect maxErr = RealVect::Zero;
          RealVect maxErr2 = RealVect::Zero;
          DataIterator dit = grids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              FluxBox& faceErr = fcError[dit];
              FluxBox& faceErr2 = fcError2[dit];
              for (int faceDir = 0; faceDir<SpaceDim; faceDir++)
                {
                  Real thisNorm = faceErr[faceDir].norm(0,dir, 1);
                  Real thisNorm2 = faceErr2[faceDir].norm(0,dir, 1);
                  if (thisNorm > maxErr[faceDir]) maxErr[faceDir] = thisNorm;
                  if (thisNorm2 > maxErr2[faceDir]) maxErr2[faceDir] = thisNorm2;

                } // end loop over face directions
            } // end loop over grids

          for (int faceDir=0; faceDir<SpaceDim; faceDir++)
            {
              maxFCGradErr[faceDir][dir] = maxErr[faceDir];
              maxFCGradErr2[faceDir][dir] = maxErr2[faceDir];
            }
        } // end loop over grad directions
      

      // write errors
      if (verbose)
        {
          pout() << "CCGradient Errors: " << endl;
          
          for (int dir=0; dir<SpaceDim; dir++)
            {
              pout() << "Dir = " << dir << ": L1 = " << L1GradErr[dir] 
                     << ",   L2 = " << L2GradErr[dir] 
                     << ",   max = " << MaxGradErr[dir] << endl;
              
              if (refinement > 0)
                {
                  Real L1rate = getConvergenceRate(L1GradErr[dir], L1GradErrSave[dir]);
                  Real L2rate = getConvergenceRate(L2GradErr[dir], L2GradErrSave[dir]);
                  Real maxrate = getConvergenceRate(MaxGradErr[dir], MaxGradErrSave[dir]);
                 
                  // compute convergence rates
                  pout() << "   rate: L1 = " << L1rate
                         << "     L2 = " << L2rate
                         << "    max = " << maxrate << endl; 

                  if (L1rate < 2.0-rateTolerance)
                    {
                      pout() << " L1 Cell-centered gradient convergence FAILS!" 
                             << endl;
                      status += 1;
                    }

                  if (L2rate < 2.0-rateTolerance)
                    {
                      pout() << " L2 Cell-centered gradient convergence FAILS!" 
                             << endl;
                      status += 1;
                    }


                  if (maxrate < 2.0-rateTolerance)
                    {
                      pout() << " MaxNorm Cell-centered gradient convergence FAILS!" 
                             << endl;
                      status += 1;
                    }
                }
              pout() << endl;
              

            } // end loop over gradient directions
          
          pout() << "FCGradient Errors: " << endl;
          for (int faceDir=0; faceDir<SpaceDim; faceDir++)
            {
              RealVect& faceErr = maxFCGradErr[faceDir];
              pout() << "   Face dir = " << faceDir << ":" 
                D_TERM(<< "  max(xErr) = " << faceErr[0],
                       << "  max(yErr) = " << faceErr[1],
                       << "  max(zErr) = " << faceErr[2]) << endl;
              if (refinement > 0)
                {
                  RealVect faceRate;
                  for (int gradDir=0; gradDir<SpaceDim; gradDir++)
                    {
                      faceRate[gradDir] = getConvergenceRate(faceErr[gradDir],
                                                             maxFCGradErrSave[faceDir][gradDir]);
                    }
                  pout() << "           rate: " 
                    D_TERM(<< " max(xErr) = " << faceRate[0],
                           << "  max(yErr) = " << faceRate[1],
                           << "  max(zErr) = " << faceRate[2]) << endl;

                  for (int gradDir=0; gradDir<SpaceDim; gradDir++)
                    {
                      if (faceRate[gradDir] < 2.0-rateTolerance)
                        {
                          pout() << "  face-centered gradient convergence test FAILS! -- gradDir = " 
                                 << gradDir << ", faceDir = " 
                                 << faceDir << endl;
                          status += 10;
                        }
                    } // end loop over grad directions

                }
            } // end loop over face directions
        

          pout() << "FCGradient2 Errors: " << endl;
          for (int faceDir=0; faceDir<SpaceDim; faceDir++)
            {
              RealVect& faceErr2 = maxFCGradErr2[faceDir];
              pout() << "   Face dir = " << faceDir << ":" 
                D_TERM(<< "  max(xErr) = " << faceErr2[0],
                       << "  max(yErr) = " << faceErr2[1],
                       << "  max(zErr) = " << faceErr2[2]) << endl;
              if (refinement > 0)
                {
                  RealVect faceRate2;
                  for (int gradDir=0; gradDir<SpaceDim; gradDir++)
                    {
                      faceRate2[gradDir] = getConvergenceRate(faceErr2[gradDir],
                                                             maxFCGradErrSave2[faceDir][gradDir]);
                    }
                  pout() << "           rate: " 
                    D_TERM(<< " max(xErr) = " << faceRate2[0],
                           << "  max(yErr) = " << faceRate2[1],
                           << "  max(zErr) = " << faceRate2[2]) << endl;

                  for (int gradDir=0; gradDir<SpaceDim; gradDir++)
                    {
                      if (faceRate2[gradDir] < 2.0-rateTolerance)
                        {
                          pout() << "  face-centered gradient convergence test FAILS! -- gradDir = " 
                                 << gradDir << ", faceDir = " 
                                 << faceDir << endl;
                          status += 10;
                        }
                    } // end loop over grad directions

                }
            }
        } // end if verbose
      
      // save errors
      L1GradErrSave = L1GradErr;
      L2GradErrSave = L2GradErr;
      MaxGradErrSave = MaxGradErr;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          maxFCGradErrSave[dir] = maxFCGradErr[dir];
          maxFCGradErrSave2[dir] = maxFCGradErr2[dir];
        }
      

#if 0
      bool writePlotFiles = true;
      if (writePlotFiles)
        {
          int nPlotComp = SpaceDim + 1;
          LevelData<FArrayBox> plotData(grids, nPlotComp);
          for (dit.begin(); dit.ok(); ++dit)
            {
              const SigmaCS& thisCS = sigmaCoords[dit];
              //const FArrayBox& DeltaFactors = thisCS.deltaFactors();
              const FArrayBox& thisPhi = ccData[dit];
              //plotData[dit].copy(DeltaFactors,0,0,SpaceDim-1);
              plotData[dit].copy(thisPhi,0,0,1);
              const FArrayBox& thisCCgrad = ccGrad[dit];
              plotData[dit].copy(thisCCgrad,0,1,SpaceDim);
            }
          
          string fileRoot("testSigmaCS.");
          // use unigrid version here
          dit.reset();
          WriteSigmaMappedUGHDF5(fileRoot, grids, plotData, sigmaCoords,
                            domainBox);
        } // end if writing plot files
#endif

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
