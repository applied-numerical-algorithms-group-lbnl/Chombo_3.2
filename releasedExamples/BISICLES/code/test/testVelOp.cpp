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
//using std::cout;
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
#include "SigmaCS.H"
#include "SSAVelocityOp.H"
#include "BCFunc.H"
#include "FABView.H"
#include "CONSTANTS.H"
#include "defineSigmaCS.H"

#include "computeNorm.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testVelOp" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
//static bool verbose = false ;
static bool verbose = true ;

static Real tolerance = 1.0e-7;
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

//int veltype = zeroVel;
//int veltype = sinusoidalxVel;
//int veltype = simpleSinusoidalVel;
//int veltype = transverseSinusoidalyVel;
//int veltype = fullSinusoidalyVel;
int veltype = sinusoidalVel;


//int thicknessType = constantThickness;
int thicknessType = doubleZb;
//int thicknessType = circle;
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
      
      vel[0] *= D_TERM(sin(2*Pi*x[0]/domainSize[0]),
                       *sin(2*Pi*x[1]/domainSize[1]),
                       *sin(2*Pi*x[2]/domainSize[2]));

      vel[1] *= D_TERM(cos(2*Pi*x[0]/domainSize[0]),
                       *cos(2*Pi*x[1]/domainSize[1]),
                       *cos(2*Pi*x[2]/domainSize[2]));
        
    }
  else
    {
      MayDay::Error("bad veltype");
    }
  return vel;
}

void initData(LevelData<FArrayBox>& a_velocity, 
              LevelData<FluxBox>& a_faceMu,
              LevelData<SigmaCS>& a_sigmaCoords,
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

  DataIterator dit=a_velocity.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const SigmaCS& thisCS = a_sigmaCoords[dit];
      const RealVect& dx = thisCS.dx();
     
      FArrayBox& velfab = a_velocity[dit];
      FluxBox& muFab = a_faceMu[dit];

      // for now, set mu = 1.0
      muFab.setVal(1.0);

      BoxIterator ccBit(velfab.box());
      RealVect ccOffset = 0.5*RealVect::Unit;
      for (ccBit.begin(); ccBit.ok(); ++ccBit)
        {
          IntVect iv = ccBit();
          RealVect mappedLoc(iv);
          mappedLoc += ccOffset;
          mappedLoc *= dx;
          RealVect realLoc = thisCS.realCoord(mappedLoc);
          RealVect thisVel = pointVelocity(realLoc, a_domainSize);

          // note that velfab has 2 components regardless of
          // dimensionality, since it's the horizontal velocity          
          velfab(iv,0) = thisVel[0];
          velfab(iv,1) = thisVel[1];
        }
    } // end loop over boxes
}

RealVect
pointExactOp(const RealVect& x, 
             const Real& H,
             const RealVect& a_domainSize)
{
  RealVect val;
  // do this in a way which is easy to extend...
  // first derivatives
  Real dudx, dudy, dvdx, dvdy, dHmuDx, dHmuDy;
  // second derivatives
  Real dduDxSqr, dduDxDy, dduDySqr, ddvDxSqr, ddvDxDy, ddvDySqr;

  if (thicknessType == constantThickness)
    {
      dHmuDx = 0;
      dHmuDy = 0;
    }
  else if (thicknessType == doubleZb)
    {
      dHmuDx = 4.0*Pi*cos(4*Pi*x[0]/a_domainSize[0])*sin(4.0*Pi*x[1]/a_domainSize[0])/a_domainSize[1];
      dHmuDy = 4.0*Pi*sin(4*Pi*x[0]/a_domainSize[0])*cos(4.0*Pi*x[1]/a_domainSize[0])/a_domainSize[1];

    }
  else
    {
      MayDay::Error("pointExactOp not implemented for this thicknessType");
    }

  
  if (veltype == zeroVel)
    {
      // first derivatives
      dudx = 0;
      dudy = 0;
      dvdx = 0;
      dvdy = 0;
      // second derivatives
      dduDxSqr = 0; 
      dduDxDy = 0; 
      dduDySqr = 0; 
      ddvDxSqr = 0; 
      ddvDxDy = 0; 
      ddvDySqr = 0;
    }
  else if (veltype == sinusoidalxVel )
    {
      // first derivatives
      dudx = 2.0*Pi*cos(2*Pi*x[0]/a_domainSize[0])/a_domainSize[0];
      dudy = 0;
      dvdx = 0;
      dvdy = 0.0;
      // second derivatives
      dduDxSqr = -4*Pi*Pi*sin(2*Pi*x[0]/a_domainSize[0])/(a_domainSize[0]*a_domainSize[0]);
      dduDxDy = 0; 
      dduDySqr = 0; 
      ddvDxSqr = 0; 
      ddvDxDy = 0; 
      ddvDySqr = 0.0;
    }
  else if (veltype == sinusoidalyVel)
    {
      MayDay::Error("pointExactOp not implemented for this velType");      
    }
  else if (veltype == simpleSinusoidalVel )
    {
      // first derivatives
      dudx = 2.0*Pi*cos(2*Pi*x[0]/a_domainSize[0])/a_domainSize[0];
      dudy = 0;
      dvdx = 0;
      dvdy = 2.0*Pi*cos(2*Pi*x[1]/a_domainSize[0])/a_domainSize[0];
      // second derivatives
      dduDxSqr = -4*Pi*Pi*sin(2*Pi*x[0]/a_domainSize[0])/(a_domainSize[0]*a_domainSize[0]);
      dduDxDy = 0; 
      dduDySqr = 0; 
      ddvDxSqr = 0; 
      ddvDxDy = 0; 
      ddvDySqr = -4*Pi*Pi*sin(2*Pi*x[1]/a_domainSize[0])/(a_domainSize[0]*a_domainSize[0]);
    }
  else if (veltype == transverseSinusoidalxVel )
    {
      MayDay::Error("pointExactOp not implemented for this velType");      
    }
  else if (veltype == transverseSinusoidalyVel )
    {
      MayDay::Error("pointExactOp not implemented for this velType");      
    }
  else if (veltype == fullSinusoidalxVel )
    {
      MayDay::Error("pointExactOp not implemented for this velType");      
    }
  else if (veltype == fullSinusoidalyVel )
    {
      MayDay::Error("pointExactOp not implemented for this velType");      
    }
  else if (veltype == sinusoidalVel )
    {
      RealVect scalX = x;
      scalX /= a_domainSize[0];
      scalX *= 2.0*Pi;
      // first derivatives
      dudx = 2.0*Pi*cos(scalX[0])*sin(scalX[1])/a_domainSize[0];
      dudy = 2.0*Pi*sin(scalX[0])*cos(scalX[1])/a_domainSize[0];;
      dvdx = -2.0*Pi*sin(scalX[0])*cos(scalX[1])/a_domainSize[0];;
      dvdy = -2.0*Pi*cos(scalX[0])*sin(scalX[1])/a_domainSize[0];;
      // second derivatives
      dduDxSqr = -4*Pi*Pi*sin(scalX[0])*sin(scalX[1])/(a_domainSize[0]*a_domainSize[0]);
      dduDxDy = 4*Pi*Pi*cos(scalX[0])*cos(scalX[1])/(a_domainSize[0]*a_domainSize[0]);
      dduDySqr = -4*Pi*Pi*sin(scalX[0])*sin(scalX[1])/(a_domainSize[0]*a_domainSize[0]);
      ddvDxSqr = -4*Pi*Pi*cos(scalX[0])*cos(scalX[1])/(a_domainSize[0]*a_domainSize[0]);
      ddvDxDy = 4*Pi*Pi*sin(scalX[0])*sin(scalX[1])/(a_domainSize[0]*a_domainSize[0]);
      ddvDySqr = -4*Pi*Pi*cos(scalX[0])*cos(scalX[1])/(a_domainSize[0]*a_domainSize[0]);
    }
  else 
    {
      MayDay::Error("pointExactOp not implemented for this velType");      
    }
  // constant mu*H part
  Real mu = 1.0;
  val[0] = H*mu*(2.0*dduDxSqr + ddvDxDy + 0.5*(dduDySqr + ddvDxDy));
  val[1] = H*mu*(0.5*(dduDxDy + ddvDxSqr) + dduDxDy + 2.0*ddvDySqr);
  if (thicknessType != constantThickness)
    {
      val[0] += dHmuDx*(2.0*dudx + dvdy) + 0.5*dHmuDy*(dudy + dvdx);
      val[1] += 0.5*dHmuDx*(dudy + dvdx) + dHmuDy*(dudx + 2.0*dvdy);
    }
  return val;
}



void
computeExactOp(LevelData<FArrayBox>& a_Error, 
               const LevelData<SigmaCS>& a_sigmaCoords, 
               const RealVect& a_domainSize)
{
  DataIterator dit = a_Error.getBoxes();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& errFab = a_Error[dit];
      const SigmaCS& thisCS = a_sigmaCoords[dit];
      const RealVect& dx = thisCS.dx();
      const FArrayBox& localH = thisCS.getH();

      //Real cellVol = D_TERM(dx[0],*dx[1],*dx[2]);
      BoxIterator bit(errFab.box());
      RealVect ccOffset = 0.5*RealVect::Unit;
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect mappedLoc(iv);
          mappedLoc += ccOffset;
          mappedLoc *= dx;
          RealVect realLoc = thisCS.realCoord(mappedLoc);
          Real H = localH(iv,0);
          RealVect thisExactVal = pointExactOp(realLoc, H, a_domainSize);
          
          // operator computed by velOp is actually cell integral 
          // multiply by cell volume to get that...
          //thisExactVal *= cellVol;
          
          // note that velfab has 2 components regardless of
          // dimensionality, since it's the horizontal velocity          
          errFab(iv,0) = thisExactVal[0];
          errFab(iv,1) = thisExactVal[1];
        }
    } // end loop over boxes        
}



void       
computeError(LevelData<FArrayBox>& a_Error, 
             const LevelData<FArrayBox>& a_LofVel, 
             const LevelData<SigmaCS>& a_sigmaCoords, 
             const RealVect& a_domainSize)
{

  const DisjointBoxLayout& grids = a_Error.getBoxes();

  // first, compute exact solution
  computeExactOp(a_Error, a_sigmaCoords, a_domainSize);

  DataIterator dit = a_Error.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      //const Box& gridBox = grids[dit];
      a_Error[dit].minus(a_LofVel[dit]);
    }

  
}


int
testVelOp();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testVelOp();

  if( status == 0 )
    pout() << indent << "VelocityOp test" << " passed." << endl ;
  else
    pout() << indent << "VelocityOp test" << " failed with return code "
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
testVelOp() 
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
  int xDir = 0;
  if (CH_SPACEDIM == 3) xDir = 1;
  for (int dir=xDir; dir<SpaceDim; dir++)
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


  // test zero-velocity, constant temperature
  {
    // isothermal ice value from Pattyn(2003)
    Real thetaVal = 238.15;
    
    // allocate storage
    DisjointBoxLayout grids(gridBoxes, procAssign, entireDomain);
    
    // note that nComp = 1 is required by LevelData API but is not used in definition of SigmaCS
    IntVect ghostVect = IntVect::Unit;
    RefCountedPtr<LevelData<SigmaCS> > sigmaCoords(new LevelData<SigmaCS>(grids,1, ghostVect));
    // now loop over grids and define coordSys on each grid
    
    DataIterator dit = grids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        Box ghostBox(grids[dit]);
        ghostBox.grow(ghostVect);
        SigmaCS& localCS = (*sigmaCoords)[dit];
        localCS.define(ghostBox, dx);
      }
    
    defineSigmaCS(*sigmaCoords, domainSize, thicknessType, basalType);

    LevelData<FArrayBox> theta(grids, 1, ghostVect);
    LevelData<FArrayBox> horizontalVel(grids, 2, ghostVect);

    for (dit.begin(); dit.ok(); ++dit)
      {
        theta[dit].setVal(thetaVal);
        horizontalVel[dit].setVal(0);
      }

    // face-centered mu:
    RefCountedPtr<LevelData<FluxBox> > faceMu(new LevelData<FluxBox>(grids,1));

    // could do this with a ConstitutiveRelation object, 
    // but for this test, it makes more sense to just set mu to one
    for (dit.begin(); dit.ok(); ++dit)
      {
        //const Box& gridBox = grids[dit];
        //const SigmaCS& localCS = sigmaCoords[dit];

        Real constMu = 1.0;
        (*faceMu)[dit].setVal(constMu);
      }
    
    // now construct SSAVelocityOp
    // single-level for now, so allocate empty DBL for coarse and fine grids
    DisjointBoxLayout emptyLevel;
    int bogusRefRatio = -1;

    // first try homogeneous Dirichlet BC's
    IntVect loBCType = IntVect::Unit;
    IntVect hiBCType = IntVect::Unit;
    RealVect loBCVal = RealVect::Zero;
    RealVect hiBCVal = RealVect::Zero;
    
    RefCountedPtr<BCFunction> bcPtr = ConstDiriNeumBC(loBCType, loBCVal,
                                                      hiBCType, hiBCVal);

    BCHolder bc(bcPtr);
    SSAVelocityOp op(grids, emptyLevel, emptyLevel, faceMu, 
                     bogusRefRatio, bogusRefRatio, entireDomain, sigmaCoords, dx,
                     dx, bc);

    LevelData<FArrayBox> L(grids, 2, IntVect::Zero);
    op.applyOp(L, horizontalVel);
    
    // compute errors
    RealVect maxErr = RealVect::Zero;
    for (dit.begin(); dit.ok(); ++dit)
      {
        const FArrayBox& thisL = L[dit];
        // exact solution here is zero
        RealVect thisMax;
        thisMax[0] = thisL.norm(0,0,1);
        thisMax[1] = thisL.norm(0,1,1);
        if (thisMax[0] > maxErr[0]) maxErr[0] = thisMax[0];
        if (thisMax[1] > maxErr[1]) maxErr[1] = thisMax[1];
      }

    
    pout() << " zero velocity test: " << endl;
    pout() << "      max err = " << maxErr[0] << ", " << maxErr[1];
    pout() << endl;
    pout() << endl;                                           

    // test for acceptable error
    if (maxErr[0] > tolerance) 
      {
        pout() << "  ZERO-VELOCITY ERR FAILS FOR x-direction !" << endl;
        status += 1;
      }

    if (maxErr[1] > tolerance) 
      {
        pout() << "  ZERO-VELOCITY ERR FAILS FOR y-direction !" << endl;
        status += 1;
      }
  }
  
  // do mesh refinement testing 

  // storage for epsSqr errors
  RealVect L1ErrSave;
  RealVect L2ErrSave;
  RealVect MaxErrSave;

  int numRefine = 4;
  for (int refinement = 0; refinement <numRefine; refinement++)
    {
      if (verbose)  pout() << "refinement " << refinement << endl;

      // allocate storage
      DisjointBoxLayout grids(gridBoxes, procAssign, entireDomain);

      IntVect ghostVect = IntVect::Unit;      
      RefCountedPtr<LevelData<SigmaCS> > sigmaCoords(new LevelData<SigmaCS>(grids, 1, ghostVect));
      // now loop over grids and define coordSys on each grid
      

      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box ghostBox(grids[dit]);
          ghostBox.grow(ghostVect);
          SigmaCS& localCS = (*sigmaCoords)[dit];
          localCS.define(ghostBox, dx);
        }
      
      defineSigmaCS(*sigmaCoords, domainSize, thicknessType, basalType);
      
      // set up data -- note that velocity is 2 components
      // regardless of dimensionality, since it's only the horizontal 
      // velocity
      LevelData<FArrayBox> horizontalVel(grids,2, ghostVect);
      RefCountedPtr<LevelData<FluxBox> > faceMu(new LevelData<FluxBox>(grids, 1, ghostVect));

      initData(horizontalVel, *faceMu, *sigmaCoords, domainSize);

      
      // now construct SSAVelocityOp
      // single-level for now, so allocate empty DBL for coarse and fine grids
      DisjointBoxLayout emptyLevel;
      int bogusRefRatio = -1;
      
      // first try homogeneous Dirichlet BC's
      IntVect loBCType = IntVect::Unit;
      IntVect hiBCType = IntVect::Unit;
      RealVect loBCVal = RealVect::Zero;
      RealVect hiBCVal = RealVect::Zero;
      
      RefCountedPtr<BCFunction> bcPtr = ConstDiriNeumBC(loBCType, loBCVal,
                                                      hiBCType, hiBCVal);
      
      BCHolder bc(bcPtr);
      SSAVelocityOp op(grids, emptyLevel, emptyLevel, faceMu, 
                       bogusRefRatio, bogusRefRatio, entireDomain, 
                       sigmaCoords, dx,
                       dx, bc);
      

      // now compute L(vel)
      LevelData<FArrayBox> LofVel(grids, 2);
      
      op.applyOp(LofVel, horizontalVel);
      
      // now evaluate error
      LevelData<FArrayBox> Error(grids, 2);
      computeError(Error, LofVel, 
                   *sigmaCoords, domainSize);

      // now compute norms of errors
      RealVect L1Err(RealVect::Zero);
      RealVect L2Err(RealVect::Zero);
      RealVect MaxErr = RealVect::Zero;
      
      DisjointBoxLayout* fakeFineGrids = NULL;
      int nRefFine = -1;
      for (int velComp=0; velComp<2; velComp++)
        {
          Interval ErrorComps(velComp, velComp);
          L1Err[velComp] = computeNorm(Error, fakeFineGrids,
                                       nRefFine, dx[0], ErrorComps,
                                       1);
          
          
          L2Err[velComp] = computeNorm(Error, fakeFineGrids,
                                       nRefFine, dx[0], ErrorComps,
                                       2);

      
          MaxErr[velComp] = computeNorm(Error, fakeFineGrids,
                                        nRefFine, dx[0], ErrorComps,
                                        0);
        }
      
      // write errors
      if (verbose)
        {
          pout() << "L(vel) Errors: " << endl;
          
          for (int velDir=0; velDir<2; velDir++)
            {
              if (velDir == 0) 
                {
                  pout() << "x-direction: " << endl;
                }
              else 
                {
                  pout() << "y-direction: " << endl;
                }
              
              pout() << "L1 = " << L1Err[velDir]
                     << ",   L2 = " << L2Err[velDir]
                     << ",   max = " << MaxErr[velDir] << endl;
              
              Real L1Rate, L2Rate, MaxRate;
              
              if (refinement > 0)
                {
                  L1Rate = getConvergenceRate(L1Err[velDir], L1ErrSave[velDir]);
                  L2Rate = getConvergenceRate(L2Err[velDir], L2ErrSave[velDir]);
                  MaxRate = getConvergenceRate(MaxErr[velDir], MaxErrSave[velDir]); 
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
            } // end loop over vel directions
          
        } // end if verbose
      
      // save errors
      L1ErrSave = L1Err;
      L2ErrSave = L2Err;
      MaxErrSave = MaxErr;

      bool writePlotFiles = true;
      if (writePlotFiles)
        {
          // 2 velocity components, then 2 L(vel) components and H
          int nPlotComp = 5;
          LevelData<FArrayBox> plotData(grids, nPlotComp);
          for (dit.begin(); dit.ok(); ++dit)
            {
              const SigmaCS& thisCS = (*sigmaCoords)[dit];
              //const FArrayBox& DeltaFactors = thisCS.deltaFactors();
              const FArrayBox& thisVel = horizontalVel[dit];
              plotData[dit].copy(thisVel,0,0,2);
              const FArrayBox& thisL = LofVel[dit];
              plotData[dit].copy(thisL,0,2,2);
              const FArrayBox& thisH = thisCS.getH();
              plotData[dit].copy(thisH,0,4,1);
            }
          
          string fileRoot("testVelOp.");
          // use unigrid version here
          dit.reset();
          WriteSigmaMappedUGHDF5(fileRoot, grids, plotData, *sigmaCoords,
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
