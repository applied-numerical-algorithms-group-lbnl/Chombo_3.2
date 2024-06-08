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
#include <vector>
#include <cassert>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "computeSum.H"
#include "computeNorm.H"
#include "FABView.H"
#include "DebugDump.H"

#ifdef CH_USE_HDF5
#include "H5FDmpio.h"
#endif

// SUNDIALS interface
#include "AmrAdaptor.H"
#include "ChomboNVector.H"

#include "UsingNamespace.H"


/// Prototypes:
int
testRegression(SUNContext ctx);

int
testConstruction(SUNContext ctx);

int
testVectorOps(SUNContext ctx);

int
testChomboConsistency(SUNContext ctx);

//int
//createSingleLevelData(SUNContext ctx, 

template <typename T> int
compareResults(const std::string&, T, T);

/// Global variables for handling output:
static const char* pgmname = "testRegression" ;

static const char* indent = "   ";
static const char* indent2 = "      " ;
static const Real tolerance = 10 * std::numeric_limits<Real>::epsilon();

int
main(int argc ,char* argv[])
{
  SUNContext ctx;

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  // int comm = MPI_COMM_WORLD;
  // SUNContext_Create(&comm, &ctx);
  SUNContext_Create(NULL, &ctx);
#else
  SUNContext_Create(NULL, &ctx);
#endif

  pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testRegression(ctx) ;
  pout() << indent << pgmname ;
  if ( status == 0 )
    {
      pout() << " passed." << endl;
    }
  else
    {
      pout() << " failed with result code " << status << endl ;
    }

  SUNContext_Free(&ctx);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
  return status ;
}

int
testRegression(SUNContext ctx)
{
  CH_TIME("testRegression");
  int status = 0;

  // TODO run once N_VNew ready
  // status += testConstruction(ctx);
  status += testVectorOps(ctx);
  // status += testChomboConsistency(ctx);

  return status;
}

int
testConstruction(SUNContext ctx)
{
  CH_TIME("testConstruction");
  int status = 0;

  // TODO use N_VNew
  N_Vector v = NULL;
  if (v == nullptr) {
    pout() << "N_VNew: expected = non-null, computed = null";
    status++;
  }

  sunindextype actualLen = N_VGetLength(v);
  compareResults("N_VGetLength", sunindextype(0), actualLen);

  N_VDestroy(v);
  if (v != nullptr) {
    pout() << "N_VNew: expected = null, computed = non-null";
    status++;
  }

  return status;
}

// Ensure the Chombo NVector acts like a traditional vector despite being
// backed by a complex AMR grid
int
testVectorOps(SUNContext ctx) {
  CH_TIME("testVectorOps");
  int status = 0;

  // Create an AMR domain, ref ratio r = 2,
  // with base 32^D, center 8^D refined once, center 2^D refined twice
  // Size of box per level = 32, 8*2, 2*2*2
  // Num pts per level = 1024, 256, 64, total = 1344
  // but total valid (32^2 - 8^2) + ((r*8)^2 - (r*2)^2) + (r*r*2)^2 = 1264
  AmrAdaptor* adaptor = new AmrAdaptor();
  Vector<LevelData<FArrayBox>*> vecAmrData;

  int nComp = 2;
  pout() << indent2 << "num comps = " << nComp << endl;
  IntVect ghostVect = IntVect::Unit;
  int maxBoxSize = 8; // will force MPI comm across levels
  const int maxLevel = 3;
  int N[maxLevel] = {32, 16, 8}; // how big the region is for each level
  int lo[maxLevel] = {0, 24, 60}; // ll corner of region for each level
  Vector<int> refRatio(2,2); // between successive levels
  Vector<DisjointBoxLayout> vecDBL;
  for (int lvl=0; lvl < maxLevel; lvl++)
  {
    IntVect loVect = lo[lvl]*IntVect::Unit;
    IntVect hiVect = loVect + (N[lvl]-1)*IntVect::Unit;
    Box b(loVect, hiVect);
    pout() << "level " << lvl << ", box " << b 
      << ", pts " << b.numPts() << endl;
    ProblemDomain domain(b);

    Vector<Box> vectBoxes;
    domainSplit(domain, vectBoxes, maxBoxSize, 1);
    Vector<int> procAssign(vectBoxes.size(), 0);
    LoadBalance(procAssign, vectBoxes);
    Box baseDomain(IntVect::Zero, IntVect::Unit*N[0]*std::pow(2,lvl)-1);
    DisjointBoxLayout dbl(vectBoxes, procAssign, baseDomain);
    // pout() << "level " << lvl << ", dbl " << dbl << endl;

    vecDBL.push_back(dbl);
    LevelData<FArrayBox>* ldf = new LevelData<FArrayBox>();
    ldf->define(dbl, nComp, ghostVect);
    vecAmrData.push_back(ldf);
  }
  adaptor->define(vecAmrData, refRatio, true);

  long numAmrPts = 0;
  long nPlvl = 0;
  for (int lvl=0; lvl < maxLevel-1; lvl++)
  {
    long nClvl = vecDBL[lvl].numCells();
    long nFlvl = vecDBL[lvl+1].numCells();
    pout() << "level " << lvl << ", pts " << nClvl << endl;
    // Subtract the number of points on the finer level
    nPlvl = nClvl - nFlvl/std::pow(refRatio[lvl],SpaceDim);
    pout() << "level " << lvl << ", uncovered pts " << nPlvl << endl;
    numAmrPts += nPlvl;
  }
  nPlvl = vecDBL[maxLevel-1].numCells(); // finest level
  pout() << "level " << (maxLevel-1) << ", uncovered pts " << nPlvl << endl;
  numAmrPts += nPlvl; // add in finest level
  pout() << "Total number of points in AMR hierarchy " << numAmrPts << endl;
  numAmrPts *= nComp;
  pout() << " w/ comps: " << numAmrPts << endl;

  // Create vector
  N_Vector v = N_VNew_Ch(ctx, adaptor);
  // The following line can be used to test a SUNDIALS MPI vector
  //N_Vector v = N_VMake_Parallel(MPI_COMM_WORLD, 1024, 1024, new sunrealtype[1024], ctx);

  // Length Test
  sunindextype actualLen = N_VGetLength(v);
  sunindextype expectedLen = (sunindextype) numAmrPts;

  status += compareResults("N_VGetLength", expectedLen, actualLen);

  // Basic Norm Tests
  N_VConst(-2, v);
  sunrealtype expectedVal = 2;
  sunrealtype actualVal = N_VMaxNorm(v);
  status += compareResults("N_VMaxNorm", expectedVal, actualVal);

  expectedVal = 2 * expectedLen;
  actualVal = N_VL1Norm(v);
  status += compareResults("N_VL1Norm", expectedVal, actualVal);

  // Clone Test
  N_Vector w = N_VClone(v);
  actualLen = N_VGetLength(w);
  status += compareResults("N_VClone length", expectedLen, actualLen);

  // Add Const Tests
  N_VConst(0, w); // Initialize w after clone
  N_VAddConst(v, 5.6, w);
  expectedVal = 3.6 * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VAddConst norm", expectedVal, actualVal);

  // Set nonuniform data in v, by level
  for (int lvl=0; lvl < maxLevel; lvl++)
  {
    LevelData<FArrayBox>& ldf = *adaptor->m_amrData[lvl];
    // Level entries are level+0.5
    for (DataIterator dit = ldf.dataIterator(); dit.ok(); ++dit)
      ldf[dit()].setVal(lvl + 0.5);
  }

  // Min Entry Test
  expectedVal = 3.6;
  actualVal = N_VMin(w);
  status += compareResults("N_VMin 1", expectedVal, actualVal);

  expectedVal = 0.5;
  actualVal = N_VMin(v);
  status += compareResults("N_VMin 2", expectedVal, actualVal);

  // Scale Test
  N_VScale(-2.5, v, w);
  // v is lvl 1=.5, lvl 2=1.5, lvl 3=2.5
  // So L1 norm is .5 * 960 + 1.5 * 240 + 2.5 * 64 = 1000 * 2 comps
  expectedVal = std::abs(-2.5*1000*nComp);
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VScale out of place", expectedVal, actualVal);

  N_VScale(0.4, w, w);
  expectedVal = 0.4 * expectedVal;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VScale in place", expectedVal, actualVal);

  // Inverse Test
  // v is lvl 1=.5, lvl 2=1.5, lvl 3=2.5
  // So 1/v L1 norm is 1/.5 * 960 + 1/1.5 * 240 + 1/2.5 * 64 = 2105.6*2 comps
  N_VInv(v, w);
  expectedVal = nComp*2105.6;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VInv out of place", expectedVal, actualVal);

  N_VInv(w, w);
  // v is lvl 1=.5, lvl 2=1.5, lvl 3=2.5
  // So L1 norm is .5 * 960 + 1.5 * 240 + 2.5 * 64 = 1000 * 2 comps
  expectedVal = nComp*1000;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VInv in place", expectedVal, actualVal);

  /* OPTIONAL - TODO, finish later
  sunbooleantype allNonzero = N_VInvTest(v, w);
  status += compareResults("N_VInvTest out of place nonzero", SUNTRUE, allNonzero);
  expectedVal = 4.0 / 3.0 * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VInvTest out of place norm", expectedVal, actualVal);
  */

  /* OPTIONAL - TODO, finish later
  N_VConst(0, w);
  allNonzero = N_VInvTest(w, w);
  status += compareResults("N_VInvTest in place nonzero", SUNFALSE, allNonzero);
  */

  // Abs Test
  N_VScale(-2, v, v); // scale so entries are -5, -3 and -1
  N_VAbs(v, w);
  expectedVal = 1;
  actualVal = N_VMin(w);
  status += compareResults("N_VAbs out of place", expectedVal, actualVal);

  N_VConst(-100.4, w);
  N_VAbs(w, w);
  expectedVal = 100.4;
  actualVal = N_VMin(w);
  status += compareResults("N_VAbs in place", expectedVal, actualVal);

  /* OPTIONAL - TODO, finish later
  // Compare Test
  N_VCompare(2, v, w);
  expectedVal = expectedLen / 2;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VCompare out of place", expectedVal, actualVal);

  N_VCompare(8.2, w, w);
  expectedVal = 0;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VCompare in place", expectedVal, actualVal);

  // Min Quotient Test
  N_VConst(6, w);
  expectedVal = -0.5;
  actualVal = N_VMinQuotient(v, w);
  status += compareResults("N_VMinQuotient 1", expectedVal, actualVal);

  expectedVal = -6;
  actualVal = N_VMinQuotient(w, v);
  status += compareResults("N_VMinQuotient 2", expectedVal, actualVal);
  */

  // Dot Product Test
  N_VConst(2, w);
  // v is lvl 1=-1, lvl 2=-3, lvl 3=-5
  // So L1 norm is -1 * 960 + -3 * 240 + -5 * 64 = -2000 * 2 comps
  expectedVal = 2 * -2000 * nComp;
  actualVal = N_VDotProd(w, v);
  status += compareResults("N_VDotProd", expectedVal, actualVal);

  expectedVal = 2 * -2000 * nComp;
  actualVal = N_VDotProd(v, w);
  status += compareResults("N_VDotProd symmetry", expectedVal, actualVal);

  // So L1 norm is 1 * 960 + 9 * 240 + 25 * 64 = 4720 * 2 comps
  expectedVal = 4720 * nComp;
  actualVal = N_VDotProd(v, v);
  status += compareResults("N_VDotProd same args", expectedVal, actualVal);

  // Division Test
  N_VConst(2, w);
  N_Vector x = N_VClone(v);
  N_VDiv(v, w, x);
  // v is lvl 1=-1, lvl 2=-3, lvl 3=-5
  // So L1 norm is -1/2 * 960 + -3/2 * 240 + -5/2 * 64 = -1000 * comps
  expectedVal = 1000 * nComp;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VDiv out of place", expectedVal, actualVal);

  // So L1 norm is 2/-1 * 960 + 2/-3 * 240 + 2/-5 * 64 = -2105.6 * comps
  N_VDiv(w, v, w);
  expectedVal = std::abs(-2105.6)*nComp;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VDiv in place 1", expectedVal, actualVal);

  N_VDiv(x, x, x);
  expectedVal = expectedLen;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VDiv in place 2", expectedVal, actualVal);

  // Product Test
  N_VConst(2, w);
  N_VProd(v, w, x);
  // v is lvl 1=-1, lvl 2=-3, lvl 3=-5
  // So L1 norm is -1*2 * 960 + -3*2 * 240 + -5*2 * 64 = -4000 * comps
  expectedVal = 4000 * nComp;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VProd out of place", expectedVal, actualVal);

  // So L1 norm is 2*-1 * 960 + 2*-3 * 240 + 2*-5 * 64 = -4000 * comps
  N_VProd(w, v, w);
  expectedVal = 4000 * nComp;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VProd in place 1", expectedVal, actualVal);

  // x is lvl 1=-2, lvl 2=-6, lvl 3=-10
  // So L1 norm is -2*-2 * 960 + -6*-6 * 240 + -10*-10 * 64 = 18880 * comps
  N_VProd(x, x, x);
  expectedVal = 18880 * nComp;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VProd in place 2", expectedVal, actualVal);

  // Weighted RMS Norm Test
  N_VConst(1, w);
  // v is lvl 1=-1, lvl 2=-3, lvl 3=-5
  // So sum sq is 1 * 960 + 9 * 240 + 25 * 64 = 4720 * 2 comps
  expectedVal = sqrt((4720.0*nComp)/expectedLen);
  actualVal = N_VWrmsNorm(v, w);
  status += compareResults("N_VWrmsNorm 1", expectedVal, actualVal);

  N_VConst(2, w);
  expectedVal = 2*expectedVal;
  actualVal = N_VWrmsNorm(w, v);
  status += compareResults("N_VWrmsNorm 2", expectedVal, actualVal);

  // Linear Sum Test
  N_VConst(3, w);
  N_VLinearSum(0, v, 0, w, x);
  expectedVal = 0;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VLinearSum out of place 1", expectedVal, actualVal);

  N_VLinearSum(2, v, 2.0 / 3.0, w, x);
  // v is lvl 1=-1, lvl 2=-3, lvl 3=-5
  // 2/3*w is 2/3 * 3 = 2
  // So 2v + 2 = x is lvl 1=0, lvl 2=-4, lvl 3=-8
  // So sum is 0 * 960 + -4 * 240 + -8 * 64 = -1472 * 2 comps
  // 2*(-1 * 960 + -3 * 240 + -5 * 64) + 2*1264 = -1472 * 2 comps
  expectedVal = abs(2*(-2000.0)*nComp + 2*expectedLen);
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VLinearSum out of place 2", expectedVal, actualVal);

  N_VLinearSum(2, v, 2.0 / 3.0, w, w);
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VLinearSum in place 1", expectedVal, actualVal);

  N_VLinearSum(-1.0 / 2.0, x, -1, v, x);
  // x is lvl 1=0, lvl 2=-4, lvl 3=-8
  // v is lvl 1=-1, lvl 2=-3, lvl 3=-5
  // So -.5*x-v sum is 
  // 1 * 960 + (-.5*(-4) - (-3)) * 240 + (-.5*(-8) - (-5)) * 64 = 2736 * 2 comps
  expectedVal = std::abs(2736)*nComp;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VLinearSum in place 2", expectedVal, actualVal);

  // TODO - maybe redo *all* the tests for comps?
  // Redo a test with different comps on max level
  AmrAdaptor* ax = dynamic_cast<AmrAdaptor*>(&NV_ADAP_CH(x));
  assert(ax != NULL);
  AmrAdaptor* av = dynamic_cast<AmrAdaptor*>(&NV_ADAP_CH(v));
  assert(av != NULL);
  Vector<LevelData<FArrayBox>*> vldfx = ax->getData();
  Vector<LevelData<FArrayBox>*> vldfv = av->getData();
  // x is lvl 1=0, lvl 2=-2, lvl 3=-4
  // v is lvl 1=-1, lvl 2=-3, lvl 3=-5 for 1st comp
  // 2nd comp of v is 2* first comp
  for (int lvl=0; lvl < maxLevel; lvl++)
  {
    LevelData<FArrayBox>& ldfx = *vldfx[lvl];
    LevelData<FArrayBox>& ldfv = *vldfv[lvl];
    for (DataIterator dit = ldfx.dataIterator(); dit.ok(); ++dit)
    {
      Real xval = -2.*lvl;
      Real vval = -1.-2.*lvl;
      ldfx[dit()].setVal(xval); // all comps
      ldfv[dit()].setVal(vval); // all comps
      ldfv[dit()].mult(2.,1,1); // 2nd comp of v is 2* first comp
    }
  }
  // So -.5*x-v sum is 
  // (0*2 - (-1)*3) * 960 + (-.5*(-2)*2 - (-3)*3) * 240 
  // + (-.5*(-4)*2 - (-5)*3) * 64 = 6736 for all comps
  N_VLinearSum(-1.0 / 2.0, x, -1, v, x);
  expectedVal = 6736.;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VLinearSum, 2 comps, in place", expectedVal, actualVal);

  // Cleanup
  N_VDestroy(v);
  N_VDestroy(w);

  return status;
}

// Ensure Chombo's representation of an AMR grid matches SUNDIALS'
// representation as a vector. Norms and other quantities should be consistent.
int
testChomboConsistency(SUNContext ctx) {
  CH_TIME("testChomboConsistency");
  // Initialization
  int status = 0;

  int numCells = 32;
  IntVect loVect = IntVect::Zero;
  IntVect hiVect = (numCells-1)*IntVect::Unit;
  Box domainBox(loVect, hiVect);
  ProblemDomain baseDomain(domainBox);

  int maxBoxSize = numCells/2;

  Vector<Box> vectBoxes;
  domainSplit(baseDomain, vectBoxes, maxBoxSize, 1);
  Vector<int> procAssign(vectBoxes.size(), 0);
  LoadBalance(procAssign, vectBoxes);
  DisjointBoxLayout level0Grids(vectBoxes, procAssign, baseDomain);
    
#if 0
  IntVect ghostVect = IntVect::Unit;
  SimpleAdaptor* adaptor = new SimpleAdaptor();
  adaptor->dp = level0Grids;
  adaptor->nComp = 1;
  adaptor->ghost = ghostVect;
  adaptor->Data = new LevelData<FArrayBox>(level0Grids, 1, ghostVect);
#endif

  // TODO consistency tests

  // Cleanup

  return status;
}

template <typename T>
int
compareResults(const std::string& funcName, const T expected, const T actual)
{
  int status = 0;

  const T error = abs(expected - actual);
  if (error > tolerance * abs(expected)) {
    status++;
    pout() << "ERROR ";
  }

  pout() << funcName << ": expected = " << expected
         << ", actual = " << actual << ", error = " << error << endl;

  return status;
}



