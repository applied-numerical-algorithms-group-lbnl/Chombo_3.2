/*
 * When using MPI, this test should be run with 1, 2, 4, or 8 MPI ranks to
 * ensure boxes are evenly distributed.
 */

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
//#include "H5FDmpio.h"

// SUNDIALS interface
#include "ChomboSundialsAdaptor.H"
#include "SimpleAdaptor.H"
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
  MPI_Finalize();
#endif
  return status ;
}

int
testRegression(SUNContext ctx)
{
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
  int status = 0;

  int numCells = 32;
  IntVect loVect = IntVect::Zero;
  IntVect hiVect = (numCells-1)*IntVect::Unit;
  Box domainBox(loVect, hiVect);
  ProblemDomain baseDomain(domainBox);

  int maxBoxSize = numCells/4;

  Vector<Box> vectBoxes;
  domainSplit(baseDomain, vectBoxes, maxBoxSize, 1);
  Vector<int> procAssign(vectBoxes.size(), 0);
  LoadBalance(procAssign, vectBoxes);
  DisjointBoxLayout level0Grids(vectBoxes, procAssign, baseDomain);

  int nComp = 2;
  IntVect ghostVect = IntVect::Unit;
  SimpleAdaptor* adaptor = new SimpleAdaptor();
  LevelData<FArrayBox>* data = new LevelData<FArrayBox>(level0Grids, nComp, ghostVect);
  adaptor->define(data, true); // will delete the LDF in destructor

  // Create vector
  N_Vector v = N_VNew_Ch(ctx, adaptor);
  // The following line can be used to test a SUNDIALS MPI vector
  //N_Vector v = N_VMake_Parallel(MPI_COMM_WORLD, 1024, 1024, new sunrealtype[1024], ctx);

  // Length Test
  sunindextype actualLen = N_VGetLength(v);
  sunindextype expectedLen = (sunindextype) nComp*pow(numCells, CH_SPACEDIM);
  status += compareResults("N_VGetLength", expectedLen, actualLen);

  // Basic Norm Tests
  N_VConst(-2, v);
  sunrealtype expectedVal = 2;
  sunrealtype actualVal = N_VMaxNorm(v);
  status += compareResults("N_VMaxNorm", expectedVal, actualVal);

  expectedVal = 2 * expectedLen;
  actualVal = N_VL1Norm(v);
  status += compareResults("N_VL1Norm", expectedVal, actualVal);

  // Set the 2nd component to be -4, check L1 norm
  for (DataIterator dit = data->dataIterator(); dit.ok(); ++dit)
    (*data)[dit].setVal(-4,1);
  expectedVal = 3 * expectedLen;
  actualVal = N_VL1Norm(v);
  status += compareResults("N_VL1Norm comps different", expectedVal, actualVal);

  // Clone Test
  N_Vector w = N_VClone(v);
  actualLen = N_VGetLength(w);
  status += compareResults("N_VClone length", expectedLen, actualLen);

  // Add Const Tests
  N_VConst(0, w); // Initialize w after clone
  N_VAddConst(v, 5.6, w);
  expectedVal = .5 * (3.6 + 1.6) * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VAddConst norm", expectedVal, actualVal);

  // Set nonuniform data in v
  DataIterator dit = data->dataIterator();
  DisjointBoxLayout dbl = data->disjointBoxLayout();
  int i = 0;
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Half of the entries are 0.5 and the other half is 1.5
    (*data)[dit].setVal((i % 2) + 0.5, dbl[dit], 0, data->nComp());
    i++;
  }

  // Min Entry Test
  expectedVal = 1.6; // in second comp
  actualVal = N_VMin(w);
  status += compareResults("N_VMin 1", expectedVal, actualVal);

  expectedVal = 0.5;
  actualVal = N_VMin(v);
  status += compareResults("N_VMin 2", expectedVal, actualVal);

  // Scale Test
  N_VScale(-2.5, v, w);
  expectedVal = 2.5 * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VScale out of place", expectedVal, actualVal);

  N_VScale(0.4, w, w);
  expectedVal = expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VScale in place", expectedVal, actualVal);

  // Inverse Test
  N_VInv(v, w);
  expectedVal = 4.0 / 3.0 * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VInv out of place", expectedVal, actualVal);

  N_VInv(w, w);
  expectedVal = expectedLen;
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
  N_VScale(-2, v, v); // scale so entries are -3 and -1
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
  expectedVal = -4 * expectedLen;
  actualVal = N_VDotProd(w, v);
  status += compareResults("N_VDotProd", expectedVal, actualVal);

  expectedVal = -4 * expectedLen;
  actualVal = N_VDotProd(v, w);
  status += compareResults("N_VDotProd symmetry", expectedVal, actualVal);

  expectedVal = 5 * expectedLen;
  actualVal = N_VDotProd(v, v);
  status += compareResults("N_VDotProd same args", expectedVal, actualVal);

  // Division Test
  N_VConst(2, w);
  N_Vector x = N_VClone(v);
  N_VDiv(v, w, x);
  expectedVal = expectedLen;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VDiv out of place", expectedVal, actualVal);

  N_VDiv(w, v, w);
  expectedVal = 4.0 / 3.0 * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VDiv in place 1", expectedVal, actualVal);

  N_VDiv(x, x, x);
  expectedVal = expectedLen;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VDiv in place 2", expectedVal, actualVal);

  // Product Test
  N_VConst(-2, w);
  N_VProd(v, w, x);
  expectedVal = 4 * expectedLen;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VProd out of place", expectedVal, actualVal);

  N_VProd(w, v, w);
  expectedVal = 4 * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VProd in place 1", expectedVal, actualVal);

  N_VProd(x, x, x);
  expectedVal = 20 * expectedLen;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VProd in place 2", expectedVal, actualVal);

  // Weighted RMS Norm Test
  N_VConst(1, w);
  expectedVal = sqrt(5);
  actualVal = N_VWrmsNorm(v, w);
  status += compareResults("N_VWrmsNorm 1", expectedVal, actualVal);

  N_VConst(2, w);
  expectedVal = 2*sqrt(5);
  actualVal = N_VWrmsNorm(w, v);
  status += compareResults("N_VWrmsNorm 2", expectedVal, actualVal);

  // Linear Sum Test
  N_VConst(3, w);
  N_VLinearSum(0, v, 0, w, x);
  expectedVal = 0;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VLinearSum out of place 1", expectedVal, actualVal);

  N_VLinearSum(2, v, 2.0 / 3.0, w, x);
  expectedVal = 2 * expectedLen;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VLinearSum out of place 2", expectedVal, actualVal);

  N_VLinearSum(2, v, 2.0 / 3.0, w, w);
  expectedVal = 2 * expectedLen;
  actualVal = N_VL1Norm(w);
  status += compareResults("N_VLinearSum in place 1", expectedVal, actualVal);

  N_VLinearSum(-1.0 / 2.0, x, -1, v, x);
  expectedVal = 3 * expectedLen;
  actualVal = N_VL1Norm(x);
  status += compareResults("N_VLinearSum in place 2", expectedVal, actualVal);

  // Cleanup
  N_VDestroy(v);
  N_VDestroy(w);

  return status;
}

// Ensure Chombo's representation of an AMR grid matches SUNDIALS'
// representation as a vector. Norms and other quantities should be consistent.
int
testChomboConsistency(SUNContext ctx) {
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
    
  IntVect ghostVect = IntVect::Unit;
  SimpleAdaptor* adaptor = new SimpleAdaptor();
  LevelData<FArrayBox>* data = new LevelData<FArrayBox>(level0Grids, 1, ghostVect);
  adaptor->define(data, true); // will delete the LDF in destructor

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


