#include <iostream>
#include <cassert>

// Chombo
#include "ParmParse.H"
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

// Local
#include "SimplePDEOpsF_F.H"

// SUNDIALS interface
#include "ChomboSundialsAdaptor.H"
#include "SimpleAdaptor.H"
#include "ChomboNVector.H"
#include "arkode/arkode_arkstep.h"

#include "UsingNamespace.H"

/* This example solves the PDE
 *
 * u_t = c * laplacian(u) + 2 * pi * cos(2 * pi * (t + x1 + ... + xd))
 *       + 4 * pi^2 * c * d * sin(2 * pi * (t + x1 + ... + xd))
 *
 * on a d-dimensional unit box with periodic boundary conditions. The exact
 * solutions is
 *
 * u(t, x1, ..., xd) = sin(2 * pi * (t + x1 + ... + xd))
 *
 * An explicit Runge-Kutta method from the ARKODE package of SUNDIALS is used
 * to integrate the problem forward in time.
 */

// A struct containing the data needed in the RHS function
struct UserData
{
  Real dx; // Mesh size
  Real diffusionCoeff;
};

// The RHS function for the PDE
int f(Real t, N_Vector y, N_Vector yDot, void *data) {
  // Extract the user data
  const UserData *ud = static_cast<UserData*>(data);
  const Real diffusionCoeff = ud->diffusionCoeff;
  const Real dx = ud->dx;

  // Extract the adaptors from the N_Vectors
  SimpleAdaptor* a_y = dynamic_cast<SimpleAdaptor*>(&NV_ADAP_CH(y));
  assert(a_y != NULL);
  SimpleAdaptor* a_yDot = dynamic_cast<SimpleAdaptor*>(&NV_ADAP_CH(yDot));
  assert(a_yDot != NULL);

  // Operate on the Chombo data
  LevelData<FArrayBox>* ly = a_y->getData();
  LevelData<FArrayBox>* lyDot = a_yDot->getData();
  DisjointBoxLayout dp = ly->disjointBoxLayout();
  assert(dp == lyDot->disjointBoxLayout());
  assert(ly->nComp() == lyDot->nComp());
  DataIterator dit = ly->dataIterator();
  a_y->exchange();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& box = dp[dit()];
    FArrayBox& yDotfab = (*lyDot)[dit()];
    const FArrayBox& yfab = (*ly)[dit()];
    FORT_SIMPLEPDERHS(CHF_CONST_FRA(yfab), CHF_FRA(yDotfab), 
          CHF_BOX(box), CHF_REAL(diffusionCoeff), CHF_REAL(dx), CHF_REAL(t));
  }

  return 0;
}

void setExactSolution(N_Vector y, const Real time, void *data) {
  const UserData *ud = static_cast<UserData*>(data);
  const Real dx = ud->dx;

  SimpleAdaptor* a_y = dynamic_cast<SimpleAdaptor*>(&NV_ADAP_CH(y));
  LevelData<FArrayBox>* ly = a_y->getData();
  DisjointBoxLayout dp = ly->disjointBoxLayout();
  DataIterator dit = ly->dataIterator();
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& box = dp[dit()];
    FArrayBox& fab = (*ly)[dit()];
    FORT_SIMPLEPDESOL(CHF_FRA(fab), CHF_BOX(box), 
        CHF_REAL(dx), CHF_REAL(time));
  }
}

void run() {
  // Pick discretization parameters
  const int numCells = 16;
  UserData ud;
  ud.dx = 1.0 / numCells;
  ud.diffusionCoeff = 1;

  // Create the mesh
  const Box domainBox(IntVect::Zero, (numCells - 1) * IntVect::Unit);
  const bool periodic[] = {true, true, true};
  const ProblemDomain baseDomain(domainBox, periodic);

  const int maxBoxSize = numCells / 2;
  Vector<Box> vectBoxes;
  domainSplit(baseDomain, vectBoxes, maxBoxSize);
  Vector<int> procAssign(vectBoxes.size(), 0);
  LoadBalance(procAssign, vectBoxes);
  DisjointBoxLayout level0(vectBoxes, procAssign, baseDomain);

  int nghost = 1;
  int ncomp = 1;
  IntVect ghostVect = nghost*IntVect::Unit;
  LevelData<FArrayBox> u(level0, ncomp, ghostVect);

  SimpleAdaptor adaptor;
  adaptor.define(&u, false);

  // Initialize SUNDIALS integrator
  const Real tspan[] = {0, 1};
  SUNContext ctx;
#ifdef CH_MPI
  void* mpi_comm_ptr = (void*)(&Chombo_MPI::comm);
  SUNContext_Create(mpi_comm_ptr, &ctx);
#else
  SUNContext_Create(nullptr, &ctx);
#endif
  N_Vector y = N_VNew_Ch(ctx, &adaptor);
  setExactSolution(y, tspan[0], &ud);

  void *integrator = ARKStepCreate(f, nullptr, tspan[0], y, ctx);
  ARKStepSetUserData(integrator, &ud);
  ARKStepSetMaxNumSteps(integrator, 1e6);
  ARKStepSStolerances(integrator, 1e-5, 1e-5);

  // Solve the PDE
  Real tOut;
  ARKStepEvolve(integrator, tspan[1], y, &tOut, ARK_NORMAL);
  
  // Get solution statistics
  long int steps;
  long int feEvals, fiEvals;
  ARKStepGetNumSteps(integrator, &steps);
  ARKStepGetNumRhsEvals(integrator, &feEvals, &fiEvals);

  N_Vector err = N_VClone(y);
  setExactSolution(err, tspan[1], &ud);
  N_VLinearSum(1, err, -1, y, err);

  // Print solution statistics
  pout() << "Start time: " << tspan[0] << std::endl;
  pout() << "Final time: " << tOut << std::endl;
  pout() << "Timesteps: " << steps << std::endl;
  pout() << "Function Evaluations: " << feEvals << std::endl;
  pout() << "Inf-Norm Error: " << N_VMaxNorm(err) << std::endl;
  pout() << "1-Norm Error: " << N_VL1Norm(err) << std::endl;
  pout() << endl << "Error values:" << endl;
  N_VPrint(err);

#ifdef CH_USE_HDF5
  adaptor.printFileHDF("solution.hdf5");
#endif

  // Free SUNDIALS data
  ARKStepFree(&integrator);
  N_VDestroy(y);
  N_VDestroy(err);
  SUNContext_Free(&ctx);
}

int main(int argc ,char* argv[]) {
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  pout() << "Starting example" << std::endl;
  run();
  pout() << "Done - exiting!" << std::endl;

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}

