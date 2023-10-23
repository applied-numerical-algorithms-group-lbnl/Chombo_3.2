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
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "memusage.H"
#include "AMRPoissonOp.H"
#include "AMRPoissonOp.H"
//#include "Proto_Helmholtz_Op.H"
#include "UsingNamespace.H"


///

//everything public and static
static int s_verbosity;


class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;
  static std::vector<int> s_bcLo, s_bcHi;
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  a_values[0]=0.;
}

void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{

  if (!a_domain.domainBox().contains(a_state.box()))
  {
    Box valid = a_valid;
    for (int i=0; i<CH_SPACEDIM; ++i)
    {
      // don't do anything if periodic
      if (!a_domain.isPeriodic(i))
      {
        Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
        Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
        if (!a_domain.domainBox().contains(ghostBoxLo))
        {
          DiriBC(a_state,
                 valid,
                 a_dx,
                 true,
                 ParseValue,
                 i,
                 Side::Lo,
                 1);
        }

        if (!a_domain.domainBox().contains(ghostBoxHi))
        {
          DiriBC(a_state,
                 valid,
                 a_dx,
                 true,
                 ParseValue,
                 i,
                 Side::Hi,
                 1);
        }
      } // end if is not periodic in ith direction
    }
  }
}

void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
            Vector<ProblemDomain>& a_amrDomains,
            Vector<int>& a_refRatios,
            Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setRHS");

  for (int lev=0; lev<=a_finestLevel; lev++)
  {
    LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
    const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();

    // rhs is cell-centered...
    RealVect ccOffset = 0.5*a_amrDx[lev]*RealVect::Unit;

    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& thisRhs = levelRhs[levelDit];
      thisRhs.setVal(1.);
    }
  }
}


///
void
setupLLCornerGrids(Vector<DisjointBoxLayout>& a_amrGrids,
                   Vector<ProblemDomain>& a_amrDomains,
                   Vector<int>& a_refRatios,
                   Vector<Real>& a_amrDx,
                   int& a_finestLevel)
{
  CH_TIME("setupGrids");

  a_finestLevel = 0;
  ParmParse ppGrids("setupGrids");

  // get grid generation parameters
  int maxLevel, maxBoxSize, blockFactor;
  Real fillRatio;

  ppGrids.get("max_level", maxLevel);

  ppGrids.get("max_box_size",maxBoxSize);

  ppGrids.get("block_factor", blockFactor);

  ppGrids.get("fillRatio", fillRatio);
  int maxNumLevels = maxLevel +1;

  a_refRatios.resize(maxNumLevels);
  ppGrids.getarr("ref_ratio", a_refRatios, 0, maxNumLevels);

  Vector<int>  is_periodic_int;
  bool is_periodic[SpaceDim];
  ppGrids.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
  {
    is_periodic[dir] = (is_periodic_int[dir] == 1);
  }

  IntVect numCells;
  Vector<int> incells(SpaceDim);
  ppGrids.getarr("num_cells", incells, 0, SpaceDim);
  numCells = IntVect(D_DECL6(incells[0],incells[1],incells[2],
                             incells[3],incells[4],incells[5]) );

  RealVect domainSize = RealVect::Unit;
  if (ppGrids.contains("domain_size"))
  {
    Vector<Real> insize(SpaceDim);
    ppGrids.getarr("domain_size", insize, 0, SpaceDim);
    domainSize = RealVect(D_DECL6(insize[0],insize[1],insize[2],
                                  insize[3],insize[4],insize[5]) );
  }

  // resize dataholders
  a_amrGrids.resize(maxNumLevels);
  a_amrDomains.resize(maxNumLevels);
  a_amrDx.resize(maxNumLevels,-1);
  a_finestLevel = 0;

  // assumes dx=dy=dz
  a_amrDx[0] = domainSize[0]/numCells[0];

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = numCells - IntVect::Unit;

  ProblemDomain baseDomain(domLo, domHi, is_periodic);
  a_amrDomains[0] = baseDomain;

  // set up refined domains, etc
  for (int lev=1; lev<= maxLevel; lev++)
  {
    a_amrDomains[lev] = a_amrDomains[lev-1];
    a_amrDomains[lev].refine(a_refRatios[lev-1]);
    a_amrDx[lev] = a_amrDx[lev-1]/a_refRatios[lev-1];
  }

  
  Vector<Vector<Box> > vectBoxes(maxLevel+1);
  {
    CH_TIME("setupGrids: BaseGridCreation");
    // generate base level grids

    domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

    Vector<int> procAssign(vectBoxes[0].size(), 0);

    LoadBalance(procAssign, vectBoxes[0]);

    DisjointBoxLayout baseGrids(vectBoxes[0], procAssign, baseDomain);

    a_amrGrids[0] = baseGrids;
  }


  if (maxLevel > 0)
  {
    CH_TIME("setupGrids: meshRefine and all that");
    int bufferSize = 0;
    BRMeshRefine meshGen(a_amrDomains[0],
                         a_refRatios,
                         fillRatio,
                         blockFactor,
                         bufferSize,
                         maxBoxSize);
      
    // to be used by MeshRefine...
    Vector<Vector<Box> > old_meshes(maxLevel+1);
    old_meshes[0] = vectBoxes[0];
    for (int lev=1; lev< old_meshes.size(); lev++)
    {
      Box whole_dom = a_amrDomains[lev].domainBox();
      old_meshes[lev] = Vector<Box>(1, whole_dom);
    }

    //just always tag IntVect::Zero;
    IntVectSet taglev(IntVect::Zero);
    Vector<IntVectSet> tags(a_finestLevel+1, taglev);
    Vector<Vector<Box> > new_meshes;

    meshGen.regrid(new_meshes, tags,  0, a_finestLevel,  old_meshes);

    for (int ilev=0; ilev < a_finestLevel; ilev++)
    {
      const auto& boxes = new_meshes[ilev];
      Vector<int> procs;
      LoadBalance(procs, boxes);
      a_amrGrids[ilev] = DisjointBoxLayout(boxes, procs);
    }
  }
} //end setupLLCornerGrids

  
void
setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver,
            LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
            const Vector<DisjointBoxLayout>& a_amrGrids,
            const Vector<ProblemDomain>& a_amrDomains,
            const Vector<int>& a_refRatios,
            const Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setupSolver");

  ParmParse ppSolver("solver");

  int numLevels = a_finestLevel+1;

  AMRPoissonOpFactory opFactory;

  // solving poisson problem here
  Real alpha =0.0;
  Real beta = 1.0;

  opFactory.define(a_amrDomains[0],
                   a_amrGrids,
                   a_refRatios,
                   a_amrDx[0],
                   &ParseBC, alpha, beta);

  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

  a_amrSolver->define(a_amrDomains[0], castFact,
                      &a_bottomSolver, numLevels);

  // multigrid solver parameters
  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ppSolver.get("num_smooth", numSmooth);
  ppSolver.get("num_mg",     numMG);
  ppSolver.get("max_iterations", maxIter);
  ppSolver.get("tolerance", eps);
  ppSolver.get("hang",      hang);

  Real normThresh = 1.0e-30;
  a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                   numMG, maxIter, eps, hang, normThresh);
  a_amrSolver->m_verbosity = s_verbosity-1;

  // optional parameters
  ppSolver.query("num_pre", a_amrSolver->m_pre);
  ppSolver.query("num_post", a_amrSolver->m_post);
  ppSolver.query("num_bottom", a_amrSolver->m_bottom);
}

int runSolver()
{
  CH_TIME("runSolver");

  int status = 0;
  ParmParse ppMain("main");

  ppMain.query("verbosity", s_verbosity);

  // set up grids&
  Vector<DisjointBoxLayout> amrGrids;
  Vector<ProblemDomain> amrDomains;
  Vector<int> refRatios;
  Vector<Real> amrDx;
  int finestLevel;

  setupLLCornerGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);
  for(int ilev = 0; ilev < amrGrids.size(); ilev++)
  {
    pout() << "ilev = " << ilev << ", grids = " << endl;
    amrGrids[ilev].print();
  }
  // initialize solver
  AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;
  amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
  BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = s_verbosity-2;
  setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
              refRatios, amrDx, finestLevel);


  // allocate solution and RHS, initialize RHS
  int numLevels = amrGrids.size();
  Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);
  // this is for convenience
  Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);

  for (int lev=0; lev<=finestLevel; lev++)
  {
    const DisjointBoxLayout& levelGrids = amrGrids[lev];
    phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
    rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
    resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
  }

  setRHS(rhs, amrDomains, refRatios, amrDx, finestLevel );

  // do solve
  int iterations = 1;
  ppMain.get("iterations", iterations);

  for (int iiter = 0; iiter < iterations; iiter++)
  {
    bool zeroInitialGuess = true;
    pout() << "about to go into solve" << endl;
    amrSolver->solve(phi, rhs, finestLevel, 0, zeroInitialGuess);
    pout() << "done solve" << endl;
  }

  // write results to file

  bool writePlots = true;
  ppMain.query("writePlotFiles", writePlots);

#ifdef CH_USE_HDF5

  if (writePlots)
  {
    int numLevels = finestLevel +1;
    Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);
       
    pout() << "Write Plots. norm=" << amrSolver->computeAMRResidual(resid,phi,rhs,finestLevel,0) << endl;
       
    for (int lev=0; lev<numLevels; lev++)
    {
      plotData[lev] = new LevelData<FArrayBox>(amrGrids[lev],
                                               3, IntVect::Zero);

      Interval phiInterval(0,0);
      phi[lev]->copyTo(phiInterval, *plotData[lev], phiInterval);
      Interval rhsInterval(1,1);
      rhs[lev]->copyTo(phiInterval, *plotData[lev], rhsInterval);
      Interval resInterval(2,2);
      resid[lev]->copyTo(phiInterval, *plotData[lev], resInterval);
    }

    string fname = "poissonOut.";

    char suffix[30];
    sprintf(suffix, "%dd.hdf5",SpaceDim);
    fname += suffix;

    Vector<string> varNames(3);
    varNames[0] = "phi";
    varNames[1] = "rhs";
    varNames[2] = "res";

    Real bogusVal = 1.0;

    WriteAMRHierarchyHDF5(fname,
                          amrGrids,
                          plotData,
                          varNames,
                          amrDomains[0].domainBox(),
                          amrDx[0],
                          bogusVal,
                          bogusVal,
                          refRatios,
                          numLevels);

    // clean up
    for (int lev=0; lev<plotData.size(); lev++)
    {
      delete plotData[lev];
    }
  } // end if writing plots
#endif // end if HDF5

  // clean up
  for (int lev=0; lev<phi.size(); lev++)
  {
    delete phi[lev];
    delete rhs[lev];
    delete resid[lev];
  }

  delete amrSolver;

  return status;
}

/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    if (argc < 2)
    {
      cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = runSolver();
    status += solverStatus;
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
