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
#include "PrChUtilities.H"
#include "DebuggingTools.H"
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

  Real alpha =4586.;
  Real beta = 4586.;
  ppSolver.get("alpha", alpha);
  ppSolver.get("beta" , beta) ;
  opFactory.define(a_amrDomains[0],
                   a_amrGrids,
                   a_refRatios,
                   a_amrDx[0],
                   &ParseBC, alpha, beta);

  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

  a_amrSolver->define(a_amrDomains[0], castFact,
                      &a_bottomSolver, numLevels);

  ///so it is doing the same relaxation as the proto version
  AMRPoissonOp::s_relaxMode = 4586;
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
  ppSolver.get("num_pre", a_amrSolver->m_pre);
  ppSolver.get("num_post", a_amrSolver->m_post);
  ppSolver.get("num_bottom", a_amrSolver->m_bottom);

  ///sets a static variable in AMRPoissonOp
  int relaxMode = 1;
  ppSolver.get("relax_mode", relaxMode);
  AMRPoissonOp::s_relaxMode =  relaxMode;
}
void
makeSmoothAndInteresting(Vector<LevelData<FArrayBox>* >  &  a_phi,
                         Vector<DisjointBoxLayout>       &  a_amrGrids,
                         Vector<ProblemDomain>           &  a_amrDomains,
                         Vector<int>                     &  a_refRatios,
                         Vector<Real>                    &  a_amrDx)
{
  double twopi   = 8.*atan(1.);
  Chombo::pout() << "makeSmoothAndInteresting: " << endl;
  for(int ilev = 0; ilev <  a_phi.size(); ilev++)
  {
    LevelData<FArrayBox>& phi_level = *a_phi[ilev];
    Real dx = a_amrDx[ilev];
    DataIterator dit = phi_level.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      FArrayBox& phi_fab = phi_level[dit[ibox]];
      for(BoxIterator bit(phi_fab.box()); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        double phi_val = 1;
        for(int idir = 0; idir < DIM; idir++)
        {
          double x_val    = dx*(double(iv[idir]) + 0.5);
          double sine_val = sin(twopi*x_val);
          phi_val *= sine_val;
        }
        phi_fab(iv, 0) =  phi_val;
      } //end loop over cells in box
    }   //end loop over boxes
    double norm = CH_XD::norm(phi_level, phi_level.interval(), 0);
    Chombo::pout() << "ilev : " << ilev << "phi_level.maxNorm = " << norm << endl;
  }     //end loop over levels
}       //end function makeSmoothAndInteresting
void
testAMROperator(Vector<      LevelData<FArrayBox>* >          & a_rhs,
                Vector<      LevelData<FArrayBox>* >          & a_phi, 
                AMRMultiGrid<LevelData<FArrayBox>  >          * a_amr_solver_ptr)
{
  Chombo::pout() << " going into testAMROperator: " << endl;
  for(int ilev = 0; ilev < a_phi.size(); ilev++)
  {
    auto& phi_level = *a_phi[ilev];
    double norm = CH_XD::norm(phi_level, phi_level.interval(), 0);
    Chombo::pout() << " ilev:           " << ilev  << endl;
    Chombo::pout() << " phi_level.norm(0) = " << norm << endl;
  }//end loop over levels

  
  int l_max  = a_rhs.size() - 1;
  int l_base = 0;
  a_amr_solver_ptr->computeAMROperator(a_rhs, a_phi, l_max, l_base, false);

  for(int ilev = 0; ilev < a_phi.size(); ilev++)
  {
    auto& rhs_level = *a_rhs[ilev];
    double norm = CH_XD::norm(rhs_level, rhs_level.interval(), 0);
    Chombo::pout() << " ilev:" << ilev  << endl;
    Chombo::pout() << " rhs_level.norm(0) = " << norm << endl;
  }//end loop over levels
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

  PrChUtilities<1>::setupLLCornerGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);
  for(int ilev = 0; ilev < amrGrids.size(); ilev++)
  {
    pout() << "ilev = " << ilev << ", grids = " << endl;
    amrGrids[ilev].print();
  }

  AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;
  amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
  BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = s_verbosity-2;
  setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
              refRatios, amrDx, finestLevel);

  int numLevels = amrGrids.size();
  Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);

  for (int lev=0; lev<=finestLevel; lev++)
  {
    const DisjointBoxLayout& levelGrids = amrGrids[lev];
    phi[lev] = new LevelData<FArrayBox>(levelGrids, 1,2*IntVect::Unit);
    rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
  }
 //because it is convenient, put in an operator test to compare with ye olde chombo
  makeSmoothAndInteresting(phi, amrGrids, amrDomains, refRatios, amrDx);
  testAMROperator(rhs, phi, amrSolver);
    
  setRHS(rhs, amrDomains, refRatios, amrDx, finestLevel );
  amrSolver->solve(phi, rhs, finestLevel, 0, true); //bool zeroInitialGuess = true;

#ifdef CH_USE_HDF5
  string fname("phi.hdf5");
  Vector<string> varNames(1, string("phi"));

  Real bogusVal = 4586.0;

  WriteAMRHierarchyHDF5(fname,
                        amrGrids,
                        phi,
                        varNames,
                        amrDomains[0].domainBox(),
                        amrDx[0],
                        bogusVal,
                        bogusVal,
                        refRatios,
                        numLevels);


#endif // end if HDF5

  // clean up
  for (int lev=0; lev<phi.size(); lev++)
  {
    delete phi[lev];
    delete rhs[lev];
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

    Chombo::pout() << "Running old_helmholtz for DIM= " <<  Chombo::SpaceDim << endl;
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
