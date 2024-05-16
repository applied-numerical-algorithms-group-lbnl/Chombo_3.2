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

#include "UsingNamespace.H"

enum muTypes {unityMu = 0,
              linearMu,
              GlensLawMu,
              numMuTypes};

enum solverTypes {multigrid = 0,
                  biCGStab = 1,
                  numSolverTypes = 2};

// define the test problem
//int muType = linearMu;
int muType = unityMu;
//int solverType = multigrid;
int solverType = biCGStab;
int thicknessType = constantThickness;
int basalType = constantZb;
int phiType = constZero;


void
setMu(LevelData<FluxBox>& a_mu,
      const PoissonParameters& a_params,
      const RealVect& a_dx)
{
  DataIterator dit = a_mu.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisMu = a_mu[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& dirFlux = thisMu[dir];
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
              D_TERM(pos[0] = a_dx[0]*(iR+offsets[0]);,
                     pos[1] = a_dx[1]*(jR+offsets[1]);,
                     pos[2] = a_dx[2]*(kR+offsets[2]));
              if (muType == unityMu)
                {
                  dirFluxR = 1.0;
                }
              else if (muType == linearMu)
                {
                  dirFluxR = D_TERM(pos[0], +pos[1], +pos[2]);
                }
              else
                {
                  MayDay::Error("bad muType");
                }
            }EndFor
       } // end loop over directions
    }
}



/******/
void getError(Vector< LevelData<FArrayBox>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const PoissonParameters&             a_params)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  RealVect domainSize = a_params.domainSize;
  a_error.resize(nlevels);
  Vector<LevelData<FArrayBox>* > phiExac(nlevels, NULL);
  Vector<LevelData<FArrayBox>* > klpExac(nlevels, NULL);
  Vector<RefCountedPtr<LevelData<FluxBox> > > mu(nlevels);
  Vector<RefCountedPtr<LevelData<SigmaCS> > > coordSys(nlevels);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      // note 2 components regardless of dimensionality
      a_error[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 2,IntVect::Unit);
      klpExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 2,IntVect::Zero);
      phiExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 2,IntVect::Unit);
      mu[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(a_grids[ilev], 1, IntVect::Zero));
      // nComp for LevelData<SigmaCS> is meaningless, but required by LevelData API
      IntVect ghostVect = IntVect::Unit;
      coordSys[ilev] = RefCountedPtr<LevelData<SigmaCS> >(new LevelData<SigmaCS>(a_grids[ilev],1, ghostVect));
      
      setMu(*mu[ilev], a_params, dxLev);


      DataIterator dit = a_grids[ilev].dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box ghostBox(a_grids[ilev][dit]);
          ghostBox.grow(ghostVect);
          SigmaCS& localCS = (*coordSys[ilev])[dit];
          localCS.define(ghostBox, dxLev);
        }
      
      defineSigmaCS(*coordSys[ilev], domainSize, thicknessType, basalType);
                                                           
      
      // this winds up being the initial guess for solve
      Real initialVal = 1.0;
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_error[ilev])[dit()].setVal(initialVal);
        }

      //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      setPhi( *phiExac[ilev], dxLev, a_params, phiType);
      setLOfPhi (*klpExac[ilev], dxLev, a_params, phiType);
      
      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  // set up solver
  if (solverType == multigrid)
    {

      AMRMultiGrid<LevelData<FArrayBox> > solver;
      BiCGStabSolver<LevelData<FArrayBox> >   bottomSolver;
      bottomSolver.m_verbosity = a_params.verbosity-3;
      
      defineMGSolver(solver, a_grids, mu, coordSys, bottomSolver, a_params);
      int lbase = 0;
      bool zeroInitialGuess = false;
      solver.solve(a_error, klpExac, a_params.numLevels-1, lbase, zeroInitialGuess);

    }
  else if (solverType == biCGStab)
    {
      
      BiCGStabSolver<Vector<LevelData<FArrayBox>* > > solver;
      // AMRMG preconditioner
      AMRMultiGrid<LevelData<FArrayBox> > mgSolver;
      BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
      defineMGSolver(mgSolver, a_grids, mu, coordSys, bottomSolver, a_params);

      Vector<ProblemDomain> vectDomain(a_params.numLevels);
      MultilevelLinearOp<FArrayBox>* mlOpPtr = NULL;
      mlOpPtr = defineBiCGStabSolver(solver, a_grids, mu, coordSys, mgSolver, a_params);

      solver.solve(a_error, klpExac);         

      delete mlOpPtr;

    }


  //create calculated data and error
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_error[ilev])[dit()] -= (*phiExac[ilev])[dit()];
        }
    }
  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete phiExac[ilev];
      delete klpExac[ilev];
    }

}
/***************/
void outputError(const Vector< LevelData<FArrayBox>* >&   a_errorFine,
                 const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
                 const Vector< DisjointBoxLayout >&       a_gridsFine,
                 const Vector< DisjointBoxLayout >&       a_gridsCoar,
                 const PoissonParameters&                 a_paramsFine,
                 const PoissonParameters&                 a_paramsCoar)
{
#if CH_SPACEDIM==2
    string fileFine("pltFineError.2d.hdf5");
    string fileCoar("pltCoarError.2d.hdf5");
#else
    string fileFine("pltFineError.3d.hdf5");
    string fileCoar("pltCoarError.3d.hdf5");
#endif
  string phiname("error");
  outputData(a_errorFine, a_gridsFine,
             a_paramsFine.coarsestDomain, a_paramsFine.refRatio,
             a_paramsFine.coarsestDx, a_paramsFine.numLevels,
             fileFine, phiname);
  outputData(a_errorCoar, a_gridsCoar,
             a_paramsCoar.coarsestDomain, a_paramsCoar.refRatio,
             a_paramsCoar.coarsestDx, a_paramsCoar.numLevels,
             fileCoar, phiname);
}
int main(int argc, char* argv[])
{
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
        pout() << "using linearSolver.inputs" << endl;
        inFile = "linearSolver.inputs";
      }

    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters paramFine, paramCoar;
    Vector<DisjointBoxLayout> gridsFine, gridsCoar;

    //read params from file
    getPoissonParameters(paramFine);
    paramCoar = paramFine;
    paramCoar.coarsen(2);
    int nlevels = paramCoar.numLevels;
    Vector<LevelData<FArrayBox>* > errorFine(nlevels, NULL);
    Vector<LevelData<FArrayBox>* > errorCoar(nlevels, NULL);

    setGrids(gridsFine,  paramFine);

    pout() << "generating fine error" << endl;
    getError(errorFine, gridsFine, paramFine);

    getCoarseLayoutsFromFine(gridsCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, paramCoar);

    int dofileout;
    pp.get("do_error_output", dofileout);
    if(dofileout == 1)
      {
        outputError(errorFine,   errorCoar,
                    gridsFine,   gridsCoar,
                    paramFine,   paramCoar);
      }

    string testname("Solution error");
    compareError(errorFine,   errorCoar,
                 gridsFine,   gridsCoar,
                 paramFine,   paramCoar,
                 testname);

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete errorFine[ilev];
        delete errorCoar[ilev];
        errorFine[ilev] = NULL;
        errorCoar[ilev] = NULL;
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
