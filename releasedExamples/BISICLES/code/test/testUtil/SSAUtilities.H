#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _SSAUTILITIES_H_
#define _SSAUTILITIES_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "MultilevelLinearOp.H"
#include "RelaxSolver.H"
#include "SSAVelocityOp.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "functionsF_F.H"
//#include "LayoutData.H"
#include "SigmaCS.H"

#include "NamespaceHeader.H"

enum phiTypes{constZero = 0,
              constOne  = 1,
              sinusoidal, 
              num_phitype};

///
/**
 */
class PoissonParameters
{
public:
  IntVect      nCells;
  int          maxGridSize;
  int          blockFactor;
  int          bufferSize;
  Real         fillRatio;
  Real         refineThresh;
  int          domBcType;
  int          maxLevel;
  int          baseLevel;
  int          numLevels;
  Vector<int>  refRatio;
  ProblemDomain coarsestDomain;
  RealVect     coarsestDx;
  RealVect     domainSize;
  RealVect     probLo;
  RealVect     probHi;
  int          coefficient_average_type;
  void coarsen(int a_factor);
  void  refine(int a_factor);
  int verbosity;
};


///
/**
 */
class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;
  static std::vector<int> s_bcLo, s_bcHi;
  static bool s_areBCsParsed, s_valueParsed;
};


extern void ParseValue(Real* pos,
                       int* dir,
                       Side::LoHiSide* side,
                       Real* a_values);

extern void ParseBC(FArrayBox& a_state,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    bool a_homogeneous);



extern void
outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain& coarsestDomain,
           const Vector<int>& vectRatio,
           RealVect dxCoarsest,
           int numlevels,
           string filename,
           string varname);


// simplest just to overload this
extern void
outputData(const Vector<RefCountedPtr<LevelData<FArrayBox> > >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain& coarsestDomain,
           const Vector<int>& vectRatio,
           RealVect dxCoarsest,
           int numlevels,
           string filename,
           string varname);

///get stuff from input file
extern void getPoissonParameters (PoissonParameters&  a_params);

extern int setGrids(Vector<DisjointBoxLayout>& vectGrids,
                    PoissonParameters&         a_params);

/// if numlevels is negative, use numlevels from params
extern int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
                  PoissonParameters&              a_params,
                  int a_phitype,
                  int numlevels = -1);


void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<RealVect>&          vectDx,
                         PoissonParameters&         a_params);

extern void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<RealVect>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels) ;

///
/**
*/
extern void compareError(const Vector< LevelData<FArrayBox>* >&   a_errorFine,
                         const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
                         const Vector< DisjointBoxLayout >&       a_gridsFine,
                         const Vector< DisjointBoxLayout >&       a_gridsCoar,
                         const PoissonParameters&                 a_paramsFine,
                         const PoissonParameters&                 a_paramsCoar,
                         const string& a_testName);


///
/**
*/
extern void compareError(const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_errorFine,
                         const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_errorCoar,
                         const Vector< DisjointBoxLayout >&       a_gridsFine,
                         const Vector< DisjointBoxLayout >&       a_gridsCoar,
                         const PoissonParameters&                 a_paramsFine,
                         const PoissonParameters&                 a_paramsCoar,
                         const string& a_testName);

///
/**
*/
extern void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                                     const Vector<DisjointBoxLayout>& a_gridsFine,
                                     const PoissonParameters&         a_params);

/********/
extern void
defineMGSolver(AMRMultiGrid<LevelData<FArrayBox> >&         a_solver,
               const Vector<DisjointBoxLayout>&             a_grids,
               Vector<RefCountedPtr<LevelData<FluxBox> > >& a_mu,
               Vector<RefCountedPtr<LevelData<SigmaCS> > >& a_coordSys,
               LinearSolver<LevelData<FArrayBox> >&         a_bottomSolver,
               const PoissonParameters&                     a_params);


extern 
MultilevelLinearOp<FArrayBox>* 
defineBiCGStabSolver(BiCGStabSolver<Vector<LevelData<FArrayBox>* > >&  a_solver,
                     const Vector<DisjointBoxLayout>&             a_grids,
                     Vector<RefCountedPtr<LevelData<FluxBox> > >& a_mu,
                     Vector<RefCountedPtr<LevelData<SigmaCS> > >& a_coordSys,
                     AMRMultiGrid<LevelData<FArrayBox> >& a_precondSolver,
                     const PoissonParameters&                     a_params);


extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory(
                      const Vector<DisjointBoxLayout>&             a_grids,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_mu,
                      Vector<RefCountedPtr<LevelData<SigmaCS> > >& a_coordSys,
                      const PoissonParameters&                     a_params);


///
/**
*/
extern void setPhi(LevelData<FArrayBox>&     a_phi,
                   const  RealVect&          a_dx,
                   const  PoissonParameters& a_params,
                   int a_phiType);

extern void setLOfPhi(LevelData<FArrayBox>&    a_kappaLOfPhi,
                      const RealVect&          a_dx,
                      const PoissonParameters& a_params,
                      int a_phiType);

#include "NamespaceFooter.H"

#endif
