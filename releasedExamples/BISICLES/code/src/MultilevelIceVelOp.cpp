#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRIO.H"
#include "MultilevelIceVelOp.H"
#include "RelaxSolver.H"

#ifdef CH_USE_PETSC
#include "PetscSolver.H"
#endif

#include "NamespaceHeader.H"

// set a couple of defaults here (do before define so that
// they can be overridden before calling define)

MultilevelIceVelOp::MultilevelIceVelOp()
{

  // set some default values
  // default preconditioner is still bicgstab
  m_bottom_solver_type = bicg;
}

/// define function
void
MultilevelIceVelOp::define(const Vector<DisjointBoxLayout>& a_vectGrids,
                           const Vector<int>& a_refRatios,
                           const Vector<ProblemDomain>& a_domains,
                           const Vector<RealVect>& a_vectDx,
                           RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >& a_opFactory,
                           int a_lBase)
{
  CH_TIME("MultilevelIceVelOp::define");
  
  int numLevels = a_vectGrids.size();
  
  // couple of sanity checks...
  CH_assert (a_lBase >= 0);
  CH_assert (a_lBase < numLevels);
  CH_assert ( a_domains.size() == numLevels);
  // since you technically only need numLevels -1 refinement ratios...
  CH_assert (a_refRatios.size() >= numLevels -1);
  
  
  m_lBase = a_lBase;
  
  m_vectGrids = a_vectGrids;
  m_refRatios = a_refRatios;
  m_domains = a_domains;
  m_vectDx = a_vectDx;
  m_vectOperators.resize(numLevels);
  
  // need to at least define the levels to lBase -1 (for coarser BC's)
  int coarsestLevel = Max(0, m_lBase-1);
  for (int level=coarsestLevel; level<numLevels; level++)
    {
      RefCountedPtr<AMRLevelOp<LevelData<FArrayBox> > > op = RefCountedPtr<AMRLevelOp<LevelData<FArrayBox> > >(a_opFactory->AMRnewOp(m_domains[level]));
      m_vectOperators[level] = op;
    }
  
  // finally, define AMRMultigrid preconditioner if required
  if (m_use_multigrid_preconditioner)
    {
      
      if (m_bottom_solver_type == bicg)
        {
          // preconditioner requires a bottom smoother
          BiCGStabSolver<LevelData<FArrayBox> >* newPtr = new BiCGStabSolver<LevelData<FArrayBox> >;
          
          m_precondBottomSolverPtr = newPtr;
          int bottomSolverVerbosity = 1;
          
          newPtr->m_verbosity = bottomSolverVerbosity;
        }
#ifdef CH_USE_PETSC
      else if (m_bottom_solver_type == PETSC)
        {
          PetscSolverViscousTensor<LevelData<FArrayBox> >* petscSolverPtr = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
          m_precondBottomSolverPtr = petscSolverPtr;
        }
#endif //if CH_USE_PETSC
      else if (m_bottom_solver_type == relax)
        {
          RelaxSolver<LevelData<FArrayBox> >* relaxSolverPtr = new RelaxSolver<LevelData<FArrayBox> >;
          m_precondBottomSolverPtr = relaxSolverPtr;
        }       
      else
        {
          MayDay::Error("bad elliptic solver type");
        }
      
      if (m_preCondSolverDepth >= 0)
        {
          m_preCondSolver.m_maxDepth = m_preCondSolverDepth;
        }
      m_preCondSolver.define(m_domains[0],
                             *a_opFactory,
                             m_precondBottomSolverPtr,
                             numLevels);
      
      // most of these have no real meaning since we're only doina a
      // single V-cycle, rather than a full solve
      Real solverEps = 1.0e-7;
      int maxIterations = 1;
      Real hang = 1.0e-7;
      Real normThresh = 1.0e-30;
      
      int num_mg = 1;
      m_preCondSolver.setSolverParameters(m_num_mg_smooth,
                                          m_num_mg_smooth,
                                          m_num_mg_smooth,
                                          num_mg,
                                          maxIterations,
                                          solverEps,
                                          hang,
                                          normThresh);
      
      // set preconditioner solver to be _really_ quiet...
      m_preCondSolver.m_verbosity = 1;
      
      // AMRMultiGrid::init has not yet been called
      // (will call later, during actual solve)
      m_isPrecondSolverInitialized = false;
      
    }
}

///
MultilevelIceVelOp::~MultilevelIceVelOp()
{
}

#include "NamespaceFooter.H"

