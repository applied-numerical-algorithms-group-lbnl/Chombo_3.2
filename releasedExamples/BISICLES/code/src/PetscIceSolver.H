#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _PETSCICESOLVER_H_
#define _PETSCICESOLVER_H_

#include "IceVelocitySolver.H"
#include "ViscousTensorOp.H"
#ifdef CH_USE_PETSC
#include "PetscSolver.H"
#endif
#include "Copier.H"


#include "NamespaceHeader.H"

/// virtual class to manage the nonlinear solve for the ice-sheet
/// momentum 
/** This is a virtual base class to facilitate the implementation of
    different solution algorithms.
*/
class PetscIceSolver : public IceVelocitySolver
{

public:

  PetscIceSolver(); 

  virtual ~PetscIceSolver();

  virtual void define(const ProblemDomain& a_coarseDomain,
                      ConstitutiveRelation* a_constRel,
		      BasalFrictionRelation* a_basalFrictionRel,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_vectRefRatio,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels);

  void setSolverType(int a_solver_type) 
  {
    pout() << "PetscIceSolver::setSolverType=" << a_solver_type << endl;
  }

  void setTolerance(Real a_tolerance) {
    m_rtol = a_tolerance; 
  }

  // set "absolute tolerance"
  /** If norm(residual) drops below the absolute tolerance (regardless 
      of the initial residual), we declare ourselves "done" -- this is a 
      simple way to avoid solver hangs.
  */
  void setAbsoluteTolerance(Real a_tolerance){ m_atol = a_tolerance; }

  // sets maximum number of iterations
  void setMaxIterations(int a_max_iter) { m_max_its = a_max_iter; }

  void setVerbosity(int a_verbosity) {m_verbosity = a_verbosity;}

  void computeAMRResidual( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_resid, 
			   const Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
			   const Vector<LevelData<FArrayBox>* >& a_rhs,
			   int a_lbase, int a_maxLevel
			   );
  void computeAMRResidualLevel( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_resid, 
				const Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
				const Vector<LevelData<FArrayBox>* >& a_rhs,
				int a_lbase, int a_maxLevel, int a_ilev
				);
  
  // solve for isothermal ice -- returns 0 if converged, 1 if not
  virtual int solve( Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		     Vector<LevelData<FArrayBox>* >& a_calvedIce,
		     Vector<LevelData<FArrayBox>* >& a_addedIce,
		     Vector<LevelData<FArrayBox>* >& a_removedIce,
		     Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		     const Real a_convergenceMetric,
		     const Vector<LevelData<FArrayBox>* >& a_rhs,
		     const Vector<LevelData<FArrayBox>* >& a_C,
		     const Vector<LevelData<FArrayBox>* >& a_C0,
		     const Vector<LevelData<FArrayBox>* >& a_A,
		     const Vector<LevelData<FArrayBox>* >& a_muCoef,
		     Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		     Real a_time,
		     int a_lbase, int a_maxLevel);

  void picardSolve_private( int a_ilev,
			    LevelData<FArrayBox> &a_horizontalVel,
			    const LevelData<FArrayBox> &a_rhs,
			    Real a_norm0, Real &a_norm, int a_numIts, int &a_it );

  void jfnkSolve_private( int a_ilev,
			  LevelData<FArrayBox> &a_horizontalVel,
			  const LevelData<FArrayBox> &a_rhs,
			  Real a_norm0, int a_numIts, int &a_it );
  
  Vector<ProblemDomain> m_domains;
  RefCountedPtr<ViscousTensorOpFactory> m_opFactoryPtr;
  Vector<RefCountedPtr<ViscousTensorOp> > m_op;
protected:
  ConstitutiveRelation* m_constRelPtr;
  BasalFrictionRelation* m_basalFrictionRelPtr;
  Vector<DisjointBoxLayout> m_grids;
  Real m_vtopSafety;
  int m_verbosity;
  bool m_isOpDefined;
  bool m_isSolverDefined;

  IceThicknessIBC* m_bc;

  int m_max_its;
  Real m_rtol;
  Real m_atol;
  int m_minPicardIterations;
  int m_amr_max_level;
  bool m_plotResidual;

  Vector<int> m_refRatios;

protected:
  // materials
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_Mu;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_Lambda;
  // beta is used to scale C below
  Vector<const LevelData<FArrayBox> *> m_Beta;
  Vector<const LevelData<FArrayBox> *> m_Beta0;
  // C here is the sliding coefficient -- the "acoeff" in terms of the viscousTensorOp
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_C;

  Real m_constThetaVal; // not used?

  /// define operator factory. In 2D, this is a ViscousTensorOp.
  void defineOpFactory(RealVect,const ProblemDomain &,int);

  /// get ViscousTensorOp scaling factors (alpha and beta)
  /** put this in a function because we need them in two separate places */
  void getOperatorScaleFactors( Real& a_alpha, Real& a_beta ) const;
public:
  // compute face-centered coefficients for tensor solver (really
  // winds up being H*mu) -- non-isothermal version...
  void computeMu( LevelData<FArrayBox> &a_horizontalVel, 
		  const LevelData<FluxBox> &a_A, 
		  const LevelData<FluxBox> &a_muCoef,
		  const RefCountedPtr<LevelSigmaCS> &a_coordSys,
		  LevelData<FArrayBox>* crseVelPtr,
		  LevelData<FArrayBox>* fineVelPtr,
		  int a_ilev,
		  Real a_time);
  // abstract methods for (PETSc) nonlinear solvers
  virtual void updateCoefs(LevelData<FArrayBox> &a_horizontalVel, int);

  // cache 
  RefCountedPtr<LevelData<FArrayBox> >m_twork1;  // work vectors
  RefCountedPtr<LevelData<FArrayBox> >m_twork2; 
  RefCountedPtr<LevelData<FluxBox> >  m_tfaceA;
  LevelData<FArrayBox>* m_tcrseVel;
  LevelData<FArrayBox> *m_tphi0; 
  RefCountedPtr<LevelData<FluxBox> > m_tmuCoef;
  RefCountedPtr<LevelSigmaCS> m_tcoordSys;
  Real m_ttime;
  // prolongation stuff
  Vector<RefCountedPtr<LevelData<FArrayBox> > >m_fineCover;
  Vector<Copier>                               m_restCopier;
  Vector<Copier>                               m_projCopier;
  void AMRProlong( LevelData<FArrayBox>&       a_fineU,
		   const LevelData<FArrayBox>& a_CrsU,
		   LevelData<FArrayBox>&       a_CrsCover,
		   Copier a_copier,
		   int a_refRatio );
  
};
#include "NamespaceFooter.H"
#endif
