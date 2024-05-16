#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "PetscIceSolver.H"
#include "ViscousTensorOp.H"
#include "ParmParse.H"
#include "CoarseAverageFace.H"
#include "IntInterpF_F.H"
#include "IceConstants.H"
#include "BisiclesF_F.H"
#include "TensorCFInterp.H"
#include "AMRIO.H"
#include "JFNKSolver.H"
#include "computeNorm.H"
#ifdef CH_USE_PETSC
#include "PetscSolver.H"
#endif

#include "NamespaceHeader.H"

////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::PetscIceSolver()
////////////////////////////////////////////////////////////////////////
PetscIceSolver::PetscIceSolver()
{
  // default constructor leaves things in an undefined state
  m_bc = NULL;
  m_constRelPtr = NULL;
  m_basalFrictionRelPtr = NULL;

  // use isothermal ice temp from Pattyn(2003)
  m_constThetaVal = 238.15;

  m_isOpDefined = false;
  m_vtopSafety = VTOP_DEFAULT_SAFETY;

  m_max_its = 20;
  m_rtol = 1.e-6;
  m_atol = 1.e-30;
  m_minPicardIterations = 3;
  m_plotResidual = false;

  ParmParse pp("petsc");
  pp.query("maxIter",m_max_its);
  pp.query("absNLTol",m_atol);
  pp.query("relNLTol",m_rtol);
  pp.query("minPicardIterations",m_minPicardIterations);
  pp.query("plotResidual",m_plotResidual);

  ParmParse pp2("amr");
  pp2.get("maxLevel",m_amr_max_level);

  // these ones don't need to be stored (at least for now), but should be set
  int mgAverageType  = CoarseAverageFace::arithmetic;
  ViscousTensorOpFactory::s_coefficientAverageType = mgAverageType;

  // set default to be linear prolongation in multigrid
  int mgProlongType = ViscousTensorOp::linearInterp;
  ViscousTensorOp::s_prolongType = mgProlongType;

#ifdef CHOMBO_TRUNK  
  /// default is "lazy gsrb"
  ParmParse pp3("solver");
  if (pp3.contains("lazyGSRB"))
    {
      bool lazyGSRB;
      pp.get("lazyGSRB", lazyGSRB);
      ViscousTensorOp::s_lazy_gsrb = lazyGSRB;
    }
#endif
  
}
////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::~PetscIceSolver() 
////////////////////////////////////////////////////////////////////////
PetscIceSolver::~PetscIceSolver() 
{
  if (m_bc != NULL)
    {
      delete m_bc;
      m_bc = NULL;
    }
}

#ifdef CH_USE_PETSC
/* ------------------------------------------------------------------- */
/* 
   FormFunction

   Input Parameters:
   user - user-defined application context
.  x - vector
.  ctx - user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  f - vector
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction( SNES snes, Vec x, Vec f, void *ctx )
{
  CH_TIME("PetscIceSolver::FormFunction");
  PetscErrorCode ierr;
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver;
  PetscIceSolver *tthis;
  PetscInt *pilev;

  PetscFunctionBegin;

  ierr = SNESGetApplicationContext(snes,(void**)&pilev); CHKERRQ(ierr);

  solver = (PetscSolverViscousTensor<LevelData<FArrayBox> >*)ctx;
  tthis = (PetscIceSolver*)solver->m_ctx;

  ierr = solver->putPetscInChombo( *tthis->m_twork2, x );     CHKERRQ(ierr);

  if ( tthis->m_tphi0 ) tthis->m_op[*pilev]->incr( *tthis->m_twork2, *tthis->m_tphi0, 1.);
  tthis->updateCoefs( *tthis->m_twork2, (int)*pilev ); 
  if ( tthis->m_tphi0 ) tthis->m_op[*pilev]->incr( *tthis->m_twork2, *tthis->m_tphi0, -1.);

  tthis->m_op[*pilev]->applyOp( *tthis->m_twork1, *tthis->m_twork2, true ); 

  ierr = solver->putChomboInPetsc( f, *tthis->m_twork1 );  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
   FormJacobian - this does not depend on PetscIceSolver so it could be moved into Mat solver.

   This can go into the PetscSolver class -- it is not BICICLES specific!!!

   Input Parameters:
.  snes - the SNES context
.  dummy - input vector
.  ctx - optional user-defined context, as set by SNESSetJacobian()

   Output Parameters:
.  jac - Jacobian matrix
.  prejac - different preconditioning matrix
.  flag - flag indicating matrix structure
*/
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"

#if PETSC_VERSION_LT(3,5,0)
PetscErrorCode FormJacobian( SNES snes,Vec x,Mat *jac,Mat *prejac,MatStructure *flag, void *ctx )
{

 CH_TIME("PetscIceSolver::FormJacobian");
  PetscErrorCode ierr;
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver;
  PetscIceSolver *tthis;
  PetscInt *pilev; // not used

  PetscFunctionBegin;
  
  ierr = SNESGetApplicationContext(snes,(void**)&pilev); CHKERRQ(ierr);
  
  solver = (PetscSolverViscousTensor<LevelData<FArrayBox> >*)ctx;
  tthis = (PetscIceSolver*)solver->m_ctx;

  // form Function was just called so do not need to update coefs
  ierr = solver->formMatrix( *prejac ); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (prejac!=jac)
    {
      ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
  *flag = SAME_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}
#else
PetscErrorCode FormJacobian( SNES snes,Vec x,Mat jac,Mat prejac, void *ctx )
{
  CH_TIME("PetscIceSolver::FormJacobian");
  PetscErrorCode ierr;
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver;
  PetscIceSolver *tthis;
  PetscInt *pilev; // not used

  PetscFunctionBegin;
  
  ierr = SNESGetApplicationContext(snes,(void**)&pilev); CHKERRQ(ierr);
  
  solver = (PetscSolverViscousTensor<LevelData<FArrayBox> >*)ctx;
  tthis = (PetscIceSolver*)solver->m_ctx;

  // form Function was just called so do not need to update coefs
  ierr = solver->formMatrix( prejac ); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(prejac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (prejac!=jac)
    {
      ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
  PetscFunctionReturn(0);
}
#endif // petsc version
#endif // petsc

////////////////////////////////////////////////////////////////////////
//  PetscIceSolver::define
////////////////////////////////////////////////////////////////////////
void 
PetscIceSolver::define(const ProblemDomain& a_coarseDomain,
		       ConstitutiveRelation* a_constRelPtr,
		       BasalFrictionRelation* a_basalFrictionRelPtr,
		       const Vector<DisjointBoxLayout>& a_vectGrids,
		       const Vector<int>& a_vectRefRatio,
		       const RealVect& a_dxCrse,
		       IceThicknessIBC* a_bc,
		       int a_numLevels)
{
  Real opAlpha, opBeta;

  getOperatorScaleFactors( opAlpha, opBeta );

  m_constRelPtr = a_constRelPtr;
  m_basalFrictionRelPtr = a_basalFrictionRelPtr;
  m_bc = a_bc->new_thicknessIBC();

  m_op.resize(a_numLevels);
  m_grids.resize(a_numLevels);
  m_domains.resize(a_numLevels);
  m_refRatios.resize(a_numLevels-1);
  m_fineCover.resize(a_numLevels-1);
  m_projCopier.resize(a_numLevels-1);
  m_restCopier.resize(a_numLevels-1);
  m_Mu.resize(a_numLevels);
  m_Lambda.resize(a_numLevels);
  m_Beta.resize(a_numLevels);
  m_Beta0.resize(a_numLevels);
  m_C.resize(a_numLevels);

  m_domains[0] = a_coarseDomain;
  for (int ilev=0;ilev<a_numLevels;ilev++)
    {
      m_grids[ilev] = a_vectGrids[ilev];
      if ( ilev < a_numLevels-1 )
	{
	  m_refRatios[ilev] = a_vectRefRatio[ilev];
	}

      m_Mu[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids[ilev], 
									     1, 
									     IntVect::Zero) );
      
      m_Lambda[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids[ilev],
										 1, 
										 IntVect::Zero) );      
      // C only has one component...
      m_C[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids[ilev], 
										1,
										IntVect::Zero));
      // not coarsest grid
      if (ilev>0)
	{
	  // make domain
	  m_domains[ilev] = m_domains[ilev-1];
	  m_domains[ilev].refine(m_refRatios[ilev-1]);

	  // make prologation stuff (index to coarse grid)
	  const int nc = 2;
	  const DisjointBoxLayout& finedbl = m_grids[ilev];
	  DisjointBoxLayout dblCoarsenedFine;
	  coarsen( dblCoarsenedFine, finedbl, m_refRatios[ilev-1]);	  
	  m_fineCover[ilev-1]=RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(dblCoarsenedFine,nc,IntVect::Unit));
	  // prolongator copier, from data to cover
	  m_projCopier[ilev-1].define(m_grids[ilev-1],dblCoarsenedFine,IntVect::Unit);
	  // restrict copier used for zero cover
	  m_restCopier[ilev-1].define(dblCoarsenedFine,m_grids[ilev-1], IntVect::Zero);
	}
    }
  
  // create op factory
  defineOpFactory( a_dxCrse, a_coarseDomain, a_numLevels );
  
  // create ops
  for (int ilev=0;ilev<a_numLevels;ilev++)
    {
      // this copies the unset data above, just needed here for dx &crdx.
      m_op[ilev] = RefCountedPtr<ViscousTensorOp>(m_opFactoryPtr->AMRnewOp(m_domains[ilev])); 
    }
}

// compute residaul on all levels, sets C-F ghosts and coarsens fine solutions
//
//
void
PetscIceSolver::computeAMRResidual( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_resid, 
				    const Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
				    const Vector<LevelData<FArrayBox>* >& a_rhs,
				    int a_lbase, int a_maxLevel
				    )
{
  for (int ilev=a_lbase;ilev<=a_maxLevel;ilev++)
    {
      computeAMRResidualLevel( a_resid, a_horizontalVel, a_rhs, 
			       a_lbase, a_maxLevel, ilev );
    }
}

void
PetscIceSolver::computeAMRResidualLevel( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_resid, 
					 const Vector<LevelData<FArrayBox>* >& a_horizontalVel, 
					 const Vector<LevelData<FArrayBox>* >& a_rhs,
					 int a_lbase, int a_maxLevel, int a_ilev
					 )
{
  if (a_ilev==a_lbase) 
    {
      if (a_ilev==a_maxLevel) // one level
	{
	  m_op[a_ilev]->residual(*a_resid[a_ilev], 
				 *a_horizontalVel[a_ilev], 
				 *a_rhs[a_ilev], true ); 
	}
      else // coarse grid
	{
	  m_op[a_ilev]->AMRResidualNC( *a_resid[a_ilev], 
				       *a_horizontalVel[a_ilev+1], 
				       *a_horizontalVel[a_ilev], 
				       *a_rhs[a_ilev], true, 
				       &(*m_op[a_ilev+1]));
	}
    }
  else if (a_ilev == a_maxLevel) // fine grid
    {
      m_op[a_ilev]->AMRResidualNF( *a_resid[a_ilev], 
				   *a_horizontalVel[a_ilev],
				   *a_horizontalVel[a_ilev-1],
				   *a_rhs[a_ilev], true );
    }
  else
    {
      m_op[a_ilev]->AMRResidual( *a_resid[a_ilev], 
				 *a_horizontalVel[a_ilev+1], 
				 *a_horizontalVel[a_ilev], 
				 *a_horizontalVel[a_ilev-1], 
				 *a_rhs[a_ilev], true, 
				 &(*m_op[a_ilev+1]));
    }
}

void PetscIceSolver::AMRProlong( LevelData<FArrayBox>&       a_fu,
				 const LevelData<FArrayBox>& a_cu,
				 LevelData<FArrayBox>&       a_CrsCover,
				 Copier a_copier,
				 int a_refRatio
				 )
{
CH_TIME("PetscIceSolver::AMRProlong");
  
  DisjointBoxLayout dbl = a_fu.disjointBoxLayout();
  DisjointBoxLayout cdbl = a_CrsCover.disjointBoxLayout();
  
  a_cu.copyTo(a_CrsCover.interval(), a_CrsCover, a_CrsCover.interval(), a_copier);
  
  for ( DataIterator dit = a_fu.dataIterator(); dit.ok(); ++dit )
    {
      FArrayBox& phi =  a_fu[dit];
      FArrayBox& coarse = a_CrsCover[dit];
      Box region = dbl[dit];

      FORT_PROLONGQUAD_ICE(CHF_FRA(phi),
			   CHF_CONST_FRA(coarse),
			   CHF_BOX(region),
			   CHF_CONST_INT(a_refRatio));
    }
}

// Picard solve in residual correction form
void PetscIceSolver::picardSolve_private( int a_ilev,
					  LevelData<FArrayBox> &a_horizontalVel,
					  const LevelData<FArrayBox> &a_rhs,
					  Real a_norm0, 
					  Real &a_norm, // out
					  int a_numIts, 
					  int &a_it   // in-out
					  )
{
#ifdef CH_USE_PETSC
  Real opAlpha, opBeta; 
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver = new PetscSolverViscousTensor<LevelData<FArrayBox> >();
  getOperatorScaleFactors( opAlpha, opBeta );
  solver->define( &(*m_op[a_ilev]), true ); // dx & crdx
  solver->define( opAlpha, opBeta,
		       &(*m_C[a_ilev]),
		       &(*m_Mu[a_ilev]),
		       &(*m_Lambda[a_ilev]) );
  solver->m_ctx = (void*)this;

  // in residual correction form (m_tphi0) and first iteration - no init guess.
  solver->setInitialGuessNonzero(true);
  
  for (int ii=0;a_it<a_numIts;a_it++,ii++)
    {
      // update coeficients
      updateCoefs( a_horizontalVel, a_ilev ); 
      
      // linear KSP solve
      solver->solve(a_horizontalVel,a_rhs);
      KSPGetResidualNorm(solver->m_ksp,&a_norm);
      
      if (m_verbosity>0)
	{
	  pout() << a_it+1 << "/" << m_max_its <<  ") Picard iteration, |r|_2=" << a_norm << ", rate=" << a_norm/a_norm0 << endl;
	}
      if (a_norm/a_norm0 < m_rtol){ a_it++; break; }
      solver->resetOperator();
    }
  delete solver;
#endif
}

void PetscIceSolver::jfnkSolve_private( int a_ilev,
					LevelData<FArrayBox> &a_horizontalVel,
					const LevelData<FArrayBox> &a_rhs,
					Real a_norm0, int a_numIts, int &a_it // in-out
					)
{  
#ifdef CH_USE_PETSC
  Real opAlpha, opBeta; PetscInt its;
  getOperatorScaleFactors( opAlpha, opBeta );
  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
  solver->define( &(*m_op[a_ilev]), true ); // dx & crdx
  solver->define( opAlpha, opBeta,
		       &(*m_C[a_ilev]),
		       &(*m_Mu[a_ilev]),
		       &(*m_Lambda[a_ilev]) );
  solver->setFunctionAndJacobian( FormFunction, FormJacobian ); // NL solve
  solver->m_ctx = (void*)this;
  
  // creates Mat and Vecs, creates SNES
  solver->setup_solver(a_horizontalVel); // creates SNES
  // set the level to index into this
  SNESSetTolerances(solver->m_snes,(PetscReal)m_atol,(PetscReal)m_rtol,PETSC_DEFAULT,a_numIts-a_it,PETSC_DEFAULT);
  SNESSetApplicationContext( solver->m_snes,(void*)&a_ilev );
  solver->solve(a_horizontalVel,a_rhs);
  SNESGetIterationNumber(solver->m_snes,&its);
  a_it += its;
  delete solver;
#endif
}

/// solve for isothermal ice
/** beta scales sliding coefficient C -- acoef in terms of the ViscousTensorOp
 */
int
PetscIceSolver::solve( Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		       Vector<LevelData<FArrayBox>* >& a_calvedIce,
		       Vector<LevelData<FArrayBox>* >& a_addedIce,
		       Vector<LevelData<FArrayBox>* >& a_removedIce,
		       Real& a_initialResidualNorm, 
		       Real& a_finalResidualNorm,
		       const Real a_convergenceMetric,
		       const Vector<LevelData<FArrayBox>* >& a_rhs,
		       const Vector<LevelData<FArrayBox>* >& a_beta,
		       const Vector<LevelData<FArrayBox>* >& a_beta0, 
		       const Vector<LevelData<FArrayBox>* >& a_A,
		       const Vector<LevelData<FArrayBox>* >& a_muCoef,
		       Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		       Real a_time, int a_lbase, int a_maxLevel )
{
  CH_assert(a_lbase==0); // todo?
  CH_assert(m_isOpDefined);
  int returnCode = 0, ilev;
  const int nc = a_horizontalVel[0]->nComp();
  Real residNorm0,levelNorm;
  IntVect ghostVect = a_horizontalVel[0]->ghostVect(); // can this be Zero???

  // copy betas and form faceA, faceMuCoef
  Vector<RefCountedPtr<LevelData<FluxBox> > > faceAs(a_maxLevel+1);
  Vector<RefCountedPtr<LevelData<FluxBox> > > faceMuCoef(a_maxLevel+1);
  for (ilev=a_lbase;ilev<=a_maxLevel;ilev++)
    {
      m_Beta[ilev] = a_beta[ilev];
      m_Beta0[ilev] = a_beta0[ilev];

      // initial implementation -- redefine solver for every iteration. 
      RefCountedPtr<LevelData<FluxBox> > faceA(new LevelData<FluxBox>(m_grids[ilev], 
								      a_A[ilev]->nComp(), 
								      IntVect::Zero));
      CellToEdge(*a_A[ilev], *faceA);
      faceAs[ilev] = faceA;

      RefCountedPtr<LevelData<FluxBox> > ref(new LevelData<FluxBox>(m_grids[ilev], 
								      a_A[ilev]->nComp(), 
								      IntVect::Zero));
      CellToEdge(*a_muCoef[ilev], *ref);
      faceMuCoef[ilev] = ref;


    }

  Vector<RefCountedPtr<LevelData<FArrayBox> > > resid(a_maxLevel+1);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > tempv(a_maxLevel+1);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > tempv2(a_maxLevel+1);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > tempv3(a_maxLevel+1);
  for (ilev=a_lbase;ilev<=a_maxLevel;ilev++)
    {
      RefCountedPtr<LevelData<FArrayBox> > v1(new LevelData<FArrayBox>(m_grids[ilev],nc,IntVect::Zero));
      resid[ilev] = v1;
      RefCountedPtr<LevelData<FArrayBox> > v2(new LevelData<FArrayBox>(m_grids[ilev],nc,ghostVect));
      tempv[ilev] = v2;
      RefCountedPtr<LevelData<FArrayBox> > v3(new LevelData<FArrayBox>(m_grids[ilev],nc,IntVect::Zero));
      tempv2[ilev] = v3;
      RefCountedPtr<LevelData<FArrayBox> > v4(new LevelData<FArrayBox>(m_grids[ilev],nc,ghostVect));
      tempv3[ilev] = v4;
    }


  // coarse grid solve 
  ilev = a_lbase;

  // update coeficients to get correct residuals, sets c-f values
  computeMu( *a_horizontalVel[ilev],
	     *faceAs[ilev],
	     *faceMuCoef[ilev],
	     a_coordSys[ilev], 
	     0,
	     0,
	     ilev,
	     a_time );
  // level 0, so no need for C/F
  CH_assert(ilev == 0);
  m_op[ilev]->residual( *resid[ilev], 
			*a_horizontalVel[ilev], 
			*a_rhs[ilev], true ); 

  // cache stuff for nonlinear solver
  m_twork1 = tempv2[ilev];
  m_twork2 = tempv3[ilev];
  m_tfaceA = faceAs[ilev];
  m_tmuCoef = faceMuCoef[ilev];
  m_tcoordSys = a_coordSys[ilev];
  m_ttime = a_time;
  m_tphi0 = 0;    // not defect correction form
  m_tcrseVel = 0;

  // start with Picard solve
  residNorm0 = m_op[ilev]->norm(*a_rhs[ilev], 2); // use |b| to compare against
  if (m_verbosity>0)
    {
      pout() << "\t\tLevel 0: |b|_2 = " << residNorm0 << endl;
    }  

  // solve
  int it = 0;
  picardSolve_private(ilev,*a_horizontalVel[ilev],*a_rhs[ilev],residNorm0,levelNorm,m_minPicardIterations,it);

  if (levelNorm/residNorm0 > m_rtol)
    { 
      // finish wit jfnk
      jfnkSolve_private(ilev,*a_horizontalVel[ilev],*a_rhs[ilev],residNorm0,m_max_its,it);
    }

  if (m_verbosity>0)
    {
      m_op[ilev]->residual( *resid[ilev], 
			    *a_horizontalVel[ilev], 
			    *a_rhs[ilev], true );      
      levelNorm = m_op[ilev]->norm(*resid[ilev],2);
      pout() << "\t\tLevel 0: " << it << " total nonlinear iterations with " << m_minPicardIterations << " Picard iterations" << " Picard |r|_2 = " << levelNorm <<endl;
    }

  for (ilev=a_lbase+1;ilev<=a_maxLevel;ilev++)
    {
      // update coeficients to get correct residuals, sets c-f values
      computeMu( *a_horizontalVel[ilev],
		 *faceAs[ilev],
		 *faceMuCoef[ilev],
		 a_coordSys[ilev], 
		 ilev==a_lbase ? 0 : a_horizontalVel[ilev-1],
		 0,
		 ilev,
		 a_time );
      // c-f just done so no need here
      m_op[ilev]->AMRResidualNF( *resid[ilev], 
                                 *a_horizontalVel[ilev], 
                                 *a_horizontalVel[ilev-1], 
                                 *a_rhs[ilev], true ); 
      
      // cashes for nonlinear and updateCoefs
      //m_twork1 = tempv2[ilev];
      //m_twork2 = tempv3[ilev];
      m_tfaceA = faceAs[ilev];
      m_tmuCoef = faceMuCoef[ilev];
      m_tcoordSys = a_coordSys[ilev];
      m_ttime = a_time;
      //m_tphi0 = a_horizontalVel[ilev]; // cache phi_0 to get correct linearization
      m_tcrseVel = a_horizontalVel[ilev-1]; 

      // prolongate to ilev - done in AmrIce.cpp, use higher order? -- need before residaul calc & c-c interp!
      // m_op[ilev]->setToZero( *a_horizontalVel[ilev] );
      // AMRProlong( *a_horizontalVel[ilev], *a_horizontalVel[ilev-1], *m_fineCover[ilev-1], 
      // 	  m_projCopier[ilev-1], m_refRatios[ilev-1] );

      residNorm0 = m_op[ilev]->norm(*a_rhs[ilev], 2); // use |b| to compare against
      if (m_verbosity>0)
	{
	  pout() << "\t\tLevel " << ilev << " |b|_2 = " << residNorm0 << endl;
	}  

      // start with Picard solve
      for( it = 0; it < m_max_its ; it++ )
	{
#ifdef CH_USE_PETSC
	  Real opAlpha, opBeta;
	  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver = new PetscSolverViscousTensor<LevelData<FArrayBox> >();
	  getOperatorScaleFactors( opAlpha, opBeta );
	  solver->define( &(*m_op[ilev]), true ); // dx & crdx
	  solver->define( opAlpha, opBeta,
			       &(*m_C[ilev]),
			       &(*m_Mu[ilev]),
			       &(*m_Lambda[ilev]) );
	  solver->m_ctx = (void*)this;
	  solver->setInitialGuessNonzero(false);
	  
	  // linear KSP solve
	  m_op[ilev]->setToZero(*tempv[ilev]);
	  solver->solve(*tempv[ilev], *resid[ilev]);
	  m_op[ilev]->incr( *a_horizontalVel[ilev], *tempv[ilev], 1.); 
	  
	  updateCoefs( *a_horizontalVel[ilev], ilev ); 
#endif // CH_USE_PETSC
	  // put in residual correction form
          if (ilev == 0)
            {
              // no coarser level
              m_op[ilev]->residual( *resid[ilev], 
                                    *a_horizontalVel[ilev], 
                                    *a_rhs[ilev], false );
            }
          else
            {
              m_op[ilev]->AMRResidualNF( *resid[ilev], 
                                         *a_horizontalVel[ilev],
                                         *a_horizontalVel[ilev-1],
                                         *a_rhs[ilev], false );
            }

	  levelNorm = m_op[ilev]->norm(*resid[ilev],2);
	  if (m_verbosity>0)
	    {
	      pout() <<it+1<<") Level "<< ilev <<" Picard |r|_2 = " << levelNorm << endl;
	    }
	  if (levelNorm/residNorm0 < m_rtol){it++; break;}
	}
      
      // for( /* void */; it < m_max_its ; it++)
      // 	{
      // 	  Real opAlpha, opBeta;
      // 	  getOperatorScaleFactors( opAlpha, opBeta );
      // 	  PetscSolverViscousTensor<LevelData<FArrayBox> > *solver = new PetscSolverViscousTensor<LevelData<FArrayBox> >;
      // 	  solver->define( &(*m_op[ilev]), true ); // dx & crdx
      // 	  solver->define( opAlpha, opBeta,
      // 			       &(*m_C[ilev]),
      // 			       &(*m_Mu[ilev]),
      // 			       &(*m_Lambda[ilev]) );
      // 	  solver->setFunctionAndJacobian( FormFunction, FormJacobian ); // NL solve
      // 	  solver->m_ctx = (void*)this;
	  
      // 	  // creates Mat and Vecs, creates SNES
      // 	  solver->setup_solver(*a_horizontalVel[ilev]); // creates SNES
      // 	  // set the level to index into this
      // 	  SNESSetTolerances(solver->m_snes,(PetscReal)m_atol,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);
      // 	  SNESSetApplicationContext( solver->m_snes, (void*)&ilev );
      // 	  m_op[ilev]->setToZero( *tempv[ilev] );
      // 	  solver->solve(*tempv[ilev], *resid[ilev]);
      // 	  delete solver;
      // 	  m_op[ilev]->incr( *a_horizontalVel[ilev], *tempv[ilev], 1.); 
	  
      // 	  updateCoefs( *a_horizontalVel[ilev], ilev ); 	  
      // 	  // put in residual correction form
      // 	  m_op[ilev]->residual( *resid[ilev], 
      // 				*a_horizontalVel[ilev], 
      // 				*a_rhs[ilev], true ); 
      // 	  levelNorm = m_op[ilev]->norm(*resid[ilev],2);
      // 	  pout() <<it+1<<") Level "<< ilev <<" SNES |r|_2 = " << levelNorm << endl;
      // 	  if (levelNorm/residNorm0 < m_rtol){it++; break;}
      // 	}
      if (m_verbosity>0)
	{
	  pout() << "\t\tLevel "<< ilev <<" done with " << it << " nonlinear iterations " << endl;      
	}
    }

  // m_op[a_lbase]->residual( *resid[a_lbase], 
  // 			   *a_horizontalVel[a_lbase], 
  // 			   *a_rhs[a_lbase], true ); 
  // levelNorm = m_op[a_lbase]->norm(*resid[a_lbase],2);
  // pout() << " base |r|_2 = " << levelNorm << endl;
  // writeLevelname(resid[a_lbase],"res0.hdf5");
  // writeLevelname(resid[a_maxLevel],"res1.hdf5");
  // if(a_lbase==a_maxLevel)
  //   writeLevelname(a_horizontalVel[a_lbase],"phi0.hdf5");
  
  if(a_maxLevel==m_amr_max_level)
    {
      int numLevels = a_maxLevel + 1;
      Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);
      int normType = 0; Real dx = m_op[0]->dx();

      // clean up with full JFNK solve
      if(a_lbase!=a_maxLevel)
	{
	  if (m_verbosity>3)
	    {
	      Vector<LevelData<FArrayBox>*> res(m_grids.size());
	      {
		MultilevelIceVelOp mlOp;
		Vector<RealVect> dxs(m_grids.size());	    
		for (ilev=a_maxLevel;ilev>=a_lbase;ilev--)
		  {
		    dxs[ilev] = m_op[ilev]->dx()*RealVect::Unit;
		    res[ilev] = &(*resid[ilev]);
		  }
		// like a cast 
		RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > 
		  opFactoryPtr = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(m_opFactoryPtr);
		
		mlOp.define( m_grids, m_refRatios, m_domains, dxs, 
			     opFactoryPtr, a_lbase);
		
		mlOp.applyOp(res, a_horizontalVel, true);
		mlOp.incr(res, a_rhs, -1);
		mlOp.scale(res, -1.0);
		m_opFactoryPtr = opFactoryPtr; // take this back
	      }
	      Real gnorm =  computeNorm(res,m_refRatios,m_op[0]->dx(),Interval(0,1),normType);
	      pout() << "AMR |r|_"<< normType << " = " << gnorm << endl;      
	    }
	  {
	    RealVect cdx(D_DECL6(dx,dx,dx,dx,dx,dx));
	    JFNKSolver* jfnkSolver;
	    jfnkSolver = new JFNKSolver();
	    jfnkSolver->define( m_domains[0],
				m_constRelPtr,
				m_basalFrictionRelPtr,
				m_grids,
				m_refRatios,
				cdx,
				m_bc,
				numLevels);
	    if (m_verbosity>0)
	      {
		pout() << "call JFNK for clean up solve" << endl;      
	      }
	    returnCode = jfnkSolver->solve( a_horizontalVel, 
					    a_calvedIce, a_addedIce, a_removedIce, 
					    a_initialResidualNorm,  a_finalResidualNorm,
					    a_convergenceMetric, a_rhs, a_beta, a_beta0, a_A, a_muCoef,
					    a_coordSys, a_time, a_lbase, a_maxLevel);
	    delete jfnkSolver;
	    if (m_verbosity>0)
	      {
		pout() << "JFNK clean up solve done" << endl;      
	      }
	  }
	}
      
      // plot residual
      if(m_plotResidual)
	{
	  Vector<LevelData<FArrayBox>*> res(m_grids.size());
	  {
	    MultilevelIceVelOp mlOp;
	    Vector<RealVect> dxs(m_grids.size());	    
	    for (ilev=a_maxLevel;ilev>=a_lbase;ilev--)
	      {
		dxs[ilev] = m_op[ilev]->dx()*RealVect::Unit;
		res[ilev] = &(*resid[ilev]);
		computeMu( *a_horizontalVel[ilev],
			   *faceAs[ilev],
			   *faceMuCoef[ilev],
			   a_coordSys[ilev], 
			   ilev==a_lbase ? 0 : a_horizontalVel[ilev-1],
			   ilev==a_maxLevel ? 0 : a_horizontalVel[ilev+1],
			   ilev,
			   a_time );
	      }
	    // like a cast 
	    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > 
	      opFactoryPtr = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(m_opFactoryPtr);

	    mlOp.define( m_grids, m_refRatios, m_domains, dxs, 
			 opFactoryPtr, a_lbase);

	    mlOp.applyOp(res, a_horizontalVel, true);
	    mlOp.incr(res, a_rhs, -1);
	    mlOp.scale(res, -1.0);
	    m_opFactoryPtr = opFactoryPtr; // take this back

	    if (m_verbosity>3)
	      {
		Real gnorm =  computeNorm(res,m_refRatios,m_op[0]->dx(),Interval(0,1),normType);
		pout() << "plot AMR resid. |r|_"<< normType << " = " << gnorm << endl;      
	      }
	  }

	  for (ilev=a_maxLevel;ilev>=a_lbase;ilev--)
	    {
	      // zero covered
	      if( ilev!=a_maxLevel)
		{
		  m_op[ilev]->zeroCovered(*resid[ilev],*m_fineCover[ilev], m_restCopier[ilev]); 
		}
	    }

	  string fname = "residual.";
	  
	  char suffix[30];
	  sprintf(suffix, "%dd.hdf5",SpaceDim);
	  fname += suffix;
	  
	  Vector<string> varNames(2);
	  varNames[0] = "phi-x";
	  varNames[1] = "phi-y";
	  
	  Real bogusVal = 1.0;
#ifdef CH_USE_HDF5
	  WriteAMRHierarchyHDF5(fname,
				m_grids,
				res,
				varNames,
				m_domains[0].domainBox(),
				m_op[0]->dx(),
				bogusVal,
				bogusVal,
				m_refRatios,
				numLevels);
#endif
	  // clean up
	  for (ilev=0; ilev<plotData.size(); ilev++)
	    {
	      delete plotData[ilev];
	    }
	}
    }
  
  return returnCode;
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::defineOpFactory()
////////////////////////////////////////////////////////////////////////
void PetscIceSolver::defineOpFactory( RealVect a_Crsdx,
				      const ProblemDomain &a_domainCoar,
				      int a_numLevels )
{
  if ( !m_isOpDefined )
    {
      if (SpaceDim == 2)
	{
	  Real alpha, beta;
	  BCHolder velSolveBC(m_bc->velocitySolveBC());
  	  
	  for (int ilev=0;ilev<a_numLevels;ilev++)
	    {
	      // so needs to be set to avoid assert failures when setting diag too early
	      DataIterator dit = m_C[ilev]->dataIterator();
	      for (dit.begin(); dit.ok(); ++dit)
		{
		  (*m_C[ilev])[dit].setVal(1.0);
		  (*m_Lambda[ilev])[dit].setVal(1.0);
		  (*m_Mu[ilev])[dit].setVal(1.0);
		}
	    }

	  // OpFactory grabs a pointer to mu,lambda,acoef,etc.
	  getOperatorScaleFactors( alpha, beta );
	  m_opFactoryPtr = RefCountedPtr<ViscousTensorOpFactory> 
	    (new ViscousTensorOpFactory( m_grids, m_Mu, m_Lambda, m_C, alpha, 
					 beta, m_refRatios, a_domainCoar, a_Crsdx[0], 
					 velSolveBC, m_vtopSafety));
	}
      else 
	{
	  MayDay::Error("PetscIceSolver::defineOpFactory not implemented for dim = SpaceDim");
	}
    }
  else 
    {
      MayDay::Error("PetscIceSolver::defineOpFactory called twice???");
    }

  m_isOpDefined = true;
}
////////////////////////////////////////////////////////////////////////
// PetscIceSolver::getOperatorScaleFactors()
////////////////////////////////////////////////////////////////////////
void
PetscIceSolver::getOperatorScaleFactors(Real& a_alpha, Real& a_beta) const
{
  a_alpha = -1.0; // sort of wrong signs, but that's what B does 
  a_beta = 1.0; 
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::updateCoefs()
////////////////////////////////////////////////////////////////////////
void
PetscIceSolver::updateCoefs( LevelData<FArrayBox> &a_horizontalVel, int a_ilev )
{
  computeMu( a_horizontalVel, 
	     *m_tfaceA, 
	     *m_tmuCoef,
	     m_tcoordSys, 
	     m_tcrseVel, 
	     0,
	     a_ilev,
	     m_ttime );
}

////////////////////////////////////////////////////////////////////////
// PetscIceSolver::computeMu()
//   side effect: sets m_Mu & m_Lambda & m_C
//   a_horizontalVel has fine grid interpoled up to it
////////////////////////////////////////////////////////////////////////
// isothermal version -- for the ViscousTensorOp, lambda = 2*mu
void 
PetscIceSolver::computeMu( LevelData<FArrayBox> &a_horizontalVel,
			   const LevelData<FluxBox> &a_faceA, 
			   const LevelData<FluxBox> &a_muCoef,
			   const RefCountedPtr<LevelSigmaCS> &a_coordSys,
			   LevelData<FArrayBox>* crseVelPtr,
			   LevelData<FArrayBox>* fineVelPtr,
			   int a_ilev,
			   Real a_time)
{
  CH_TIME("PetscIceSolver::computeMu");
  ProblemDomain levelDomain = m_domains[a_ilev];
  const DisjointBoxLayout& levelGrids = m_grids[a_ilev];
  const LevelSigmaCS& levelCS = *a_coordSys;
  LevelData<FArrayBox>& levelVel = a_horizontalVel;
  LevelData<FluxBox>& levelMu = *m_Mu[a_ilev];

  const LevelData<FluxBox>& levelA = a_faceA;
  LevelData<FluxBox>& levelLambda = *m_Lambda[a_ilev];
  const LevelData<FArrayBox>& levelBeta = *m_Beta[a_ilev];
  const LevelData<FArrayBox>& levelBeta0 = *m_Beta0[a_ilev];
  LevelData<FArrayBox>& levelC = *m_C[a_ilev];
  DataIterator dit = levelGrids.dataIterator();

  // // first thing, if there is a finer level, average-down
  // // the current velocity field
  if ( fineVelPtr )
    {
      CoarseAverage averager(fineVelPtr->getBoxes(),
			     levelGrids,
			     fineVelPtr->nComp(),
			     m_refRatios[a_ilev]);
      
      averager.averageToCoarse(levelVel, *fineVelPtr);
    }

  // this is needed
  levelVel.exchange();

  // first set BC's on vel
  m_bc->velocityGhostBC(levelVel,
			levelCS,
			levelDomain, a_time);

  //slc : qcfi.coarseFineInterp fills the edges of lev > 0 cells
  //but not the corners. We need them filled to compute the
  //rate-of-strain invariant, so here is a bodge for now
  // if (SpaceDim == 2)
  //   {
  //     for (dit.begin(); dit.ok(); ++dit)
  // 	{
  // 	  Box sbox = levelVel[dit].box();
  // 	  sbox.grow(-1);
  // 	  FORT_EXTRAPCORNER2D(CHF_FRA(levelVel[dit]),
  // 			      CHF_BOX(sbox));
  // 	}

  //   }
  
  // // actually need to use a cornerCopier, too...
  // CornerCopier cornerCopier(levelGrids, levelGrids, 
  // 			    levelDomain,levelVel.ghostVect(),
  // 			    true);
  // levelVel.exchange(cornerCopier);

  int refToCrs = crseVelPtr ? m_refRatios[a_ilev-1] : -1;
  IntVect muGhost = IntVect::Zero;
  m_constRelPtr->computeFaceMu( levelMu,
				levelVel, 1.0,
				crseVelPtr,
				refToCrs,
				levelA,
				levelCS,
				levelDomain,
				muGhost);

  // now multiply by ice thickness H
  const LevelData<FluxBox>& faceH = levelCS.getFaceH();
  // Real muMax = 1.23456789e+300;
  // Real muMin = 0.0;
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  FArrayBox& thisMu = levelMu[dit][dir];
	  const Box& box = thisMu.box();
	  	  
	  // FORT_MAXFAB1(CHF_FRA(thisMu),
	  // 	       CHF_CONST_REAL(muMin),
	  // 	       CHF_BOX(box));
	  
	  thisMu.mult(faceH[dit][dir],box,0,0,1);
	  thisMu.mult((a_muCoef)[dit][dir],box,0,0,1);
	    
	  // FORT_MINFAB1(CHF_FRA(thisMu),
	  // 	       CHF_CONST_REAL(muMax),
	  // 	       CHF_BOX(box));
	}    
      
      // also update alpha (or C)
      const Box& gridBox = levelGrids[dit];
      m_basalFrictionRelPtr->computeAlpha
	(levelC[dit], levelVel[dit],levelBeta[dit], 1.0, 
	 levelCS, dit, a_ilev ,gridBox);
      
      levelC[dit] += levelBeta0[dit];

// #if CH_SPACEDIM==2
//       {
// 	Real mu0 = 1.0;
// 	Real C0 = 1.0;
	
// 	FORT_ENFORCEWELLPOSEDCELL
// 	  (CHF_FRA1(levelC[dit],0),
// 	   CHF_FRA1(levelMu[dit][0],0),
// 	   CHF_FRA1(levelMu[dit][1],0),
// 	   CHF_CONST_REAL(mu0),
// 	   CHF_CONST_REAL(C0),
// 	   CHF_BOX(levelGrids[dit]));
	
//       }
// #endif

      // lambda = 2*mu
      FluxBox& lambda = levelLambda[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  lambda[dir].copy(levelMu[dit][dir]);
	  lambda[dir] *= 2.0;
	}

    }
}
#include "NamespaceFooter.H"
