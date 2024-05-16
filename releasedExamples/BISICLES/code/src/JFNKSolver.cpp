#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "JFNKSolver.H"
#include "QuadCFInterp.H"
#include "CornerCopier.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"
#include "ParmParse.H"
#include "BisiclesF_F.H"
#include "LinearSolver.H"
#include "RelaxSolver.H"
#include "BiCGStabSolver.H"
#include "GMRESSolver.H"
#include "CGSolver.H"
#include "IceUtility.H"
#include <sstream>
#include "NamespaceHeader.H"

// static member initialization
int LinearizedVTOp::m_bottom_solver_type = bicg;
int LinearizedVTOp::m_MG_solver_depth = -1;

void JFNKSolver::imposeMaxAbs(Vector<LevelData<FArrayBox>*>& a_u,
		 Real a_limit)
{

  for (int lev = 0; lev < a_u.size(); ++lev)
    {
      LevelData<FArrayBox>& levelU = *a_u[lev];
      for (DataIterator dit (levelU.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  FArrayBox& thisU = levelU[dit];
	  FORT_ABSLIMITFAB(CHF_FRA(thisU), 
			   CHF_CONST_REAL(a_limit), 
			   CHF_BOX(thisU.box()));
	}
    }
  

}

int LinearizedVTOp::m_residualID(0);

LinearizedVTOp::LinearizedVTOp(NonlinearViscousTensor* a_currentState , 
	       Vector<LevelData<FArrayBox>*>& a_u,
	       Real a_h,  Real a_err, Real a_umin, bool a_hAdaptive, 
	       Vector<DisjointBoxLayout>& a_grids,
	       Vector<int>& a_refRatio,
	       Vector<ProblemDomain>& a_domains,
	       Vector<RealVect>& a_dxs,
	       int a_lBase, 
	       int a_numMGSmooth,
	       int a_numMGIter,
	       LinearizationMode a_mode) 
  : m_u(a_currentState),
    m_h(a_h), m_err(a_err), m_umin(a_umin), m_hAdaptive(a_hAdaptive),
    m_grids(a_grids), m_refRatio(a_refRatio),
    m_domains(a_domains), m_dxs(a_dxs), m_lBase(a_lBase),
    m_mode(a_mode),m_writeResiduals(false)
    
{
  
  Vector<LevelData<FArrayBox>*> localU( m_grids.size());
  for (int lev = 0; lev < m_grids.size(); ++lev)
    {
      localU[lev] = a_u[lev];
    }

  RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
    opFactoryPtr = m_u->opFactoryPtr();
  
  m_mlOp.m_bottom_solver_type = m_bottom_solver_type;
  m_mlOp.m_preCondSolverDepth = m_MG_solver_depth;
  m_mlOp.m_num_mg_smooth =  a_numMGSmooth;
  m_mlOp.m_num_mg_iterations = a_numMGIter;
  m_mlOp.define(m_grids , m_refRatio, m_domains, m_dxs, 
		opFactoryPtr ,a_lBase);
  
  create(m_fu,localU);
  setU(localU);

  if (a_mode == JFNK_LINEARIZATION_MODE)
    {
      
      m_uPerturbed = m_u->newNonlinearViscousTensor();
      RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
	pOpFactoryPtr = m_uPerturbed->opFactoryPtr();
      m_perturbedMlOp.define(m_grids , m_refRatio, m_domains, m_dxs, 
			 pOpFactoryPtr , a_lBase);
  
      create(m_uplushv,localU);
    }
  else
    {
      m_uPerturbed = NULL;
    }
  
    
  }

void  LinearizedVTOp::outerResidual(Vector<LevelData<FArrayBox>*>& a_lhs, 
			     const Vector<LevelData<FArrayBox>*>& a_u, 
			     const Vector<LevelData<FArrayBox>*>& a_rhs, 
			     bool a_homogeneous )
  {
    
    m_mlOp.applyOp(a_lhs, a_u, a_homogeneous);
    incr(a_lhs, a_rhs, -1);
    scale(a_lhs, -1.0);
    if (m_writeResiduals)
    {
      writeResidual(a_u,a_lhs);
    }
  }
  

void  LinearizedVTOp::residual(Vector<LevelData<FArrayBox>*>& a_lhs, 
		       const Vector<LevelData<FArrayBox>*>& a_v, 
		       const Vector<LevelData<FArrayBox>*>& a_rhs, 
		       bool a_homogeneous)
{    
  applyOp(a_lhs, a_v, a_homogeneous);
  incr(a_lhs, a_rhs, -1.0);
  scale(a_lhs, -1.0);
  if (m_writeResiduals)
    {
      if (m_mode == JFNK_LINEARIZATION_MODE)
	writeResidual(m_uplushv,a_lhs);
      else
	writeResidual(a_v,a_lhs);
    }

}


void LinearizedVTOp::setU(Vector<LevelData<FArrayBox>*>& a_u)
{
  CH_TIME("LinearizedVTOp::setU");
  m_u->setState(a_u);
  m_mlOp.applyOp(m_fu, a_u);
  
}
Real LinearizedVTOp::finiteh(const Vector<LevelData<FArrayBox>*>& a_v)
{
  Real h = m_h;

  if (m_hAdaptive)
    {
      
      //calculation of h from PETSc SNES (see license notice at the end of this file)
      // h = err*u'v/||v||^2 if  |u'v| > umin*||v||_{1}
      //  = err*umin*sign(u'v)*||v||_{1}/||v||^2   otherwise
      Real utv = dotProduct(m_u->getState(), a_v); 
      Real scale = D_TERM(m_dxs[0][0], *m_dxs[0][1], *m_dxs[0][2]);
      utv /= scale;
      Real vL1norm = norm(a_v, 1);
      Real vL2norm = norm(a_v, 2);
      if (vL1norm > 0.0)
	{	  
	  if (utv >= 0.0 && utv <  m_umin*vL1norm) utv = m_umin*vL1norm;
	  else if (utv < 0.0 && utv > - m_umin*vL1norm) utv = -m_umin*vL1norm;
	  h = m_err * utv / (vL2norm *  vL2norm);
	}
      else
	{
	  h =  m_err;
	}
    }
  return h;

}


void LinearizedVTOp::applyOp(Vector<LevelData<FArrayBox>*>& a_lhs, 
		       const Vector<LevelData<FArrayBox>*>& a_v, 
		       bool a_homogeneous )
{

  CH_TIME("LinearizedVTOp::applyOp");
 
  

  if (m_mode == JFNK_LINEARIZATION_MODE)
    {
      //"Newton mode"
      setToZero(m_uplushv);
      Real h = finiteh(a_v);
      axby(m_uplushv, m_u->getState(), a_v, 1.0, h);
      CH_assert(norm(m_uplushv,0) < HUGE_NORM);
      m_uPerturbed->setState(m_uplushv);
      m_perturbedMlOp.applyOp(a_lhs, m_uplushv , a_homogeneous);

      incr(a_lhs, m_fu, -1.0);
      scale(a_lhs, 1.0 / h);

    }
  else if (m_mode == PICARD_LINEARIZATION_MODE)
    {
      //"Picard mode"
      m_mlOp.applyOp(a_lhs, a_v, a_homogeneous);
    }
  else 
    {
      MayDay::Error("Unknown LinearizationMode in LinearizedVTOp::applyOp");
    }

}


void LinearizedVTOp::writeResidual  
(const Vector<LevelData<FArrayBox> *>& a_u,
 const Vector<LevelData<FArrayBox> *>& a_residual)
{

  Vector<std::string> names;
  names.push_back("xU");
  names.push_back("yU");
  names.push_back("xRes");
  names.push_back("yRes");
  names.push_back("C");
  names.push_back("muSum");
  

  Vector<LevelData<FArrayBox>* > data(a_u.size(),NULL);
  for (int lev = 0; lev < a_u.size(); lev++)
    {
      data[lev] = new LevelData<FArrayBox>(m_grids[lev],int(names.size()),IntVect::Zero);
      int j = 0;
      a_u[lev]->copyTo(Interval(0,SpaceDim-1),*data[lev],Interval(j,j+SpaceDim-1)); j+=SpaceDim;
      a_residual[lev]->copyTo(Interval(0,SpaceDim-1),*data[lev],Interval(j,j+SpaceDim-1)); j+=SpaceDim;
      m_u->alpha()[lev]->copyTo(Interval(0,0),*data[lev],Interval(j,j)); j+=1;

      const LevelData<FluxBox>& mu =  *m_u->mu()[lev];
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  FArrayBox muSum; 
	  muSum.define(  Interval(j,j),(*data[lev])[dit]);
	  for (BoxIterator bit( m_grids[lev][dit]);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
              D_TERM(
              const IntVect ive = iv + BASISV(0);,
              const IntVect ivn = iv + BASISV(1);, )
              muSum(iv,0) = D_TERM(mu[dit][0](iv) + mu[dit][0](ive), 
                                   + mu[dit][1](iv) +  mu[dit][1](ivn), );
	    }
	}
    }

  // for (int lev = a_u.size() -1; lev > 0; lev--)
  //   {
  //     CoarseAverage avg(m_grids[lev],m_grids[lev-1],data[lev]->nComp(),m_refRatio[lev-1]);
  //     avg.averageToCoarse(*data[lev-1], *data[lev]);
  //   }
  
  //decide on the file name
  char file[32];
  sprintf(file,"jfnkopres.%06d.2d.hdf5",m_residualID);
  Real dt(1.0); Real time(m_residualID);
  pout() << "writing " << file << std::endl;
#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5(file ,m_grids, data ,names, m_domains[0].domainBox(),
			m_dxs[0][0], dt, time, m_refRatio, data.size());
#endif

  for (int lev = 0; lev < a_u.size(); lev++)
    {
      delete data[lev];
    }
  m_residualID++;
}

void LinearizedVTOp::preCond(Vector<LevelData<FArrayBox>*>& a_cor,
		     const Vector<LevelData<FArrayBox>*>& a_residual)
{
  //CH_assert(norm(a_residual,0) < HUGE_NORM);
  CH_TIME("LinearizedVTOp::precond");
  m_mlOp.preCond( a_cor, a_residual);
  //CH_assert(norm(a_cor,0) < HUGE_NORM);
}

JFNKSolver::Configuration::~Configuration()
{

}

JFNKSolver::Configuration::Configuration()
{
  
  //sensible defaults
  m_residualOnly = false;
  m_linearSolverType=Relax; 
  m_maxIter = 15;
  m_absTol = 1.0e-10;
  m_relTol = 1.0e-10;
  m_BiCGStabRelTol = 1.0e-3;
  m_maxBiCGStabIter = 10;
  m_CGRelTol = 1.0e-3;
  m_maxCGIter = 10;
  m_RelaxRelTol = 1.0e-3;
  m_maxRelaxIter = 10;
  m_RelaxHang = 0.25;
  m_GMRESRelTol = 1.0e-3;
  m_maxGMRESIter = 10;
  m_normType = 0;
  m_verbosity = 5;
  m_vtopSafety = 0.5;
  m_vtopRelaxMinIter = 8;
  m_vtopRelaxTol = 1.0e-2;
  m_numMGSmooth = 8;
  m_numMGIter = 1;
  m_h = 0.025;
  m_hAdaptive=false;
  m_err = 1.0e-6;
  m_umin = 1.0e-6;
  m_switchRate = 1.5;
  m_minPicardIter = 0;
  m_muMax = 1.23456789e+300;
  m_muMin = 0.0;
  m_writeResiduals = false;
  m_minStepFactor = 1.0;
  m_eliminateFastIce = false;
  m_eliminateFastIceSpeed = 1.0e+5;
  m_eliminateFastIceEdgeOnly = false;
  m_eliminateRemoteIceTol = 1.0;
  m_eliminateRemoteIceMaxIter = 10;

  //// artificial drag that applies everywhere
    /* if used, need to ensure that 
       (m_artificial_drag_coef * |u|)^m_artifical_drag_power << rho * g * h * grad(s) ~ 10^3-4
       when |u| ~ 10^4. m_artificial_drag_coef =1.0e-4 is about right.
    */   
  m_artificialDragCoef = 0.0;
  m_artificialDragPower = 8.0;
  // these ones don't need to be stored (at least for now), but should be set
  int mgAverageType  = CoarseAverageFace::arithmetic;
  ViscousTensorOpFactory::s_coefficientAverageType = mgAverageType;

  // set default to be linear prolongation in multigrid
  int mgProlongType = ViscousTensorOp::linearInterp;
  ViscousTensorOp::s_prolongType = mgProlongType;

  ViscousTensorOp::s_lazy_gsrb = false; // true default 
  
}


void JFNKSolver::Configuration::parse(const char* a_prefix)
{
  // set parameters based on parmParse
  ParmParse pp(a_prefix);
  pp.query("maxIter",m_maxIter);
  pp.query("absTol", m_absTol);
  pp.query("relTol",m_relTol);
  pp.query("normType",m_normType);
  pp.query("verbosity",m_verbosity);
  pp.query("vtopSafety",m_vtopSafety);
  pp.query("vtopRelaxTol",m_vtopRelaxTol);
  pp.query("vtopRelaxMinIter",m_vtopRelaxMinIter);
  pp.query("numMGSmooth",m_numMGSmooth);
  pp.query("numMGIter",m_numMGIter);
  pp.query("h",m_h);
  pp.query("hAdaptive",m_hAdaptive);
  pp.query("err",m_err);
  pp.query("umin",m_umin);
  pp.query("switchRate",m_switchRate);
  pp.query("minPicardIterations", m_minPicardIter);
  pp.query("muMax", m_muMax);
  pp.query("muMin", m_muMin);
  pp.query("writeResiduals", m_writeResiduals);
  pp.query("minStepFactor", m_minStepFactor);
  pp.query("eliminateFastIce",m_eliminateFastIce);
  pp.query("eliminateFastIceSpeed",m_eliminateFastIceSpeed);
  pp.query("eliminateFastIceEdgeOnly",m_eliminateFastIceEdgeOnly);

  {
    ParmParse ppAmr("amr"); // ugly, but we generally want to inherit these if the are not specified
    ppAmr.query("eliminate_remote_ice_tol",m_eliminateRemoteIceTol);
    ppAmr.query("eliminate_remote_ice_max_iter",m_eliminateRemoteIceMaxIter);


    pp.query("eliminateRemoteIceTol",m_eliminateRemoteIceTol);
    pp.query("eliminateRemoteIceMaxIter",m_eliminateRemoteIceMaxIter);
  }

  {
    /// since we tuned everthign to m/a, scaling to m/a by default 
    Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
    ParmParse ppc("constants");
    ppc.query("seconds_per_unit_time",seconds_per_unit_time);
    m_scale = SECONDS_PER_TROPICAL_YEAR/ seconds_per_unit_time;
  }

  
  
  pp.query("scale", m_scale);

  pp.query("artificialDragCoef",m_artificialDragCoef);
  pp.query("artificialDragPower",m_artificialDragPower);
  
  if (pp.contains("solverType") )
    {
      int solverIntType = m_linearSolverType;
      pp.query("solverType", solverIntType);
      if (solverIntType == BiCGStab)
        {
          m_linearSolverType = BiCGStab;
	  pp.query("BiCGStabRelTol",m_BiCGStabRelTol);
	  pp.query("maxBiCGStabIter",m_maxBiCGStabIter);
        }
      else if (solverIntType == Relax)
        {
          m_linearSolverType = Relax;
	  pp.query("RelaxTol", m_RelaxRelTol);
	  pp.query("RelaxRelTol", m_RelaxRelTol);
	  pp.query("maxRelaxIter", m_maxRelaxIter);
	  pp.query("RelaxHang", m_RelaxHang);
	  m_numMGIter = 1; // m_numMGIter > 1 doesn't help
        }
      else if (solverIntType == GMRES)
        {
          m_linearSolverType = GMRES;
	  MayDay::Error("JFNKSolver -- GMRES is on the blink");
        }
      else if (solverIntType == CG)
        {
	  m_linearSolverType = CG;
	  pp.query("CGRelTol",m_CGRelTol);
	  pp.query("maxCGIter",m_maxCGIter);
	}
      else if (solverIntType == petsc)
        {
#ifdef CH_USE_PETSC
          m_linearSolverType = petsc;
#else
          // default back to relax if petsc isn't compiled in
          // (this is just to simplify comparisons btwn petsc and MG)
          m_linearSolverType = Relax;
	  pp.query("RelaxTol", m_RelaxRelTol);
	  pp.query("RelaxRelTol", m_RelaxRelTol);
	  pp.query("maxRelaxIter", m_maxRelaxIter);
	  pp.query("RelaxHang", m_RelaxHang);
	  m_numMGIter = 1; // m_numMGIter > 1 doesn't help
#endif          
        }      
      else 
        {
          MayDay::Error("JFNKSolver -- bad linear solver type");
        }
    }
  
  int bs_type = LinearizedVTOp::m_bottom_solver_type;
  pp.query("bottom_solver_type", bs_type);
  // if petsc not compiled in but petsc specified, then fall back to default
#ifndef CH_USE_PETSC  
  if (bs_type == PETSC)
    {
      bs_type = LinearizedVTOp::m_bottom_solver_type;
    }
#endif


  LinearizedVTOp::m_bottom_solver_type = bs_type;

  int mg_depth = LinearizedVTOp::m_MG_solver_depth;
  pp.query("mg_solver_depth", mg_depth);
  LinearizedVTOp::m_MG_solver_depth = mg_depth;

  int mgAverageType  = CoarseAverageFace::arithmetic;
  pp.query("mgCoefficientAverageType", mgAverageType);
  ViscousTensorOpFactory::s_coefficientAverageType = mgAverageType;

  // set default to be linear prolongation in multigrid
  int mgProlongType = ViscousTensorOp::linearInterp;
  pp.query("mgProlongType", mgProlongType);

  ViscousTensorOp::s_prolongType = mgProlongType;

#ifdef CHOMBO_TRUNK
  /// default is "lazy gsrb"
  if (pp.contains("lazyGSRB"))
    {
      bool lazyGSRB;
      pp.get("lazyGSRB", lazyGSRB);
      ViscousTensorOp::s_lazy_gsrb = lazyGSRB;
    }
#endif
  
}

void JFNKSolver::define(const ProblemDomain& a_coarseDomain,
			ConstitutiveRelation* a_constRel,
			BasalFrictionRelation* a_basalFrictionRel,
			const Vector<DisjointBoxLayout>& a_grids,
			const Vector<int>& a_refRatios,
			const RealVect& a_dxCrse,
			IceThicknessIBC* a_bc,
			int a_numLevels)
{
  m_config.parse("JFNKSolver");

  m_constRelPtr = a_constRel;
  m_basalFrictionRelPtr = a_basalFrictionRel;
  m_bcPtr = a_bc;
  
  m_grids.resize(a_numLevels); 
  m_grids[0] = a_grids[0];

  m_refRatios.resize(a_numLevels);
  m_refRatios = a_refRatios;

  m_dxs.resize(a_numLevels);
  m_dxs[0] = a_dxCrse;

  m_domains.resize(a_numLevels);
  m_domains[0] = a_coarseDomain;
  
  
  for (int lev = 1; lev < a_numLevels; ++lev)
    {
      m_dxs[lev] = m_dxs[lev-1] / m_refRatios[lev-1];
      m_domains[lev] = m_domains[lev-1];
      m_domains[lev].refine(m_refRatios[lev-1]);
      m_grids[lev] = a_grids[lev];
      
    }

}

//IceVelocitySolver full solve
inline
int JFNKSolver::solve(Vector<LevelData<FArrayBox>* >& a_u,
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
		      Real a_time,  int a_lbase, int a_maxLevel)
{

  //Vector<LevelData<FluxBox>* > muCoef(a_maxLevel,NULL);
 
  int rc =  solve(a_u, a_calvedIce, a_addedIce, a_removedIce,
		  a_initialResidualNorm, a_finalResidualNorm,
		  a_convergenceMetric, false, a_rhs, a_C, a_C0, a_A, 
		  a_muCoef, a_coordSys,
		  a_time, a_lbase,  a_maxLevel);


  return rc;

}			    


#ifdef CH_USE_PETSC
#undef __FUNCT__
#define __FUNCT__ "solve"
#endif

void JFNKSolver::eliminateFastIce(Vector<LevelData<FArrayBox>* >& a_velocity,
				  Vector<LevelData<FArrayBox>* >& a_calvedIce,
				  Vector<LevelData<FArrayBox>* >& a_addedIce,
				  Vector<LevelData<FArrayBox>* >& a_removedIce,
				  Vector<LevelData<FArrayBox>* >& a_rhs,
				  Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
				  IceNonlinearViscousTensor&  a_current)
				  
				  
{
  CH_TIME("JFNKSolver::eliminateFastIce");
  // optionally get rid of ice with excessive |u|
  if ( m_config.m_eliminateFastIce)
    {
      int eliminated = IceUtility::eliminateFastIce
	(a_coordSys, a_velocity, a_calvedIce, a_addedIce, a_removedIce,
	 m_grids , m_domains, 
	 m_refRatios, m_dxs[0][0], a_velocity.size() -1, 
	 m_config.m_eliminateRemoteIceMaxIter,  m_config.m_eliminateRemoteIceTol, 
	 m_config.m_eliminateFastIceSpeed,  m_config.m_eliminateFastIceEdgeOnly, m_config.m_verbosity);
      
      //need to redfine RHS
      if (eliminated > 0)
	{
	  IceUtility::defineRHS(a_rhs, a_coordSys, m_grids, m_dxs);
	  //a_current.setState(a_velocity);
	  //NB the viscosity and residual will also need computing, we are relying on a subsequent call
	}
    }
}


int JFNKSolver::solve(Vector<LevelData<FArrayBox>* >& a_u,
		      Vector<LevelData<FArrayBox>* >& a_calvedIce,
		      Vector<LevelData<FArrayBox>* >& a_addedIce,
		      Vector<LevelData<FArrayBox>* >& a_removedIce,
		      Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		      const Real a_convergenceMetric,
		      const bool a_linear,
		      const Vector<LevelData<FArrayBox>* >& a_rhs,
		      const Vector<LevelData<FArrayBox>* >& a_C,
		      const Vector<LevelData<FArrayBox>* >& a_C0,
		      const Vector<LevelData<FArrayBox>* >& a_A,
		      const Vector<LevelData<FArrayBox>* >& a_muCoef,
		      Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		      Real a_time,  int a_lbase, int a_maxLevel)
{
  CH_TIME("JFNKSolver::solve");
  CH_assert(a_lbase == 0); //\todo : support lBase != 0
  int returnCode = 0;
  
  Vector<LevelData<FArrayBox>* > localU(a_maxLevel+1);
  Vector<LevelData<FArrayBox>* > localRhs(a_maxLevel+1);
  Vector<LevelData<FArrayBox>* > localC(a_maxLevel+1);
  Vector<LevelData<FArrayBox>* > localC0(a_maxLevel+1);
 
 for (int lev = 0; lev < a_maxLevel + 1; ++lev)
   {
     localU[lev] = a_u[lev];
     localRhs[lev] = a_rhs[lev];
     localC[lev] = a_C[lev]; 
     localC0[lev] = a_C0[lev];

     /// scale u (e.g if u is in m/s, scaling to m/a)
     for (DataIterator dit(a_u[lev]->disjointBoxLayout()); dit.ok(); ++ dit)
       {
	 (*a_u[lev])[dit] *= m_config.m_scale;
       }
   }
 
 
 //cell face A
 Vector<LevelData<FluxBox>* > faceA(a_maxLevel+1, NULL);
 Vector<LevelData<FluxBox>* > faceMuCoef(a_maxLevel+1, NULL);
 for (int lev = 0; lev < a_maxLevel + 1; ++lev)
   {
     faceA[lev] = new LevelData<FluxBox>
       (m_grids[lev], a_A[lev]->nComp(), IntVect::Zero);
     CellToEdge(*a_A[lev] , *faceA[lev]);
     faceMuCoef[lev] = new LevelData<FluxBox>
       (m_grids[lev],  a_muCoef[lev]->nComp(),  IntVect::Zero);
      CellToEdge(*a_muCoef[lev] , *faceMuCoef[lev]);

   }
 
 //An IceNonlinearViscousTensor object is used to compute the viscosity and drag coefficents
 //(which depend on the velocity u) for the current outer iteration
 IceNonlinearViscousTensor current
   (m_grids , m_refRatios, m_domains,m_dxs , a_coordSys, localU ,
    localC, localC0, a_maxLevel, *m_constRelPtr, *m_basalFrictionRelPtr, *m_bcPtr, 
    a_A, faceA, a_time, m_config.m_vtopSafety, m_config.m_vtopRelaxMinIter, m_config.m_vtopRelaxTol, 
    m_config.m_muMin,  m_config.m_muMax, m_config.m_scale, m_config.m_artificialDragCoef, m_config.m_artificialDragPower);

 // eliminate fast ice if required:
 eliminateFastIce(localU, a_calvedIce, a_addedIce, a_removedIce, localRhs, a_coordSys, current);
 current.setState(localU);
 current.setFaceViscCoef(faceMuCoef);

 // LinearizedVTOp J0 provides several methods useful for initializtion, 
 LinearizedVTOp J0
   (&current, localU, m_config.m_h, m_config.m_err, m_config.m_umin, m_config.m_hAdaptive,  m_grids, m_refRatios, 
    m_domains, m_dxs, a_lbase, m_config.m_numMGSmooth, m_config.m_numMGIter, PICARD_LINEARIZATION_MODE);
 
 
 
 //storage for residual and correction data
 Vector<LevelData<FArrayBox>*> residual;
 Vector<LevelData<FArrayBox>*> du;
 J0.create(residual,localU);  
 J0.create(du,localU);
 



 //calculate the initial residual and its norm
 J0.m_writeResiduals =  m_config.m_writeResiduals;

 if (a_linear)
   {
     J0.setToZero(localU);
   }
 J0.outerResidual(residual,localU,localRhs);
 Real oldResNorm;
 Real& resNorm = a_finalResidualNorm;
 resNorm = a_initialResidualNorm = J0.norm(residual, m_config.m_normType);
 

 if (m_config.m_verbosity > 0)
   {
     pout() << "JFNK initial residual norm = " << resNorm << std::endl;
   }
 Real convergenceMetric = std::max(a_convergenceMetric,a_initialResidualNorm);
 if (m_config.m_verbosity > 0)
    {
      pout() << "JFNK convergence metric = " <<  convergenceMetric << std::endl;
    }
  if (m_config.m_residualOnly)
    return 0;

  bool done =  (resNorm  < m_config.m_absTol) |  (resNorm < a_convergenceMetric * m_config.m_relTol);
  
  if (a_linear)
    {
     
      linearSolve(J0, localU, localRhs, PICARD_LINEARIZATION_MODE);
      J0.residual(residual,localU,localRhs);
      resNorm = J0.norm(residual, m_config.m_normType);
      done = true;
    }
  else
    {     
      LinearizationMode mode =  (m_config.m_minPicardIter < 0)?JFNK_LINEARIZATION_MODE:PICARD_LINEARIZATION_MODE;
      int iter = 0;
      
      while (!done && iter < m_config.m_maxIter)
	{
	  oldResNorm = resNorm;
	  
	  
	  //create a linearization (either the Jacobian of f or an approximation to it) around the current a_u
	  LinearizedVTOp J
	    (&current, localU, m_config.m_h, m_config.m_err, m_config.m_umin, m_config.m_hAdaptive, m_grids, 
	     m_refRatios, m_domains, m_dxs, a_lbase, m_config.m_numMGSmooth, m_config.m_numMGIter, mode);
	  J.m_writeResiduals =  m_config.m_writeResiduals;
	  // we don't *always* need to re-evalue the residual, but we do if we e.g removed any ice. We could test for those cases,
	  // but this is not much more expensive  
	  J.outerResidual(residual, localU, localRhs);
	  
	  //solve the linear system J du = r (with r = -f(u))
	  J.setToZero(du);
	  linearSolve(J, du, residual, mode);
	  
	  //update u <- u + w*du for some w with minW < w <= 1.0.
	  //When du is a JFNK step, allow w < 1 and w = 0 if ||f(u + minW * du)|| > ||f(u)||
	  Real minW = (mode == JFNK_LINEARIZATION_MODE)? m_config.m_minStepFactor:1.0;
	  bool resetOnFail = (mode == JFNK_LINEARIZATION_MODE);

	  resNorm = lineSearch(localU, residual, localRhs,  du, J, 
			       current, minW, resetOnFail);
	  
	  //test the state and decide whether to switch modes
	  if (mode == PICARD_LINEARIZATION_MODE)
	    {
	      Real rate = oldResNorm / resNorm;
	      if (m_config.m_verbosity > 0)
		{
		  pout() << " Picard iteration " << iter 
			 << " residual norm = " << resNorm 
			 << " rate = " <<  rate
			 << std::endl;
		}
	      
	      if (iter >= m_config.m_minPicardIter && rate > 1.0 && rate < m_config.m_switchRate){
		//since the iterations are progressing slowly
		//but positively, try JFNK mode
		mode = JFNK_LINEARIZATION_MODE;
	      }
	    }
	  else if (mode == JFNK_LINEARIZATION_MODE)
	    {
	      Real rate_switch_back = 1.0 + 0.25*(m_config.m_switchRate - 1.0);
	      if (resNorm >= oldResNorm/rate_switch_back)
		{
		  if (m_config.m_verbosity > 0){
		    pout() << "JFNK iteration " << iter  <<  " did not reduce residual (much)" << std::endl;
		    if (m_config.m_verbosity > 1)
		      {
			pout() << " -- old residual = " << oldResNorm  
			       << " -- new residual = " << resNorm << std::endl;
		      }
		  }
		  
		  //The JFNK step has failed, so switch back to Picard mode 
		  mode = PICARD_LINEARIZATION_MODE;
		  resNorm = oldResNorm;
		}
	      else
		{
		  if (m_config.m_verbosity > 0)
		    {
		      pout() << " JFNK iteration " << iter 
			     << " residual norm = " << resNorm 
			 << " rate = " <<  oldResNorm / resNorm 
			     << std::endl;
		      
		    }
		}
	    }
	  
	  done = resNorm < m_config.m_absTol || resNorm < convergenceMetric * m_config.m_relTol;
	  
	  //eliminate fast ice if needed
	  eliminateFastIce(localU, a_calvedIce, a_addedIce, a_removedIce, localRhs, a_coordSys, current);
	  
	  ++iter;
	}
       pout() << "JFNK final residual norm = " << resNorm << std::endl;

       // check to see if we actually solved the problem
       if (done)
	 {
	   returnCode = 0;
	   if (m_config.m_verbosity > 0)
	     {
	       pout() << "JFNKSolver converged -- final norm(resid) = "
		      << resNorm << " after " << iter << " iterations"
		      << endl;
	     }
	 }
       else 
	 {
	   returnCode = 1;
	   if (m_config.m_verbosity > 0)
	     {
	       pout() << "JFNKSolver NOT CONVERGED -- final norm(resid) = "
		      << resNorm << " after " << iter << " iterations"
		      << endl;          
	     }
	 }
       
    }// end !a_linear


  //unscale u
  for (int lev = 0; lev <= a_maxLevel; lev++)
    {
      for (DataIterator dit(a_u[lev]->disjointBoxLayout()); dit.ok(); ++ dit)
	{
	  (*a_u[lev])[dit] /= m_config.m_scale;
	}
    }
  
  // clean up storage
  for (int lev=0; lev<residual.size(); lev++)
    {
      if (residual[lev] != NULL)
	{
	  delete residual[lev];
	  residual[lev] = NULL;
	}
      
      if (du[lev] != NULL)
	{
	  delete du[lev];
	  du[lev] = NULL;
	}

      if (faceA[lev] != NULL)
	{
	  delete faceA[lev];
	  faceA[lev] = NULL;
	}

      if (faceMuCoef[lev] != NULL)
	{
	  delete faceMuCoef[lev];
	  faceMuCoef[lev] = NULL;
	}

    }
  
  return returnCode;

} 



LinearSolver<Vector<LevelData<FArrayBox>* > >* 
JFNKSolver::newLinearSolver (LinearizedVTOp& a_op,  LinearizationMode a_mode)
{

  LinearSolver<Vector<LevelData<FArrayBox>* > >* solver = NULL;
  if (m_config.m_linearSolverType == Configuration::BiCGStab)
    {
      BiCGStabSolver<Vector<LevelData<FArrayBox>* > >* biCGStabSolver 
	= new BiCGStabSolver<Vector<LevelData<FArrayBox>* > >;
      
      biCGStabSolver->define(&a_op , false);
      biCGStabSolver->m_verbosity = m_config.m_verbosity - 1;
      biCGStabSolver->m_reps =  m_config.m_BiCGStabRelTol;
      biCGStabSolver->m_imax =  m_config.m_maxBiCGStabIter; 
      biCGStabSolver->m_normType = m_config.m_normType;
      //JFNK mode will only work well if 
      //the linear system is solved quickly.
      //By being intolerant here, we revert to 
      //the cheaper Picard mode sooner rather than later
      if (a_mode == JFNK_LINEARIZATION_MODE)
	{
	  biCGStabSolver->m_numRestarts = 0;
	  biCGStabSolver->m_hang =  m_config.m_RelaxHang;
	}
      solver = biCGStabSolver;
    }
  else if (m_config.m_linearSolverType == Configuration::CG)
    {
      CGSolver<Vector<LevelData<FArrayBox>* > >* cGSolver 
	= new CGSolver<Vector<LevelData<FArrayBox>* > >;
      
      cGSolver->define(&a_op , false);
      cGSolver->m_verbosity = m_config.m_verbosity - 1;
      cGSolver->m_eps =  m_config.m_CGRelTol;
      cGSolver->m_imax =  m_config.m_maxCGIter;
      cGSolver->m_normType =  m_config.m_normType;
      
      cGSolver->m_numRestarts = 0;
      cGSolver->m_hang = 1.0;
      //JFNK mode will only work well if 
      //the linear system is solved quickly.
      //By being intolerant here, we revert to 
      //the cheaper Picard mode sooner rather than later
      if (a_mode == JFNK_LINEARIZATION_MODE)
	{
	  cGSolver->m_numRestarts = 0;
	  cGSolver->m_hang = 1.0;
	}
      solver = cGSolver;
    }
  else if (m_config.m_linearSolverType == Configuration::GMRES)
    {
      GMRESSolver<Vector<LevelData<FArrayBox>* > >* gmresSolver 
	= new GMRESSolver<Vector<LevelData<FArrayBox>* > >();
      
      gmresSolver->define(&a_op , false);
      gmresSolver->m_verbosity = m_config.m_verbosity - 1;
      gmresSolver->m_reps =  m_config.m_GMRESRelTol;
      gmresSolver->m_imax = m_config.m_maxGMRESIter;
      gmresSolver->m_normType =  m_config.m_normType;
      solver = gmresSolver;
    }
  else if (m_config.m_linearSolverType == Configuration::Relax)
    {
      
      RelaxSolver<Vector<LevelData<FArrayBox>* > >* relaxSolver
	= new RelaxSolver<Vector<LevelData<FArrayBox>* > >();
      
      relaxSolver->define(&a_op,false);
      relaxSolver->m_verbosity = m_config.m_verbosity - 1;
      relaxSolver->m_normType =  m_config.m_normType;
      relaxSolver->m_eps =  m_config.m_RelaxRelTol;
      relaxSolver->m_imax =  m_config.m_maxRelaxIter;
      relaxSolver->m_hang =  m_config.m_RelaxHang;
      solver = relaxSolver;
    }
  else
    {
      CH_assert(solver != NULL);
      MayDay::Error("JFNKSolver::newLinearSolver : unknown linear solver type");
    }
  return solver;
}


Real JFNKSolver::lineSearch(Vector<LevelData<FArrayBox>* >& a_u, 
			    Vector<LevelData<FArrayBox>* >& a_residual,
			    const Vector<LevelData<FArrayBox>* >& a_rhs,
			    const Vector<LevelData<FArrayBox>* >& a_du, 
			    LinearizedVTOp& a_op,
			    NonlinearViscousTensor& a_nvt, Real a_minW, 
			    bool a_resetOnFail)
{
  CH_TIME("JFNKSolver::lineSearch");
  //try u + du, and if the residual is not reduced try a sequence
  //of smaller steps u + w*du, halving w until either the residual
  // is reduced or w < a_minW
  Real w = 1.0; 
  Real oldResNorm = a_op.norm(a_residual, m_config.m_normType);
  Real resNorm = oldResNorm;
  a_op.incr(a_u,a_du,w);


  do {
    
    LinearizedVTOp testOp (&a_nvt, a_u, m_config.m_h, m_config.m_err, m_config.m_umin, m_config.m_hAdaptive, m_grids, 
			   m_refRatios, m_domains, m_dxs, 0, 
			   m_config.m_numMGSmooth, m_config.m_numMGIter,PICARD_LINEARIZATION_MODE );
    testOp.m_writeResiduals = m_config.m_writeResiduals;
    testOp.outerResidual(a_residual,a_u,a_rhs);
    resNorm = testOp.norm(a_residual, m_config.m_normType);
    if (resNorm >= oldResNorm )
      {
	if (m_config.m_verbosity > 0)
	  pout()  << "JFNKSolver::lineSearch residual norm = " << resNorm << std::endl;
	w *= 0.5;
	
	if (w >= a_minW)
	  {
	    if (m_config.m_verbosity > 0)
	      pout()  << "JFNKSolver::lineSearch halving step length, w =   " << w << std::endl; 
	    //step halfway back to the last U
	    testOp.incr(a_u,a_du,-1.0*w);	
	    a_nvt.setState(a_u); 
	  }
	else if (a_resetOnFail)
	  {
	    //give up, return to the last U
	    testOp.incr(a_u,a_du,-2.0*w); 
	    a_nvt.setState(a_u);  
	    testOp.outerResidual(a_residual,a_u,a_rhs);
	    resNorm = testOp.norm(a_residual, m_config.m_normType);
	  } 
      }
  } while ( resNorm >= oldResNorm && w >= a_minW);
  
return resNorm;
}




void 
JFNKSolver::linearSolve(LinearizedVTOp& a_op, 
			Vector<LevelData<FArrayBox>* >& a_u,
			const Vector<LevelData<FArrayBox>* >& a_rhs, 
			LinearizationMode a_mode)
{
  if (m_config.m_linearSolverType == Configuration::petsc)
    {
#ifdef CH_USE_PETSC
      // petsc solver setup is expensive and it's easy 
      // to swap the coefficients, so keep a pre-defined 
      // solver around and just reset the coefficients
      if (m_petscSolver == NULL)
	{
	  m_petscSolver = new PetscAMRSolver(m_config.m_verbosity);
	  m_petscSolver->m_petscCompMat.setVerbose(m_config.m_verbosity-1);
	}
      RefCountedPtr<CompBC> bcfunc = RefCountedPtr<CompBC>(m_bcPtr->velocitySolveBC());
      BCHolder bc(bcfunc);
      Real opAlpha, opBeta;
      opAlpha = -1.0;
      opBeta = 1.0;
      m_petscSolver->m_petscCompMat.define(m_domains[0],m_grids,m_refRatios,bc,m_dxs[0]); // generic AMR setup
      m_petscSolver->m_petscCompMat.defineCoefs(opAlpha,opBeta,a_op.current().mu(),
						a_op.current().lambda(),a_op.current().alpha());

      m_petscSolver->solve_mfree(a_u,a_rhs,&a_op);
      //CHKERRQ(ierr);
#else
      MayDay::Error("linearSolverType is petsc, but code compiled w/o PETSC=TRUE");
#endif
    }
  else 
    {  
      // Chombo solver setup is cheap and it's hard 
      // to swap the coefficients, so create a new solver
      LinearSolver<Vector<LevelData<FArrayBox>* > >* krylovSolver = 
	newLinearSolver(a_op,a_mode);
      krylovSolver->solve(a_u, a_rhs);
      delete krylovSolver;	      
    }
}

#include "NamespaceFooter.H"


/*
  PETSc Licensing Notification
  
  Permission to use, reproduce, prepare derivative works, and to redistribute to others this software, derivatives of this software, and future versions of this software as well as its documentation is hereby granted, provided that this notice is retained thereon and on all copies or modifications. This permission is perpetual, world-wide, and provided on a royalty-free basis. UChicago Argonne, LLC and all other contributors make no representations as to the suitability and operability of this software for any purpose. It is provided "as is" without express or implied warranty.
  Principal Software authors
  
  
  
  Mathematics and Computer Science Division
  Argonne National Laboratory,
  Argonne IL 60439
  Any questions or comments on the software may be directed to petsc-maint@mcs.anl.gov.
  
  Portions of this software are copyright by UChicago Argonne, LLC. Argonne National Laboratory with facilities in the state of Illinois, is owned by The United States Government, and operated by UChicago Argonne, LLC under provision of a contract with the Department of Energy.
  DISCLAIMER
  
  PORTIONS OF THIS SOFTWARE WERE PREPARED AS AN ACCOUNT OF WORK SPONSORED BY AN AGENCY OF THE UNITED STATES GOVERNMENT. NEITHER THE UNITED STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR THE UNIVERSITY OF CHICAGO, NOR ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. REFERENCE HEREIN TO ANY SPECIFIC COMMERCIAL PRODUCT, PROCESS, OR SERVICE BY TRADE NAME, TRADEMARK, MANUFACTURER, OR OTHERWISE, DOES NOT NECESSARILY CONSTITUTE OR IMPLY ITS ENDORSEMENT, RECOMMENDATION, OR FAVORING BY THE UNITED STATES GOVERNMENT OR ANY AGENCY THEREOF. THE VIEW AND OPINIONS OF AUTHORS EXPRESSED HEREIN DO NOT NECESSARILY STATE OR REFLECT THOSE OF THE UNITED STATES GOVERNMENT OR ANY AGENCY THEREOF. */
