#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INVERSE_VI_VEL_SOLVER_H_
#define _INVERSE_VI_VEL_SOLVER_H_

#include "JFNKSolver.H"
#include "LevelDataOps.H"
#include "SurfaceFlux.H"
#include "BasalFriction.H"
#include "MuCoefficient.H"
#include "IceConstants.H"
#include "LevelDataBasalFriction.H"
#include "CGOptimize.H"
#include "NamespaceHeader.H"


/// Abstract subclass of IceVelocitySolver  that supports inverse velocity problems
/** 
    In this context, inverse velocity problems compute basal friction (\f$C\f$) and stiffness (\f$\phi\f$)
    parameters for a stress balance problem such that its solution (\f$u\f$) is a fit to the a 
    known (observed) velocity (\f$u_o\f$). Just as any IceVelocitySolver, they must compute
    a velocity given geometry etc, but must also expose BasalFriction and MuCoefficient 
    data to any user that understands the idea of the inverse velocity problem

    Includes static members that could be useful to both implementors and users
*/
class InverseIceVelocitySolver :  public IceVelocitySolver
{

  static void copyAMR ( Vector<RefCountedPtr<LevelData<FArrayBox> > >&a_copy, 
			const Vector<LevelData<FArrayBox>*>& a_orig )
  {
    a_copy.resize(0);

    for (int lev = 0; lev < a_orig.size(); lev++)
      {
	if (a_orig[lev])
	  {
	    const DisjointBoxLayout& dbl = a_orig[lev]->disjointBoxLayout();
	    LevelData<FArrayBox>* lptr = new LevelData<FArrayBox>
	      (dbl, a_orig[lev]->nComp(),a_orig[lev]->ghostVect());
	    for (DataIterator dit(dbl); dit.ok(); ++dit)
	      {
		(*lptr)[dit].copy( (*a_orig[lev])[dit]);
	      }
	    a_copy.push_back(RefCountedPtr<LevelData<FArrayBox> >( lptr));
	  }	
      }
  }


public:
  
  virtual void setPreviousTime(Real a_time) =0;
  
  /**
     Create a BasalFriction object given an AMR hierarchy of data

     @param  a_C     AMR hierarchy of friction coefficient\f$C\f$
     @param  a_ratio  refinement ratios
     @param  a_dxCrse mesh spacing in coarsest level
     @return a pointer to a MultiLevelDataBasalFriction object
  */ 
  static BasalFriction* basalFriction
  (const Vector<LevelData<FArrayBox>*>& a_C, const Vector<int>& a_ratio, const RealVect& a_dxCrse )
  {
    Vector<RefCountedPtr<LevelData<FArrayBox> > > C;
    copyAMR(C, a_C);
    MultiLevelDataBasalFriction* ptr = new  MultiLevelDataBasalFriction(C,a_dxCrse,a_ratio);
    return ptr;
  }

  /**
     Create a MuCoefficient object given an AMR hierarchy of data

     @param  a_mu     AMR hierarchy of stiffness coefficient \f$\phi\f$
     @param  a_ratio  refinement ratios
     @param  a_dxCrse mesh spacing in coarsest level
     @return a pointer to a MultiLevelDataMuCoefficient object
  */
  static MuCoefficient* muCoefficient
  (const Vector<LevelData<FArrayBox>*>& a_mu, const Vector<int>& a_ratio, const RealVect& a_dxCrse )
  {
    Vector<RefCountedPtr<LevelData<FArrayBox> > > mu;
    copyAMR(mu, a_mu);
    MultiLevelDataMuCoefficient* ptr = new  MultiLevelDataMuCoefficient(mu,a_dxCrse,a_ratio);
    return ptr;
  }

  /// Create a BasalFriction object containing the solution \f$C\f$ 
  virtual BasalFriction* basalFriction() = 0;
  /// is the basal friction an improvement on the initial guess?
  
  /// Create a MuCoefficient object containing the solution \f$\phi\f$ 
  virtual MuCoefficient* muCoefficient() = 0;


  
  
};


///Implementation of InverseIceVelocitySolver for the vertically integrated (SSA,SSA*,L1L2) stress models
/**

   Computes a velocity \f$u\f$, a basal friction coefficient \f$C\f$ and a stiffening coefficient
   \f$\phi\f$ which mimimises an objective function

   \f$
   J = || f(u,u_o) ||^2_2 + || g(\nabla \cdot ( u h)  , \nabla \cdot ( u_o h_o) ) ||^2_2 
   +  \alpha _{0} || \nabla C ||^2_2  + \alpha _{1} || \nabla \phi ||^2_2  
   +  \alpha _{2} || \nabla X_0 ||^2_2  + \alpha _{3} || \nabla X_1 ||^2_2
   +  \alpha _{4} || X_0 ||^2_2  + \alpha _{5} || X_1 ||^2_2
   \f$

   This is essentially the inverse problem described in e.g Cornford et al, 2015, The Cryopshere. 

   Because \f$C\f$ and \f$\phi\f$ are positive, they are expressed as \f$C = C_0 \exp (X_0) \f$
   and \f$\phi = \phi_0 \exp (X_1) \f$.

   The functions f,g are composed from a few parts to: 
   - weight observations according to their accuracy. For example, f(u,u_o)  can be zero in regions where   u_o is unknown
   - choose a particular objective function, e.g, f (u,u_o) = |u| - |u_o|  can be selected to match speeds rather than velocities.
   
   The first part of J is minimal when the model and observed velocities match, given the form of f.

   The second part is minimal when the model  and observed flux divergence  match. As this is often a noisy field, there is 
   an option to smooth it 
   
   The remaining parts of \f$J\f$ are used to regularize the solution. Any or all
   the \f$\alpha _{n}\f$ can be zero. Typically, in time-independent problems 
   only \f$\alpha _{0}\f$ and \f$\alpha _{1}\f$ are non-zero - this is Tikhonov regularization. 
   
   New to the code, and experimental, is some support for time independent problems. For now
   this is limited to recomputing \f$C(t)\f$ and  \f$\phi(t)\f$ to minimise a \f$ J(t) \f$.
   The factors \f$\alpha _{4} , \alpha _{5} \f$ can be made to depend on the time-step, so that
   \f$C(t + \Delta t)\f$ will tend to be closer to \f$C(t)\f$ as  \f$ \Delta t \rightarrow 0\f$.

   J is reduced in a series of nonlinear conjugate gradient (CG) iterations. To that end,
   the class implements the interface implied by the function template CGOptimize. 
   The major method is  computeObjectiveAndGradient, which calculates J and
   its directional derivatives with respect to X_0 and X_1. This function also
   writes the current state of the problem to hdf5 files, either every 'outer' CG iteration, 
   or even more frequently, on every call, which will correspond to CG's 'inner' (secant) 
   iterations
   
   Uses JFNKSolver to solve a number of internal problems. 
*/
class InverseVerticallyIntegratedVelocitySolver : public InverseIceVelocitySolver, 
						  public ObjectiveWithGradient<Vector<LevelData<FArrayBox>*> >
{

public:

  /// Data read through ParmParse
  class Configuration
  {
  public: 

    Configuration();

    ~Configuration();

    /// parse parameters from <prefix>.* section of ParmParse table
    /** 
	@param a_prefix prefix string in the ParmParse table
     */
    void parse(const char* a_prefix = "control");

    /// Do not attempt to optimize unless the mesh has at least m_minLevel levels
    /**
       Do not atttempt optimizion until the mesh has been 
       refined to some desired level. Instead, just fill \f$ C \f$ and \f$ \phi \f$
       with the initial guess, and return the observed velocity
       as the solution. This really only makes sense if the 
       we want to use the observed velocity in a refinemnt criterion.
       That approach was sucessful in the older AMRIceControl class
       (which did its own refinement), so this parameter exists
       to support that
       
     */
    int m_minLevelForOptimization;

    /// time delay between optmizations. Allows time-dependent inverse problems
    /// without re-doing the optimization at every time step.
    Real m_minTimeBetweenOptimizations;

    /// estimate of typical timestep. Just here to turn times into integers for file names
    Real m_dtTypical;
    
    enum BoundMethod {none, projection, MAX_BOUND_METHOD};
    /// Method to enforce constraints: none, or projection. default projection
    BoundMethod m_boundMethod;
    
    enum VelocityMisfitType {speed, velocity, log_speed, MAX_VELOCITY_MISFIT_TYPE};
    /// Type of velocity misfit : speed, velocity, or log_speed. Set with <prefix>.vel_misfit_type, default speed
    VelocityMisfitType m_velMisfitType;

    /// attempt to optimize w.r.t X0?
    bool m_optimizeX0;
    /// attempt to optimize w.r.t X1?
    bool m_optimizeX1;
    /// evaluate unregluarized gradient of objective with respect to mucoef (or rather, X1) in the shelf only?
    bool m_gradMuCoefShelfOnly;
    
    /// Lower bound on \f$X_0 = \ln C - \ln C_0\f$. Set with with <prefix>.lowerX0
    Real m_lowerX0;

    /// Upper bound on \f$X_0 = \ln C - \ln C_0\f$. Set with with <prefix>.upperX0
    Real m_upperX0; 

    /// Lower bound on \f$X_1 = \ln \phi - \ln \phi_0\f$. Set with with <prefix>.lowerX1
    Real m_lowerX1;

    /// Upper bound on \f$X_1 = \ln \phi - \ln \phi_0\f$. Set with with <prefix>.upperX1
    Real m_upperX1; 

    /// Lower bound on \f$\phi_0\f$. Set with <prefix>.initialLowerMuCoef
    Real m_initialLowerMuCoef;

    /// Upper bound on \f$\phi_0\f$. Set with <prefix>.initialUpperMuCoef
    Real m_initialUpperMuCoef;  

    /// Lower bound on \f$C_0\f$. Set with <prefix>.initialLowerC
    Real m_initialLowerC;

    /// Upper bound on \f$C_0\f$. Set with <prefix>.initialUpperC
    Real m_initialUpperC;

    /// Coefficient of \f$\int _\Omega \alpha (u - u_0) d \Omega \f$ in \f$J\f$. Set with <prefix>.velMisfitCoefficient, default 1.0
    Real m_velMisfitCoefficient;

    /// Coefficient of \f$\int _\Omega \beta (\nabla \cdot. (uh) - a^{\prime} ) d \Omega \f$ in \f$J\f$. Set with <prefix>.divuhMisfitCoefficient default 1.0
    Real m_divuhMisfitCoefficient;

    /// If \f$> 0\f$, smooth \f$(\nabla \cdot. (uh) - a^{\prime} )\f$. Set with <prefix>.divuhMisfitSmooth, default 0.0
    Real m_divuhMisfitSmooth;

    /// Coefficient of \f$\int _\Omega (\nabla C )^2 d \Omega \f$ in \f$J\f$. Set with <prefix>.gradCsqRegularization
    Real m_gradCsqRegularization; 

    /// Coefficient of \f$\int _\Omega (\nabla \phi )^2 d \Omega \f$ in \f$J\f$. Set with <prefix>.gradMuCoefsqRegularization
    Real m_gradMuCoefsqRegularization; 

    

    
    /// Coefficient of \f$\int _\Omega (\nabla X_0 )^2 d \Omega \f$ in \f$J\f$. Set with <prefix>.gradX0sqRegularization
    Real m_gradX0sqRegularization; 

    /// Coefficient of \f$\int _\Omega (\nabla X_1 )^2 d \Omega \f$ in \f$J\f$. Set with <prefix>.gradX1sqRegularization
    Real m_gradX1sqRegularization;

    /// Coefficient of \f$\int _\Omega ( X_0 )^2 d \Omega \f$ in \f$J\f$. Set with <prefix>.X0Regularization
    Real m_X0Regularization;

    /// Coefficient of \f$\int _\Omega ( X_1 )^2 d \Omega \f$ in \f$J\f$. Set with <prefix>.X1Regularization
    Real m_X1Regularization;

    /// Coefficient of \f$\int _\Omega ( X_0 )^2 / \Delta t d \Omega \f$ in \f$J\f$. Set with <prefix>.X0TimeRegularization
    Real m_X0TimeRegularization;

    /// Coefficient of \f$\int _\Omega ( X_1 )^2 / \Delta t d \Omega \f$ in \f$J\f$. Set with <prefix>.X1TimeRegularization
    Real m_X1TimeRegularization;

    /// Limit domain of integration \f$\Omega\f$ to regions where ice thickness exceeds this value. Set with <prefix>.thicknessThreshold, default 100.0 m
    Real m_thicknessThreshold;

    /// Maximum number of 'outer' nonlinear conjugate gradient iterations. Set with <prefix>.CGmaxIter, default 16
    Real m_CGmaxIter;

    /// Tolerance for the 'outer' nonlinear conjugate gradient method. Set with <prefix>.CGtol, default 1.0e-3
    Real m_CGtol;

    /// Size of the initial 'inner' secant step. Set with <prefix>.CGsecantParameter, default 1.0e-7
    Real m_CGsecantParameter;

    /// \f$a\f$ such that secant step \f$n+1\f$ has a maximum length \f$l_{n+1} \leq a l_{n}\f$. Set with <prefix>.CGsecantStepMaxGrow, default 2.0
    Real m_CGsecantStepMaxGrow;

    /// Maximum number of 'outer' secant iterations. Set with <prefix>.CGsecantMaxIter, default 20
    int m_CGsecantMaxIter;

    /// Tolerance for the 'inner' secant method. Set with <prefix>.CGsecantTol, default 1.0e-1
    Real m_CGsecantTol;

    /// Give up when succesive CG steps don't progress quickly enough . Set with <prefix>.CGhang, default 0.999
    Real m_CGhang;

    /// Restart CG at a fixed interval. Set with <prefix>.CGrestartInterval, default 9999
    Real m_CGrestartInterval;

    ///observed velocity data coefficient
    SurfaceData* m_velObs_c;
    
    ///observed x-velocity data
    SurfaceData* m_velObs_x;
    
#if CH_SPACEDIM > 1
    ///observed y-velocity data
    SurfaceData* m_velObs_y;
#endif

    /// write inner (secant) steps to plot files
    bool m_writeInnerSteps;

    /// base of filename for 'inner' plots
    std::string m_innerStepFileNameBase;

    /// base of filename for 'outer' plots
    std::string m_outerStepFileNameBase;

    /// observed flux divergence data coefficient
    SurfaceData* m_divuhObs_c;

    /// observed flux divergence data
    SurfaceData* m_divuhObs_a;    

    /// used to amplify (damp) the gradient non-uniformly
    SurfaceData* m_gradientFactor;

    int m_outerCounter;

    

  };

private:

  Configuration m_config;

  ///standard operations for LevelData<FArrayBox> objects
  Vector<LevelDataOps<FArrayBox> > m_vectOps;

  /// create an AMR Hierarchy worth of data
  template <class T>
  void create(Vector<LevelData<T>*>& a_data, 
	      const int a_ncomp, const IntVect& a_ghost)
  {
    if (a_data.size() <= m_finest_level)
      a_data.resize(m_finest_level+1);

    for (int lev =0; lev <= m_finest_level; lev++)
      {
	if (a_data[lev] != NULL)
	  delete a_data[lev];
	a_data[lev] = new LevelData<T>(m_grids[lev], a_ncomp, a_ghost);
      }
  }

  /// free (any) storage in AMR Hierarchy a_a
  template <class T>
  void free(Vector<LevelData<T>*>& a_data)
  {
    for (int lev =0; lev < a_data.size(); lev++)
      {
	if (a_data[lev] != NULL)
	  {
	    delete a_data[lev];
	    a_data[lev] = NULL;
	  }
      }
  }

  // set a_x = 0.0
  void setToZero(Vector<LevelData<FArrayBox>* >& a_a);


  void  mapX(const Vector<LevelData<FArrayBox>* >& a_x);

  void writeState
  (const std::string& a_file, int a_counter,
   const Vector<LevelData<FArrayBox>* >& a_x,
   const Vector<LevelData<FArrayBox>* >& a_g) const;


  ///compute cell-centered div(UH) 
  void computeDivUH();

  //calculate the rhs to the adjoint equation, and the misfits
  void computeAdjointRhs();

  // add directional derivatives of unregularized problem to a_g
  void computeGradient(Vector<LevelData<FArrayBox>* >& a_g, 
		       const  Vector<LevelData<FArrayBox>* >& a_x);

  /// add Tikhonov regularization to a_g
  void regularizeGradient(Vector<LevelData<FArrayBox>* >& a_g, 
			  const  Vector<LevelData<FArrayBox>* >& a_x);

  /// modify a_g at constraint boundaries
  void applyProjection(Vector<LevelData<FArrayBox>* >& a_g, 
		       const  Vector<LevelData<FArrayBox>* >& a_x);
    
  ///solve the forward or adjoint stress balance equation
  void solveStressEqn (Vector<LevelData<FArrayBox>* >& a_u,
		       const bool a_adjoint,
		       const Vector<LevelData<FArrayBox>* >& a_rhs,
		       const Vector<LevelData<FArrayBox>* >& a_C,
		       const Vector<LevelData<FArrayBox>* >& a_C0,
		       const Vector<LevelData<FArrayBox>* >& a_A,
		       const Vector<LevelData<FArrayBox>* >& a_muCoef);
  /// Current time
  Real m_time; 

  /// Time of the previous call to solve() 
  Real m_prev_time;
  
  /// effectively, 1/Delta t
  Real XnTimeCoef()
  {
    return  (m_time > m_prev_time + TINY_NORM)?( 1.0 / (m_time-m_prev_time)):0.0;
  }

  /// combines timestep-dependent and timestep-independent coefficients of \f$\int _\Omega ( X_0 )^2 d \Omega \f$ in \f$J\f$
  Real X0Regularization()
  {
    return m_config.m_X0Regularization +  m_config.m_X0TimeRegularization * XnTimeCoef();
  }

  /// combines timestep-dependent and timestep-independent coefficients of \f$\int _\Omega ( X_1 )^2 d \Omega \f$ in \f$J\f$
  Real X1Regularization()
  {
    return m_config.m_X1Regularization +  m_config.m_X1TimeRegularization * XnTimeCoef();
  }

  int m_innerCounter, m_outerCounter;
  
  bool m_velocityInitialised;

  ///finest level 
  int m_finest_level;

  /// grids for each level
  Vector<DisjointBoxLayout> m_grids;
  
  /// refinement ratios  for each level
  Vector<int> m_refRatio;
  
  /// mesh spacings for each level
  Vector<RealVect> m_dx;

  /// problem domains for each level
  Vector<ProblemDomain> m_domain;

  /// write state to outerStateFile() on 'outer' iterations
  std::string outerStateFile() const;

  /// write state to innerStateFile() on 'inner' iterations
  std::string innerStateFile() const;

  /// co-ordinate system, stores thickness, topography, surfaces, mask
  Vector<RefCountedPtr<LevelSigmaCS > > m_coordSys;
  
  /// rate factor, and 'c0'
  Vector<LevelData<FArrayBox>*> m_A;
  
  /// 'C_0 \todo rename this'
  Vector<LevelData<FArrayBox>*> m_C0;

  // coefficients of basal friction.
  // m_C = exp(X[0]) * m_COrigin;
  Vector<LevelData<FArrayBox>*> m_COrigin; // original value of C
  Vector<LevelData<FArrayBox>*> m_C      ; // current value of C
  Vector<LevelData<FArrayBox>*> m_Cmasked ; // C with shelf values set 
  Vector<LevelData<FArrayBox>*> m_lapC;
  Vector<LevelData<FArrayBox>*> m_gradCSq;

  /// coefficients of mu.
  /// m_muCoef = exp(X[1]) * m_referenceMuCoefOrigin;
  Vector<LevelData<FArrayBox>*> m_muCoefOrigin; // original value of muCoef
  Vector<LevelData<FArrayBox>*> m_muCoef      ; // current value of muCoef
  Vector<LevelData<FArrayBox>*> m_lapMuCoef;
  Vector<LevelData<FArrayBox>*> m_gradMuCoefSq;
  //Vector<LevelData<FluxBox>*> m_faceMuCoef;
  
  /// X-dependents
  Vector<LevelData<FArrayBox>*> m_lapX;
  Vector<LevelData<FArrayBox>*> m_gradXSq;
  
  ///right-hand side of the velocity equation
  Vector<LevelData<FArrayBox>*> m_rhs;

  ///right-hand side of the adjoint velocity equation
  Vector<LevelData<FArrayBox>*> m_adjRhs;

  ///lagrange multipler (at base of ice)
  Vector<LevelData<FArrayBox>*> m_adjVel;
  
  ///velocity observations on AMR hierarchy
  Vector<LevelData<FArrayBox>*> m_velObs;

  ///local velocity observation confidence (1/sigma^2) on AMR hierarchy
  Vector<LevelData<FArrayBox>*> m_velCoef;

  ///flux divergence observations on AMR hierarchy
  Vector<LevelData<FArrayBox>*> m_divuhObs;

  ///flux divergence observation confidence (1/sigma^2) on AMR hierarchy
  Vector<LevelData<FArrayBox>*> m_divuhCoef;

  ///velocity (at base of ice)
  Vector<LevelData<FArrayBox>*> m_velb;

  ///velocity (at surface of ice)
  Vector<LevelData<FArrayBox>*> m_vels;

  ///flux divergence
  Vector<LevelData<FArrayBox>*> m_divuh;
  
  /// misfit term  \f$|| f(u,u_o) ||^2_2\f$
  Vector<LevelData<FArrayBox>*> m_velocityMisfit;
  
  /// misfit term \f$|| g(\nabla \cdot ( u h)  , \nabla \cdot ( u_o h_o) ) ||^2_2\f$
  Vector<LevelData<FArrayBox>*> m_divuhMisfit;

  //physics pointers
  BasalFrictionRelation* m_basalFrictionRelation;
  ConstitutiveRelation* m_constitutiveRelation;
  IceThicknessIBC* m_thicknessIBC;

  Real m_bestMisfit;
  Vector<LevelData<FArrayBox>*> m_bestVel;
  Vector<LevelData<FArrayBox>*> m_bestC;
  Vector<LevelData<FArrayBox>*> m_bestMuCoef;
  bool m_optimization_done; // don't want to provide the dat aabove unless it is meaningful
  
  const AmrIceBase* m_amrIce;

  // need to pass these to functions like JFNKSolver::solve
  Vector<LevelData<FArrayBox>* > m_calvedIce, m_addedIce, m_removedIce;
  
public:

  InverseVerticallyIntegratedVelocitySolver();

  ~InverseVerticallyIntegratedVelocitySolver();

  void setPreviousTime(Real a_time)
  {
    m_prev_time = a_time;
  }
  
  void define(const ProblemDomain& a_coarseDomain,
	      ConstitutiveRelation* a_constRel,
	      BasalFrictionRelation* a_FrictionRel,
	      const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<int>& a_refRatio,
	      const RealVect& a_dxCrse,
	      IceThicknessIBC* a_bc,
	      int a_numLevels)
  {
    MayDay::Error("not sure why we got here...");
  }

  void define(const AmrIceBase& a_amrIce,
	      const ProblemDomain& a_coarseDomain,
	      ConstitutiveRelation* a_constRel,
	      BasalFrictionRelation* a_FrictionRel,
	      const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<int>& a_refRatio,
	      const RealVect& a_dxCrse,
	      IceThicknessIBC* a_bc,
	      int a_numLevels);
  
 
  int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
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




  //implementation of  ObjectiveWithGradient
  void computeObjectiveAndGradient
  (Real& a_fm, Real& a_fp,
   Vector<LevelData<FArrayBox>* >& a_g, 
   const  Vector<LevelData<FArrayBox>* >& a_x, 
   bool a_inner);
  
  void restart(){;}

  //allocate storage for a_ca, copy from a_b
  void create(Vector<LevelData<FArrayBox>* >& a_a,  const  Vector<LevelData<FArrayBox>* >& a_b);
  
  //apply preconditioner s = M^{-1}r
  void preCond(Vector<LevelData<FArrayBox>* >& a_s,  const Vector<LevelData<FArrayBox>* >& a_r);


  ///need to make this explicit
  void free(Vector<LevelData<FArrayBox>*>& a_data)
  {
    free<FArrayBox>(a_data);
  }


  // add a scalar to a_a
  void plus(Vector<LevelData<FArrayBox>* >& a_a, Real a_r );

  // set a_x = s * a_x
  void scale(Vector<LevelData<FArrayBox>* >& a_x, 
	     const  Real a_s);

  // set a_y = a_x
  void assign(Vector<LevelData<FArrayBox>* >& a_y, 
	      const Vector<LevelData<FArrayBox>* >& a_x);
  
  // set a_y = a_y + a_s * a_x
  void incr(Vector<LevelData<FArrayBox>* >& y, 
	    const Vector<LevelData<FArrayBox>* >& x, Real s);
  
  // return a_y.a_x
  Real dotProduct(Vector<LevelData<FArrayBox>* >& a_y, 
		  const Vector<LevelData<FArrayBox>* >& a_x);

  // no of degrees of freedom
  int nDoF(const Vector<LevelData<FArrayBox>* >& x);

  BasalFriction* basalFriction()
  {
    BasalFriction* ptr = NULL;
    if (m_optimization_done)
      ptr = InverseIceVelocitySolver::basalFriction(m_bestC, m_refRatio, m_dx[0]);
    return ptr;
  }
  
  MuCoefficient* muCoefficient()
  {
    MuCoefficient* ptr = NULL;
    if (m_optimization_done)
      ptr =  InverseIceVelocitySolver::muCoefficient(m_bestMuCoef, m_refRatio, m_dx[0]);
    return ptr;
  }
};



#include "NamespaceFooter.H"
#endif
