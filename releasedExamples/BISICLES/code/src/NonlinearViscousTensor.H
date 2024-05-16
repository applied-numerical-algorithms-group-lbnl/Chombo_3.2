#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NONLINEARVISCOUSTENSOR_H_
#define _NONLINEARVISCOUSTENSOR_H_


#include "IceVelocitySolver.H"
#include "AMRMultiGrid.H"
#include "MultilevelIceVelOp.H"
#include "ViscousTensorOp.H"
#ifdef CH_USE_PETSC
#include "PetscAMRSolver.H"

#endif

//#31556926.0;
#include "NamespaceHeader.H"

/** 
    Virtual base for classes which encapsulate a nonlinear viscosity tensor. 
   
    States are defined by the velocity, and the important dependents are the 
    coefficients in the linear ViscousTensorOp: alpha, mu, lambda. 
*/ 
class NonlinearViscousTensor
{
public:
  virtual ~NonlinearViscousTensor(){};
  ///factory method
  virtual NonlinearViscousTensor* newNonlinearViscousTensor()=0;
  /// set the state (velocity)
  /**
     Modify state, and thus the  coefficients in the linear ViscousTensorOp: alpha, mu, lambda. 

     @param[in] a_state an amr hierarchy of velocity data
   */
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_state)=0;

  /// fill ghost cells of a_state as per the boundary conditions
  virtual void applyBC(Vector<LevelData<FArrayBox>*>& a_phi)=0;

  /// retrieve the state
  /**
     @return a reference to the current state (velocity)
   */
  virtual const Vector<LevelData<FArrayBox>*>& getState() = 0;

  ///Retrieve the drag coefficient,  alpha
  /**
     @return an amr hierarchy representing the coefficient \f$ \alpha(x,y) \f$, 
     that appears in ViscousTensorOp as the coeffcient of the velocity \f$ u \f$.
   */
  virtual Vector<RefCountedPtr< LevelData<FArrayBox> > >& alpha() =0;

  ///Retrieve the viscosity coefficient, mu
  /**
     @return an amr hierarchy representing the coefficient \f$ \mu(x,y) \f$, 
     that appears in ViscousTensorOp as the coeffcient of the strain rate 
     \f$ \nabla \vec{u}  + ( \nabla \vec{u})^T \f$.
   */  
  virtual Vector<RefCountedPtr< LevelData<FluxBox> > >& mu()  =0;

  ///Retrieve the compressibility coefficient, lambda
  /**
     @return an amr hierarchy representing the coefficient \f$ \lambda(x,y) \f$, 
     that appears in ViscousTensorOp as the coeffcient of the compression rate 
     \f$ \nabla. \vec{u} I \f$.
   */
  virtual Vector<RefCountedPtr< LevelData<FluxBox> > >& lambda() =0;

  /// provide a pointer to an AMRLevelOpFactory<LevelData<FArrayBox> > defined at the current state
  virtual RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > opFactoryPtr()=0;
};

/// Implement NonlinearViscousTensor for typical (velocity-depenendent) ice rheology
/** 
    Provides the means to compute the coefficients \f$ \mu(u) \f$, 
    \f$ \lambda(u)  = 2*\mu(u) \f$ and \f$ \alpha(u) \f$ 
    in \f$ L[u] = \nabla \cdot ( \mu(h,A(T),u)  (\nabla u + (\nabla u)^T + \lambda  \nabla. u I)  
    - \alpha(u) u \f$
    and provide the appropriate ViscousTensorOpFactory for a given \f$u\f$. Also provides
    the means to compute the stress tensor components 
    
    The velocity dependencies of \f$\alpha(u)\f$ and \f$\mu(u)\f$, \f$\lambda(u)\f$ are determined 
    by BasalFrictionRelation and ConstitutiveRelation  objects, while the velocity 
    independent parts are set seperately, in 
    the contructor and can be modified using setState, and setFaceViscCoef.

    The effective basal drag coefficient is computed as \f$alpha(u) = C  f_b(u)  + C_0\f$ , 
    with fb defined by a BasalFrictionRelation, ie \f$ \tau _b = [ C  f_b(u)  + C_0 ]  u \f$

    The effective viscosity is computed as \f$ \mu(u) = \phi  h  f_v(A(T),u)\f$,
    with \f$f_v \f$
    defined by a ConstitutiveRelation, and \f$ \lambda(u) = 2  \mu(u)\f$
  
*/
class IceNonlinearViscousTensor : public NonlinearViscousTensor
{
 
  Vector<LevelData<FArrayBox>*> m_u;
  Vector<LevelData<FArrayBox>*> m_C;
  Vector<LevelData<FArrayBox>*> m_C0;
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > m_opFactoryPtr;
  
  //mesh parameters
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_refRatio;
  Vector<ProblemDomain> m_domains;
  Vector<RealVect> m_dxs;
  int m_finestLevel;
  Vector<RefCountedPtr<LevelSigmaCS > > m_coordSys;
  
  //constitutive relations
  const ConstitutiveRelation* m_constRelPtr;
  const BasalFrictionRelation* m_basalFrictionRelPtr;
  
  //boundary conditions
  IceThicknessIBC* m_bcPtr;

  //pointers to temperature fields
  Vector<LevelData<FArrayBox>*> m_A;
  Vector<LevelData<FluxBox>*> m_faceA;

  //additional face-centred coefficient of H*mu, usually NULL
  Vector<LevelData<FluxBox>*> m_muCoef;

  // time : this class doesn't care, but m_bcPtr does
  Real m_time;

  // safety factor for VicousTensorOp
  Real m_vtopSafety;
  //  minimum number of smnooths in viscous tensor op mg relax
  int m_vtopRelaxMinIter;
   // relax tolerance for  viscous tensor op
  Real m_vtopRelaxTol;

  //lower and upper bounds on the viscosity
  Real m_muMin, m_muMax;

  Real m_scale; // scale mu,alpha by 1/m_scale

  
  // artificial drag that applies everywhere
    /* if used, need to ensure that 
       (m_artificial_drag * |u|)^m_artifical_drag_power << rho * g * h * grad(s) ~ 10^4
       when |u| ~ 10^4
    */
    Real m_artificialDragCoef, m_artificialDragPower;
    
  
  // BCHolder
  RefCountedPtr<BCFunction> m_velSolveBC;

  //coefficient storage. 
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_mu;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_lambda;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_alpha;

  void setupCoeffs();

public:

  ~IceNonlinearViscousTensor();

  virtual NonlinearViscousTensor* newNonlinearViscousTensor();

  IceNonlinearViscousTensor(const IceNonlinearViscousTensor&);

  /**
     The usual constructor.
     @param a_grids hierarchy of AMR layouts
     @param a_refRatio refinement ratio from levels l to levels l + 1
     @param a_domains domain boundaries for each level
     @param a_dxs     mesh spacing for each level
     @param a_coordSys the ice sheet geometry
     @param a_u ice  velocity
     @param a_C velocity independent part of the basal friction 
     @param a_C0 additional basal friction, such that Tb = C * fb(u)*u + C0 * u
     @param a_finestLevel finest AMR level
     @param a_constRel velocity dependent part of the viscosity, fv(A, u)
     @param a_basalFrictionRel velocity dependent part of the basal friction, fb(u)
     @param a_bc Ice thickness and velocity boundary conditions
     @param a_A velocity independent part of the viscosity, A
     @param a_faceFlowLawCoefficient another velocity independent part of the viscosity, phi
     @param a_time 
     @param a_vtopSafety [DEPRECATE?]
     

   */
  IceNonlinearViscousTensor
  (const Vector<DisjointBoxLayout>& a_grids,
   const Vector<int>& a_refRatio,
   const Vector<ProblemDomain>& a_domains,
   const Vector<RealVect>& a_dxs,
   const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   const Vector<LevelData<FArrayBox>*>& a_u,
   const Vector<LevelData<FArrayBox>*>& a_C,
   const Vector<LevelData<FArrayBox>*>& a_C0,
   const int a_finestLevel,
   const ConstitutiveRelation& a_constRel,
   const BasalFrictionRelation& a_basalFrictionRel,
   IceThicknessIBC& a_bc,
   const Vector<LevelData<FArrayBox>*>& a_A,
   const Vector<LevelData<FluxBox>*>& a_faceFlowLawCoefficient,
   Real a_time, 
   Real a_vtopSafety,
   int a_vtopRelaxMinIter,
   Real a_vtopRelaxTol,
   Real a_muMin = 0.0,
   Real a_muMax = 1.234567e+300,
   Real a_scale = 1.0,
   Real a_artificialDragCoef = 0.0, Real a_artificialDragPower=8.0);

  //implementation of NonlinearViscousTensor interface 
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_u);

  virtual void applyBC(Vector<LevelData<FArrayBox>*>& a_phi)
  {
    for (int lev =0; lev < a_phi.size(); lev++)
      {
	for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	  {
	    (*m_velSolveBC)((*a_phi[lev])[dit], m_grids[lev][dit], m_domains[lev], 
                            m_dxs[lev][0], dit(), false);
	  }
      }
  }


  virtual Vector<LevelData<FArrayBox>*>& getState() 
  {
    return m_u;
  }
 
  virtual Vector<RefCountedPtr<LevelData<FluxBox> > >& mu()
  {
    return m_mu;
  }
  virtual Vector<RefCountedPtr<LevelData<FluxBox> > >& lambda()
  {
    return m_lambda;
  }
  virtual Vector<RefCountedPtr<LevelData<FArrayBox> > >& alpha() 
  {
    return m_alpha;
  }

  virtual RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
  opFactoryPtr() {
    return m_opFactoryPtr;
  } 

  /** 
      Retrieve the basal friction coeffcient (not the effective drag coeffient, alpha() ) 
   */
  virtual Vector<LevelData<FArrayBox>*>& getC() 
  {
    return m_C;
  }

  /** 
      An overloaded setState, because beta changes with time
   */
  virtual void setState(const Vector<LevelData<FArrayBox>*>& a_u, 
			const Vector<LevelData<FArrayBox>*>& a_basalFrictionCoefficent)
  {
    m_C = a_basalFrictionCoefficent;
    setState(a_u);
  }

  /** 
      Allows an additional coefficient f of the face-centered viscosity to be specified,
      so that H*mu*(grad(u) ... ) is replaced with  f*H*mu*(grad(u) ... ).
  */
  virtual void setFaceViscCoef(const Vector<LevelData<FluxBox>*>& a_muCoef)
  {
    int nlev = min(a_muCoef.size(),m_muCoef.size());

    for (int lev=0; lev < nlev; lev++)
      {
	m_muCoef[lev] = a_muCoef[lev];
      }
  }

  /** 
      Reset to default (NULL) the coefficient set in setFaceViscCoef
  */
  virtual void unsetFaceViscCoef()
  {
    for (int lev=0; lev < m_muCoef.size(); lev++)
      {
	m_muCoef[lev] = NULL;
      }
  }

  /** 
      Provide the components of the Viscous Tensor G[ij] at cell faces, such that
      j-faces we have G[ij] for all i and so on
  */
  virtual void computeViscousTensorFace(const Vector<LevelData<FluxBox>*>& a_viscousTensor);
    
};

#include "NamespaceFooter.H"
#endif
