#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _L1L2CONSTITUTIVERELATION_H_
#define _L1L2CONSTITUTIVERELATION_H_

#include "ConstitutiveRelation.H"
#include "IceThicknessIBC.H"
#include "NamespaceHeader.H"


#if BISICLES_Z == BISICLES_LAYERED
/// Implemnets ConstitutiveRelation to give the effective viscosity in the L1L2 stress model.
/** 
    The L1L2 stress model is described by Schoof and Hindmarsh 2010.

    Like GlensFlowRelation, this class computes \f$\mu\f$ such that (stress =
    \f$\mu\f$*strainrate ) assuming that \f$\mu = \mu(\nabla(\vec{u}),A(T))\f$.
    But \f$\mu\f$\ is now found by solving a nonlinear problem 
    \f$1 - 2 \mu A (4 \mu ^2 \epsilon^2 + \tau_d^2)^((n-1)/2) = 0 \f$
    numerically for mu at each point and integrating vertically.
    \f$\tau_d = \rho g h \nabla{s} \f$ is the gravitaional driving stress,
    and its appearance means that vertical shear strains computed 
    accoding to the shallow ice approximation are included in the 
    effective viscosity (so that ice softens in the region of steep slopes).
    Although an analytic solution does exist when n = 3, 
    it is expensive/unreliable to compute
*/
class L1L2ConstitutiveRelation : public ConstitutiveRelation
{

  GlensFlowRelation glensFlowRelation;

public:
  ///
  L1L2ConstitutiveRelation() 
    :m_solverTol(1.0e-6),
     m_additionalVelocitySIAGradSLimit(1.0e-2),
     m_effectiveViscositySIAGradSLimit(1.0e+10),
     m_additionalVelocitySIAOnly(true),
     m_startFromAnalyticMu(false),
     m_layerCoarsening(0)
  {;}
  ///
  
  void parseParameters();

  virtual ~L1L2ConstitutiveRelation() {;}

  GlensFlowRelation* getGlensFlowRelationPtr()
  { return &glensFlowRelation;}
  GlensFlowRelation& getGlensFlowRelation()
  { return glensFlowRelation;}

  virtual Real power() const 
  { return glensFlowRelation.power(); }

  /// creates a new copy of this ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;

  /// computes cell-centered \f$\mu_{AS}\f$ based on the cell-centered velocity
  /**
     \param a_mu -- mu_{AS} based on the local velocity field.
     \param a_vel -- Cell-centered velocity field.
     \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     \param a_nRefCrse -- refinement ratio to next coarser level     
     \param a_A: Cell-centered flow law coefficient (Glenn's A) field
     \param a_coordSys: LevelSigmaCS object containing the geometry.
     \param a_ghostVect: how the boxes on which we want to compute mu compare
                         with the DisjointBoxLayout. (can be negative)
  */
  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                         const LevelData<FArrayBox>* a_crseVelPtr,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const;
  
 
  /// Compute a cell centred bulk dissipation \f$\Phi = 4\mu\epsilon^2\f$ 
  /// (heat source) at the cell centres. This ought to have the same dependencies 
  /// as cell centred mu
  /**
     a_dissipation -- \f$\Phi\f$ based on the local velocity field.
     a_vel -- Cell-centered velocity field.
     a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
     (NULL if no coarser level)
     a_nRefCrse -- refinement ratio to next coarser level
     a_A: Cell-centered flow law coefficient (Glenn's A) field
     a_coordSys:  SigmaCS object containing the geometry of this patch.
     a_box: cell-centered box over which to do this computation
  **/
  virtual void computeDissipation(LevelData<FArrayBox>& a_dissipation,
				  const LevelData<FArrayBox>& a_vel, 
				  const LevelData<FArrayBox>* a_crseVel,
				  int a_nRefCrse,
				  const LevelData<FArrayBox>& a_A,
				  const LevelSigmaCS& a_coordSys,
				  const ProblemDomain& a_domain,
				  const IntVect& a_ghostVect = IntVect::Zero) const;
  
  /// computes face-centered \f$\mu_{AS}\f$ based on cell-centered velocity
  /** \param a_mu: face-centered mu_{AS} based on the local velocity field.
      \param a_vel: Cell-centered velocity field.
      \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                          (NULL if no coarser level)
      \param a_nRefCrse -- refinement ratio to next coarser level     
      \param a_A: Cell-centered flow law coefficient (Glenn's A) field
      \param a_coordSys: LevelSigmaCS object containing the geometry.
      \param a_ghostVect: how the boxes on which we want to compute mu compare
                          with the DisjointBoxLayout. (can be negative)
  */  
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                             const LevelData<FArrayBox>* a_crseVelPtr,
                             int nRefCrse,
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const;



  /// modify 2D thickness advection velocity and diffusion coefficients, also 3D ice velocity
  /**
     \param a_cellVel : cell-centered basal velocity field \f$u(z=b)\f$
     \param a_crseVelPtr :  coarser level cell-centered velocity field (NULL if no coarser level)
     \param a_crseVelPtr :  coarser level cell-centered diffusivity field (NULL if no coarser level)
     \param a_nRefCrse : refinement ratio to next coarser level 
     \param a_coordSys: LevelSigmaCS object containing the geometry.
     \param a_grids : DisjointBoxLayout of grids for the level
     \param a_domain : ProblemDomain for the level
     \param a_A: (x,y,sigma) cell-centered flow law coefficient (Glenn's A) field
     \param a_sA: (x,y) cell-centered A field at the ice surface (\f$\sigma = 0\f$)
     \param a_bA: (x,y) cell-centered A field at the ice base ( \f$\sigma = 1\f$)
     \param faceVelAdvection : face-centered vertically averaged advection velocity. Should contain the velocity at the ice base on entry. 
     \param a_faceVelTotal : face-centered vertically averaged total velocity (including diffusion). Should contain the velocity at the ice base on entry
     \param a_faceDiffusivity: face-centered vertically thickness diffusivity. 
     \parama_cellDiffusivity: cell-centered vertically thickness diffusivity.
     \param a_layerXYFaceXYVel : (x,y) face-centered, sigma cell-centered 3D x and y velocity
     \param a_layerSFaceXYVel : (x,y) cell-centered, sigma face centered 3D x and y velocity
   */
  virtual void modifyTransportCoefficients
  (const LevelData<FArrayBox>& a_cellVel,
   const LevelData<FArrayBox>* a_crseVelPtr,
   const LevelData<FArrayBox>* a_crseDiffusivityPtr,
   int a_nRefCrse,
   const LevelSigmaCS& a_coordSys,
   const DisjointBoxLayout& a_grids,
   const ProblemDomain& a_domain,
   const LevelData<FArrayBox>& a_A,
   const LevelData<FArrayBox>& a_sA,
   const LevelData<FArrayBox>& a_bA,
   LevelData<FluxBox>& a_faceVelAdvection,
   LevelData<FluxBox>& a_faceVelTotal,
   LevelData<FluxBox>& a_faceDiffusivity,
   LevelData<FArrayBox>& a_cellDiffusivity,
   LevelData<FluxBox>& a_layerXYFaceXYVel,
   LevelData<FArrayBox>& a_layerSFaceXYVel) const;
    

  /// sets tolerance for secant solve
  void solverTolerance(Real a_solverTol) 
  { m_solverTol = a_solverTol; }


private:
  ///compute mu for given sigma at a_box cell centres,
  /// plus velocity gradient & rate-of-strain invariant
  /// at base
  ///a_mu has one component for each layer
  virtual void computeMuZ(LevelData<FArrayBox>& a_mu,
			  LevelData<FArrayBox>& a_gradVel,
			  LevelData<FArrayBox>& a_epsSqr,
			  const Vector<Real>& a_sigma,
			  const LevelData<FArrayBox>& a_vel, const Real& a_scale,
                          const LevelData<FArrayBox>* a_crseVelPtr,
                          int a_nRefCrse,
			  const LevelData<FArrayBox>& a_A,
			  const LevelSigmaCS& a_coordSys,
			  const ProblemDomain& a_domain,
			  const IntVect& a_ghostVect = IntVect::Unit) const;
  
  

  //compute mu for given sigma at a_faceDir cell faces
  /// plus velocity gradient & rate-of-strain invariant
  /// at base
  ///a_mu has one component for each layer
  virtual void computeFaceMuZ(LevelData<FluxBox>& a_mu,
			      LevelData<FluxBox>& a_gradVel,
			      LevelData<FluxBox>& a_epsSqr,
			      const Vector<Real>& a_sigma,
			      LevelData<FArrayBox>& a_vel, const Real& a_scale,
                              const LevelData<FArrayBox>* a_crseVelPtr,
                              int a_nRefCrse,
			      const LevelData<FluxBox>& a_A,
			      const LevelSigmaCS& a_coordSys,
			      const ProblemDomain& a_domain,
			      const IntVect& a_ghostVect) const;
  

  /// calculate mu for given sigma at a cell faces or centers
  /// plus velocity gradient & rate-of-strain invariant
  /// at base
  ///a_mu has one component for each layer
  virtual void computeEitherMuZ(FArrayBox& a_mu,
				const Vector<Real>& a_sigma,
				const FArrayBox& a_grads,
				const FArrayBox& a_epsSqr,
				const FArrayBox& a_A,const Real& a_scale,
				const FArrayBox& a_H,
				const Real& a_rhog,
				const Box& a_box,
				int a_layerCoarsen) const;

  /// computes the vertically averaged, face-centred, face direction
  /// component of the velocity field
  /** 
      \param a_vel: the velocity field at the ice base at cell centers
      \param a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                          (NULL if no coarser level)
      \param a_nRefCrse -- refinement ratio to next coarser level     \
      \param a_faceVel: the vertically averaged velocity field at cell faces
      \param a_faceDir : face direction
      \param a_cellBox : box of cell centres: computes face flux for each of these
      The L1L2 constutive relation implies a z-varying
      velocity field, which is constructed from a 2D velocity
      field computed at the ice base. 
  */ 
  

  virtual void computeFaceFluxVelocity(LevelData<FArrayBox>& a_vel,
                                       const LevelData<FArrayBox>* a_crseVelPtr,
                                       int a_nRefCrse,
                                       const LevelData<FArrayBox>& a_A,
                                       const LevelData<FArrayBox>& a_thickness,
                                       const RealVect& a_dx,
                                       LevelData<FluxBox>& a_fluxVel,
  				       LevelData<FluxBox>& a_layerXYFaceXYVel,
  				       LevelData<FArrayBox>& a_layerSFaceXYVel,
                                       const LevelSigmaCS& a_coordSys,
  				       const ProblemDomain& a_domain,
                                       const IntVect& a_cellGhost = IntVect::Zero) const;

 
  // tolerance for the secant solve when computing mu
  Real m_solverTol;
  // if > 0, limit the magnitude of grad(s) when computing the SIA part of u(z) - u(z=base)
  Real m_additionalVelocitySIAGradSLimit; 
  // if > 0, limit the magnitude of grad(s) when computing the effective viscosity
  Real m_effectiveViscositySIAGradSLimit;
  // neglect non SIA contribution to u(z) - u(z=base)  
  bool m_additionalVelocitySIAOnly; 

  // for n = 3, there is an analytic solution to f(mu) = 0; sometimes
  // it is preferable to start from that value
  bool m_startFromAnalyticMu;

  // Set m_layerCoarsening > 0 to skip layers in expensive calculations and interpolate instead
  int m_layerCoarsening;
  

};
#endif 
#include "NamespaceFooter.H"
#endif
