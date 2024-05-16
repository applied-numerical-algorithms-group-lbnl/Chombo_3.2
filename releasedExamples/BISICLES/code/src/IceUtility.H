#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _ICEUTILITY_H_
#define _ICEUTILITY_H_
#include "Vector.H"
#include "LevelSigmaCS.H"
#include "ConstitutiveRelation.H"
#include "IceThicknessIBC.H"
#include "NamespaceHeader.H"

/// IceUtility namespace : general purpose ice sheet functions common to the 'big' classes (AmrIce, JFNKSolver, etc)
namespace IceUtility
{
  
  /// Compute a C0 such that basal traction Tb = C0 u + f(C,u)
  /** 
      Basal traction Tb = C0 u + f(C,u)
      
      f(C,u) is the conventional part of the basal traction.
      
      the C0 u term is used to add a linear drag, and for now is used
      only for 'wall drag', which allows for a crude representation
      of ice friction along rocky walls. Two Parmparse parameters, 
      wall_drag.basic and wall_drag.extra, determine the calculation.
      
      if (wall_drag.basic) is true, then C0 = (C + wall_drag.extra) u
      
      For backward compatibility, wall_drag.basic == amr.WallDrag
      and wall_drag.extra = amr.wallDragExtra
      
  */
  void computeC0(Vector<LevelData<FArrayBox>* >& a_vectC0,
		 const Vector<LevelData<FArrayBox>* >& a_vectC,
		 const Vector<DisjointBoxLayout>& a_grids,
		 const Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys,
		 const Vector<Real> a_dx, int a_finest_level);
  
  
  void addWallDrag(FArrayBox& a_drag, 
		   const BaseFab<int>& a_mask,
		   const FArrayBox& a_usrf,
		   const FArrayBox& a_thk,
		   const FArrayBox& a_topg,
		   const FArrayBox& a_beta,
		   const Real& a_extra, 
		   const RealVect& a_dx,
		   const Box& a_box);


  void addThinIceDrag(FArrayBox& a_drag, 
		      const BaseFab<int>& a_mask,
		      const FArrayBox& a_thk,
		      const Real& a_extra,
		      const Real& a_thin,
		      const Box& a_box);
  
  
  void defineRHS(Vector<LevelData<FArrayBox>* >& a_rhs,
		 const Vector<RefCountedPtr<LevelSigmaCS > >& a_CS,
		 const Vector<DisjointBoxLayout>& a_grids,
		 const Vector<RealVect>& a_dx);


  /// compute the operator L phi = a * phi + b * grad^2 (phi) 
  void applyHelmOp
  (LevelData<FArrayBox>& a_lapPhi,
   const LevelData<FArrayBox>& a_phi, 
   const Real& a_a, const Real& a_b,
   const DisjointBoxLayout& a_grids,
   const RealVect& a_dx);

  ///apply cell-centred  operator div(phi) to face-centred u. 
  ///Assumes that ghost cells have been set
  void applyDiv
  (LevelData<FArrayBox>& a_divU,
   const LevelData<FluxBox>& a_u, 
   const DisjointBoxLayout& a_grids,
   const RealVect& a_dx);

  ///apply cell-centred  operator [grad(phi)]^2 to face-centred u. 
  ///Assumes that ghost cells have been set
  void applyGradSq
  (LevelData<FArrayBox>& a_gradPhiSq,
   const LevelData<FArrayBox>& a_phi, 
   const DisjointBoxLayout& a_grids, 
   const RealVect& a_dx);
 
  ///compute face-centered us given cell-centered
  ///vector u and scalar s 
  ///Assumes that boundary values have been set.
  ///USES CENTERED SCHEME
  void computeFaceFluxCentered
  (LevelData<FluxBox>& a_us,
   const LevelData<FArrayBox>& a_u, 
   const LevelData<FArrayBox>& a_s, 
   const DisjointBoxLayout& a_grids);

  ///compute face-centered us given cell-centered
  ///vector u and scalar s 
  ///Assumes that boundary values have been set.
  ///USES  FIRST ORDER UPWIND SCHEME
  void computeFaceFluxUpwind
  (LevelData<FluxBox>& a_us,
   const LevelData<FluxBox>& a_u, 
   const LevelData<FArrayBox>& a_s, 
   const DisjointBoxLayout& a_grids);

  ///compute cell centered rate factor A from the temperature 
  void computeA
  (LevelData<FArrayBox>& a_A,
   const Vector<Real>& a_sigma,
   const LevelSigmaCS& a_coordSys,
   const RateFactor* a_rateFactor,
   const LevelData<FArrayBox>& a_internalEnergy);

  /// extrapolate face centered velocity field (usually derived by cell-to-face average) to the margins
  void extrapVelocityToMargin(LevelData<FluxBox>& a_faceVel, 
					  const LevelData<FArrayBox>& a_cellVel, 
					  const LevelSigmaCS& a_coordSys);

  ///compute face-centered velocity and thickness diffusion coefficients from cell-centered velocity
  void computeFaceVelocity
  (LevelData<FluxBox>& a_faceVelAdvection,
   LevelData<FluxBox>& a_faceVelTotal,
   LevelData<FluxBox>& a_faceDiffusivity,
   LevelData<FArrayBox>& a_cellDiffusivity,
#if BISICLES_Z == BISICLES_LAYERED
   LevelData<FluxBox>& a_layerXYFaceXYVel,
   LevelData<FArrayBox>& a_layerSFaceXYVel,
#endif
   const LevelData<FArrayBox>& a_velocity,
   const LevelSigmaCS& a_coordSys,
   const IceThicknessIBC* a_bc,
   const LevelData<FArrayBox>& a_A,
#if BISICLES_Z == BISICLES_LAYERED
   const LevelData<FArrayBox>& a_sA,
   const LevelData<FArrayBox>& a_bA,
#endif			 
   const LevelData<FArrayBox>* a_crseVelocity,
   const LevelData<FArrayBox>* a_crseDiffusivity,
   int a_nRefCrse,
   const ConstitutiveRelation* a_constitutiveRelation,
   bool a_additionalVelocity, bool a_implicitDiffusion); 


  void computeSigmaVelocity
  (LevelData<FArrayBox>& a_uSigma,
   const LevelData<FluxBox>& a_layerThicknessFlux,
   const LevelData<FArrayBox>& a_layerSFaceXYVel,
   const LevelData<FArrayBox>& a_dHdt,
   const DisjointBoxLayout a_grid,
   const LevelData<FArrayBox>& a_surfaceThicknessSource,
   const LevelData<FArrayBox>& a_basalThicknessSource,
   const Vector<Real>& a_dSigma,
   const RealVect& a_dx,
   const Real& a_dt );

  
  ///Identify regions of fast ice  and eliminate them. return the total number of cells eliminated
  int eliminateFastIce
  (Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   Vector<LevelData<FArrayBox>* >& a_vel,
   Vector<LevelData<FArrayBox>* >& a_calvedIce,
   Vector<LevelData<FArrayBox>* >& a_addedIce,
   Vector<LevelData<FArrayBox>* >& a_removedIce,
   const Vector<DisjointBoxLayout>& a_grids,
   const Vector<ProblemDomain>& a_domain,
   const Vector<int>& a_refRatio, Real a_crseDx,
   int a_finestLevel, int a_maxIter, Real a_thinIceTol, Real a_fastIceTol,
   bool a_edgeOnly, int a_verbosity = 0);

  ///Identify regions of floating ice that are remote
  ///from grounded ice and eliminate them.
  void eliminateRemoteIce
  (Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   Vector<LevelData<FArrayBox>* >& a_vel,
   Vector<LevelData<FArrayBox>* >& a_calvedIce,
   Vector<LevelData<FArrayBox>* >& a_addedIce,
   Vector<LevelData<FArrayBox>* >& a_removedIce,
   const Vector<DisjointBoxLayout>& a_grids,
   const Vector<ProblemDomain>& a_domain,
   const Vector<int>& a_refRatio, Real a_crseDx,
   int a_finestLevel, int a_maxIter, Real a_tol, 
   int a_verbosity = 0);

  /// subgrid grounding line interpolation : multiply a_u by the grounded portion of each cell
  void multiplyByGroundedFraction
  (LevelData<FArrayBox>& a_u, const LevelSigmaCS& a_coords,
   const DisjointBoxLayout& a_grids, int a_subdivision);

  /// set C = 0 in floating region, potentially using subgrid interpolation 
  void setFloatingBasalFriction
  (LevelData<FArrayBox>& a_C, const LevelSigmaCS& a_coords,
   const DisjointBoxLayout& a_grids);

}

#include "NamespaceFooter.H"
#endif
