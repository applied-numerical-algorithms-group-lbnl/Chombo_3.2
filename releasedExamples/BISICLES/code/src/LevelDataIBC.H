#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  LevelDataIBC.H
// ============
//
// PhysIBC-derived class which stores initial topography and thickness data
// and imposes either periodic or reflection boundary conditions

#ifndef _LEVELDATAIBC_H_
#define _LEVELDATAIBC_H_

#include "IceThicknessIBC.H"

#include "NamespaceHeader.H"

/// class LevelDataIBC stores initial topography and thickness data on a single levels and imposes either periodic or reflection boundary conditions
/**
   LevelDataIBC constructs initial thickness and topography 
   data from uniformly gridded DEMs such as ALBMAP (usually loaded from 
   a file). It also imposes reflection or periodic boundary conditions on the 
   velocity field. The "source" level domain need not coincide with any
   of the "destination" level domains stored by  (say) 
   AmrIce, but do need to be compatible with them such that the usual 
   coarsening and refining operations can be used to construct the destination
   data from the source.
 */
class LevelDataIBC : public IceThicknessIBC
{
public:
  /// Constructor
  /**
   */
  LevelDataIBC(RefCountedPtr<LevelData<FArrayBox> > a_thck, 
	       RefCountedPtr<LevelData<FArrayBox> > a_topg,
	       const RealVect& a_dx,
               Real a_defaultThickness = 0.0,
               Real a_defaultTopography = -10000,
               bool a_setDefaultValues = false);

  LevelDataIBC();

  /// Destructor
  /**
   */
  virtual ~LevelDataIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new IceThicknessIBC object with 
     m_isDefined = false (i.e., its define() must be called before it
     is used). 
   */
  virtual IceThicknessIBC* new_thicknessIBC();


  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U);
 
  ///modify the ice sheet geometry after regrid
  bool regridIceGeometry(LevelSigmaCS& a_coords,
			 const RealVect& a_dx,
			 const RealVect& a_domainSize,
			 const Real& a_time,
			 const LevelSigmaCS* a_crseCoords,
			 const int a_refRatio);

  ///set up the initial ice sheet geometry
  /** Thickness and topography data are determined, in order of preference, by
      1. copying from one of the stored levels, if possible, or
      2. coarse averaging from a stored level, if possible but 1. is not, or
      3. interpolating from a_crseCoords
      4. interpolating from a coarser stored level
   */ 
  virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time,
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio);


  /// Set boundary fluxes
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  virtual
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  /// return boundary condition for Ice velocity solve
  /** 
   */
  virtual RefCountedPtr<CompGridVTOBC> velocitySolveBC();

  /// if appropriate, modify velocity solve RHS in a problem-dependent way. 
  virtual void modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
                                 LevelSigmaCS& a_coords,
                                 const ProblemDomain& a_domain,
                                 Real a_time, Real a_dt)
  {;}


  /// set non-periodic ghost cells for surface height z_s. 
  virtual void setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                   LevelSigmaCS& a_coords,
                                   const ProblemDomain& a_domain,
                                   const RealVect& a_dx, 
                                   Real a_time, Real a_dt);
  /// set non-periodic geometry (thickness, topography) ghost cells
  virtual void setGeometryBCs(LevelSigmaCS& a_coords,
			      const ProblemDomain& a_domain,
			      const RealVect& a_dx,
			      Real a_time, Real a_dt);

private:

  RefCountedPtr<LevelData<FArrayBox> > m_thck;
  RefCountedPtr<LevelData<FArrayBox> > m_topg;
  RealVect m_dx;
  RefCountedPtr<CompGridVTOBC> m_velBCs;
  bool m_isBCsetUp;
  bool m_verbose;
  /// default value for thickness
  Real m_default_thickness;
  // default value for topogaphy
  Real m_default_topography;
  /// if true, set thickness and topography to default values as a hedge
  /** Useful if we're filling a grown domain from a smaller dataset
   */
  bool m_set_default_values;
  void setupBCs();

private:
  // Disallowed for all the usual reasons
  void operator=(const LevelDataIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelDataIBC(const LevelDataIBC& a_input)
    :m_thck(a_input.m_thck),m_topg(a_input.m_topg),m_dx(a_input.m_dx)
  {
    MayDay::Error("invalid operator");
  }

  void groundCoarse(LevelSigmaCS& a_coords, int a_refRatio);

};

#include "NamespaceFooter.H"
#endif
