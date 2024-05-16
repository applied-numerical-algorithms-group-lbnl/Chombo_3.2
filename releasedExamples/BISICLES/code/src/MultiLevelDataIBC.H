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
// MultiLevelDataIBC.H
// ============
//
// PhysIBC-derived class which stores initial topography and thickness data
// and imposes either periodic or reflection boundary conditions
#ifndef _MULTILEVELDATAIBC_H_
#define _MULTILEVELDATAIBC_H_

#include "IceThicknessIBC.H"
#include "NamespaceHeader.H"

/// class MultiLevelDataIBC stores initial topography and thickness data on multiple levels and imposes either periodic or reflection boundary conditions
/**
   MultiLevelDataIBC is similar to LevelDataIBC (and might eventually subsume
   it entirely): it imposes either reflection or periodic boundary conditions 
   on the velocity field, and provides an initial thickness and topgraphy 
   derived from  discrete data (typically loaded from a file). Unlike LevelDataIBC, 
   it stores thickness and topography data in an AMR hierarchy, that is
   as Vector<LevelData<FArrayBox>*>. The "source" level domains that MultiLevelDataIBC 
   stores need not be the same as the "destination" level domains stored bu  (say) 
   AmrIce, but do need to be compatible with them such that the usual 
   coarsening and refining operations can be used to construct the destination
   data from the source.

*/
class MultiLevelDataIBC : public IceThicknessIBC
{
public:
  MultiLevelDataIBC(const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_thck,
		    const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_topg, 
		    const RealVect& a_dxCrse, const Vector<int> & a_refRatio);
		    
  virtual ~MultiLevelDataIBC();

  
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  virtual IceThicknessIBC* new_thicknessIBC();

  virtual void initialize(LevelData<FArrayBox>& a_U);

  ///set up the initial ice sheet geometry
  /** Thickness and topography data are determined, in order of preference, by
      1. copying from one of the stored levels, if possible, or
      2. coarse averaging from a finer stored level, if possible but 1. is not, or
      3. interpolating from a_crseCoords
      4. interpolating from a coarser stored level
   */    
  virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio);
   
  ///modify the ice sheet geometry after regrid
  bool regridIceGeometry(LevelSigmaCS& a_coords,
			 const RealVect& a_dx,
			 const RealVect& a_domainSize,
			 const Real& a_time,
			 const LevelSigmaCS* a_crseCoords,
			 const int a_refRatio);

  /// Set boundary thickness fluxes
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
                     const Real&      a_time)
  {;}

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time)
  {;}

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
  
  // Disallowed for all the usual reasons
  void operator=(const MultiLevelDataIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  MultiLevelDataIBC(const MultiLevelDataIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

private:
  RefCountedPtr<CompGridVTOBC> m_velBCs;
  bool m_isBCsetUp, m_verbose;
  //levels of thickness data
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_thck;
  //levels of topography data
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_topg;
  // count of levels, should be == m_thck.size() == m_topg.size()
  int m_nLevel;
  //spacing of the coarsest mesh
  const RealVect m_dxCrse;
  //refinement ratios of the usual kind
  const Vector<int>  m_refRatio;

};
#include "NamespaceFooter.H"
#endif
