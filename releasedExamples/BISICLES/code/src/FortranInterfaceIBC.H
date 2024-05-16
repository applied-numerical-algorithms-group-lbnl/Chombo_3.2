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
//  FortranInterfaceIBC.H
// ============
//
// PhysIBC-derived class for interfacing with fortran -- IBC fields like 
// initial ice thickness and topography are passed in as fortran arrays, which 
// are then used to set values in the BISICLES code
//

#ifndef _FORTRANTHICKNESSIBC_H_
#define _FORTRANTHICKNESSIBC_H_

#include "IceThicknessIBC.H"
#include "RefCountedPtr.H"
#include "NamespaceHeader.H"

class AmrIce;

enum FIBC_types {FIBC_Dirichlet = 0,
                 FIBC_Neumann,
                 FIBC_iceDivide,
                 FIBC_iceStreamXBC,
                 NUM_FIBC_TYPES};


/// Physical/domain initial and boundary conditions
/**
   
 */
class FortranInterfaceIBC : public IceThicknessIBC
{

 
public:
  /// Constructor
  /**
   */
  FortranInterfaceIBC();

  /// Destructor
  /**
   */
  virtual ~FortranInterfaceIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// static function to set reference fabs
  /** if a_nodal is false, then data is cell-centered and both a_fab and 
      a_ccFab are set to be aliases which point to a_data_ptr. If a_nodal is
      true, then a_fab is set to be an aliase which points to a_data_ptr, 
      while a_ccFab is allocated and filled with cell-averaged data.
  */
  static void  setFAB(Real* a_data_ptr,
                      const int* a_dimInfo,
                      const int* boxlo, const int* boxhi, 
                      const Real* a_dew, const Real* a_dns,
                      const IntVect& a_offset,
                      const IntVect& a_nGhost,
                      FArrayBox& a_fab,
                      FArrayBox& a_ccFab,
                      const bool a_nodal,
                      const bool a_verbose = false);

  /// sets reference fabs for (horizontal) velocity
  /** if a_nodal is false, then data is cell-centered and both a_fab and 
      a_ccFab are set to be aliases which point to a_data_ptr. If a_nodal is
      true, then a_fab is set to be an aliase which points to a_data_ptr, 
      while a_ccFab is allocated and filled with cell-averaged data.
  */
  virtual void  setVelFAB(Real* a_uVelPtr, Real* a_vVelPtr,
                          const int* a_dimInfo,
                          const int* boxlo, const int* boxhi, 
                          const Real* a_dew, const Real* a_dns,
                          const IntVect& a_offset,
                          const IntVect& a_nGhost,
                          FArrayBox& a_uFab,
                          FArrayBox& a_vFab,
                          FArrayBox& a_uccFab,
                          FArrayBox& a_vccFab,
                          const bool a_nodal = false,
                          const bool a_fillData = false);


  /// nGhost is the number of ghost cells in in the input data (which are ignored)
  virtual void setThickness(Real* a_data_ptr,
			    const int* a_dimInfo,
                            const int* boxlo, const int* boxhi, 
			    const Real* a_dew, const Real* a_dns,
			    /// amount to offset indices (due to C-Fortran numbering, for example)
			    const IntVect& a_offset, 
                            const IntVect& a_nGhost = IntVect::Zero, 
			    const bool a_nodal = false);

  virtual void setTopography(Real* a_data_ptr,
			     const int* a_dimInfo,
                             const int* boxlo, const int* boxhi, 
			     const Real* a_dew, const Real* a_dns,
			     /// amount to offset indices (due to C-Fortran numbering, for example)
			     const IntVect& a_offset,
                             const IntVect& a_nGhost = IntVect::Zero, 
			     const bool a_nodal = false);

  /// sets storage for upper surface. 
  /** Currently only needed in order to return BISICLES-computed upper 
      surface to CISM */
  virtual void setSurface(Real* a_data_ptr,
                          const int* a_dimInfo,
                          const int* boxlo, const int* boxhi, 
                          const Real* a_dew, const Real* a_dns,
                          /// amount to offset indices (due to C-Fortran numbering, for example)
                          const IntVect& a_offset,
                          const IntVect& a_nGhost = IntVect::Zero, 
                          const bool a_nodal = false);
  

  /// regions where we artificially set thickness to zero
  /** this is done in setThickness, for lack of a better place, so 
      this needs to be set before setThickness is called. 
      a_clearRegions defines logically-rectangular regions where 
      the thickness is artificially set to zero. These regions are 
      defined relative to the original input data (i.e. before any shifting 
      due to ghost cells)
  */
  virtual void setThicknessClearRegions(const Vector<Box>& a_clearRegions);

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

  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
   */
  virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio) ;

  /// update the topgraphy in LevelSigmaCS but *not* the thickness
  virtual bool regridIceGeometry(LevelSigmaCS& a_coords,
				 const RealVect& a_dx,
				 const RealVect& a_domainSize,
				 const Real& a_time, 
				 const LevelSigmaCS* a_crseCoords,
				 const int a_refRatio) ;

  /// flatten thickness and basal topography back to input FArrayBoxes
  /** 
      re-store in m_inputThickness and m_inputTopography -- not const 
      because we modify the values in the data holders
  */
  virtual void flattenIceGeometry(const Vector<RefCountedPtr<LevelSigmaCS> > & a_amrGeometry);

  /// flatten an arbitrary dataset back to input FArrayBoxes
  virtual void flattenData(Real* a_data_ptr,
                           const int* a_dimInfo,
                           const int* a_boxlo, const int* a_boxhi,
                           const Real* a_dew, const Real* a_dns,
                           const IntVect& a_offset,
                           const Vector<LevelData<FArrayBox>* >& a_amrData,
                           const Vector<int>& a_vectRefRatio,
                           const Vector<Real>& a_amrDx,
                           int a_srcComp, 
                           int a_destComp,
                           int a_nComp,
                           const IntVect& a_nGhost,
                           const bool a_nodal);
               

  /// flatten velocity field to the 3D velocity data holder
  void flattenVelocity(Real* a_uVelPtr, Real* a_vVelPtr,
                       const int* a_dimInfo,
                       const int* a_boxlo, const int* a_boxhi,
                       const Real* a_dew, const Real* a_dns,
                       const IntVect& a_offset,
                       const Vector<LevelData<FArrayBox>* >& a_amrVel,
                       const Vector<int>& a_vectRefRatio,
                       const Vector<Real>& a_amrDx,
                       const IntVect& a_nGhost,
                       const bool a_nodal);
  

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

  void setVerbose(bool a_verbose) {m_verbose = a_verbose;}

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
  
  /// set grids using Boxes passed in from CISM
  /** creates a DisjointBoxLayout using the grid boxes and 
      processor distribution used by CISM. */
  static void setGrids(DisjointBoxLayout& a_grids,
                       const Box& a_gridBox,
                       const ProblemDomain& a_domain,
                       bool a_verbose);

  /// has setGrids been called?
  bool gridsSet() const {return m_gridsSet;}


  /// utility function to fill in holes in topography
  /** looks for isolated values of holeVal and replaces then
      with average of neighbors */
  virtual void fillTopographyHoles(Real holeVal);

  virtual void checkOK() const;

  // accessor
  const RefCountedPtr<LevelData<FArrayBox> >& inputThicknessPtr() const
  {return  m_inputThicknessLDF;}

protected:

  // ice sheet thickness at domain boundaries
  Real          m_boundaryThickness;

  RealVect m_domainSize;

  bool m_verbose;

  // have bc's been set up?
  bool m_isBCsetUp;

  /// set up bc's based on parmparse inputs
  void setupBCs();

  bool m_extrapBoundary; //extrapolate the boundary geoemtry? or reflect it?

  RefCountedPtr<CompGridVTOBC> m_velBCs;

  DisjointBoxLayout m_grids;
  bool m_gridsSet;

  // ice sheet thickness from glimmer
  /// alias to data array passed in from glimmer-CISM
  FArrayBox m_inputThickness;
  /// true if m_inputThickness is nodal
  bool m_nodalThickness;
  /// cell-centered input thickness 
  /** if nodalThickness is false, this is identical to m_inputThickness. 
      If nodalThickness is true, then this contains the cell-averaged version 
      of the inputThickness
  */
  FArrayBox m_ccInputThickness;
      
  RealVect m_inputThicknessDx;
  IntVect m_thicknessGhost; 
  /// distributed version of thickness
  RefCountedPtr<LevelData<FArrayBox> > m_inputThicknessLDF;

  /// regions where we will artificially set ice thickness to zero
  Vector<Box> m_thicknessClearRegions;

  /// alias to data array passed in from glimmer-CISM
  FArrayBox m_inputTopography;
  // true if m_inputTopography is nodal
  bool m_nodalTopography;
  /// cell-centered input topography
  /**
     if nodalTopography is false, this is identical to m_inputTopography.
     if nodalTopography is true, this contains the cell-averaged version
     of the inputTopography
  */
  FArrayBox m_ccInputTopography;

  RealVect m_inputTopographyDx;
  IntVect m_topographyGhost;
  /// distributed version of topography
  RefCountedPtr<LevelData<FArrayBox> > m_inputTopographyLDF;


  /// alias to data array passed in from glimmer-CISM
  FArrayBox m_inputSurface;
  // true if m_inputSurface is nodal
  bool m_nodalSurface;
  /// cell-centered input upper ice surface
  /**
     if nodalSurface is false, this is identical to m_inputSurface;
     if nodalSurface is true, this contains the cell-averaged version
     of the inputSurface
  */
  FArrayBox m_ccInputSurface;

  RealVect m_inputSurfaceDx;
  IntVect m_surfaceGhost;
  /// distributed version of topography
  RefCountedPtr<LevelData<FArrayBox> > m_inputSurfaceLDF;

private:
  // Disallowed for all the usual reasons
  void operator=(const FortranInterfaceIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  FortranInterfaceIBC(const FortranInterfaceIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
