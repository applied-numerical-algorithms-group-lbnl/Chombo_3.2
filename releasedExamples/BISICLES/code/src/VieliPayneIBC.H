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
//  PayneVieliIBC.H
// ============
//
// PhysIBC-derived class for grounding-line testing
//

#ifndef _VIELIPAYNEIBC_H_
#define _VIELIPAYNEIBC_H_

#include "IceThicknessIBC.H"

#include "NamespaceHeader.H"


/// BCFunction-derived class to handle velocity-solve BC's 
/** in particular, this BCFunction contains the 
    Vector<LevelSigmaCS* > which we need to define the marine boundary 
    condition */
class VieliPayneBCFunction : public IceVelBCFunction
{
public:
  VieliPayneBCFunction() {;}
  
  VieliPayneBCFunction(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                       Vector<ProblemDomain>& a_vectDomain);
  
  
  ~VieliPayneBCFunction(); 

  virtual void define(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                      Vector<ProblemDomain>& a_vectDomain);
  
  
  /// this function throws an error since we need the DataIndex
  /** since we need the DataIndex to do the marine boundary condition, 
      this function throws a MayDay::Error...
  */
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous);

  /// compute ghost-cell values for the velocity solve
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          const DataIndex&     a_index,
                          bool                 a_homogeneous);


protected:
  
  Vector<RefCountedPtr<LevelSigmaCS > > m_vectCS;
  Vector<ProblemDomain> m_vectDomain;
};


/// Physical/domain initial and boundary conditions
/**
 */
class VieliPayneIBC : public IceThicknessIBC
{
public:
  /// Constructor
  /**
   */
  VieliPayneIBC();

  /// Destructor
  /**
   */
  virtual ~VieliPayneIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// set parameters
  void setParameters(const Real& a_thickness,
                     const RealVect& a_slope,
                     const Real& a_originElevation,
                     const Real& a_seaLevel);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new IceThicknessIBC object with 
     m_isDefined = false (i.e., its define() must be called before it
     is used). 
   */
  virtual IceThicknessIBC* new_thicknessIBC();

#if 0
  /// set boundary ice-sheet thickness
  void setBoundaryThickness(Real a_boundaryThickness)
  { m_boundaryThickness = a_boundaryThickness; }
#endif

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

  /// fill ghost cells on velocity
  virtual void velocityGhostBC(LevelData<FArrayBox>& a_velocity,
                               LevelSigmaCS& a_coords,
                               const ProblemDomain& a_domain,
                               Real a_time);

  /// if appropriate, modify velocity solve RHS in a problem-dependent way. 
  /** 
      In this case, what's done is to set the boundary thickness for 
      use in computing the inhomogeneous form of the velocity BC.
      (actual RHS is left unmodified, which suggests that maybe this
      function should really be called something else...)
   */
  virtual void modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
                                 LevelSigmaCS& a_coords,
                                 const ProblemDomain& a_domain,
                                 Real a_time, Real a_dt);


  /// set non-periodic ghost cells for surface height z_s. 
  /** use linear extrapolation
   */
  virtual void setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                   LevelSigmaCS& a_coords,
                                   const ProblemDomain& a_domain,
                                   const RealVect& a_dx, 
                                   Real a_time, Real a_dt);

  /// set AMR grid hierarchy (for the BC's which need this)
  virtual void setGridHierarchy(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                                Vector<ProblemDomain>& a_vectDomain);
    

  /// ice thickness at edge of domain for use in computing velocity Bc's
  /** this is an initial hack intended to get things up and 
      running, and assumes that thickness is constant along the
      front... This is public so that the BCFunc can access it.
  */
  static RealVect s_edgeThickness;

protected:

  RealVect m_slope;

  ///
  Real m_originElevation;

  /// ice sheet thickness (initially constant)
  Real m_thickness;

  RealVect m_domainSize;

  Real m_seaLevel;

  // have bc's been set up?
  bool m_isBCsetUp;

  // have parameters been set
  bool m_paramsSet;

  /// set up bc's based on parmparse inputs
  void setupBCs();

  RefCountedPtr<CompGridVTOBC> m_velBCs;
  RefCountedPtr<VieliPayneBCFunction> m_BCfunction;


private:
  // Disallowed for all the usual reasons
  void operator=(const VieliPayneIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  VieliPayneIBC(const VieliPayneIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};




#include "NamespaceFooter.H"
#endif
