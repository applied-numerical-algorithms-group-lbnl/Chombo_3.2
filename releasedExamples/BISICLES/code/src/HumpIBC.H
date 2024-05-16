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
//  HumpIBC.H
// ============
//
// PhysIBC-derived class for "hump" test case
// (a ellipsoidal hump of ice on a flat plain)


#ifndef _HUMPIBC_H_
#define _HUMPIBC_H_

#include "IceThicknessIBC.H"

#include "NamespaceHeader.H"


/// Physical/domain initial and boundary conditions
/**
   
 */
class HumpIBC : public IceThicknessIBC
{
public:
  /// Constructor
  /**
   */
  HumpIBC();

  /// Destructor
  /**
   */
  virtual ~HumpIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// sets parameters for the hump problem
  void setParameters(Real a_maxThickness,
                     Real a_radSqr,
                     Real a_baseElevation,
                     Real a_minThickness,
                     RealVect a_center,
		     const Real& a_seaLevel,
                     const RealVect& a_widthScale = RealVect::Unit);

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

protected:

  ///
  Real m_baseElevation;

  /// max ice sheet thickness
  Real m_maxThickness;

  // it's easier to work in terms of radius squared...
  Real m_radSqr;

  RealVect m_center;

  // width factor (0 means no variation in that direction)
  /** useful for testing 1D versions of the hump
   */
  RealVect m_widthScale;

  RealVect m_domainSize;

  // thickness of ice outside of hump (normally zero)
  Real m_minThickness;

  Real m_seaLevel;

  // have bc's been set up?
  bool m_isBCsetUp;

  // have parameters been set
  bool m_paramsSet;

  /// set up bc's based on parmparse inputs
  void setupBCs();

  RefCountedPtr<CompGridVTOBC> m_velBCs;

  // do higher-order averaged initial condition (for convergence testing)
  bool m_doHOinit;

private:
  // Disallowed for all the usual reasons
  void operator=(const HumpIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  HumpIBC(const HumpIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
