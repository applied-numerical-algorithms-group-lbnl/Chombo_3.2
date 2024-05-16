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
//  BasicThicknessIBC.H
// ============
//
// PhysIBC-derived class for basic (simple) IBC for ice-sheet thickness.
// this case sets the thickness at domain boundaries to a constant 
// (set by the setBoundaryThickness function)
//

#ifndef _BASICTHICKNESSIBC_H_
#define _BASICTHICKNESSIBC_H_

#include "IceThicknessIBC.H"

#include "NamespaceHeader.H"


/// Physical/domain initial and boundary conditions
/**
   
 */
class BasicThicknessIBC : public IceThicknessIBC
{
public:
  /// Constructor
  /**
   */
  BasicThicknessIBC();

  /// Destructor
  /**
   */
  virtual ~BasicThicknessIBC();

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
  /** reads info from ParmParse and sets up ice sheet geometry,
      often by calling the defineLevelSigmaCS function in the util directory
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

protected:

  // ice sheet thickness at domain boundaries
  Real          m_boundaryThickness;

  RealVect m_domainSize;

  // have bc's been set up?
  bool m_isBCsetUp;

  /// set up bc's based on parmparse inputs
  void setupBCs();

  RefCountedPtr<CompGridVTOBC> m_velBCs;

private:
  // Disallowed for all the usual reasons
  void operator=(const BasicThicknessIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  BasicThicknessIBC(const BasicThicknessIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};



#include "NamespaceFooter.H"
#endif
