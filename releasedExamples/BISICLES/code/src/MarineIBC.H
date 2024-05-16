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
//  MarineIBC.H
// ============
//
// PhysIBC-derived class for grounding-line testing
// An alternative to Dan's VieliPayneIBC class which
// imposes the boundary condition through a source 
// term s = 1/Dx * 0.5 *  ri/rw * (ri - rw) * g * H^2  
// (rather than setting grad(u).n)

#ifndef _MARINEIBC_H_
#define _MARINEIBC_H_

#include "IceThicknessIBC.H"
#include "RealFunction.H"
#include "BCFunc.H"
#include "NamespaceHeader.H"


class MISOMIPBedrockElevation : public RealFunction<RealVect>
{
  //shallow parameters
  Real B0; 
  Real B2;
  Real B4;
  Real B6;

public:
  MISOMIPBedrockElevation() :  B0(-150.0), B2(-728.8), B4(343.91), B6(-50.57)
  {;}

  Real operator()(const RealVect& x)
  {
    Real y = x[1] - 40e3;
    Real topography = Bx(x[0]) + By(y,4.0e+3,5.0e+2,24.0e+3);
    topography = max(topography, -720.0);
    return topography;
  }

protected:
  Real Bx(const Real x)
  {
    Real xx = x/300.0e3;
    Real xx2 = xx*xx;
    Real xx4 = xx2*xx2;
    Real xx6 = xx4*xx2;
    return B0 + B2*xx2 + B4*xx4 + B6*xx6;
  }

  Real By(const Real y, const Real fc, const Real dc, const Real wc)
  {
    return dc/(1.0+exp(-2.0*(y-wc)/fc)) + dc/(1.0+exp(2.0*(y+wc)/fc));
  }
        


};


class SchoofBedrockElevation : public RealFunction<RealVect>
{
  Real m_lengthScale; 
  Real m_originElevation;
  Real m_coeff2, m_coeff4, m_coeff6;
public:

  SchoofBedrockElevation(Real a_lengthScale, Real a_originElevation, 
                         Real a_coeff2, Real a_coeff4, Real a_coeff6)
  : m_lengthScale(a_lengthScale),m_originElevation(a_originElevation),
    m_coeff2(a_coeff2),m_coeff4(a_coeff4),m_coeff6(a_coeff6)
  {;}

  Real operator()(const RealVect& x)
  {
    Real x2 = x[0]*x[0]/ (0.25*m_lengthScale*m_lengthScale);
    Real x4 = x2 * x2;
    Real x6 = x4 * x2;

    return m_originElevation + m_coeff2 * x2 + m_coeff4 * x4 + m_coeff6*x6;

  }

};

class KatzBedrockElevation : public RealFunction<RealVect>
{
  Real m_domainLength, m_domainWidth;
  // channel width, background slope parameters
  Real m_originElevation,  m_alpha, m_sigma, m_s; 
  // profile in the channel
  SchoofBedrockElevation m_channelProfile;  

public:

  KatzBedrockElevation(Real a_domainLength, Real a_domainWidth, 
		       Real a_originElevation, Real a_alpha, Real a_sigma,
		       Real a_lengthScaleFactor,
		       Real a_coeff2, Real a_coeff4, Real a_coeff6)
    : m_domainLength(a_domainLength),  m_domainWidth(a_domainWidth),
      m_originElevation(a_originElevation), m_alpha(a_alpha), m_sigma(a_sigma),
      m_channelProfile(a_lengthScaleFactor*a_domainLength,  
		       a_originElevation, a_coeff2, 
		       a_coeff4, a_coeff6)
  {
    m_s =1.0 / ( 2.0 * m_sigma * m_sigma);
  }

  Real operator()(const RealVect& x)
  {

    Real y = x[1]-m_domainWidth/2.0;
    Real G = std::exp(-y*y * m_s);

    return (m_alpha*x[0] + m_originElevation)*(1.0 - G)
      + m_channelProfile(x) * G;

  }


};


/// Physical/domain initial and boundary conditions
/**
   
 */
class MarineIBC : public IceThicknessIBC
{
public:
  /// Constructor
  /**
   */
  MarineIBC();

  /// Destructor
  /**
   */
  virtual ~MarineIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  void setParameters(RefCountedPtr<RealFunction<RealVect > > a_thicknessFunction,
		     RefCountedPtr<RealFunction<RealVect > > a_bedrockFunction,
		     const Real& a_seaLevel);


  // Vector of bedrockFunctions allows for simple composition to create geometries
  void setParameters(RefCountedPtr<RealFunction<RealVect > > a_thicknessFunction,
		     Vector<RefCountedPtr<RealFunction<RealVect > > > a_bedrockFunction,
		     const Real& a_seaLevel);

  // /// set parameters : inclined plane geomety
  // void setParameters(const Real& a_thickness,
  // 		     const Real& a_originElevation,
  // 		     const RealVect& a_slope,
  //                    const Real& a_seaLevel);

  // //set parameters : katz geometry
  // void setParametersKatz(const Real& a_thickness,
  // 			 const Real& a_domainLength, 
  // 			 const Real& a_domainWidth, 
  // 			 const Real& a_originElevation, 
  // 			 const Real& a_alpha, 
  // 			 const Real& a_sigma,
  // 			 const Real& a_coeff2, 
  // 			 const Real& a_coeff4, 
  // 			 const Real& a_coeff6,
  // 			 const Real& a_seaLevel);


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

  bool regridIceGeometry(LevelSigmaCS& a_coords,
			 const RealVect& a_dx,
			 const RealVect& a_domainSize,
			 const Real& a_time,
			 const LevelSigmaCS* a_crseCoords,
			 const int a_refRatio);

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


  /// ice thickness at edge of domain for use in computing velocity Bc's
  /** this is an initial hack intended to get things up and 
      running, and assumes that thickness is constant along the
      front... This is public so that the BCFunc can access it.
  */
  static RealVect s_edgeThickness;

protected:

  Vector<RefCountedPtr<RealFunction<RealVect> > > m_bedrockFunction;
  RefCountedPtr<RealFunction<RealVect> > m_thicknessFunction;

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


private:
  // Disallowed for all the usual reasons
  void operator=(const MarineIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  MarineIBC(const MarineIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
