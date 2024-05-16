#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

///@file Complex surface fluxes, which contain one or more SurfaceFlux objects

#ifndef _COMPLEXSURFACEFLUX_H_
#define _COMPLEXSURFACEFLUX_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "SurfaceFlux.H"
#include "NamespaceHeader.H"


/// class which takes an input fortran array 
/** averages or interpolates as necessary to fill the flux
 */
class fortranInterfaceFlux : public SurfaceFlux
{
public:
 
  fortranInterfaceFlux();

  virtual ~fortranInterfaceFlux() {}

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

  /// set fortran array-valued surface flux
  void setFluxVal(Real* a_data_ptr,
                  const int* a_dimInfo,
                  const int* a_boxlo, const int* a_boxhi, 
                  const Real* a_dew, const Real* a_dns,
                  const IntVect& a_offset,
                  const IntVect& a_nGhost,
                  const ProblemDomain& a_domain,
                  const bool a_nodal);
    

  /// has setGrids been called?
  bool gridsSet() const {return m_gridsSet;}


protected:
  IntVect m_fluxGhost;
  bool m_verbose;

  ProblemDomain m_domain;
  DisjointBoxLayout m_grids;
  bool m_gridsSet;

  // input flux from CISM
  // alias to data array passed in from CISM
  FArrayBox m_inputFlux;
  /// true if inputFlux is nodeal
  bool m_nodalFlux;

  /// cell-centered input flux
  /**
     if nodalFlux is false, this is identical to m_inputFlux.
     if nodalFlux is true, this contains the cell-averaged version
     of the inputFlux
  */
  FArrayBox m_ccInputFlux;

  RealVect m_inputFluxDx;

  /// distributed version of flux
  RefCountedPtr<LevelData<FArrayBox> > m_inputFluxLDF;

  bool m_isValSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const fortranInterfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  fortranInterfaceFlux(const fortranInterfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

/// class which multiplies two surface fluxes together
/** 
    Useful for using an externally defined mask, among other things.
    
    Run-time configuration is carried out in SurfaceFlux::parse. 
    If the string `A.B.C.type = productFlux` is found, then the expected form is:
    
     - `A.B.C.flux1 =  <valid SurfaceFlux type>`
     - `A.B.C.flux2 =  <valid SurfaceFlux type>`

 */
class ProductSurfaceFlux : public SurfaceFlux
{

public:
  
  ProductSurfaceFlux(SurfaceFlux* a_Flux1, SurfaceFlux* a_flux2);

  virtual ~ProductSurfaceFlux();

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);


protected:
  SurfaceFlux* m_flux1;
  SurfaceFlux* m_flux2;

  
private:
  // Disallowed for all the usual reasons
  void operator=(const ProductSurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  ProductSurfaceFlux(const ProductSurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  ProductSurfaceFlux()
  {
    MayDay::Error("invalid operator");
  }

};
  
/// Applies a different surface flux depending on the mask
/**
   Applies a different surface fluxes for grounded and floating ice
   and for open sea and land.

   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = maskedFlux` is found, then the expected form is:
   
    - `A.B.C.grounded = <valid SurfaceFlux type (optional, default ZeroFlux) >`
    - `A.B.C.floating = <valid SurfaceFlux type (optional, default ZeroFlux) >`
    - `A.B.C.openland = <valid SurfaceFlux type (optional, default ZeroFlux) >`
    - `A.B.C.opensea  = <valid SurfaceFlux type (optional, default ZeroFlux) >`

 */
class MaskedFlux : public SurfaceFlux
{

public:
  /// Construct a MaskedFlux given its component fluxes
  /**
     @param a_groundedIceFlux flux over grounded ice
     @param a_floatingIceFlux flux over floating ice
     @param a_openSeaFlux  flux over ice free areas below sea level
     @param a_openLandFlux flux over ice free areas above sea level
   */
  MaskedFlux(SurfaceFlux* a_groundedIceFlux, SurfaceFlux* a_floatingIceFlux,
	     SurfaceFlux* a_openSeaFlux, SurfaceFlux* a_openLandFlux);

  virtual ~MaskedFlux(){};

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);


protected:
  SurfaceFlux* m_groundedIceFlux;
  SurfaceFlux* m_floatingIceFlux;
  SurfaceFlux* m_openSeaFlux;
  SurfaceFlux* m_openLandFlux;
  
private:
  // Disallowed for all the usual reasons
  void operator=(const MaskedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  MaskedFlux(const MaskedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  MaskedFlux()
  {
    MayDay::Error("invalid operator");
  }

};
  

///Applies a flux inside a box and (optionally inside a time interval) and zero outside
/**

   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = boxBoundedFlux` is found, then the expected form is:

   - `A.B.C.lo = <float>,<float>` (minimum x,y values)
   - `A.B.C.hi = <float>,<float>`  (maximum x,y values)
   - `A.B.C.time = <float>,<float>` (minimum,maximum time values , default all times)
   - `A.B.C.flux =  <valid SurfaceFlux type>`

 */
class BoxBoundedFlux : public SurfaceFlux
{
  RealVect m_lo,m_hi;
  Real m_startTime, m_endTime;
  SurfaceFlux* m_fluxPtr;
public:
  BoxBoundedFlux
(const RealVect& a_lo,
 const RealVect& a_hi,
 const Real& a_startTime,
 const Real& a_endTime,
 SurfaceFlux* a_fluxPtr)
    : m_lo(a_lo),m_hi(a_hi),m_startTime(a_startTime),m_endTime(a_endTime)
  {
    m_fluxPtr = a_fluxPtr->new_surfaceFlux();
    if (!(m_hi > m_lo))
      {
	pout() <<  "!(m_hi > m_lo), m_lo = " << m_lo << " m_hi = " << m_hi << endl;
	CH_assert( (m_hi > m_lo) );
	MayDay::Error("BoxBoundedFlux::  !(m_hi > m_lo)");
      }
  }
  

   /// destructor
  virtual ~BoxBoundedFlux(){};

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				  const AmrIceBase& a_amrIce, 
				  int a_level, Real a_dt);

private:
  // Disallowed for all the usual reasons
  void operator=(const BoxBoundedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  BoxBoundedFlux(const BoxBoundedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  BoxBoundedFlux()
  {
    MayDay::Error("invalid operator");
  }

};

///Provide the flux a*x + b*y, where x,y are fluxes and a,b are scalars
/**
   Run-time configuration is carried out in SurfaceFlux::parse. If
   the string `A.B.C.type = axbyFlux` is found, then the
   expected form is:
    - `A.B.C.a = <float>`
    - `A.B.C.x.type = <valid SurfaceFlux type>`
    - `A.B.C.x.<etc> = <the remaining A.B.C.x configuration>`
    - `A.B.C.b = <float>`
    - `A.B.C.y.type = <valid SurfaceFlux type> `
    -` A.B.C.y.<etc> = <the remaining A.B.C.y configuration>`
 */
class AxbyFlux : public SurfaceFlux
{
  SurfaceFlux* m_x;
  SurfaceFlux* m_y;
  Real m_a,m_b;
  
public:
  /// Construct an AxbyFlux  a*x + b*y given two fluxes x,y and two scalars a,b.
  /**
     @param a_a scalar a
     @param a_x flux x
     @param a_b scalar b
     @param a_y flux y
   */
  AxbyFlux(const Real& a_a, SurfaceFlux* a_x, 
	   const Real& a_b, SurfaceFlux* a_y);

  virtual ~AxbyFlux();

  virtual SurfaceFlux* new_surfaceFlux();
  
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

private:
  
  // Disallowed for all the usual reasons
  void operator=(const AxbyFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AxbyFlux(const AxbyFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  AxbyFlux()
  {
    MayDay::Error("invalid operator");
  }

};


///Adds together  a list of fluxes
/**
   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = compositeFlux` is found, then the expected form is:

   - `A.B.C.nElements = <int>` (number of fluxes)
   - `A.B.C.element0 =  <valid SurfaceFlux type>`
   - `A.B.C.element1 =  <valid SurfaceFlux type>`
   - ...
   - `A.B.C.elementn =  <valid SurfaceFlux type>`
 */
class CompositeFlux : public SurfaceFlux
{
  
  Vector<SurfaceFlux*> m_fluxes;

public:
  CompositeFlux(const Vector<SurfaceFlux*>& a_fluxes);

  virtual ~CompositeFlux();

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);
  
private:
  // Disallowed for all the usual reasons
  void operator=(const CompositeFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  CompositeFlux(const CompositeFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  CompositeFlux()
  {
    MayDay::Error("invalid operator");
  }


};

/// SurfaceFlux which has a piecewise linear dependence on ice thickness 
/**

   Usually used in conjunction with MaskedFlux to impose a depth 
   dependent ice shelf melt rate. Provided that the cavity water
   thickness \f$ d = s-h-b > d_0\f$, the surface flux is computed
   from  \f$ n \f$ pairs of abissiciae and ordinates \f$ ( h_i, f_i ) \f$ 

   \f$ f = f_0 \f$  for  \f$ h < h_0 \f$ 

   \f$ f = f_{i}(h-h_{i+1})/(h{i}-h_{i+1})  + f_{i+1}(h-h_{i})/(h_{i+1}-h_{i})  + \f$ for \f$ h_i \leq h < h_{i+1} \f$

   \f$ f = f_{n-1} \f$ for \f$ h > h_{n-1} \f$

   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = piecewiseLinearFlux` is found, then the expected form is:
    
   - `A.B.C.n = <int>` (number of abscissae and ordinates)
   - `A.B.C.abscissae = <float> <float> ...` (list of thickness absissae)
   - `A.B.C.ordinates = <float> <float> ...` (list of flux ordinates) 
   - `A.B.C.minWaterDepth = <float (optional, default -1.0)>

*/
class PiecewiseLinearFlux : public SurfaceFlux
{
public:
  /// 
  /**
     @param a_abscissae list of thickness abscissae \f$h_i\f$
     @param a_abscissae list of flux ordinates \f$f_i\f$
     @param a_minWaterDepth set flux to zero in sub-shelf cavities thinner than a_minWaterDepth
   */
  PiecewiseLinearFlux(const Vector<Real>& a_abscissae, 
		      const Vector<Real>& a_ordinates,
		      const Real a_minWaterDepth = -1.0);

  ~PiecewiseLinearFlux(){};

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);
  

protected:
  Vector<Real> m_abscissae;
  Vector<Real> m_ordinates;
  Real m_minWaterDepth;

private:
   // Disallowed for all the usual reasons
  void operator=(const PiecewiseLinearFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PiecewiseLinearFlux(const PiecewiseLinearFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  PiecewiseLinearFlux()
  {
    MayDay::Error("invalid operator");
  }

};


/// Normalizes any flux so that the integrated flux has a specified value
/**
   A surface flux that satisfies
   
   \f$ \int _\Omega  f^2(x,y) d \Omega = a^2  \int _\Omega   d \Omega \f$

   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = normalizedFlux` is found, then the expected form is:

   - `A.B.C.direction = <valid SurfaceFlux>` (spatial pattern of the flux)
   - `A.B.C.amplitude = <Real>` (amplitude, \f$ a \f$)

 */
class NormalizedFlux : public SurfaceFlux
{
  
  SurfaceFlux* m_direction;
  Real m_amplitude;

public:
  /// Construct a normalized flux from a flux and an amplitude
  /**
     @param a_direction : defines the spatial pattern
     @param a_amplitude : defines the amplitude
   */
  NormalizedFlux(SurfaceFlux* a_direction, const Real& a_amplitude);

  virtual ~NormalizedFlux();

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

  
private:
  // Disallowed for all the usual reasons
  void operator=(const NormalizedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  NormalizedFlux(const NormalizedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  NormalizedFlux()
  {
    MayDay::Error("invalid operator");
  }


};

typedef NormalizedFlux NormalizedData;


/// Redistrubutes a flux over a given area to peak around regions of peak negative ice flux divergence in floating ice 
/**

   This surfaceFlux is intended for coupling with ocean models and
   is probably not suitable for any other purpose. The assumption is that 
   ocean models running on a coarser mesh result in a stripe of strong
   melting downstream from the 'correct' location and we want to chose an 
   interpolation scheme that place this melt closer to steep slopes
   present in the ice shelf provided that can be done within a distance
   representing the ocesn model mesh.

   Computes a flux \f$ f \f$, given a length \f$ L \f$ and a flux \f$ g \f$. 
   The domain is divided into a uniform mesh with spacing  \f$ L \f$. 
   Within each cell, a constant
   
   \f$ Q = \int _\Omega g -  w min(0,\nabla . (\vec{u} h)) d \Omega \f$
   
   is computed, while 

   \f$ f = min(0,\nabla . (\vec{u} h)) +  Q / \int _\Omega w d \Omega $

   The total flux is conserved, and applied to 
   *the floating portion* of each L x L cell, and 
 
   \f$ J = \int _\Omega (f -  w \nabla . (\vec{u} h) d \Omega \f)^2$ \f$
   
   is minimal.  \f$ w \f$ is 1 in floating ice and zero elsewhere.

   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = divUHLocalizedFlux` is found, then the expected form is:

   - `A.B.C.flux = <valid SurfaceFlux>` (spatial pattern of the flux)
   - `A.B.C.mesh_spacing = <valid SurfaceFlux>` (mesh spacing)

 */
class FloatingDivUHLocalizedFlux : public SurfaceFlux
{
  
  SurfaceFlux* m_flux;
  Real m_mesh_spacing;
 
public:
  /// Construct a divuhLocalized flux from a flux and a mesh spacing
  /**
     @param a_flux : the flux to be redistributed
     @param a_mesh_spacing : defines the size of the cells over which flux is redistrubuted
   */
  FloatingDivUHLocalizedFlux(SurfaceFlux* a_flux, const Real& a_mesh_spacing);

  virtual ~FloatingDivUHLocalizedFlux();

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

  
private:
  // Disallowed for all the usual reasons
  void operator=(const FloatingDivUHLocalizedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  FloatingDivUHLocalizedFlux(const FloatingDivUHLocalizedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  FloatingDivUHLocalizedFlux()
  {
    MayDay::Error("invalid operator");
  }


};



/// A flux proportional to the difference between current and target thickness
/**
    
   A flux proportional to the difference between current and target thickness.
   f = 1/timesscale * ( thickness_target - thickness)

   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = targetThicknessFlux` is found, then the expected form is:

   - `A.B.C.target = <valid SurfaceFlux>` (spatial pattern of the flux)
   - `A.B.C.timescale = <Real>` (timescale, \f$ a \f$)

 */
class TargetThicknessFlux : public SurfaceFlux
{
  
  SurfaceFlux* m_target;
  Real m_timescale;

public:
  /// Construct a normalized flux from a flux and an amplitude
  /**
     @param a_direction : defines the spatial pattern
     @param a_amplitude : defines the amplitude
   */
  TargetThicknessFlux(SurfaceFlux* a_target, const Real& a_timescale);

  virtual ~TargetThicknessFlux();

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

  
private:
  // Disallowed for all the usual reasons
  void operator=(const TargetThicknessFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  TargetThicknessFlux(const TargetThicknessFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  TargetThicknessFlux()
  {
    MayDay::Error("invalid operator");
  }


};


#include "NamespaceFooter.H"
#endif // multiple-include preventer

