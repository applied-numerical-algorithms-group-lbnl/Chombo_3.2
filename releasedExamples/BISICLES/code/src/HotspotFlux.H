#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _HOTSPOTFLUX_H_
#define _HOTSPOTFLUX_H_
#include "SurfaceFlux.H"

#include "NamespaceHeader.H"

// simple version which sets surface flux to a hot spot in space and time
class HotspotFlux : public SurfaceFlux
{
public:
  /// constructor
  HotspotFlux();

  /// destructor
  virtual ~HotspotFlux() {}
  
  /// factory method
  /** return a pointerto a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep. flux should be defined in meters/second in the current 
      implementation. 
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

  /// set flux value in meters/year
  void setFluxVal(const Real& a_fluxVal);
  
  /// set location of (circular) hot spot
  void setSpotLoc(Real a_radius, RealVect a_center);

  /// set location of (ellipsoid) hot spot
  void setSpotLoc(RealVect a_radius, RealVect a_center);

  /// set start and stop times
  void setSpotTimes(Real a_startTime, Real a_stopTime);

protected:
  /// amount of flux in hot spot
  Real m_fluxVal;

  /// radius of hot spot
  RealVect m_radius;

  /// center of hot spot
  RealVect m_center;

  /// time at which hot spot turns on
  Real m_startTime;

  /// time at which hot spot turns off
  Real m_stopTime;

  /// checks
  bool m_isValSet;
  bool m_isLocSet;
  
private:
  // Disallowed for all the usual reasons
  void operator=(const constantFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  HotspotFlux(const constantFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


#include "NamespaceFooter.H"
#endif
