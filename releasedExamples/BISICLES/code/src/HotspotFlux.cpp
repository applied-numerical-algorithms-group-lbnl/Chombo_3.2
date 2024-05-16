#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "BisiclesF_F.H"
#include "HotspotFlux.H"
#include "AmrIce.H"

#include "NamespaceHeader.H"

/// constructor
HotspotFlux::HotspotFlux() 
  : m_fluxVal(0.0), m_radius(D_DECL(-10.0,-10.0,-10.0)), m_center(D_DECL(0,0,0)), m_startTime(-1.23456e300), m_stopTime(1.23456e300), m_isValSet(false), m_isLocSet(false)
{
}

/// factory method
/** return a pointerto a new SurfaceFlux object
 */
SurfaceFlux* 
HotspotFlux::new_surfaceFlux()
{
  HotspotFlux* newPtr = new HotspotFlux;
  newPtr->m_fluxVal = m_fluxVal;

  newPtr->m_radius = m_radius;
  newPtr->m_center = m_center; 
  newPtr->m_startTime = m_startTime;
  newPtr->m_stopTime = m_stopTime;

  newPtr->m_isValSet = m_isValSet;
  newPtr->m_isLocSet = m_isLocSet;

  return static_cast<SurfaceFlux*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep. flux should be defined in meters/second in the current 
      implementation. 
  */
void 
HotspotFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
                                  const AmrIceBase& a_amrIce, 
                                  int a_level, Real a_dt)
{
  Real time = a_amrIce.time();

  // get geometry info
  RealVect dxLevel = a_amrIce.dx(a_level);
  Real distSqr;
  //Real radSqr = m_radius[0]*m_radius[1];
  Real radSqr = 1.0;

  // initialize to zero
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisFlux = a_flux[dit];
      thisFlux.setVal(0.0);

      if ((time > m_startTime) && (time < m_stopTime))
        {
          BoxIterator bit(thisFlux.box());
          // compute distance from spot center
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              // location of iv cell center 
              RealVect loc(iv);
              loc += 0.5*RealVect::Unit;
              loc *= dxLevel;

              // subtract center
              loc -= m_center;
              // divide by radii
              loc /= m_radius;
              distSqr = D_TERM(loc[0]*loc[0],+loc[1]*loc[1],+loc[2]*loc[2]);
              if (distSqr < radSqr)
                {
                  thisFlux(iv,0) = m_fluxVal;
                } // end if in hot spot
            } // end loop over cells in this box
        } // end if time is active
    } // end loop over boxes
}

/// set flux value in meters/year
void
HotspotFlux::setFluxVal(const Real& a_fluxVal)
{
  m_fluxVal = a_fluxVal;
  m_isValSet = true;
}
  

/// set location of (circular) hot spot
void
HotspotFlux::setSpotLoc(Real a_radius, RealVect a_center)
{ 
  m_radius = a_radius*RealVect::Unit; 
  m_center = a_center;
  m_isLocSet = true;
}


/// set location of (ellipsoid) hot spot
void
HotspotFlux::setSpotLoc(RealVect a_radius, RealVect a_center)
{ 
  m_radius = a_radius; 
  m_center = a_center;
  m_isLocSet = true;
}

/// set start and stop times
void
HotspotFlux::setSpotTimes(Real a_startTime, Real a_stopTime)
{
  m_startTime = a_startTime;
  m_stopTime = a_stopTime;
}



#include "NamespaceFooter.H"
