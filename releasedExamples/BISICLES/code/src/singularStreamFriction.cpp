#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "singularStreamFriction.H"
#include "IceConstants.H"
#include "CONSTANTS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

// -----------------------------------------------------------------
//   singularStreamFriction
// -----------------------------------------------------------------


/// constructor

singularStreamFriction::singularStreamFriction
(const Real& a_slippyC,
 const Real& a_stickyC,
 const Real& a_width,
 const Real& a_twistNumber,
 const Real& a_twistAmplitude,
 const RealVect& a_domainSize)
  :m_slippyC(a_slippyC),m_stickyC(a_stickyC),m_width(a_width),
   m_twistNumber(a_twistNumber),m_twistAmplitude(a_twistAmplitude),
   m_domainSize(a_domainSize)
{
}
 

/// destructor
singularStreamFriction::~singularStreamFriction() 
{
}

/// factory method
/** return a pointer to a new BasalFriction object
 */
BasalFriction* 
singularStreamFriction::new_basalFriction() const
{
  singularStreamFriction* newPtr = new singularStreamFriction
    (m_slippyC,m_stickyC,m_width,m_twistNumber,m_twistAmplitude,m_domainSize);
  return static_cast<BasalFriction*>(newPtr);

}

/// define basal friction coefficient C and place in a_C
/** time and dt are included in case this is time-dependent. Units 
    should be Pa*a/m (any conversion to mks units is internal to the
    AmrIce code)       
*/
void
singularStreamFriction::setBasalFriction(LevelData<FArrayBox>& a_C,
					 LevelSigmaCS& a_coordSys,
					 Real a_time,
					 Real a_dt)
{
  // not quite correct in 3d
  CH_assert(SpaceDim != 3);
  const int xDir = 0; const int yDir = 1; 
  const Real omega = 2.0 * Pi * m_twistNumber / m_domainSize[xDir];
  const Real A = m_twistAmplitude * m_domainSize[yDir]; 
  RealVect dx = a_coordSys.dx();  
  for (DataIterator dit = a_C.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& C = a_C[dit];
      C.setVal(m_stickyC);
      BoxIterator bit(C.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= dx;

	  Real y0 = A* sin(omega*x[xDir]) + 0.5* m_domainSize[yDir];
	  if (std::abs(x[yDir] - y0) < 0.5*m_width)
	    {
	      C(iv) = m_slippyC;
	    }
          
        }                         
    }
}



#include "NamespaceFooter.H"
