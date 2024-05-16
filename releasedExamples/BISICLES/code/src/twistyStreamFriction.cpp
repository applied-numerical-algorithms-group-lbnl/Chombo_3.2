#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "twistyStreamFriction.H"
#include "IceConstants.H"
#include "CONSTANTS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

// -----------------------------------------------------------------
//   twistyStreamFriction
// -----------------------------------------------------------------


/// constructor
twistyStreamFriction::twistyStreamFriction()
{
}

twistyStreamFriction::twistyStreamFriction(const Real& a_betaVal,
                                           const RealVect& a_omega,
                                           const Real& a_magOffset,
                                           const Real& a_eps,
                                           const RealVect& a_domainSize)
  : m_betaVal(a_betaVal), m_omega(a_omega), m_magOffset(a_magOffset), m_eps(a_eps), m_domainSize(a_domainSize)
{

}
 

/// destructor
twistyStreamFriction::~twistyStreamFriction() 
{
}

/// factory method
/** return a pointer to a new BasalFriction object
 */
BasalFriction* 
twistyStreamFriction::new_basalFriction() const
{
  twistyStreamFriction* newPtr = new twistyStreamFriction;
  newPtr->m_betaVal = m_betaVal;
  newPtr->m_omega = m_omega;
  newPtr->m_magOffset = m_magOffset;
  newPtr->m_eps = m_eps;
  newPtr->m_domainSize = m_domainSize;
  return static_cast<BasalFriction*>(newPtr);

}

/// define basal friction coefficient beta^2 and place in a_betaSqr
/** time and dt are included in case this is time-dependent. Units 
    should be Pa*a/m (any conversion to mks units is internal to the
    AmrIce code)       
*/
void
twistyStreamFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                       LevelSigmaCS& a_coordSys,
                                       Real a_time,
                                       Real a_dt)
{
  // not quite correct in 3d
  CH_assert(SpaceDim != 3);

  RealVect twoPiOmega = m_omega;
  twoPiOmega /= m_domainSize;
  twoPiOmega *= 2.0*Pi;

  RealVect dx = a_coordSys.dx();  
  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisBeta = a_betaSqr[dit];
      
      BoxIterator bit(thisBeta.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= dx;
          
          Real offset = m_magOffset*sin(twoPiOmega[0]*x[0]);
          Real betaVal = m_betaVal*(1.0 +m_eps +D_TERM(1.0,
                                                       *sin(twoPiOmega[1]*x[1]+offset),
                                                       *sin(twoPiOmega[2]*x[2])));

          thisBeta(iv,0) =betaVal;
        }                         
    }
}

void 
twistyStreamFriction::setParameters(const Real& a_betaVal, 
                                    const RealVect& a_omega,
                                    const Real& a_magOffset,
                                    const Real& a_eps,
                                    const RealVect& a_domainSize)
{
  m_betaVal = a_betaVal;
  m_omega = a_omega;
  m_magOffset= a_magOffset;
  m_eps = a_eps;
  m_domainSize = a_domainSize;
}
   



#include "NamespaceFooter.H"
