#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Virtual base class encapsulating coefficient of basal friction 
// for ice-sheet model

#ifndef _SINGULARSTREAMFRICTION_H_
#define _SINGULARSTREAMFRICTION_H_

#include "BasalFriction.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"

/// A straight or twisting ice stream of width W with a discontinous C
/** 

    Domain dimensions are 0 < x < X , 0 < y < Y

    C = C_slippy if |y - a(x)| < W/2, C_sticky otherwise
    a(x) = Y/2 + A * sin(2 * pi * n * x / X)

*/
class singularStreamFriction : public BasalFriction
{
public:
  /// constructor
  singularStreamFriction(const Real& a_slippyC,
			 const Real& a_stickyC,
			 const Real& a_width,
			 const Real& a_twistNumber,
			 const Real& a_twistAmplitude,
			 const RealVect& a_domainSize);
    
  /// destructor
  virtual ~singularStreamFriction();

  /// factory method
  /** return a pointer to a new BasalFriction object
   */
  virtual BasalFriction* new_basalFriction() const;

  /// define basal friction coefficient beta^2 and place in a_betaSqr
  /** time and dt are included in case this is time-dependent. Units 
      should be Pa*a/m (any conversion to mks units is internal to the
      AmrIce code)       
  */
  virtual void setBasalFriction(LevelData<FArrayBox>& a_C,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);

  /// set friction value in Pa*a/m)
  void setParameters(const Real& a_slippyC,
		     const Real& a_stickyC,
		     const Real& a_width,
		     const Real& a_twistNumber,
		     const Real& a_twistAmplitude,
		     const RealVect& a_domainSize);

  Real m_slippyC, m_stickyC, m_width, m_twistNumber, m_twistAmplitude;
  RealVect m_domainSize;

private:
  // Disallowed for all the usual reasons
  void operator=(const singularStreamFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  singularStreamFriction(const constantFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }
};
  


#include "NamespaceFooter.H"
#endif
