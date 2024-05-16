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

#ifndef _TWISTYSTREAMFRICTION_H_
#define _TWISTYSTREAMFRICTION_H_

#include "BasalFriction.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"

/// simple basal friction resulting in a twisting ice stream
/** 
    beta^2 = betaVal*(1 + eps +
    sin(2*Pi*Omega[1]*y + offset))

    where offset = magOffset*sin(2*Pi*Omega[0]*x)
*/
class twistyStreamFriction : public BasalFriction
{
public:
  /// constructor
  twistyStreamFriction();

  twistyStreamFriction(const Real& a_betaVal,
                       const RealVect& a_omega,
                       const Real& m_magOffset,
                       const Real& a_eps,
                       const RealVect& a_domainSize);
    
  /// destructor
  virtual ~twistyStreamFriction();

  /// factory method
  /** return a pointer to a new BasalFriction object
   */
  virtual BasalFriction* new_basalFriction() const;

  /// define basal friction coefficient beta^2 and place in a_betaSqr
  /** time and dt are included in case this is time-dependent. Units 
      should be Pa*a/m (any conversion to mks units is internal to the
      AmrIce code)       
  */
  virtual void setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);

  /// set friction value in Pa*a/m)
  void setParameters(const Real& a_betaVal, 
                     const RealVect& a_omega,
                     const Real& a_magOffset,
                     const Real& a_eps,
                     const RealVect& a_domainSize);
   
  Real m_betaVal;
  RealVect m_omega;
  Real m_magOffset;
  Real m_eps;
  RealVect m_domainSize;

private:
  // Disallowed for all the usual reasons
  void operator=(const twistyStreamFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  twistyStreamFriction(const constantFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }
};
  

#include "NamespaceFooter.H"


#endif // multiple-include preventer
