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

#ifndef _BASALFRICTION_H_
#define _BASALFRICTION_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"

/// interface class encapsulating coefficient of basal friction
/// (beta^2) for ice-sheet model
/** This is a virtual base class -- it is expected that specific cases
    will be derived from this
*/
class BasalFriction
{
public:
  /// virtual destructor
  virtual ~BasalFriction() {};

  /// factory method
  /** 
      @return a pointer to a new BasalFriction object
  */
  virtual BasalFriction* new_basalFriction() const = 0;

  /// compute basal friction coefficient beta^2 and place in a_betaSqr
  /** 
      time and dt are included in case this is time-dependent. Units
      should usually be chosen so that the basal traction \f$\tau_b(\beta^2, u) \f$
      is measured in Pascals when \f$u \f$ is measured in meters per year.

      @param a_betaSqr basal friction coefficient beta^2 data
      @param a_coordSys ice sheet geometry
      @param a_time current time
      @param a_dt current time step
      
  */
  virtual void setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt) = 0;

  static BasalFriction *parse(const char* a_prefix, const RealVect& a_domainSize);

private:
  // Disallowed for all the usual reasons
  void operator=(const BasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // // Disallowed for all the usual reasons
  // BasalFriction(const BasalFriction& a_input)
  // {
  //   MayDay::Error("invalid operator");
  // }
};


/// BasalFriction which sets basal friction to zero. Probably completely useless.
class zeroFriction : public BasalFriction
{
public:
  zeroFriction() {};

  ~zeroFriction() {};

  BasalFriction* new_basalFriction() const;

  void setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);


private:
  // Disallowed for all the usual reasons
  void operator=(const zeroFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  zeroFriction(const zeroFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


/// BasalFriction where \f$ \beta^2 \f$ is constant in space and time
/**
   Run time configuration is carried out in BasalFriction::parse(A, ... )
   If the option `A.beta_type == constantBeta` is found, then the 
   expected form is:
   
   - `A.betaValue = <float>` 

 */
class constantFriction : public BasalFriction
{
public:
  constantFriction();

  /// Construct a constant friction with a known value
  /**
     @param a_value basal friction coefficient
   */
  constantFriction(Real a_betaSqr);

  ~constantFriction() {}

  BasalFriction* new_basalFriction() const;

  void setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);

  /// Change the constanst basal friction coefficient
  /**
     @param a_betaSqr basal friction coefficient
   */
  void setFrictionVal(const Real& a_betaSqr);

protected:

  Real m_frictionVal;
  bool m_isValSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const constantFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  constantFriction(const constantFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


/// simple version which sets basal friction to a sinusoidal profile
/** 
    beta^2 = betaVal*(1 + eps +
    sin(2*Pi*Omega[0]*x)*sin(2*Pi*Omega[1]*y))
*/
class sinusoidalFriction : public BasalFriction
{
public:
  /// constructor
  sinusoidalFriction();

  sinusoidalFriction(const Real& a_betaVal,
                     const RealVect& a_omega,
                     const Real& a_eps,
                     const RealVect& a_domainSize);
    
  /// destructor
  virtual ~sinusoidalFriction();

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
  void setSinParameters(const Real& a_betaVal, 
                        const RealVect& a_omega,
                        const Real& a_eps,
                        const RealVect& a_domainSize);
   
  Real m_betaVal;
  RealVect m_omega;
  Real m_eps;
  RealVect m_domainSize;

private:
  // Disallowed for all the usual reasons
  void operator=(const sinusoidalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  sinusoidalFriction(const constantFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


#include "NamespaceFooter.H"


#endif // multiple-include preventer
