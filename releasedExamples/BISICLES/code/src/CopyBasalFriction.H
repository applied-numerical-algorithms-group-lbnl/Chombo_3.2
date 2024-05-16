#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// coefficient of basal friction stored in a CopyIBC
#ifndef _COPYBASALFRICTION_H_
#define _COPYBASALFRICTION_H_

#include "BasalFriction.H"
#include "RefCountedPtr.H"
#include "NamespaceHeader.H"

/// Allows pre-defined basal friction coefficients to be passed in and copied
/** use of the refCountedPtr means that we can reference a single location 
    which can be set and modified elsewhere. This is primarily intended for 
    use with the solverBenchmark code, but may find uses within the main 
    BISICLES code as well. Will obviously need to be redefined at 
    regridding time*/
class CopyBasalFriction : public BasalFriction
{
public:
  /// null constructor 
  CopyBasalFriction();

  /// Set reference value for friction coefficients -- calls associated define
  CopyBasalFriction(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_vectBeta,
                    Vector<ProblemDomain>& a_vectDomain);
  
  /// destructor
  virtual ~CopyBasalFriction()
  {;}
  
  /// Set reference value for friction coefficients
  void define(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_vectBeta,
              Vector<ProblemDomain>& a_vectDomain);

  /// factory method
  /** return a pointer to a new BasalFriction object
   */
  virtual BasalFriction* new_basalFriction() const;

  /// define basal friction coefficient C and place in a_betaSqr
  /** time and dt are included in case this is time-dependent. Units 
      should be Pa*a/m 
  */
  virtual void setBasalFriction(LevelData<FArrayBox>& a_C,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
				Real a_dt);
  
protected:
  /// frictional coefficient values
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_vectBeta;
  /// Problem domains are used to index levels
  Vector<ProblemDomain> m_vectDomains;
  bool m_verbose;

private:
  void operator=(const CopyBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  CopyBasalFriction(const CopyBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

 

};

#include "NamespaceFooter.H"

#endif
