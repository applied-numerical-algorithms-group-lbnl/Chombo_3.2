#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// coefficient of basal friction stored in a FortranInterfaceIBC
#ifndef _FORTRANBASALFRICTION_H_
#define _FORTRANBASALFRICTION_H_

#include "BasalFriction.H"
#include "NamespaceHeader.H"


class FortranInterfaceBasalFriction : public BasalFriction
{
public:
  /// constructor
  FortranInterfaceBasalFriction();
  
  /// destructor
  virtual ~FortranInterfaceBasalFriction()
  {;}
  
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
  
  virtual void setReferenceFAB(Real* a_data_ptr,
			       const int* a_dimInfo,
			       const RealVect& a_dx,
			       const IntVect& a_ghost,
			       const bool a_nodal = false);
protected:
  FArrayBox m_fab;
  RealVect m_dx;
  IntVect m_ghost;
  bool m_verbose;

private:
  void operator=(const FortranInterfaceBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  FortranInterfaceBasalFriction(const FortranInterfaceBasalFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

 

};

#include "NamespaceFooter.H"

#endif
