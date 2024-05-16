#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _GAUSSIANBUMPFRICTION_H_
#define _GAUSSIANBUMPFRICTION_H_

#include "BasalFriction.H"

#include "NamespaceHeader.H"

/// basal friction with a gaussian bump centred on b
/* 
C = C0 * (1 - a exp( - (x-bx)^2/(2 * cx^2) - (y-by)^2/(2 * cy^2))

C0, a, b, c can be time dependent, such that C0 = C0[i] etc where t < t[i]

*/

class GaussianBumpFriction : public BasalFriction
{

public:
  GaussianBumpFriction()
  {;}

  GaussianBumpFriction(const Vector<Real>& a_t,
		       const Vector<Real>& a_C0, 
		       const Vector<Real>& a_a,
		       const Vector<RealVect>& a_b,
		       const Vector<RealVect>& a_c);

  virtual ~GaussianBumpFriction()
  {;}

  virtual BasalFriction* new_basalFriction() const;

 virtual void setBasalFriction(LevelData<FArrayBox>& a_C,
			       LevelSigmaCS& a_coordSys,
			       Real a_time,
			       Real a_dt);
 
private:
  
  Vector<Real> m_t;
  Vector<Real> m_C0;
  Vector<Real> m_a;
  Vector<RealVect> m_b;
  Vector<RealVect> m_c;

   // Disallowed for all the usual reasons
  GaussianBumpFriction(const GaussianBumpFriction& a_input)
  {
    MayDay::Error("invalid operator");
  }

};
#include "NamespaceFooter.H"
#endif
