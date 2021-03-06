#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SMOOTHINTERSECTION_H_
#define _SMOOTHINTERSECTION_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
   This implicit function is the maximum of two (of more) implicit
   functions with which it is constructed.  This results in an implicit
   function whose interior is the intersection of the interiors of the given
   implicit functions.
   For the maximum, it uses the smooth absolute value Ae(a,b) \approx |a-b|
   max(a,b) = 0.5*(a+b + | a-b |)
*/
class SmoothIntersection: public BaseIF
{
public:

  ///
  /**
     Constructor specifying any number of implicit functions as inputs.
  */
  SmoothIntersection(const Vector<BaseIF*>& a_impFuncs,
                     const Real           & a_delta);

  /// Destructor
  virtual ~SmoothIntersection();

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  ///
  virtual Real derivative(const  IntVect& a_deriv,
                          const RealVect& a_point) const;


  virtual BaseIF* newImplicitFunction() const;

  ///max(a,b) = 0.5*(a+b + | a-b |)
  /**
     loops through all the functions
   */
  virtual Real smoothMax(const  IntVect & a_deriv,
                         const RealVect & a_point,
                         const int      & a_closestIF, 
                         const int      & a_nextClosestIF
                         ) const;

  void findClosest(int            & a_closestIF, 
                   int            & a_nextClosestIF,
                   int            & a_numWithinDelta,
                   const RealVect & a_point) const;

protected:
  int             m_numFuncs; // number of implicit functions
  Vector<BaseIF*> m_impFuncs; // implicit functions
  Real            m_delta   ; //smoothing length scale


  //debug stuff
   IntVect  m_ivDebug;
   Real     m_dxDebug;
  RealVect  m_rvDebug;

  

private:
  SmoothIntersection()
  {
    MayDay::Abort("SmoothIntersection uses strong construction");
  }

  SmoothIntersection(const SmoothIntersection& a_inputIF)
  {
    MayDay::Abort("SmoothIntersection does not have copy construction");
  }

  void operator=(const SmoothIntersection& a_inputIF)
  {
    MayDay::Abort("SmoothIntersection doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
