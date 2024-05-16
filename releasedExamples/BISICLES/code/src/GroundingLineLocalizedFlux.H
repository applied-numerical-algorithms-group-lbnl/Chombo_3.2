#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _GROUNDING_LINE_LOCALIZED_FLUX_H_
#define _GROUNDING_LINE_LOCALIZED_FLUX_H_

#include "SurfaceFlux.H"
#include "NamespaceHeader.H"
/// surface flux which has a piecwise linear dependence
/// on thickness H. 
class GroundingLineLocalizedFlux : public SurfaceFlux
{

  SurfaceFlux* m_groundingLineFlux; 
  SurfaceFlux* m_ambientFlux;
  Real m_powerOfThickness;
  

public:
  //constructor
  GroundingLineLocalizedFlux
  (SurfaceFlux* a_groundingLineFlux, 
   SurfaceFlux* a_ambientFlux,
   const Real& a_powerOfThickness
   );

  /// destructor
  ~GroundingLineLocalizedFlux(){};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep. flux should be defined in meters/second in the current 
      implementation. 
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);
  


protected:
  Vector<Real> m_abscissae;
  Vector<Real> m_ordinates;

private:
   // Disallowed for all the usual reasons
  void operator=(const GroundingLineLocalizedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  GroundingLineLocalizedFlux(const GroundingLineLocalizedFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  GroundingLineLocalizedFlux()
  {
    MayDay::Error("invalid operator");
  }

};
#include "NamespaceFooter.H"
#endif
