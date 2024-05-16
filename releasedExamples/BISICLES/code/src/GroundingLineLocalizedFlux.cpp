
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "GroundingLineLocalizedFlux.H"
#include "AmrIceBase.H"
#include "BisiclesF_F.H"
#include "NamespaceHeader.H"

GroundingLineLocalizedFlux::GroundingLineLocalizedFlux
(SurfaceFlux* a_groundingLineFlux, 
 SurfaceFlux* a_ambientFlux, 
 const Real& a_powerOfThickness)
  :m_groundingLineFlux(a_groundingLineFlux),
   m_ambientFlux(a_ambientFlux),
   m_powerOfThickness(a_powerOfThickness)
{
  CH_assert(a_ambientFlux != NULL);
  CH_assert(a_groundingLineFlux != NULL);
}

SurfaceFlux* GroundingLineLocalizedFlux::new_surfaceFlux()
{
  SurfaceFlux* g = m_groundingLineFlux->new_surfaceFlux();
  SurfaceFlux* a = m_ambientFlux->new_surfaceFlux();
  return static_cast<SurfaceFlux*>(new GroundingLineLocalizedFlux(g,a,m_powerOfThickness));
}

void GroundingLineLocalizedFlux::surfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{

  const LevelData<FArrayBox>& proximity = 
    *a_amrIce.groundingLineProximity(a_level);

  if (m_groundingLineFlux)
    m_groundingLineFlux->surfaceThicknessFlux(a_flux,a_amrIce,a_level,a_dt);

  LevelData<FArrayBox> ambientFlux;
  ambientFlux.define(a_flux);
      
  if (m_ambientFlux)
    m_ambientFlux->surfaceThicknessFlux(ambientFlux,a_amrIce,a_level,a_dt);

  for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit)
    {

      Real pprox = 1.0;
      Real pthck = m_powerOfThickness;
      const FArrayBox& thck = a_amrIce.geometry(a_level)->getH()[dit];
    
      FORT_PROXIMITYFILL(CHF_FRA1(a_flux[dit],0),
			 CHF_CONST_FRA1(ambientFlux[dit],0),
			 CHF_CONST_FRA1(proximity[dit],0),
			 CHF_CONST_REAL(pprox),
			 CHF_CONST_FRA1(thck,0),
			 CHF_CONST_REAL(pthck),
			 CHF_BOX(proximity[dit].box()));

      
      
    }
}
#include "NamespaceFooter.H"
