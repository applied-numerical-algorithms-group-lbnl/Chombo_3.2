#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
// VerticalConductionInternalEnergyIBC.cpp
// ============
//
// PhysIBC-derived class which computes an initial temperature by 
// solving a vertical conduction problem given a surface temperature
// and a basal heat flux,  and imposes either periodic or reflection boundary conditions

#include "VerticalConductionInternalEnergyIBC.H"
#include "IceConstants.H"
#include "IceThermodynamics.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "ReadLevelData.H"
#include "AmrIceBase.H"
#include "amrIceF_F.H"
#include "NamespaceHeader.H"


VerticalConductionInternalEnergyIBC* 
VerticalConductionInternalEnergyIBC::parse
(ParmParse& a_pp)
{

  SurfaceFlux* s = SurfaceFlux::parse("VerticalConductionInternalEnergyIBC.dissipation");
  return new VerticalConductionInternalEnergyIBC(s);
}

VerticalConductionInternalEnergyIBC::VerticalConductionInternalEnergyIBC(SurfaceFlux* a_basalDissipation)
  :m_basalDissipation(NULL)
{
  if (a_basalDissipation != NULL)
    {
      m_basalDissipation = a_basalDissipation->new_surfaceFlux();
    }
  else
    {
      m_basalDissipation = static_cast<SurfaceFlux*>(new zeroFlux());
    }
}

VerticalConductionInternalEnergyIBC::~VerticalConductionInternalEnergyIBC()
{
   if (m_basalDissipation != NULL)
    {
      delete m_basalDissipation;
      m_basalDissipation = NULL;
    }
}

void VerticalConductionInternalEnergyIBC::define(const ProblemDomain& a_domain,
						 const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

void VerticalConductionInternalEnergyIBC::basalHeatFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  CH_TIME("VerticalConductionInternalEnergyIBC::basalHeatFlux");
  const DisjointBoxLayout& dbl = a_flux.disjointBoxLayout();
  for (DataIterator dit(dbl); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(0.0);
    }
}

void VerticalConductionInternalEnergyIBC::initializeIceInternalEnergy
(LevelData<FArrayBox>& a_E,
 LevelData<FArrayBox>& a_tillWaterDepth,
 LevelData<FArrayBox>& a_surfaceE, 
 LevelData<FArrayBox>& a_basalE,
 const AmrIceBase& a_amrIce,
 int a_level, Real a_dt)
{
  CH_TIME("VerticalConductionInternalEnergyIBC::initializeIceInternalEnergy");

  const LevelSigmaCS& coordSys = *a_amrIce.geometry(a_level); 
  const DisjointBoxLayout dbl = coordSys.grids();
  
  Real fakeDt = 1.0e+6;
 
  if (!a_amrIce.surfaceHeatBoundaryDirichlett())
    {
      CH_assert(a_amrIce.surfaceHeatBoundaryDirichlett());
      MayDay::Error("VerticalConductionInternalEnergyIBC requires a dirichlett boundary at the upper surface");
    }
  a_amrIce.surfaceHeatBoundaryData().evaluate(a_surfaceE, a_amrIce, a_level, fakeDt);
  
  //basal heat flux is geothermal flux + dissipation at base
  LevelData<FArrayBox> geoflux(dbl,1,IntVect::Unit);
  a_amrIce.basalHeatBoundaryData().evaluate(geoflux, a_amrIce, a_level, fakeDt);
  LevelData<FArrayBox> dissipation(dbl,1,IntVect::Unit);
  basalDissipation().evaluate(dissipation, a_amrIce, a_level, fakeDt);
  
  for (DataIterator dit(dbl); dit.ok(); ++dit)
    {
      const Box& box = dbl[dit];

      // should not be needed, so don't set
      FArrayBox scaledSurfaceHeatFlux(box,1);

      int surfaceTempDirichlett = 1;
      a_basalE[dit].copy(a_surfaceE[dit]);

      FArrayBox& scaledBasalHeatFlux = geoflux[dit];
      scaledBasalHeatFlux += dissipation[dit];
      scaledBasalHeatFlux /= coordSys.iceDensity();
      
      int nLayers = a_E.nComp();
      FArrayBox rhs(box,nLayers); rhs.setVal(0.0);
      FArrayBox usig(box,nLayers+1); usig.setVal(0.0);

      for (int layer =0; layer < nLayers; layer++)
	{
	  a_E[dit].copy(a_surfaceE[dit],0,layer,1);
	}

      Real time = 0.0;
      

      const Real& rhoi = coordSys.iceDensity();
      const Real& rhoo = coordSys.waterDensity();
      const Real& gravity = coordSys.gravity();
      // ensure these are set
      MayDay::Error("VerticalConductionInternalEnergyIBC not functional");
      //columnThermodynamicsSetConstants(rhoi, rhoo, gravity);
      // FORT_UPDATEINTERNALENERGY
      // 	(CHF_FRA(a_E[dit]),
      // 	 CHF_FRA1(a_tillWaterDepth[dit],0),
      // 	 CHF_FRA1(a_surfaceE[dit],0), 
      // 	 CHF_FRA1(a_basalE[dit],0),
      // 	 CHF_CONST_FRA1(scaledSurfaceHeatFlux,0),
      // 	 CHF_CONST_FRA1(scaledBasalHeatFlux,0),
      // 	 CHF_CONST_FIA1(coordSys.getFloatingMask()[dit],0),
      // 	 CHF_CONST_FIA1(coordSys.getFloatingMask()[dit],0),
      // 	 CHF_CONST_FRA(rhs),
      // 	 CHF_CONST_FRA1(coordSys.getH()[dit],0),
      // 	 CHF_CONST_FRA1(coordSys.getH()[dit],0),
      // 	 CHF_CONST_FRA(usig),
      // 	 CHF_CONST_VR(coordSys.getFaceSigma()),
      // 	 CHF_CONST_VR(coordSys.getDSigma()),
      // 	 CHF_CONST_REAL(time), 
      // 	 CHF_CONST_REAL(fakeDt),
      // 	 CHF_CONST_INT(nLayers),
      // 	 CHF_CONST_INT(surfaceTempDirichlett),
      // 	 CHF_BOX(box));
    }
    
  int dbg=0;dbg++;

}

VerticalConductionInternalEnergyIBC* 
VerticalConductionInternalEnergyIBC::new_internalEnergyIBC()
{
  return new VerticalConductionInternalEnergyIBC(m_basalDissipation);
}


#include "NamespaceFooter.H"
