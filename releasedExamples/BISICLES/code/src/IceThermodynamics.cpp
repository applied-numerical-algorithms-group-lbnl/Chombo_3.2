#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "IceThermodynamics.H"
#include "IceThermodynamicsF_F.H"
#include "IceConstants.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

#if defined(CH_USE_FLOAT) || defined(CH_LANG_CC)
#define iceconductivity 2.1e0
#define ICECONDUCTIVITY 2.1e0
#define moistureconductivity 0.2e0
#define MOISTURECONDUCTIVITY 0.2e0
#define iceheatcapacity 2009.0e0
#define ICEHEATCAPACITY 2009.0e0
#define icelatentheat  334.0e+3
#define ICELATENTHEAT  334.0e+3
#define icepmeltfactor 9.7456e-8
#define ICEPMELTFACTOR 9.7456e-8
#define triplepoint 273.15e0
#define TRIPLEPOINT 273.15e0
#define MIN_TEMPERATURE 200.0e0
#else
#define iceconductivity 2.1d0
#define ICECONDUCTIVITY 2.1d0
#define moistureconductivity 0.2d0
#define MOSITURECONDUCTIVITY 0.2d0
#define iceheatcapacity 2009.0d0
#define ICEHEATCAPACITY 2009.0d0
#define icelatentheat  334.0d+3
#define ICELATENTHEAT  334.0d+3
#define icepmeltfactor 9.7456d-8
#define ICEPMELTFACTOR 9.7456d-8
#define triplepoint 273.15d0
#define TRIPLEPOINT 273.15d0
#define MIN_TEMPERATURE 200.0d0
#endif


//default defs for static data members
// common physical constants
Real IceThermodynamics::m_ice_conductivity(ICECONDUCTIVITY);
Real IceThermodynamics::m_ice_heat_capacity(ICEHEATCAPACITY);
Real IceThermodynamics::m_ice_latent_heat_fusion(ICELATENTHEAT);
Real IceThermodynamics::m_moisture_conductivity(MOISTURECONDUCTIVITY);
Real IceThermodynamics::m_ice_pressure_melt_factor(ICEPMELTFACTOR);
Real IceThermodynamics::m_triple_point(TRIPLEPOINT);
Real IceThermodynamics::m_ice_density(ICE_DENSITY);
Real IceThermodynamics::m_water_density(SEA_WATER_DENSITY);
Real IceThermodynamics::m_gravity(GRAVITY);
Real IceThermodynamics::m_seconds_per_unit_time(SECONDS_PER_TROPICAL_YEAR);
// idiosyncratic parameters 
Real IceThermodynamics::m_water_fraction_drain(0.01); 
Real IceThermodynamics::m_water_fraction_max(0.05);
Real IceThermodynamics::m_water_drain_factor(0.02);
Real IceThermodynamics::m_till_water_drain_factor(0.001);
Real IceThermodynamics::m_till_water_max(4.0);
Real IceThermodynamics::m_floating_base_max_heat_flux(1.2345678e+300);

///compose internal energy F(T,w) from temperature T and water fraction  
void IceThermodynamics::composeInternalEnergy
(FArrayBox& a_F, const FArrayBox& a_T, const FArrayBox& a_W, 
 const Box& a_box, bool a_test)
{
  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_T.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_W.nComp() == a_T.nComp());
  CH_assert(a_F.nComp() == a_T.nComp());

  FORT_COMPOSEINTERNALENERGYICE(CHF_FRA(a_F), 
				CHF_CONST_FRA(a_T), 
				CHF_CONST_FRA(a_W), 
				CHF_CONST_REAL(m_ice_heat_capacity),
				CHF_CONST_REAL(m_ice_latent_heat_fusion),
				CHF_BOX(a_box));

  if (a_test)
    {
      CH_assert(a_F.norm(a_box, 0, 0, a_F.nComp() ) <= 2.0* m_triple_point * m_ice_heat_capacity);
    } 

}

///decompose internal energy F into temperarure T and water fraction W, given pressure P
void IceThermodynamics::decomposeInternalEnergy
(FArrayBox& a_T, FArrayBox& a_W, 
 const FArrayBox& a_F, const FArrayBox& a_P, const Box& a_box)
{
  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_T.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_P.box().contains(a_box));
  CH_assert(a_W.nComp() == a_T.nComp());
  CH_assert(a_F.nComp() == a_T.nComp());
  CH_assert(a_P.nComp() == a_T.nComp());

  Real Tmin(MIN_TEMPERATURE);
  
  FORT_DECOMPOSEINTERNALENERGYICE(CHF_FRA(a_T), 
				  CHF_FRA(a_W), 
				  CHF_CONST_FRA(a_F),
				  CHF_CONST_FRA(a_P),
				  CHF_CONST_REAL(m_ice_heat_capacity),
				  CHF_CONST_REAL(m_ice_latent_heat_fusion),
				  CHF_CONST_REAL(m_ice_pressure_melt_factor),
				  CHF_CONST_REAL(m_triple_point),
				  CHF_CONST_REAL(Tmin),  
				  CHF_BOX(a_box));


}

/// Advance column thermodyamics through one time step
void IceThermodynamics::timestep(FArrayBox& a_internalEnergy,
				 FArrayBox& a_tillWaterDepth,
				 FArrayBox& a_surfaceInternalEnergy,
				 FArrayBox& a_basalInternalEnergy,
			         const FArrayBox& a_scaledSurfaceHeatFlux,
				 const FArrayBox& a_scaledBasalHeatFlux,
				 const FArrayBox& a_tillWaterDrainFactor,
				 const BaseFab<int>& a_oldMask,
				 const BaseFab<int>& a_newMask,
				 const FArrayBox& a_oldH,
				 const FArrayBox& a_newH,
				 const FArrayBox& a_uSigma,
				  const FArrayBox& a_rhs,
				 const Vector<Real>& a_faceSigma,
				 const Vector<Real>& a_dSigma,
				 Real a_halftime,
				 Real a_dt,
				 int a_nLayers,
				 bool a_surfaceTempDirichlett,
				 const Box& a_box)
{

  int surfaceTempDirichlett = a_surfaceTempDirichlett?1:0;
  
  FORT_UPDATEINTERNALENERGY
    (CHF_FRA(a_internalEnergy),
     CHF_FRA1(a_tillWaterDepth,0), 
     CHF_FRA1(a_surfaceInternalEnergy,0), 
     CHF_FRA1(a_basalInternalEnergy,0),
     CHF_CONST_FRA1(a_scaledSurfaceHeatFlux,0),
     CHF_CONST_FRA1(a_scaledBasalHeatFlux,0),
     CHF_CONST_FRA1(a_tillWaterDrainFactor,0),
     CHF_CONST_FIA1(a_oldMask,0),
     CHF_CONST_FIA1(a_newMask,0),
     CHF_CONST_FRA(a_rhs),
     CHF_CONST_FRA1(a_oldH,0),
     CHF_CONST_FRA1(a_newH,0),
     CHF_CONST_FRA(a_uSigma),
     CHF_CONST_VR(a_faceSigma),
     CHF_CONST_VR(a_dSigma),
     CHF_CONST_REAL(a_halftime), 
     CHF_CONST_REAL(a_dt),
     CHF_CONST_INT(a_nLayers),
     CHF_CONST_INT(surfaceTempDirichlett),
     CHF_BOX(a_box));
 
}



/// Send some args, constant, and ParmParse data to f90 via ChF77
void IceThermodynamics::setConstants(Real a_rhoi, Real a_rhow, Real a_gravity, Real a_seconds_per_unit_time)
{

  m_ice_density = a_rhoi;
  m_water_density = a_rhow;
  m_gravity = a_gravity;
  
 
  ParmParse pp("ColumnThermodynamics");

  m_ice_conductivity = iceconductivity;
  m_moisture_conductivity = moistureconductivity;
  m_ice_heat_capacity = iceheatcapacity;
  m_ice_latent_heat_fusion = icelatentheat;
  m_ice_pressure_melt_factor = icepmeltfactor;
  m_triple_point = triplepoint;
  m_seconds_per_unit_time = a_seconds_per_unit_time;

  // maximum water fraction for no drainage 
  m_water_fraction_drain = 0.01 ;
  pp.query("water_fraction_drain", m_water_fraction_drain);
 
  // maximum permitted water fraction. Should -> infinity
  m_water_fraction_max = 0.05;
  pp.query("water_fraction_max", m_water_fraction_max);
  
  // Rate of drainage from temperate ice.
  // Roughly equivalent to Aschwanden 2010 [ doi: 10.3189/2012JoG11J088]
  m_water_drain_factor = 0.02;
  pp.query("water_drain_factor", m_water_drain_factor);

  // Rate of drainage from the till.
  // Default 1 mm a^-1, from Aschwanden et all 2016 [doi:10.1038/ncomms10524]
  m_till_water_drain_factor = 0.001;
  
  pp.query("till_water_drain_factor", m_till_water_drain_factor);


  // max till water depth.
  // max value 4.0. 
  m_till_water_max = 4.0;
  
  pp.query("till_water_max", m_till_water_max);

  // limit heat flux across the ice shelf base (useful when e.g ice shelf geometry is fixed)
  pp.query("floating_base_max_heat_flux",m_floating_base_max_heat_flux);
  
  
  FORT_COLUMNTHERMODYAMICSSETCONSTANTS
    (CHF_CONST_REAL(m_seconds_per_unit_time),
     CHF_CONST_REAL(m_ice_density),
     CHF_CONST_REAL(m_water_density),
     CHF_CONST_REAL(m_gravity),
     CHF_CONST_REAL(m_ice_heat_capacity),
     CHF_CONST_REAL(m_ice_latent_heat_fusion),
     CHF_CONST_REAL(m_ice_conductivity),
     CHF_CONST_REAL(m_moisture_conductivity),
     CHF_CONST_REAL(m_ice_pressure_melt_factor),
     CHF_CONST_REAL(m_triple_point),
     CHF_CONST_REAL(m_water_fraction_drain),
     CHF_CONST_REAL(m_water_fraction_max),
     CHF_CONST_REAL(m_water_drain_factor),
     CHF_CONST_REAL(m_till_water_drain_factor),
     CHF_CONST_REAL(m_till_water_max),
     CHF_CONST_REAL(m_floating_base_max_heat_flux));			       
   
}


#include "NamespaceFooter.H"
