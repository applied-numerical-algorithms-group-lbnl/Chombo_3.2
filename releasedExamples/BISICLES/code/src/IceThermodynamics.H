#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef ICE_THERMODYNAMICS_H
#define ICE_THERMODYNAMICS_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "NamespaceHeader.H"

/// Common function for thermodynamics, including transport within columns
class IceThermodynamics
{



public:

  static Real m_ice_conductivity,
    m_ice_heat_capacity, 
    m_ice_latent_heat_fusion,
    m_ice_pressure_melt_factor,
    m_triple_point,
    m_seconds_per_unit_time,
    m_moisture_conductivity,
    m_water_fraction_drain,
    m_water_fraction_max,
    m_water_drain_factor,
    m_till_water_drain_factor,
    m_till_water_max,
    m_ice_density,
    m_water_density,
    m_gravity,
    m_floating_base_max_heat_flux;
  
  static Real icepmeltfactor() {return m_ice_pressure_melt_factor;}
  static Real IceConductivity() {return m_ice_conductivity;}
  static Real IceHeatCapacity() {return m_ice_heat_capacity;}
  
  /// Send some args, constant, and ParmParse data to f90 via ChF77
  static void setConstants(Real a_rhoi, Real a_rhow, Real a_gravity, Real a_seconds_per_unit_time);
  
  ///compose internal energy F(T,w) from temperature T and water fraction  
  static void composeInternalEnergy(FArrayBox& a_F, const FArrayBox& a_T, const FArrayBox& a_w, 
				    const Box& a_box, bool a_test = true);
    
  ///compose internal energy F(T,w) from temperature T assuming water fraction is zero  
  static void composeInternalEnergy(FArrayBox& a_F, const FArrayBox& a_T,  
				    const Box& a_box, bool a_test = true)
  {
    FArrayBox w(a_box,a_T.nComp()); w.setVal(0.0);
    composeInternalEnergy(a_F, a_T, w, a_box, a_test);
  }

  ///decompose internal energy F into temperarure T and water fraction W, given pressure P
  static void decomposeInternalEnergy(FArrayBox& a_T,  FArrayBox& a_w, 
			       const FArrayBox& a_F, const FArrayBox& a_P, 
			       const Box& a_box);

  ///Advance column thermodyamics through one time step
  static void timestep(FArrayBox& a_internalEnergy,
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
		       const Box& a_box);
  

};


#include "NamespaceFooter.H"
#endif 
