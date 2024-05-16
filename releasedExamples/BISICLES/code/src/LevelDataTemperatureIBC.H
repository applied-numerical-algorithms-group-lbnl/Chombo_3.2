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
//  LevelDataTemperatureIBC.H
// ============

#ifndef _LEVELDATATEMPERATUREIBC_H_
#define _LEVELDATATEMPERATUREIBC_H_

#include "IceInternalEnergyIBC.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   PhysIBC-derived class which stores initial temperature data
   and imposes either periodic or reflection boundary conditions
 */
class LevelDataTemperatureIBC : public ReflectionIceInternalEnergyIBC
{  
  RefCountedPtr<LevelData<FArrayBox> > m_bulkData; // a temperature or internal energy
  Real m_defaultTemperature, m_defaultE;
  bool m_dataIsInternalEnergy;
#if BISICLES_Z == BISICLES_LAYERED
  RefCountedPtr<LevelData<FArrayBox> > m_surfaceData ; // a temperature or internal energy
  RefCountedPtr<LevelData<FArrayBox> > m_basalHeatFlux ; // usually a heat flux
#endif
  RealVect m_dx;
public:
  /// Constructor
  /**
   */
  LevelDataTemperatureIBC(RefCountedPtr<LevelData<FArrayBox> > a_bulkData, 
#if BISICLES_Z == BISICLES_LAYERED
			  RefCountedPtr<LevelData<FArrayBox> > a_surfaceData,
			  RefCountedPtr<LevelData<FArrayBox> > a_basalHeatFlux,
#endif
			  const RealVect& a_dx,const Real& a_defaultTemperature,
			  const bool& a_dataIsInternalEnergy );

  /// Destructor
  /**
   */
  virtual ~LevelDataTemperatureIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new IceInternalEnergyIBC object with 
     m_isDefined = false (i.e., its define() must be called before it
     is used). 
   */
  virtual LevelDataTemperatureIBC* new_internalEnergyIBC();

  /// set a basal heat flux. units are Joules / Year
  virtual void basalHeatFlux(LevelData<FArrayBox>& a_flux,
			     const AmrIceBase& a_amrIce, 
			     int a_level, Real a_dt);

#if BISICLES_Z == BISICLES_LAYERED
  virtual void initializeIceInternalEnergy(LevelData<FArrayBox>& a_E, 
					   LevelData<FArrayBox>& a_tillWaterDepth,
					   LevelData<FArrayBox>& a_surfaceE, 
					   LevelData<FArrayBox>& a_basalE, 
					   const AmrIceBase& a_amrIce, 
					   int a_level, Real a_dt);

#elif BISICLES_Z == BISICLES_FULLZ
#error BISICLES_FULLZ not implemented
#endif


  static LevelDataTemperatureIBC* parse(ParmParse& a_pp);

private:
  


  // Disallowed for all the usual reasons
  void operator=(const LevelDataTemperatureIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelDataTemperatureIBC(const LevelDataTemperatureIBC& a_input)
    :m_bulkData(a_input.m_bulkData)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
