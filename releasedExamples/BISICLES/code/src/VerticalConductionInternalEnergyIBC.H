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
//  VerticalConductionInternalEnergyIBC.H
// ============

#ifndef _VERTICALCONDUCTIONINTERNALENERGYIBC_H_
#define _VERTICALCONDUCTIONINTERNALENERGYIBC_H_

#include "IceInternalEnergyIBC.H"
#include "ParmParse.H"
#include "SurfaceFlux.H"
#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   PhysIBC-derived class which computes an initial temperature by 
   solving a vertical conduction problem given a surface temperature
   and a basal heat flux, and imposes either periodic or reflection 
   boundary conditions
*/
class VerticalConductionInternalEnergyIBC : public ReflectionIceInternalEnergyIBC
{  
  
public:
  /// Constructor
  /**
   */
  VerticalConductionInternalEnergyIBC(SurfaceFlux* a_basalDissipation);
  /// Destructor
  /**
   */
  virtual ~VerticalConductionInternalEnergyIBC();

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
  virtual VerticalConductionInternalEnergyIBC* new_internalEnergyIBC();

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

  static VerticalConductionInternalEnergyIBC* parse(ParmParse&);

private:
  
  SurfaceFlux* m_basalDissipation;

  SurfaceFlux& basalDissipation() const
  {
    if (m_basalDissipation == NULL)
      {
	CH_assert(m_basalDissipation != NULL);
	MayDay::Error("no basal dissipation");
      }
    return *m_basalDissipation;
  }


  // Disallowed for all the usual reasons
  void operator=(const VerticalConductionInternalEnergyIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  VerticalConductionInternalEnergyIBC(const VerticalConductionInternalEnergyIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
