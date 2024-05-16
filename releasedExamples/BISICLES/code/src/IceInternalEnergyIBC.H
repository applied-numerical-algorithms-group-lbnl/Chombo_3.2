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
//  IceInternalEnergyIBC.H
// ============
//
// Virtual bases that extends PhysIBC for internal energy transport. Also provides boundary conditions
// needed in the treatment of internal energy diffusion, which need to 
// be handled quite differently in the layered 2D and full 3D
// models.

#ifndef _ICEINTERNALENERGY_H_
#define _ICEINTERNALENERGY_H_

#include "PhysIBC.H"
#include "LevelSigmaCS.H"
#include "NamespaceHeader.H"

class AmrIceBase;
/// Common virtual base class for  internal energy transport IBC
class IceInternalEnergyIBC : public PhysIBC
{
public:
  virtual ~IceInternalEnergyIBC() {}
  //Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object. Default is to just call
     PhysIBC::define
  */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx)
  {PhysIBC::define(a_domain, a_dx);}

    
  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual PhysIBC* new_physIBC() 
  {return static_cast<PhysIBC*>(new_internalEnergyIBC());}
  
  /// same as new_physIBC, except pointer cast to an IceInternalEnergyIBC
  virtual IceInternalEnergyIBC* new_internalEnergyIBC() = 0;
  

  virtual void initialize(LevelData<FArrayBox>& a_U);

  /// set a basal heat flux : by default set to zero. units are Joules / Year
  /// default implementation sets zero flux.
  virtual void basalHeatFlux(LevelData<FArrayBox>& a_flux,
			     const AmrIceBase& a_amrIce, 
			     int a_level, Real a_dt)=0;


#if BISICLES_Z == BISICLES_LAYERED
  virtual void initializeIceInternalEnergy(LevelData<FArrayBox>& a_E,
					   LevelData<FArrayBox>& a_tillWaterDepth,
					   LevelData<FArrayBox>& a_surfaceE, 
					   LevelData<FArrayBox>& a_basalE,
					   const AmrIceBase& a_amrIce, 
					   int a_level, Real a_dt)=0;

  virtual void setIceInternalEnergyBC(LevelData<FArrayBox>& a_E,
				      LevelData<FArrayBox>& a_tillWaterDepth,
				      LevelData<FArrayBox>& a_surfaceE, 
				      LevelData<FArrayBox>& a_basalE,
				      const LevelSigmaCS& a_coordSys)=0;


#elif BISICLES_Z == BISICLES_FULLZ
#error BISICLES_FULLZ not implemented
#endif

  /// Set boundary fluxes
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time);
  
  /// Set boundary slopes
  virtual void setBdrySlopes(FArrayBox&       a_dW,
			     const FArrayBox& a_W,
			     const int&       a_dir,
			     const Real&      a_time)
  {
    // one sided differences
  }
  
  
  virtual void artViscBC(FArrayBox&       a_F,
			 const FArrayBox& a_U,
			 const FArrayBox& a_divVel,
			 const int&       a_dir,
			 const Real&      a_time)
  {
    //shouldn't get to here, so
    MayDay::Error("IceInternalEnergyIBC::artViscBC not implemented");
  }


};

/// Concrete IceInternalEnergyIBC that sets ice temperature to be constant
class ConstantIceTemperatureIBC : public IceInternalEnergyIBC
{
  Real m_T;
  
public:
  virtual ~ConstantIceTemperatureIBC(){};
  ConstantIceTemperatureIBC(Real a_T) : m_T(a_T) {}
     
  virtual ConstantIceTemperatureIBC* new_internalEnergyIBC();

  /// set a basal heat flux to zero. units are Joules / Year
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

  
  virtual void setIceInternalEnergyBC(LevelData<FArrayBox>& a_E, 
				      LevelData<FArrayBox>& a_tillWaterDepth,
				     LevelData<FArrayBox>& a_surfaceT, 
				     LevelData<FArrayBox>& a_basalT,
				     const LevelSigmaCS& a_coordSys);

#elif BISICLES_Z == BISICLES_FULLZ
#error BISICLES_FULLZ not implemented
#endif
  void set(Real a_T)
  {
    m_T = a_T;
  }
  Real T() 
  {
    return (m_T);
  }

private:
  ConstantIceTemperatureIBC(){}

};

/// Partial implementation of IceInternalEnergyIBC, provides reflection (or periodic) boundary conditions
class ReflectionIceInternalEnergyIBC : public IceInternalEnergyIBC
{  
  
public:
  
  virtual void setIceInternalEnergyBC
  (LevelData<FArrayBox>& a_E, 
   LevelData<FArrayBox>& a_tillWaterDepth,
   LevelData<FArrayBox>& a_surfaceE, 
   LevelData<FArrayBox>& a_basalE,
   const LevelSigmaCS& a_coordSys);
  
};




#include "NamespaceFooter.H"
#endif
