#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Virtual base class encapsulating surface fluxes for ice-sheet model

#ifndef _SURFACEFLUX_H_
#define _SURFACEFLUX_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "NamespaceHeader.H"

class AmrIceBase;

class SurfaceFlux;

typedef SurfaceFlux SurfaceData;

/// abstract class defining the surface flux interface
/** 
    The ice sheet model include several inputs natually described
    as surface flux densities: ice accumulation and ablation,
    energy influx, bedrock uplift rates. 

    \todo In fact, since all that is really needed here is a methods to fill 2D fields,
    a very much larger set of inputs could be described in the
    same way. For now, a simple 'typedef SurfaceFlux SurfaceData'
    is our nod to this, but a refactoring exercise where, instead
    SurfaceFlux is a minor specialization of SurfaceData is on the cards. 
*/
class SurfaceFlux
{
  
public:
  
  /// virtual destructor
  virtual ~SurfaceFlux() {};

  /// factory method: return a pointer to a new SurfaceFlux object
  virtual SurfaceFlux* new_surfaceFlux() = 0;

  /// define source term for thickness evolution and place it in flux
  /** 
      @param a_flux output flux data 
      @param     a_amrIce reference to the ice sheet state
      @param     a_level mesh level of a_flux
      @param     a_dt current timestep

      a_dt is included in case one needs integrals or averages over a
      timestep.  flux should be defined in meters per year in the current 
      implementation. 
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt) = 0;
  
  /// Sub-classes of this class turned out to be useful for various data apart from thickness fluxes...
  /** 
      @param     a_data  output data 
      @param     a_amrIce reference to the ice sheet state
      @param     a_level mesh level of a_flux
      @param     a_dt current timestep 
  */
  virtual void evaluate(LevelData<FArrayBox>& a_data,
			const AmrIceBase& a_amrIce, 
			int a_level, Real a_dt)
  {
    surfaceThicknessFlux(a_data, a_amrIce, a_level,a_dt);
  }


  ///assemble a SurfaceFlux* object from ParmParse inputs
  /**
     
     @param  a_prefix prefix string in the ParmParse table

     Run-time configuration of SurfaceFlux objects is complicated by 
     the 'complex fluxes', e.g MaskedFlux that include a number of subsidiary SurfaceFlux members. 
     We support this by parsing strings of the form recursively. For example,
     say that the inputs.example file contained the strings 'surfaceFlux.type = maskedFlux',
     'surfaceFlux.grounded.type = axbyFlux', 'surfaceFlux.floating.type = productFlux'.
     SurfaceFlux::parse would first (call 1) be called with the argument a_prefix = 'surfaceFlux',
     and would itself call SurfaceFlux::parse twice more, once with the argument
     a_prefix = 'surfaceFlux.grounded' (call 1.1), and one with a_prefix = 'surfaceFlux.floating' (call 1.2).
     Call 1.1 would parse the string 'surfaceFlux.grounded.type = axbyFlux', and in turn
     call SurfaceFlux::parse twice with the arguments  a_prefix = 'surfaceFlux.grounded.x' (call 1.1.1) 
     and a_prefix = 'surfaceFlux.grounded.y' (call 1.1.2).
     
     

    
   */
  static SurfaceFlux* parse(const char* a_prefix);

  // /// re-implement to set something up once 
  // virtual void init(const AmrIceBase& a_amrIce){};
  // /// re-implement to set something up once per time step, rather
  // /// than per-level per-timestep
  // virtual void initTimeStep(const AmrIceBase& a_amrIce){};

private:
  // Disallowed for all the usual reasons
  void operator=(const SurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // // Disallowed for all the usual reasons
  // SurfaceFlux(const SurfaceFlux& a_input)
  // {
  //   MayDay::Error("invalid operator");
  // }
};


/// simple implementation of SurfaceFlux which sets surface flux to zero
/**
   Run-time configuration is carried out in SurfaceFlux::parse. A
   Zeroflux object will be created if the string `A.B.C.type = zeroFlux` is found
 */
class zeroFlux : public SurfaceFlux
{
public:
  zeroFlux() {};

  virtual ~zeroFlux() {};

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

private:
  // Disallowed for all the usual reasons
  void operator=(const zeroFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  zeroFlux(const zeroFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


typedef zeroFlux ZeroData;

/// simple implementation of SurfaceFlux which sets surface flux to a constant in space and time
/**

   Run-time configuration is carried out in SurfaceFlux::parse. 
   If the string `A.B.C.type = constantFlux` is found, then the expected form is:

   - `A.B.C.flux_value = <float>` 

*/
class constantFlux : public SurfaceFlux
{
public:
  constantFlux();

  /// Construct a constantFlux with a given constant value
  constantFlux(const Real& a_fluxVal):
    m_fluxVal(a_fluxVal), m_isValSet(true){;};

  virtual ~constantFlux() {}

  virtual SurfaceFlux* new_surfaceFlux();

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);

  /// set flux value in meters/year
  void setFluxVal(const Real& a_fluxVal);

protected:
  Real m_fluxVal;

  bool m_isValSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const constantFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  constantFlux(const constantFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

typedef constantFlux ConstantData;


#include "NamespaceFooter.H"


#endif // multiple-include preventer
