#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRDAMAGE_H_
#define _AMRDAMAGE_H_

#include "AmrIce.H"
#include "LevelSigmaCS.H"
#include "CalvingModel.H"
#include "IceConstants.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

class DamageConstitutiveRelation;
class DamageCalvingModel;
#define DAMAGE_MAX_VAL 0.9
#define DAMAGE_N_GHOST 4
#define DAMAGE_N_COMP 1
#define WATER_N_GHOST 4
#define WATER_N_COMP 1

/// Continuum damage model. 
class AMRDamage
{

  friend class DamageConstitutiveRelation;
  friend class DamageCalvingModel;
  int m_time_step;
  Real m_time;

  /// cell-centered damage 
  //Vector<LevelData<FArrayBox>*> m_damage;
  /// AMR mesh hierarchy
  int m_finestLevel;
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_ratio;
  Vector<RealVect> m_dx;
  
  //storage for the face-centered minimum damage. This
  //actually gets written to by an external object (a DamageConstitutiveRelation), 
  //so, maybe a be of redesign is in order... for now:
  Vector<LevelData<FluxBox>*> m_faceMinDamage;


public:
  ~AMRDamage();
  AMRDamage();

  const LevelData<FArrayBox>* damage(int a_level) const;

  const LevelData<FArrayBox>* water(int a_level) const;
  
  int finestLevel() const
  {
    return m_finestLevel;
  }
  void define(const Vector<DisjointBoxLayout>& a_grids, 
	      const Vector<int>& a_ratio,
	      int a_finestLevel,
	      const RealVect& a_crseDx);

  /// advance by timestep a_dt. a_dt must satisfy CFL
  void timestep(Real a_dt, AmrIce& a_amrIce);

#ifdef CH_USE_HDF5 
  //plot file methods

  /// add the names of variables to add to the plot file to a_var, in the order they will be written
  void addPlotVars(Vector<std::string>& a_vars);

  /// copy level a_level plot data to   LevelData<FArrayBox>& a_data
  void writePlotData(LevelData<FArrayBox>& a_data, int a_level);

  //check point methods
  
  /// fill a_var with the names of variables to add to the checkpoint file 
  void addCheckVars(Vector<std::string>& a_vars);
  
  /// write level a_level checkpoint data
  void writeCheckData(HDF5Handle& a_handle, int a_level);

  
  /// read level a_level checkpoint data
  void readCheckData(HDF5Handle& a_handle,  HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids);
#endif

protected:
  //cell-centered damage
  mutable Vector<LevelData<FArrayBox>* > m_damage;
  //cell-centered 'water depth' parameter, which affects crevasse depth calculations.
  mutable Vector<LevelData<FArrayBox>* > m_water;
  
private:
  
  void computeDamageSource(Vector<LevelData<FArrayBox>* >& a_source, 
		     const Vector<LevelData<FArrayBox>* >& a_damage,
		     const Vector<RefCountedPtr<LevelSigmaCS> >& geometry,
		     const Vector<LevelData<FArrayBox>const *>& a_visTensor,
		     const Vector<LevelData<FArrayBox>const * >& a_surfThickSource,
		     const Vector<LevelData<FArrayBox>const * >& a_baseThickSource,
		     Real a_dt);

  void computeWaterSource(Vector<LevelData<FArrayBox>* >& a_source,  Real a_dt, AmrIce& a_amrIce);

  void updateWater(Vector<LevelData<FArrayBox>* >& a_water,
		   const Vector<LevelData<FluxBox>* >& a_faceFlux,
		   const Vector<LevelData<FluxBox>* >& a_faceVel,
		   const Vector<LevelData<FArrayBox>* >& a_source,
		   Real a_dt);

  void computeFlux( Vector<LevelData<FluxBox>* >& a_faceFlux,
		    const Vector<LevelData<FArrayBox>* >& a_damage,
		    const Vector<LevelData<FArrayBox>* >& a_source,
		    const Vector<LevelData<FArrayBox>* >& a_cellVel, 
		    const Vector<LevelData<FluxBox>* >& a_faceVel,
		    Real a_dt);

  void updateDamage(Vector<LevelData<FArrayBox>* >& a_damage,
		    const Vector<LevelData<FluxBox>* >& a_faceFlux,
		    const Vector<LevelData<FluxBox>* >& a_faceMinDamage,
		    const Vector<LevelData<FluxBox>* >& a_faceVel,
		    const Vector<LevelData<FArrayBox>* >& a_source,
		    const Vector<RefCountedPtr<LevelSigmaCS> >& geometry,
		    const Vector<LevelData<FArrayBox> const * >& a_vt,
		    Real a_dt);
  
  void initField(Vector<LevelData<FArrayBox>* >& a_field,
		 int a_nComp,
		 const IntVect& a_ghost,
		 bool a_setToScalar = false,
		 const Real& a_scalar = 1.12345678e+300);

  void regridField(Vector<LevelData<FArrayBox>* >& a_field);

   
  /// read checkpoint data for a single LevelData<FArrayBox> 
  LevelData<FArrayBox>* readCheckData(
				      const std::string& a_varLabel,
				      const std::string& a_dataLabel,
				      HDF5Handle& a_handle,
				      HDF5HeaderData&  a_header,
				      const DisjointBoxLayout& a_grids); 
};
  
  

/// Contains an AMRDamage, responds to changes in the ice sheet state
class DamageIceObserver : public AmrIce::Observer 
{
  AMRDamage* m_damagePtr;

public:

  DamageIceObserver();
  
  ~DamageIceObserver();

  AMRDamage& damage() const;
  
  /// General purpose notify method: do something in response to a signal from an AmrIce
  void notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce);

  // plot file methods 
  
  /// add the names of variables to add to the plot file to a_var, in the order they will be written
  void addPlotVars(Vector<std::string>& a_vars);

  /// copy level a_level plot data to   LevelData<FArrayBox>& a_data
  void writePlotData(LevelData<FArrayBox>& a_data, int a_level);

  //check point methods
  
  /// fill a_var with the names of variables to add to the checkpoint file 
  void addCheckVars(Vector<std::string>& a_vars);
  
  /// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
  void writeCheckData(HDF5Handle& a_handle, int a_level);
  
  /// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
  void readCheckData(HDF5Handle& a_handle,  HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids);

};


/// PhysIBC for the damage model
class DamagePhysIBC : public PhysIBC
{

public:

  DamagePhysIBC()
  {
  }

  void define(const ProblemDomain& a_domain,
	      const Real&          a_dx);
  
  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  PhysIBC* new_physIBC();

  /// Set up initial conditions
  void initialize(LevelData<FArrayBox>& a_U);

  /// Set boundary primitive values
  void primBC(FArrayBox&            a_WGdnv,
	      const FArrayBox&      a_Wextrap,
	      const FArrayBox&      a_W,
	      const int&            a_dir,
	      const Side::LoHiSide& a_side,
	      const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one-sided difference
     approximations.  If this function doesn't change them, they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
		     const Real&      a_time);
  
  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  void artViscBC(FArrayBox&       a_F,
		 const FArrayBox& a_U,
		 const FArrayBox& a_divVel,
		 const int&       a_dir,
		 const Real&      a_time);   
};


class DamageCalvingModel : public CalvingModel
{

  CalvingModel* m_calvingModelPtr;
  AMRDamage* m_damageModelPtr;
public:

  DamageCalvingModel(CalvingModel* a_calvingModelPtr, AMRDamage* a_damageModelPtr):
    m_calvingModelPtr(a_calvingModelPtr),  m_damageModelPtr(a_damageModelPtr)
  {
    
  }

  DamageCalvingModel(const DamageCalvingModel& a)
  {
    m_calvingModelPtr = a.m_calvingModelPtr->new_CalvingModel();
    m_damageModelPtr = a.m_damageModelPtr;
  }
  
  virtual CalvingModel* new_CalvingModel()
  {
    return new DamageCalvingModel( *this);
  }
  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      CalvingModel::Stage a_stage);
  

  
};


#include "NamespaceFooter.H"

#endif
