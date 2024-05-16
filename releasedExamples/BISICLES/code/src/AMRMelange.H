#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRMELANGE_H_
#define _AMRMELANGE_H_

#include "AmrIce.H"
#include "LevelSigmaCS.H"
#include "SurfaceFlux.H"
#include "NamespaceHeader.H"

#define MELANGE_N_GHOST 1
#define MELANGE_N_COMP 1

/// Continuum melange transport model.
/**
   
   The melange transport model is inended for coupling
   with ocean and/or iceberg model. Ice removed from the ice sheet 
   through calving, either by flowing across the front or
   by larger scale ice shelf collapse is transferred to 
   a melange thickness field defined on the same domain.
   To enable the model, set the ParmParse option

   - `melange_model = true`

   The ocean / iceberg model is expected to provide
   a melange source as a SurfaceFlux (i.e, much as
   an oceanic melt rate would be supplied to the
   ice sheet). Provided that this field balances the
   source of melange from the ice sheet, the melange
   thickness will not grow over time.

   To specify the surface flux, set the components of 
   melange_model.external_source. For example

   - `melange_model.external_source.type = LevelData`
   - `melange_model.external_source.fileFormat = melange_source.hdf5`
   - `melange_model.external_source.name = melange_source`
   
   One potential problem is that the ocean model domain
   may not entirely cover BISICLES open ocean regions. 
   BISICLES wil normally be run at a finer resolution and
   so will resolve narrower fjords. By setting

   - `melange_model.diffusion_factor > 0.0`

   melange thickness will be trasnported diffusively, 
   with a diffusion coeffcient \f$K = f * \max(1-s, 0) * \max (s-h-b,0)\f$  
   so that there is zero diffusion inside the ice sheet
   and on open land.


 */
class AMRMelange
{

  int m_time_step;
  Real m_time;
  Real m_diffusion_factor;
  
  /// cell-centered melange 
  //Vector<LevelData<FArrayBox>*> m_melange;
  /// AMR mesh hierarchy
  int m_finestLevel;
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_ratio;
  Vector<RealVect> m_dx;
  SurfaceFlux* m_external_source;

public:
  ~AMRMelange();
  AMRMelange();

  const LevelData<FArrayBox>* melange(int a_level) const;
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

  /// increment the melange thickness by a_scale * ice thickness
  void increment(AmrIce& a_amrIce, Real a_scale);

  
#ifdef CH_USE_HDF5 
  //plot file methods

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
#endif

  const LevelData<FArrayBox>* melangeThickness(int lev)
  {
    return m_melange[lev];
  }
  
protected:
  //cell-centered melange
  mutable Vector<LevelData<FArrayBox>* > m_melange;

private:
  
  void computeSource(Vector<LevelData<FArrayBox>* >& a_source,
		     AmrIce& a_amrIce,   Real a_dt);

  void updateMelange(Vector<LevelData<FArrayBox>* >& a_melange,
		     const Vector<LevelData<FArrayBox>* >& a_source,
		     const Vector<RefCountedPtr<LevelSigmaCS> >& geometry,
		     Real a_dt);

};
  
  

/// Contains an AMRMelange, responds to changes in the ice sheet state
class MelangeIceObserver : public AmrIce::Observer 
{
  AMRMelange* m_melangePtr;
  bool  m_next_increment_positive; //used to make sure that m_melangePtr->increment calls occur in pairs
  
public:

  MelangeIceObserver();
  ~MelangeIceObserver();

  AMRMelange& melange() const;
  
  /// General purpose notify method: do something in response to a signal from an AmrIce
  void notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce);

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


/// PhysIBC for the melange model
class MelangePhysIBC : public PhysIBC
{

public:

  MelangePhysIBC()
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



#include "NamespaceFooter.H"

#endif
