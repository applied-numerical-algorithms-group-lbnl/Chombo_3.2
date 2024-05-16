#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// concrete class encapsulating surface fluxes determined  
// by copying, coarsening, or interpolating a LevelData<FArrayBox>
// covering an entire domain. Objects can be defined either
// by specifying a RefCountedPtr<LevelData<FArrayBox> > , or by specifying
// a std::map<Real,string> mapping time ti to an hdf5 file f. In the
// latter case, the flux at time t is found by linear interploation
// between time max(ti <= t) and min(ti > t) 

// \todo replace the std::map<Real,string> mechanism with
// a suitable abstraction

#ifndef _ISMIP6OCEANFORCING_H_
#define _ISMIP6OCEANFORCING_H_

#include "SurfaceFlux.H"
#include "ParmParse.H"
#include <map>
#include "NamespaceHeader.H"

class ISMIP6OceanForcing : public SurfaceFlux
{
public:

  
  ISMIP6OceanForcing(ParmParse& a_pp);

  /// destructor
  virtual ~ISMIP6OceanForcing() {};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** 
  */

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);


protected:

  std::map<int,std::string> m_year_file;
  int m_start_year, m_end_year;
  std::string m_name;
  std::string m_basin_file, m_basin_var_name;
  std::string m_deltaT_file,  m_deltaT_var_name;
  bool m_anomaly, m_local;
  Real m_factor, m_basin_mean_min_thickness;
  bool m_tf_is_source;
  

  //- 1.0 * melt rate on a uniform mesh
  RefCountedPtr<LevelData<FArrayBox> > m_uniform_source;
  // mesh spacing of m_uniform_source
  RealVect m_dx;
  // timestamp for
  int m_uniform_source_year;

  // layer spacing of m_TF (between components)
  Real m_dz;
  int m_n_layer;

  // basin masks: mask(x,y,n) = basin_number(x,y) == n ? 1 : 0  
  RefCountedPtr<LevelData<FArrayBox> > m_basin_mask;
  int m_n_basin; // number of basins
  
  // deltaT(x,y)
  RefCountedPtr<LevelData<FArrayBox> > m_deltaT; // no-need to keep reading this
  
  /// compute source from thermal forcing
  void computeSource(LevelData<FArrayBox>& a_source,
		     LevelData<FArrayBox>& a_TFb,
		     LevelData<FArrayBox>& a_TFb_basin_mean,
		     LevelData<FArrayBox>& a_deltaT,
		     Real a_gamma0);

  /// compute thermal forcing TFb at the ice shelf base from the 3D 'layered' field TF
  void computeTFb(LevelData<FArrayBox>& a_TFb,
		  const LevelData<FArrayBox>& a_TF,
		  const AmrIceBase& a_amrIce);


  /// compute the basin means of TFb, store in TFb_mean
  void computeBasinMeans(LevelData<FArrayBox>&a_TFb_basin_mean,
			 LevelData<FArrayBox>&a_TFb,
			 LevelData<FArrayBox>&a_basin_mask,
			 const AmrIceBase& a_amrIce);
  
  /// update the uniform source
  void updateUniformSource(Real a_time, const AmrIceBase& a_amrIce);

  ///  Read a field TF(x,y,z) and compute source(x,y) 
  void readUniformSource ( RefCountedPtr<LevelData<FArrayBox> >& a_source,
			   RealVect& a_dx,
			   int a_year,
			   const AmrIceBase& a_amrIce);
    
   
};


#include "NamespaceFooter.H"

#endif
