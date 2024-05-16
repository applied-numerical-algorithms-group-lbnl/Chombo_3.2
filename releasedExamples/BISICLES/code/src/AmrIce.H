#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef AMR_ICE_H
#define AMR_ICE_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "IceVelocitySolver.H"
#include "RealVect.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "SurfaceFlux.H"
#include "BasalFriction.H"
#include "IceThicknessIBC.H"
#include "IceInternalEnergyIBC.H"
#include "CalvingModel.H"
#include "MuCoefficient.H"
#include "CH_HDF5.H"
#include "DomainDiagnosticData.H"
#include "AmrIceBase.H"
#include "ParmParse.H"
#include <sstream>
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "NamespaceHeader.H"

/// class to manage non-subcycled AMR ice-sheet model
/**
 */
class AmrIce : public AmrIceBase
{

public:
  /// Default constructor
  /** At the moment, there is only one constructor, which defines itself
      based on an inputs file.  (This may change, of course)
  */
  AmrIce();

  /// destructor
  /** destructor
   */
  virtual ~AmrIce();

  /// set default values before initialization
  void setDefaults();

  /// initializes object based on inputs data passed in through parmParse
  /** initializes new object based on data in ParmParse database.  This 
      may change to a more explicit initialization if it's thought to be
      a good idea.
  */
  void initialize();

  /// set domain size (in meters)
  void setDomainSize(const RealVect& a_domainSize)
  { m_domainSize = a_domainSize;}
    

  /// set constitutive relation
  /** needs to be done before calling initialize
   */
  void setConstitutiveRelation(ConstitutiveRelation* a_constRelPtr) 
  {m_constitutiveRelation = a_constRelPtr->getNewConstitutiveRelation(); }

  //const ConstitutiveRelation* constitutiveRelation() const
  //{return m_constitutiveRelation;}

   /// set rate factor
  /** needs to be done before calling initialize
   */
  void setRateFactor(RateFactor* a_rateFactorPtr)
  {m_rateFactor = a_rateFactorPtr->getNewRateFactor(); }


   /// set rate factor at base
  /** needs to be done before calling initialize
   */
  void setBasalRateFactor(RateFactor* a_rateFactorPtr)
  {m_basalRateFactor = a_rateFactorPtr->getNewRateFactor(); }
  

  /// set basal Friction relation
  /** needs to be done before calling initialize
   */
  void setBasalFrictionRelation(BasalFrictionRelation* a_basalFrictionRelPtr) 
  {m_basalFrictionRelation = a_basalFrictionRelPtr->getNewBasalFrictionRelation(); }


  /// set BC for thickness advection
  void setThicknessBC( IceThicknessIBC* a_thicknessIBC);
  
  /// set BC for internalEnergy advection
  void setInternalEnergyBC( IceInternalEnergyIBC* a_internalEnergyIBC);

  /// set surfact flux for ice sheet
  void setSurfaceFlux(SurfaceFlux* a_surfaceFluxPtr)
  {
    if (m_surfaceFluxPtr != NULL)
      {
	delete(m_surfaceFluxPtr);
	m_surfaceFluxPtr = NULL;
      }
    m_surfaceFluxPtr = a_surfaceFluxPtr->new_surfaceFlux(); 
  }
  
  /// set basal flux for ice sheet
  void setBasalFlux(SurfaceFlux* a_basalFluxPtr)
  {
    if (m_basalFluxPtr != NULL)
      {
	delete(m_basalFluxPtr);
	m_basalFluxPtr = NULL;
      }
    m_basalFluxPtr = a_basalFluxPtr->new_surfaceFlux(); 
  }

  /// set surface Heat flux for ice sheet
  void setSurfaceHeatBoundaryData(SurfaceFlux* a_surfaceHeatBoundaryDataPtr, bool a_dirichlett, bool a_temperature)
  {
    if (m_surfaceHeatBoundaryDataPtr != NULL)
      {
	delete(m_surfaceHeatBoundaryDataPtr);
	m_surfaceHeatBoundaryDataPtr = NULL;
      }
   
    m_surfaceBoundaryHeatDataDirichlett = a_dirichlett;
    m_surfaceBoundaryHeatDataTemperature = a_temperature;

    if (a_surfaceHeatBoundaryDataPtr != NULL)
      {
	m_surfaceHeatBoundaryDataPtr = a_surfaceHeatBoundaryDataPtr->new_surfaceFlux(); 
      }
  }
  


  /// set basal Heat flux for ice sheet. Does not include the heat flux due to sliding friction
  void setBasalHeatBoundaryData(SurfaceFlux* a_basalHeatBoundaryDataPtr)
  {
    if (m_basalHeatBoundaryDataPtr != NULL)
      {
	delete(m_basalHeatBoundaryDataPtr);
	m_basalHeatBoundaryDataPtr = NULL;
      }
    if (a_basalHeatBoundaryDataPtr != NULL)
      {
	m_basalHeatBoundaryDataPtr = a_basalHeatBoundaryDataPtr->new_surfaceFlux(); 
      }
  }
  
  


  /// set bedrock flux for ice sheet
  void setTopographyFlux(SurfaceFlux* a_topographyFluxPtr)
  {
    if (m_topographyFluxPtr != NULL)
      {
	delete(m_topographyFluxPtr);
	m_topographyFluxPtr = NULL;
      }
    m_topographyFluxPtr = a_topographyFluxPtr->new_surfaceFlux(); 
  }

  /// set basal friction
  void setBasalFriction(const BasalFriction* a_basalFrictionPtr)
  {m_basalFrictionPtr = a_basalFrictionPtr->new_basalFriction();}
   
  void setMuCoefficient(const MuCoefficient* a_muCoefficientPtr)
  {m_muCoefficientPtr = a_muCoefficientPtr->new_muCoefficient();}

  void setCalvingModel(CalvingModel* a_calvingModelPtr)
  {
    if (m_calvingModelPtr) delete m_calvingModelPtr;
    m_calvingModelPtr =  a_calvingModelPtr->new_CalvingModel();
  }

  
  /// advance solution until either max_time or max_step are reached
  void run(Real a_max_time, int a_max_step);

  /// compute one timestep
  void timeStep(Real a_dt);

  // compute face-centered half-time ice thickness using ppm
  void computeH_half(Vector<LevelData<FluxBox>* >& a_H_half, Real a_dt);

  // compute thickness fluxes
  void computeThicknessFluxes(Vector<LevelData<FluxBox>* >& a_vectFluxes,
                              const Vector<LevelData<FluxBox>* >& a_H_half,
                              const Vector<LevelData<FluxBox>* >& a_faceVelAdvection);

  // update ice thickness *and* bedrock elevation
  void updateGeometry(Vector<RefCountedPtr<LevelSigmaCS> >& a_vect_coordSys_new, 
                       Vector<RefCountedPtr<LevelSigmaCS> >& a_vectCoords_old, 
                       const Vector<LevelData<FluxBox>* >& a_vectFluxes, 
                       Real a_dt);
  
  

  /// do regridding
  void regrid();

  /// diagnostic function -- integrates thickness over domain
  Real computeTotalIce() const;


  /// fill flattened Fortran array of data with ice thickness
  void getIceThickness(Real* a_data_ptr, int* a_dim_info, 
		       Real* a_dew, Real* a_dns) const;

  /// is this object defined and initialized?
  /**
   */
  bool isDefined() const;

#ifdef CH_USE_HDF5

private:
  /// General purpose storage for (small) metadata.
  /**  
       Anything placed here will get written to the plot and 
       checkpoint files on write and read from checkpoint file on restart 
   */
  HDF5HeaderData m_headerData;

  /// write BISICLES specific metadata to HDF5 attributes
  void writeMetaDataHDF5(HDF5Handle& a_handle) const;

public:
  void setHeader(const std::string& a_key, const std::string& a_val)
  {
    m_headerData.m_string[a_key] = a_val;
  }
  void getHeader(const std::string& a_key, std::string& a_val)
  {
    a_val =  m_headerData.m_string[a_key];
  }
  void setHeader(const std::string& a_key, const int& a_val)
  {
    m_headerData.m_int[a_key] = a_val;
  }
  void getHeader(const std::string& a_key, int& a_val)
  {
    a_val = m_headerData.m_int[a_key];
  }
  void setHeader(const std::string& a_key, const Real& a_val)
  {
    m_headerData.m_real[a_key] = a_val;
  }
  
  void getHeader(const std::string& a_key, Real& a_val)
  {
    a_val =  m_headerData.m_real[a_key];
  }

  void writeAMRHierarchyHDF5(HDF5Handle& a_handle, 
				     const Vector<DisjointBoxLayout>& a_grids,
				     const Vector<LevelData<FArrayBox>* > & a_data,
				     const Vector<string>& a_name,
				     const Box& a_domain,
				     const Real& a_dx,
				     const Real& a_dt,
				     const Real& a_time,
				     const Vector<int>& a_ratio,
			     const int& a_numLevels) const;
  
  void writeAMRHierarchyHDF5(const string& filename,
			     const Vector<DisjointBoxLayout>& a_grids,
			     const Vector<LevelData<FArrayBox>* > & a_data,
			     const Vector<string>& a_name,
			     const Box& a_domain,
			     const Real& a_dx,
			     const Real& a_dt,
			     const Real& a_time,
			     const Vector<int>& a_ratio,
			     const int& a_numLevels) const;
 
  /// write hdf5 plotfile to the standard location
  /** writes plotfile to the standard location; filename is
      <plot_root>.<step>.<DIM>d.hdf5
   */
  void writePlotFile();
  void writeAMRPlotFile();
  /// write checkpoint file out for later restarting
  /** writes checkpoint data to the standard location; filename is
      <check_root>.<step>.<DIM>d.hdf5
   */
  void writeCheckpointFile() ;

  /// write checkpoint file out for later restarting
  /** writes checkpoint data to the location specified by a_file
   */
  void writeCheckpointFile(const string& a_file) ;

  /// read checkpoint file for restart 
  /** read checkpoint data from file pointed to by a_handle
   */
  void readCheckpointFile(HDF5Handle& a_handle);

  /// set up for restart
  void restart(const string& a_restart_file);

  



private:

  /// store CF output data on a uniform mesh 
  LevelData<FArrayBox> m_uniform_cf_data;
  /// mesh level for cf data
  int m_cf_level;
  /// (CF) names for fields to be written 
  Vector<std::string> m_uniform_cf_data_name;
  /// (CF) standard names for fields to be written (some fields share cf name) 
  Vector<std::string> m_uniform_cf_standard_name;
   /// Long names for fields to be written 
  Vector<std::string> m_uniform_cf_long_name;
  /// Units for fields to be written 
  Vector<std::string> m_uniform_cf_units;
  /// function/interval pair (maybe should store as std::pairs) used to evaluate cf data on a given AMR level
  Vector<std::function<const LevelData<FArrayBox>*(int, LevelData<FArrayBox>&) > >m_cf_field_function;
  // the function prototype include a LevelData<FArrayBox>& that can be used for storage by the implementation
  
  LevelData<FArrayBox> m_cf_tmp_data;
  Vector<Interval> m_cf_field_interval;
  Real m_cf_dt;
  Real m_cf_start_time;
  bool m_plot_style_cf, m_plot_style_amr;
  enum OutputFileNumbering { time_step = 0, time_seconds = 1, time_years = 2, time_yyyymmdd_360 = 3, NUM_PLOT_NUMBERING};
  OutputFileNumbering m_output_file_numbering;

  /// a number that represents the current time
  long outputNumbering(); 

  friend class DomainDiagnosticData;
  DomainDiagnosticData m_cf_domain_diagnostic_data;

  /// set options for plotfile,checkpoint IO

  void setOutputOptions(ParmParse& a_pp);
  /// write hdf5 plotfile with CF names and units (see http://cfconventions.org/)
  /** 
      flush CF data to the standard location; filename is
      <cfplot_root>.<time>.<DIM>d.hdf5. Data is taken from the 
      accumulated (time-mean or time sum data) managed by calls, 
      and is reset once written
      
  */
  void flushCFData();
  /// accumlate CF (http://cfconventions.org/) data to be written by writeCFPlotFile()
  /**
     
     Allows a time-mean (or sum) of CF (http://cfconventions.org/) data to 
     be written.

     set m_uniform_CF_Data += dt / (m_plot_time) * uniform_CF_Data 

     call a_reset - true to set m_uniform_CF_Data <- 0
     
  */
  void accumulateCFData(Real a_dt, bool a_reset = false);

  /// set up various run-long data for CF Output
  void initCFData();
public:
  
#endif

  /// compute total ice and total grounded ice (volume, at this point)
  //Real computeTotalGroundedIce() const;

  /// compute volume over flotation
  //Real computeVolumeAboveFlotation() const;


  /// fill a_area with per-cell area fraction given by a_maskval
  /**
     e.g when a_maskVal == GROUNDEDMASKVAL, a_area will be
     filled with values from 0 (not grounded ice) to 1 (entirely grounded ice)
   */
  void computeAreaFraction(LevelData<FArrayBox>& a_area,
			   int a_maskVal, int a_level) const;

  /// compute area of ice with grounded/floating status given by a_maskval
  // Real computeArea(int a_maskVal) const;
  
  /// compute area of grounded ice.
  //Real computeGroundedArea() const;

  /// compute area of floating ice.
  //Real computeFloatingArea() const;

  /// compute flux over ice
  //Real computeFluxOverIce(const Vector<LevelData<FArrayBox>*> a_flux);

  /// compute change in volume over ice
  //Real computeDeltaVolumeOverIce() const;

  /// compute flux domain
  //Real computeTotalFlux(const Vector<LevelData<FArrayBox>*> a_flux);


  enum velSolverTypes { Picard = 0,
			JFNK = 1,
			KnownVelocity = 2,
			PetscNLSolver = 3,
			FASMGAMR = 4,
			Python = 5,
			InverseVerticallyIntegrated = 6,
                        NUM_SOLVER_TYPES};

  void setTime(Real a_time)
  {
    m_time = a_time - m_offsetTime;
  }

  /// read only access members implementing AMRIceBase
  Real time() const
  {
    return m_time + m_offsetTime;
  };

  Real dt() const
  {
    return m_dt;
  }

  int finestLevel() const
  {
    return m_finest_level;
  }

  const RealVect dx(int a_level) const
  {
    return m_vect_coordSys[a_level]->dx();
  };
  
  const RefCountedPtr<LevelSigmaCS> geometry(int a_level) const 
  {
    return m_vect_coordSys[a_level];
  };
  const DisjointBoxLayout& grids(int a_level) const
  {
    return m_amrGrids[a_level];
  }
  const Vector<DisjointBoxLayout>& grids() const
  {
    return m_amrGrids;
  }
  const Vector<RefCountedPtr<LevelSigmaCS> >& amrGeometry() const
  {
    return m_vect_coordSys;
  }

  /// nearly direct access to the ice thickness 
  void incrementIceThickness(Vector<LevelData<FArrayBox>*> a_thk);
  
  const Vector<int>& refRatios() const { return m_refinement_ratios;}

  const Vector<Real>& amrDx() const {return m_amrDx;}

  const LevelData<FArrayBox>* velocity(int a_level) const
  {
    return m_velocity[a_level];
  };

  const Vector<LevelData<FArrayBox>* >& amrVelocity() 
  {
    return m_velocity;
  }

  SurfaceFlux& surfaceHeatBoundaryData() const
  {
    if (m_surfaceHeatBoundaryDataPtr == NULL)
      {
	CH_assert(m_surfaceHeatBoundaryDataPtr != NULL);
	MayDay::Error("no surface heat boundary data");
      }
    return *m_surfaceHeatBoundaryDataPtr;
  }

  bool surfaceHeatBoundaryDirichlett() const
  {
    return m_surfaceBoundaryHeatDataDirichlett;
  }

  /// return pointer to real-valued ice fraction on level a_level
  /** (0 = no ice, 1 = full of ice, 0<phi<1 -- partially filled 
      (useful for calving model and for communication with POP  **/
  LevelData<FArrayBox>* iceFrac(int a_level) 
  {
    return m_iceFrac[a_level];
  };

  bool surfaceHeatBoundaryTemperature() const
  {
    return m_surfaceBoundaryHeatDataTemperature;
  }

  SurfaceFlux& basalHeatBoundaryData() const
  {
    if (m_basalHeatBoundaryDataPtr == NULL)
      {
	CH_assert(m_basalHeatBoundaryDataPtr != NULL);
	MayDay::Error("no surface heat boundary data");
      }
    return *m_basalHeatBoundaryDataPtr;
  }

  /// return pointer to real-valued ice fraction on level a_level
  /** (0 = no ice, 1 = full of ice, 0<phi<1 -- partially filled 
      (useful for calving model and for communication with POP  **/
  const LevelData<FArrayBox>* iceFrac(int a_level) const
  {
    return m_iceFrac[a_level];
  };
  /// return const ref to real-valued 'mu coefficient' on level a_level
  //** useful for calving models **/
  const LevelData<FArrayBox>& muCoef(int a_level) const
  {
    return *m_cellMuCoef[a_level];
  }
  
  /// return reference to full AMR ice fraction
  Vector<LevelData<FArrayBox>* >& amrIceFrac() 
  {
    return m_iceFrac;
  }

  /// return const reference to full AMR ice fraction
  const Vector<LevelData<FArrayBox>* >& amrIceFrac() const
  {
    return m_iceFrac;
  }

  const LevelData<FArrayBox>* groundingLineProximity(int a_level) const;

  Real groundingLineProximityScale() const
  {
    return m_groundingLineProximityScale;
  }
  

  const LevelData<FArrayBox>* surfaceThicknessSource(int a_level) const;
  const LevelData<FArrayBox>* basalThicknessSource(int a_level) const;
  const LevelData<FArrayBox>* viscousTensor(int a_level) const;
  const LevelData<FArrayBox>* dragCoefficient(int a_level) const;
  const LevelData<FArrayBox>* viscosityCoefficient(int a_level) const;

  /// access function for face-centered velocities
  const Vector<LevelData<FluxBox>* >& faceVelocities() {return m_faceVelTotal;}

  /// access function for layered face-centered velocities
  const Vector<LevelData<FluxBox>* >& faceVelocitiesLayered() const {return m_layerXYFaceXYVel;} 

  /// access function for calved ice data
  const  Vector<LevelData<FArrayBox>* >& calvedIceThickness() const {return m_calvedIceThickness; }
  const  Vector<LevelData<FArrayBox>* >& removedIceThickness() const {return m_removedIceThickness; }
  const  Vector<LevelData<FArrayBox>* >& addedIceThickness() const {return m_addedIceThickness; }
  
  /// access function for surface temperature
  const  Vector<LevelData<FArrayBox>* >& surfaceInternalEnergy() const {return m_sInternalEnergy; }

  /// access function for surface heat flux
  const  Vector<LevelData<FArrayBox>* >& surfaceHeatFlux() const {return m_sHeatFlux; }

protected:

  



  /// compute tags for regridding
  void tagCells(Vector<IntVectSet>& a_tags);

  /// compute tags for the level a_level
  void tagCellsLevel(IntVectSet& a_tags, int a_level);
 
  /// compute tags at initial time
  void tagCellsInit(Vector<IntVectSet>& a_tags);

  /// initialize grids at initial time
  void initGrids(int a_finest_level);
  
  /// set up grids from grids file
  void setupFixedGrids(const std::string& a_gridFile);

  /// set up storage, etc on an AMR level
  void levelSetup(int a_level, const DisjointBoxLayout& a_grids);
                       
  /// initialize data on hierarchy
  void initData(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCoordSys,
                Vector<LevelData<FArrayBox>* >& a_initialVelocity);

  /// initialize diagnostics
  //void initDiagnostics();

  /// define nonlinear ellipic solver for computing velocity field
  virtual void defineSolver();

  /// solve for velocity field (or just set up some auxilliary quantities)
  virtual void solveVelocityField(bool a_forceSolve = false, Real a_convergenceMetric = -1.0);

  /// compute RHS for velocity field solve
  virtual void defineVelRHS(Vector<LevelData<FArrayBox>* >& a_vectRhs);

  /// set basal friction coefficients C,C0 prior to velocity solve
  virtual void setBasalFriction(Vector<LevelData<FArrayBox>* >& a_C, Vector<LevelData<FArrayBox>* >& a_C0);

  /// set mu coefficient (phi) prior to velocity solve
  virtual void setMuCoefficient(Vector<LevelData<FArrayBox>* >& a_cellMuCoef);

  /// given the current cell centred velocity field, compute a face centred velocity field
  void computeFaceVelocity(Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
			   Vector<LevelData<FluxBox>* >& a_faceVelTotal,
			   Vector<LevelData<FluxBox>* >& a_diffusivity,
			   Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
			   Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel);


  // /// compute div(vel*H) at a given time
  // void computeDivThicknessFlux(Vector<LevelData<FArrayBox>* >& a_divFlux,
  //                              Vector<LevelData<FluxBox>* >& a_flux,
  //                              Vector<LevelData<FArrayBox>* >& a_thickness,
  //                              Real a_time, Real a_dt);

  // increment phi := phi + dt*dphi
  void incrementWithDivFlux(Vector<LevelData<FArrayBox>* >& a_phi,
                            const Vector<LevelData<FArrayBox>* >& a_dphi,
                            Real a_dt);

  // incremente coordSys with new thickness
  void updateCoordSysWithNewThickness(const Vector<LevelData<FArrayBox>* >& a_thickness);

  /// initialize real-valued ice fraction  
  /** this is entirely based on ice thickness h, as opposed
      to the advective update
   */
  void setIceFrac(const LevelData<FArrayBox>& a_thickness, int a_level);


  /// update real-valued ice fraction after thickness update 
  /** this is different from setIceFrac in that it preserves 
      fractional values, but still sets values to zero in cells 
      with no ice. Thickness is non-const because it may set very small
      values of H (less than epsilon) to zero
   */
  void updateIceFrac(LevelData<FArrayBox>& a_thickness, int a_level);

  void setIceFracBC(LevelData<FArrayBox>& a_iceFrac, const ProblemDomain& a_domain);
  
  /// update real-valued ice fraction through advection from neighboring cells
  void advectIceFrac(Vector<LevelData<FArrayBox>* >& a_iceFrac,
                     const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,  
                     Real a_dt);


  // /// update real-valued ice fraction through advection from neighboring cells
  // /**
  //    Also increment a_thicknessSource based on the calving rate
  // */
  // // void advectIceFrac(Vector<LevelData<FArrayBox>* >& a_iceFrac,
  // // 		     Vector<LevelData<FArrayBox>* >& a_thicknessSource,
  // //                    const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
  // // 		     const Vector<LevelData<FluxBox>* >& a_thicknessFlux,
  // //                    Real a_dt);
  
  // /// update real-valued ice fraction through application of the calving rate
  // /**
  //    Also increment a_thicknessSource based on the calving rate
  // */
  // // void advectIceFracUpstream(Vector<LevelData<FArrayBox>* >& a_iceFrac,
  // // 			     Vector<LevelData<FArrayBox>* >& a_thicknessSource,
  // // 			     const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
  // // 			     const Vector<LevelData<FArrayBox> *>& a_calvingRate,
  // // 			     Real a_subDt, Real a_dt);

  // // /// update real-valued ice fraction through advection from neighboring cell
  // // void advectIceFracDownstream(Vector<LevelData<FArrayBox>* >& a_iceFrac,
  // // 			       Vector<LevelData<FArrayBox>* >& a_thicknessSource,
  // // 			       const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
  // // 			       const Vector<LevelData<FluxBox>* >& a_thicknessFlux,
  // // 			       const Vector<LevelData<FArrayBox> *>& a_calvingRate,
  // // 			       Real a_subDt, Real a_dt);

  /// update covered cells, ghost cells etc, of a_iceFrac
  void updateInvalidIceFrac(Vector<LevelData<FArrayBox> *> a_iceFrac);
  
  /// compute timestep
  Real computeDt();

  /// compute timestep at initial time
  Real computeInitialDt();

  //set an entire amr hierarchy to zero
  void setToZero(Vector<LevelData<FArrayBox>*>& a_data);

  /// implicit solve for diffusive fluxes
  void implicitThicknessCorrection(Real a_dt,
				   const Vector<LevelData<FArrayBox>* >& a_sts,
				   const Vector<LevelData<FArrayBox>* >& a_bts,
				   const Vector<LevelData<FArrayBox>* >& a_vts
				   );
  
  void applyCalvingCriterion(CalvingModel::Stage a_stage);
  
  void eliminateRemoteIce();

  /// solve a_alpha * phi - a_beta * laplacian(phi) = a_rhs
  /// and grad(phi).n = 0 on the domain boundaries
  void helmholtzSolve(Vector<LevelData<FArrayBox>* >& a_phi,
		      const Vector<LevelData<FArrayBox>* >& a_rhs, 
		      Real a_alpha, Real a_beta) const;

  /// solve a_alpha* phi - a_beta * laplacian(phi)  = phi'
  /// and grad(phi).n = 0 on the domain boundaries
  void helmholtzSolve(Vector<LevelData<FArrayBox>* >& a_phi,
		      Real a_alpha, Real a_beta) const;
 

#if BISICLES_Z == BISICLES_LAYERED
 

  /// compute the flow rate coefficient A given bulk internal energy & geometry
  void computeA(Vector<LevelData<FArrayBox>* >& a_A, 
			Vector<LevelData<FArrayBox>* >& a_sA,
			Vector<LevelData<FArrayBox>* >& a_bA,
			const Vector<LevelData<FArrayBox>* >& a_internalEnergy, 
			const Vector<LevelData<FArrayBox>* >& a_sInternalEnergy,
			const Vector<LevelData<FArrayBox>* >& a_bInternalEnergy,
			 const Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys) const;

  //compute the face- and layer- centered internal energy (a_layerEH_half)
  //and thickness (a_layerH_half) at time a_time + 1/2 * a_dt
  void computeInternalEnergyHalf(Vector<LevelData<FluxBox>* >& a_layerEH_half, 
				 Vector<LevelData<FluxBox>* >& a_layerH_half,
				 const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel, 
				 const Real a_dt, const Real a_time);

 
  
  /// Update the internalEnergy and tillWaterDepth fields
  void updateInternalEnergy(Vector<LevelData<FluxBox>* >& a_layerTH_half, 
			    Vector<LevelData<FluxBox>* >& a_layerH_half,
			    const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
			    const Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel,
			    const Real a_dt, const Real a_time, 
			    Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysNew,
			    Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysOld,
			    const Vector<LevelData<FArrayBox>*>& a_surfaceThicknessSource,
			    const Vector<LevelData<FArrayBox>*>& a_basalThicknessSource,
			    const Vector<LevelData<FArrayBox>*>& a_volumeThicknessSource);
  /// update the temperature (a derived field)
  void updateTemperature();
#endif


  static void postInterpolationReFloat(LevelData<FArrayBox>& a_H,
				       const LevelData<FArrayBox>& a_coarseH,
				       const LevelData<FArrayBox>& a_coarseBed,
				       const DisjointBoxLayout a_newDBL,
				       const ProblemDomain& a_domain,
				       int a_refRatio,
				       Real a_seaLevel, 
				       Real a_waterDensity, 
				       Real a_iceDensity);

  // diagnostic routine -- compute discharge
  void computeDischarge(const Vector<LevelData<FluxBox>* >& a_vectFluxes);

  // diagnostic routine 
  //void endTimestepDiagnostics();


  /// temporal accuracy 
  /** 1 = 1st-order (use old-time velocity for flux computation)
      2 = 2nd-order (recompute velocity at half-time for flux computation)
      3 = 3rd-order Runge-Kutta update (if eventually implemented)
      4 = 4th-order Runge-Kutta update
  */
  int m_temporalAccuracy;
  
  /// number of ghost cells we need for m_old_thickness
  /** this is dependent on the advection scheme we're using
      AMRGodunov-based -- 4
      RK4 -- 1
   */
  int m_num_thickness_ghost;

  /// what type of nonlinear solver to use for velocity field
  int m_solverType;

  /// max number of iterations for velocity solver
  int m_maxSolverIterations;

  /// solver tolerance for velSolver.
  Real  m_velocity_solver_tolerance;

  /// solve the velocity problem every m_velocitySolveInterval timesteps
  int m_velocity_solve_interval;

  /// Residual norms at the start and end of the last velocity solve
  Real  m_velocitySolveInitialResidualNorm;
  Real  m_velocitySolveFinalResidualNorm;
  /// nonlinear elliptic solver for velocity field
  IceVelocitySolver* m_velSolver;
  

  /// constitutive relation
  ConstitutiveRelation* m_constitutiveRelation;
  /// rate factor
  RateFactor* m_rateFactor;
  /// basal Friction relation
  BasalFrictionRelation* m_basalFrictionRelation;
  /// rate factor at base
  RateFactor* m_basalRateFactor;

  /// if false, don't reset basal friction for floating ice to zero (default is true)
  /** This option exists primarly so that we can run the CISM circular 
      shelf test problem
  */
  bool  m_reset_floating_friction_to_zero;

  /// pointer to a GodunovPhysics for thickness advection
  GodunovPhysics* m_thicknessPhysPtr;

  /// PatchGodunov used for thickness advection -- one for each level
  Vector<PatchGodunov*> m_thicknessPatchGodVect;

  /// IBC for thickness advection
  IceThicknessIBC* m_thicknessIBCPtr;
 
  /// IBC for internalEnergy advection
  IceInternalEnergyIBC* m_internalEnergyIBCPtr;

  /// Calving model
  CalvingModel* m_calvingModelPtr;

  /// sets surface flux for ice sheet
  SurfaceFlux* m_surfaceFluxPtr;
  /// sets basal flux for ice sheet
  SurfaceFlux* m_basalFluxPtr;
  
   /// sets surface heat flux or internalEnergy for ice sheet
  SurfaceFlux* m_surfaceHeatBoundaryDataPtr;
  /// sets basal heat flux or temperature for ice sheet
  SurfaceFlux* m_basalHeatBoundaryDataPtr;
  
  /// sets topography flux (which raises/lowers bedrock)
  SurfaceFlux* m_topographyFluxPtr;
  /// stores the accumulated bedrock difference. This is needed because when 
  /// m_interpolate_zb == false, the regrid() function re-reads the bedrock, which
  /// in many cases is a DEM (e.g, bedmap2), so that the evolution of topography would be obliterated
  Vector<LevelData<FArrayBox>*> m_deltaTopography;

  /// sets basal friction coefficient (beta^2)
  BasalFriction* m_basalFrictionPtr;

  /// sets mu coefficient (phi)
  MuCoefficient* m_muCoefficientPtr;
  
  /// max number ofSflux levels
  int m_max_level;
  
  /// current finest level
  int m_finest_level;
  
  /// if 0 <  m_finest_timestep_level < m_finest_level and >= 0, 
  /// timestep only levels 0 to m_finest_timestep_level, and interpolate
  /// any levels above that. The idea is to avoid being hammered too hard by 
  /// diminishing stability limits on dt
  int m_finest_timestep_level;
  int finestTimestepLevel() const
  {
    return (0 <= m_finest_timestep_level)?
      std::min(m_finest_timestep_level,m_finest_level):m_finest_level;
  }


  /// maximum level at which refinement tags are computed
  int m_tag_cap;

  /// blocking factor
  int m_block_factor;

  /// only allow tagging within a subset of the coarse domain
  IntVectSet m_tag_subset;

  Vector<IntVectSet> m_vectTagSubset;

  /// grid efficiency
  Real m_fill_ratio;

  /// proper nesting radius
  int m_nesting_radius;

  /// max box size
  int m_max_box_size;

  /// max box size for level 0 (defaults to max_box_size)
  /** this allows a different max box size for the level 0 domain,
      which can be useful since the number of boxes on level 0 is often 
      driven by load balancing concerns, while max_box_size on other levels 
      is generally driven by issues of cache size, etc) */
  int m_max_base_grid_size;

  /// regrid interval
  int m_regrid_interval;
  int m_n_regrids;

  /// if false, don't do initial velocity solve during problem setup. 
  /// (default is true)
  bool m_doInitialVelSolve;

  /// if true, use nonzero initial guess for velocity field (default is false)
  /** 
      if true, make semi-informed initial guess for first pass through 
      ice velocity solver
  */
  bool m_doInitialVelGuess;
  // types of initial guess
  enum initialGuessTypes {SlidingLaw=0, // sliding law with the form Cu = rhs 
			  ConstMu=1, // solve the elliptic equations with a supplied mu
			  Function=2, //use a RealFunction<RealVect> 
			  NUM_INITIAL_GUESS_TYPES};
  int m_initialGuessType;

  //paramamters which specify the linear solver in the case m_initialGuessType = ConstMu
  Real m_initialGuessConstMu; // domain wide viscosity
  int m_initialGuessSolverType; // solver type, defaults to the Picard solver
  RealVect m_initialGuessConstVel; // initial velocity, useful when the basal friction relation is nonlinear

  /// if true, interpolate things like base topography at regrid time 
  /** if false, use IceThicknesIBC::regridIceGeometry function to set values
   */
  bool m_interpolate_zb;
  /// thickness interpolation method to use on regrid, see LevelSigmaCS::ThicknessInterpolationMethod
  int m_regrid_thickness_interpolation_method;

  // if > 0, then the coefficient of friction at the base
  // of the ice and the RHS will be smoothed.
  Real m_basalLengthScale;
 
  /// stores tags from previous regrids
  Vector<IntVectSet> m_vectTags; 

  /// tag on grad(velocity)?
  bool m_tagOnGradVel;
  
  /// tagging value (undivided gradient(vel) threshold for regridding)
  Real m_tagging_val;

  /// tag on laplacian(velocity)
  bool m_tagOnLapVel;

  /// tag on laplacian(velocity) in grounded regions only
  bool m_tagOnGroundedLapVel;

  /// tagging threshold value for undivided lap(vel)
  Real m_laplacian_tagging_val;
  /// when tagging on undivided lap(vel), ignore cells where the basal friction coeffcient exceeds 
  Real m_laplacian_tagging_max_basal_friction_coef;
  /// tag on strain rate invariant?
  bool m_tagOnEpsSqr;

  /// tagging value for strain rate invariant
  /** tag if epsSqr*dx > than this value
   */
  Real m_epsSqr_tagVal;

  /// tag on velocity solve RHS
  bool m_tagOnVelRHS;

  /// tagging value for velocity RHS
  Real m_velRHS_tagVal;

  /// tag on div(H grad (vel)) 
  bool m_tagOndivHgradVel;

  /// tagging threshold value for undivided div(H grad(vel))
  Real m_divHGradVel_tagVal;

  /// tag the grounding line
  bool m_tagGroundingLine;
  /// threshold velocity at the grounding line : tag if |u| > value
  Real m_groundingLineTaggingMinVel;
  /// threshold Basal friction coefficient at the grounding line : tag if C < value
  Real m_groundingLineTaggingMaxBasalFrictionCoef;
  
  ///tag  where |vel| * dx > m_maxVelDxVal
  bool m_tagVelDx;
  Real m_velDx_tagVal;
  int m_velDx_tagVal_finestLevelGrounded;
  int m_velDx_tagVal_finestLevelFloating;

  // tag thin ice shelf cavities. Refinining thin cavities 
  // sometimes reveals pinning points
  bool m_tag_thin_cavity;
  Real m_tag_thin_cavity_thickness;

  /// tag at the ice margin
  bool m_tagMargin;
  int m_margin_tagVal_finestLevel;

  /// tag any cell with ice in it
  bool m_tagAllIce;

  /// tag any cell with ice in it on level 0 
  /** (useful when running a case like Antarctica or Greenland where
      the useful domain is subset of the computational domain) */
  bool m_tagAllIceOnLevel0;
  
  /// tag all cells for refinement (refine entire domain)
  bool m_tagEntireDomain;


  //attempt to eliminate floating ice that is unconnected to grounded ice 
  bool m_eliminate_remote_ice;
  int m_eliminate_remote_ice_max_iter;
  Real m_eliminate_remote_ice_tol;

  //attempt to eliminate unconnected floating ice after regridding
  bool m_eliminate_remote_ice_after_regrid;
  
#ifdef HAVE_PYTHON
  /// tag cells according to a python function
  bool m_tagPython;
  PyObject* m_tagPythonModule;
  PyObject* m_tagPythonFunction;
#endif
  /// amount to buffer tags used in regridding
  // grow tags in all directions by m_tags_grow
  int m_tags_grow; 
  //grow tags in all direction dir by max(0,m_tags_grow_dir[dir] -  m_tags_grow)
  IntVect m_tags_grow_dir; 

  /// refinement ratios
  Vector<int>  m_refinement_ratios;

  /// cell spacing at each level
  Vector<Real> m_amrDx;

  /// problem domains at each level
  Vector<ProblemDomain> m_amrDomains;
  
  // problem domain size
  RealVect m_domainSize;

  /// current grids
  Vector<DisjointBoxLayout> m_amrGrids;

  /// keeps track of which levels are completely covered
  Vector<int> m_covered_level;

  /// book-keeping; keeps track of number of cells per level
  Vector<int> m_num_cells;

  /// current time
  Real m_time; Real m_offsetTime;

  /// most recent timestep
  Real m_dt;

  /// most recent stable timestep (> m_dt if we took a step just to hit a plot interval)
  Real m_stable_dt;

  /// timestep scaling
  Real m_cfl;

  /// cfl number for initial timestep (useful if initial data needs small cfl)
  Real m_initial_cfl;

  /// maximum amount cfl number may grow in one timestep
  Real m_max_dt_grow;

  /// maximum acceptable velocity: if max(u) > m_vel_exit , then dump a plot file and exit
  Real m_velocity_exit;

  // current step
  int m_cur_step;

  // isothermal model ?
  bool m_isothermal;
 
  //depth of water in surface crevasses
  Real m_waterDepth;

  /// should a Dirichlett internal energy condition be imposed at the upper surface?
  bool m_surfaceBoundaryHeatDataDirichlett;
  /// is the  Dirichlett data a temperature (rather than an internal energy)?
  bool m_surfaceBoundaryHeatDataTemperature;

  //ice density
  Real m_iceDensity;

  //sea water density
  Real m_seaWaterDensity;

  // acceleration due to gravity
  Real m_gravity;

  // size of unit time (for velocity & time step *only*, other quantities are SI) in seconds
  Real m_seconds_per_unit_time;
  
  /// current old-time data
  Vector<LevelData<FArrayBox>*> m_old_thickness;

 

  /// horizontal velocity field (found by an IceVelocitySolver)
  Vector<LevelData<FArrayBox>*> m_velocity;

  /// velocity solve RHS (in case we want to save it)
  Vector<LevelData<FArrayBox>* > m_velRHS;

  // coeffient of |u|^m-1 u in the basal traction
  Vector<LevelData<FArrayBox>* > m_velBasalC; 

  //coefficient of mu in the stress tensor
  Vector<LevelData<FArrayBox>* > m_cellMuCoef;
  //Vector<LevelData<FluxBox>* > m_faceMuCoef;

  // thickness diffusion coefficient
  Vector<LevelData<FluxBox>* > m_diffusivity;
  Real m_additionalDiffusivity;
  
  //face velocities
  Vector<LevelData<FluxBox>* > m_faceVelAdvection;
  Vector<LevelData<FluxBox>* > m_faceVelTotal; // m_faceVelTotal = m_faceVelAdvection - m_diffusivity * grad(H)
#if BISICLES_Z == BISICLES_LAYERED
  Vector<LevelData<FArrayBox>* > m_layerSFaceXYVel;
  Vector<LevelData<FArrayBox>* > m_layerSFaceSVel;
  Vector<LevelData<FluxBox>* > m_layerXYFaceXYVel;
#endif

 

  // SigmaCS mapping
  Vector<RefCountedPtr<LevelSigmaCS > > m_vect_coordSys;
  
  /// real-valued ice fraction
  ///  real-valued ice fraction on level a_level
  /** (0 = no ice, 1 = full of ice, 0<phi<1 -- partially filled 
      (useful for calving model and for communication with POP  **/
  Vector<LevelData<FArrayBox>* > m_iceFrac;

  // Cell-centered internalEnergy field
  Vector<LevelData<FArrayBox>* > m_internalEnergy, m_temperature;
  //cell-centered till water depth (accompanies the internal energy)
  Vector<LevelData<FArrayBox>* > m_tillWaterDepth;
  // Cell-centered A field
  Vector<LevelData<FArrayBox>* > m_A;
  // used to decide whether A needs updating (after a regrid, or a internalEnergy update)
  mutable bool m_A_valid;
#if BISICLES_Z == BISICLES_LAYERED
  // InternalEnergy and A at the top and bottom surfaces of 
  // the ice, horizontally cell centred (but at vertical faces)
  Vector<LevelData<FArrayBox>* > m_sInternalEnergy, m_sTemperature ;
  Vector<LevelData<FArrayBox>* > m_bInternalEnergy, m_bTemperature;
  Vector<LevelData<FArrayBox>* > m_sHeatFlux;
  Vector<LevelData<FArrayBox>* > m_bHeatFlux;
  Vector<LevelData<FArrayBox>* > m_sA;
  Vector<LevelData<FArrayBox>* > m_bA;
#elif BISICLES_Z == BISICLES_FULL
  // possibly won't even need surface and basal temperature
  // field, but if we do they will need to be 2D
#endif
  //the thickness source will have three components, from the surface and base
  //and distributed throught the volume.
  //Although they are combined for thickness evolution, they must
  //be kept seperate for the internalEnergy calculation. 
  mutable Vector<LevelData<FArrayBox>* > m_surfaceThicknessSource;
  mutable Vector<LevelData<FArrayBox>* > m_basalThicknessSource;
  mutable Vector<LevelData<FArrayBox>* > m_volumeThicknessSource;

  /// modify the source terms (surface,basal,volume) to keep all or part of the ice sheet/shelf steady
  void setStableSources(FArrayBox& a_sts,FArrayBox& a_bts, FArrayBox& a_vts,
			const FArrayBox& a_divuh, const BaseFab<int>& a_mask,
			const Box& a_box) const;
  
  Vector<LevelData<FArrayBox>* > m_divThicknessFlux; // set to -dH/dt
  //thickness of calved ice 
  Vector<LevelData<FArrayBox>* > m_calvedIceThickness;
  Vector<LevelData<FArrayBox>* > m_removedIceThickness;
  Vector<LevelData<FArrayBox>* > m_addedIceThickness;
  Vector<LevelData<FArrayBox>* > m_calvedThicknessSource;


  //cell centered storage for the stress-balance equation coefficients
  mutable Vector<LevelData<FArrayBox>* > m_dragCoef;
  mutable Vector<LevelData<FArrayBox>* > m_viscosityCoefCell;

  //cell centered storage for the viscous tensor components
  mutable Vector<LevelData<FArrayBox>* > m_viscousTensorCell;
  //face centered storage for the viscous tensor components
  mutable Vector<LevelData<FluxBox>* > m_viscousTensorFace;
  // cache validity flag
  mutable bool m_viscousTensor_valid;
  void updateViscousTensor() const;

 
  

  //Scalar field that varies in the ice shelf from 0 at the grounding line
  mutable Vector<LevelData<FArrayBox>* > m_groundingLineProximity;
  // cache validity flag;
  mutable bool m_groundingLineProximity_valid;
  // scale length for grounding line proximity
  Real m_groundingLineProximityScale;
  int m_groundingLineProximityCalcType;
  void updateGroundingLineProximity() const;
  
  /// if true, report discharge of solid ice from the ice edge
  bool m_report_discharge;

  /// report diagnostic every m_report_time_interval years (rather than every time step)
  Real m_report_time_interval;
  Real m_next_report_time;

  /// type of basal friction distribution
  int m_beta_type;

  /// background value for beta
  Real m_betaVal;

  /// small parameter for use in computing beta, if needed
  Real m_betaEps;

  /// background slope of bottom topography
  /** for various reasons (mostly having to do with periodicity)
      it's simpler to just keep this as a parameter here in AmrIce, 
      where it's added to the grad(z_surf) when it comes time to 
      compute the RHS for the velocity solve, rather than include it 
      in the sigmaCS definition.

      The way this is implemented, the basal height at a realVect
      location X is defined as 
      
      z_b = sigmaCS.getTopography() +basalSlope[0]*X[0] +basalSlope[1]*X[1]
  */
  RealVect m_basalSlope;


  /// how verbose should we be?
  static int s_verbosity;

  /// is this object initialized?
  bool m_is_defined;
  
  /// 
  bool m_do_restart;

  /// if starting from a restart, timestep of restart
  int m_restart_step;

  /// 
  string m_plot_prefix;
  
  ///
  string m_check_prefix;

  bool m_check_overwrite; // overwrite check point files
  bool m_check_exit; // exit immediately after writing checkpoint?


  /// 
  int m_plot_interval;  
  /// write a plot file every m_plot_time_interval years (rather than steps)
  Real m_plot_time_interval;

  /// if true, reduce the number of fields written to plot files
  bool m_reduced_plot;

  /// if true, include mask in plotfiles
  bool m_write_mask;

  /// if true, include dH/dt in plotfiles
  bool m_write_dHDt;

  // if true, include (cell-averaged) flux velocity in plotfiles
  bool m_write_fluxVel;

  // if true, include (cell-averaged) viscous tensor components in plotfiles
  bool m_write_viscousTensor;

  // if true, include base velocity in plotfiles
  bool m_write_baseVel;

  /// if true, include velocity solver RHS in plotfiles
  bool m_write_solver_rhs;

  //if true, include internalEnergy in plot files
  bool m_write_internal_energy;

  //if true, include thickness source terms in plot files
  bool m_write_thickness_sources;

  //if true, include fields required for ISMIP6 in plot files
  bool m_write_ismip6;

  /// if true, write out .map.hdf5 file
  bool m_write_map_file;
  
  /// if true, write out plotfile before doing velocity solves
  /** this can be convenient when working with the solverBenchmark
   */
  bool m_write_presolve_plotfiles;

  /// 
  int m_check_interval;

  // if true, allow thickness to change in time
  bool m_evolve_thickness; 

  // if true, allow ice frac to change in time
  bool m_evolve_ice_frac;
  
  //if true, allow velocity field to change in time
  bool m_evolve_velocity;

  // if true, bedrock topogography evolves rather than surface in grounded ice
  bool m_evolve_topography_fix_surface;

  //if true, keep floating or grounded ice stable;
  bool m_grounded_ice_stable ;
  bool m_floating_ice_stable ;
  //if true, use the basal flux to set dhdt
  bool m_floating_ice_basal_flux_is_dhdt;
  bool m_floating_ice_basal_flux_is_min_dhdt;
  bool m_grounded_ice_basal_flux_is_dhdt;

  /** if true, scale thickness sources (surface flux and basal flux) by 
      iceFrac values (allows for effective partial cells)
  */
  bool m_frac_sources;


  enum DiffusionTreatment {NONE,IMPLICIT,EXPLICIT};
  DiffusionTreatment m_diffusionTreatment;
 
  // if true, time step are always an integer power of two
  bool m_timeStepTicks;

  //ised to store a fixed time step > 0 , if desired
  Real m_fixed_dt;

  // if true, report sum of grounded ice (volume) -- default is false
  bool m_reportGroundedIce;




#if BISICLES_Z == BISICLES_LAYERED
  //additional members used in poor-mans multidim mode
public:
  const Vector<Real>& getFaceSigma() const
  {
    return m_faceSigma;
  }
protected:
  void setLayers(const Vector<Real>& a_sigma);


  bool m_sigmaSet;
  
  // Number layers within the ice sheet 
  unsigned int m_nLayers;

  //arrangment of layers within the bedrock
  Vector<Real> m_faceSigma;
  
  //if true, include layer velocities in plot files
  bool m_write_layer_velocities;
  bool m_additionalVelocity;

#endif
private: 
  template <class T>
  void levelAllocate(LevelData<T>** a_ptrPtr, const DisjointBoxLayout& a_grids, 
		     int a_nComp, const IntVect& a_ghost)
  {
    CH_assert(a_ptrPtr != NULL);
    if (*a_ptrPtr != NULL)
      delete *a_ptrPtr;
    *a_ptrPtr = new LevelData<T>(a_grids, a_nComp, a_ghost);
  }


public:
  /// part of Observer-Observee implementation
  /** 
      Similar to the implementation in Chombo/lib/src/AMREllitic/Multigrid.H
   */
  class Observer {

    friend class AmrIce;

    AmrIce* m_observee;

    // setObservee is private so that only friends (AmrIce) can access
    void setObservee(AmrIce* a_amrIce)
    {
      m_observee = a_amrIce;
    }

    // clearObservee is private so that only friends (AmrIce) can access
    void clearObservee()
    {
      m_observee = NULL;
    }

  public:

    enum Notification {PostInitGrids,PreVelocitySolve,PostVelocitySolve,PostGeometryUpdate,PreCalving,PostCalving};
    
    virtual ~Observer()
    {
      if (m_observee != NULL)
	{
	  m_observee->removeObserver(this);
	}
    }

    /// General purpose notify method
    virtual void notify(const Notification, AmrIce&)=0;


    /// Fat interface members.
    /** For now, these allow observers to add data
	to the plot and checkpoint files, and read checkpoint data. This
	might not be the most graceful design, but it fits with the 
	way plot and checkpoint files have been done up till now 
    */
#ifdef CH_USE_HDF5
    /// add the names of variables to add to the plot file to a_var, in the order they will be written
    virtual void addPlotVars(Vector<std::string>& a_vars){;}
    
    /// copy level a_level plot data to   LevelData<FArrayBox>& a_data
    virtual void writePlotData(LevelData<FArrayBox>& a_data, int a_level){;}

    /// fill a_var with the names of variables to add to the checkpoint file 
    virtual void addCheckVars(Vector<std::string>& a_vars){;}
    
    /// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
    virtual void writeCheckData(HDF5Handle& a_handle, int a_level){;}

    /// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
    virtual void readCheckData(HDF5Handle& a_handle, HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids){;}
#endif
  };

  void addObserver(Observer* a_observer)
  {
    // add a_observer to the list so long as it isn't there already
    if (std::find(m_observers.begin(), m_observers.end(), a_observer) == m_observers.end())
      {
	m_observers.push_back(a_observer);
	a_observer->setObservee(this);
      }
  }

  void removeObserver(Observer* a_observer)
  {
    // remove observer from the list
    std::vector<Observer*>::iterator it = std::find(m_observers.begin(), m_observers.end(), a_observer);
    if ( it != m_observers.end())
      {
	a_observer->clearObservee();
	m_observers.erase(it);
      }
  }


private:

  std::vector<Observer*> m_observers; //needs to be a std::vector for begin(), end()

  /// notify observers that something or other has occured
  void notifyObservers(const Observer::Notification a_n)
  {
    for (int i = 0; i < m_observers.size(); i++)
      m_observers[i]->notify(a_n, *this);
  }

};  
  
#include "NamespaceFooter.H"

#endif


