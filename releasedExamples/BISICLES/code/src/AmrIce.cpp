#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>

using std::ifstream; 
using std::ios;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::string;
#include "BISICLES_VERSION.H"
#include "Box.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "ParmParse.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "parstream.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "FineInterp.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "MayDay.H"
#include "AmrIce.H"
#include "computeNorm.H" 
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "PiecewiseLinearFillPatch.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "DerivativesF_F.H"
#include "DivergenceF_F.H"
#include "computeSum.H"
#include "CONSTANTS.H"
#include "IceConstants.H"
#include "ExtrapBCF_F.H"
#include "amrIceF_F.H"
#include "BisiclesF_F.H"
#include "IceThermodynamics.H"
#include "JFNKSolver.H"
#include "InverseVerticallyIntegratedVelocitySolver.H"
#include "PetscIceSolver.H"
#include "RelaxSolver.H"
#ifdef CH_USE_FAS
#include "FASIceSolver.H"
#endif
#include "KnownVelocitySolver.H"
#include "VCAMRPoissonOp2.H"
#include "AMRPoissonOpF_F.H"
#include "CH_HDF5.H"
#include "IceUtility.H"
#include "LevelMappedDerivatives.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif

#include "NamespaceHeader.H"

// small parameter defining when times are equal
#define TIME_EPS 1.0e-12

int AmrIce::s_verbosity = 1;
      
/// fill flattened Fortran array of data with ice thickness
void
AmrIce::getIceThickness(Real* a_data_ptr, int* a_dim_info, 
			Real* a_dew, Real* a_dns) const
{
  // dimInfo is (SPACEDIM, nz, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering
  IntVect hiVect(D_DECL(a_dim_info[2]-1,a_dim_info[3]-1, a_dim_info[1]-1));
  Box fabBox(IntVect::Zero, hiVect);

  FArrayBox exportHfab(fabBox, 1, a_data_ptr); 

  // now pack this into a LevelData  -- this will need to be expanded in MPI
  Vector<Box> exportBoxes(1,fabBox);
  Vector<int> procAssign(1,0);
  // ignore periodicity, since we don't have ghost cells
  DisjointBoxLayout exportGrids(exportBoxes, procAssign);
  LevelData<FArrayBox> exportLDF(exportGrids, 1);

  // this isn't correct in 3d...
  CH_assert(SpaceDim != 3);
  RealVect exportDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // assume that dx = dy, at least for now
  CH_assert (exportDx[0] == exportDx[1]);

  // start at level 0, then work our way up to finest level, 
  // over-writing as we go. An optimzation would be to check 
  // to see if finer levels cover the entire domain...
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      const LevelSigmaCS& levelCS = *(m_vect_coordSys[lev]);
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      const LevelData<FArrayBox>& levelThickness = levelCS.getH();

      // refinement ratio
      Real refRatio = exportDx[0]/m_amrDx[lev];

      Real tolerance = 1.0e-6;

      if (refRatio > 1.0 + tolerance)
        {
          // current level finer than export level -- average solution
          // onto output
          int nRef = (int)(refRatio + tolerance);

          CoarseAverage averager(levelGrids,exportGrids,
                                 1, nRef);
          averager.averageToCoarse(exportLDF, levelThickness);
        }
      else if (refRatio < 1.0-tolerance)
        {
          // current level coarser than export level -- interpolate solution
          int nRef = (int)(1.0/refRatio + tolerance);
          
          // FineInterp needs a problem domain
          ProblemDomain exportDomain(m_amrDomains[lev]);
          exportDomain.refine(nRef);
                                     

          FineInterp interpolator(exportGrids, 1, nRef, exportDomain);
          interpolator.interpToFine(exportLDF, levelThickness);
        }
      else
        {
          // same resolution -- just do a copyTo
          levelThickness.copyTo(exportLDF);
        }

    } // end loop over levels

  // now copy to input fab
  DataIterator exportDit = exportLDF.dataIterator();
  for (exportDit.begin(); exportDit.ok(); ++exportDit)
    {
      // in parallel, this should only be on proc 0
      exportHfab.copy(exportLDF[exportDit]);
    }
  
}

bool 
AmrIce::isDefined() const
{
  return m_is_defined;
}


AmrIce::AmrIce() : m_velSolver(NULL),
                   m_constitutiveRelation(NULL),
		   m_rateFactor(NULL),
		   m_basalFrictionRelation(NULL),
		   m_basalRateFactor(NULL),
                   m_thicknessPhysPtr(NULL),
                   m_thicknessIBCPtr(NULL), 
                   m_surfaceFluxPtr(NULL),
		   m_basalFluxPtr(NULL),
		   m_surfaceHeatBoundaryDataPtr(NULL),
		   m_basalHeatBoundaryDataPtr(NULL),
		   m_topographyFluxPtr(NULL),
                   m_basalFrictionPtr(NULL),
		   m_calvingModelPtr(NULL)
		   
{
  setDefaults();
}

void
AmrIce::setDefaults()
{
  m_sigmaSet = false;
  // set some bogus values as defaults 
  m_is_defined = false;
  m_max_level = -1;
  m_finest_level = -1;
  m_finest_timestep_level = -1;
  m_tag_cap = 100;
  m_block_factor = -1;
  m_fill_ratio = -1;
  m_do_restart = false;
  m_restart_step = -1;
  //  m_constitutiveRelation = NULL;
  m_solverType = JFNK;
  // at the moment, 1 is the only one which works
  m_temporalAccuracy = 1;
  m_num_thickness_ghost = 4;
  // default is -1, which means use the solver's own defaults
  m_maxSolverIterations = -1;
  
  m_velocity_solver_tolerance = 1e-10;

  //by default, solve the full velocity problem on every timestep
  m_velocity_solve_interval = 1;
  //m_velSolver = NULL;
  m_domainSize = -1*RealVect::Unit;

  //m_beta_type = constantBeta;
  m_betaVal = 1000.0;
  m_betaEps = 0.0;
  m_basalSlope = RealVect::Zero;
  
  m_interpolate_zb = true;
  m_regrid_thickness_interpolation_method = 0;

  // set the rest of these to reasonable defaults
  m_nesting_radius = 1;
#ifdef CH_USE_PETSC
  m_nesting_radius = 3;
#endif
  m_tagOnGradVel = false;
  m_tagging_val = 0.1;
  m_tagOnLapVel = false;
  m_tagOnGroundedLapVel = false;
  m_laplacian_tagging_val = 1.0;
  m_laplacian_tagging_max_basal_friction_coef = 1.2345678e+300;
  m_tagOnEpsSqr = false;  
  m_epsSqr_tagVal =0.1;
  m_tagOnVelRHS = false;
  m_velRHS_tagVal = 1.0;
  m_tagOndivHgradVel = false;
  m_divHGradVel_tagVal = 1.0;
  m_tagGroundingLine  = false;
  m_tagVelDx = false;
  m_velDx_tagVal = 0.0;
  m_velDx_tagVal_finestLevelGrounded = -1;
  m_velDx_tagVal_finestLevelFloating = -1;
  m_tagMargin  = false;
  m_margin_tagVal_finestLevel = -1;
  m_tagAllIce  = false;
  m_tagAllIceOnLevel0 = false;  
  m_tagEntireDomain = false;
  m_groundingLineTaggingMinVel = 200.0;
  m_groundingLineTaggingMaxBasalFrictionCoef = 1.2345678e+300;
  m_tag_thin_cavity = false;
  m_tag_thin_cavity_thickness = TINY_THICKNESS;
#ifdef HAVE_PYTHON
  m_tagPython = false;
  m_tagPythonModule = NULL;
  m_tagPythonFunction = NULL;
#endif

  m_tags_grow = 1;
  m_tags_grow_dir = IntVect::Zero;
  m_cfl = 0.25;
  m_max_dt_grow = 1.5;
  m_dt = 1.0e20;
  m_stable_dt = m_dt;
  m_velocity_exit = HUGE_VEL;
  m_max_box_size = 64;
  m_max_base_grid_size = -1;
  m_isothermal = true;
  m_waterDepth = 0.0;
  m_surfaceBoundaryHeatDataDirichlett = true;
  m_surfaceBoundaryHeatDataTemperature = true;
  m_iceDensity = ICE_DENSITY; 
  m_seaWaterDensity = SEA_WATER_DENSITY;
  m_gravity = GRAVITY;
  m_seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
  
  m_report_discharge = false;
  m_report_time_interval = 0.01;
  m_eliminate_remote_ice = false;
  m_eliminate_remote_ice_max_iter = 10;
  m_eliminate_remote_ice_tol = 1.0;
  m_eliminate_remote_ice_after_regrid = false;

  m_diffusionTreatment = NONE;
  m_additionalDiffusivity = 0.0;
  m_additionalVelocity = false;
  m_timeStepTicks = false;
  m_fixed_dt  = 0.0;

  m_velocitySolveInitialResidualNorm = -1.0;
  m_doInitialVelSolve = true; 
  m_doInitialVelGuess = false; 
  m_initialGuessType = SlidingLaw;
  m_initialGuessConstMu = 1.0e+7; // a number that seems plausible, only needed
                                  // if m_initialGuessType = ConstMu
  m_initialGuessSolverType = JFNK; 
  m_initialGuessConstVel = RealVect::Zero; // only needed if m_initialGuessType = ConstMu *and* the basal traction relation is nonlinear
  m_reset_floating_friction_to_zero = true; // set basal friction to zero where ice is floating
 
  m_basalLengthScale = 0.0; // don't mess about with the basal friction / rhs by default
 
  m_evolve_thickness = true;
  m_evolve_ice_frac = true;
  m_evolve_velocity = true;
  m_evolve_topography_fix_surface = false;
  m_grounded_ice_stable = false;
  m_floating_ice_stable = false;
  m_floating_ice_basal_flux_is_dhdt = false;
  m_floating_ice_basal_flux_is_min_dhdt = false;
  m_grounded_ice_basal_flux_is_dhdt = false;
  m_frac_sources = false;

  m_groundingLineProximityScale = 1.0e+4;
  m_groundingLineProximityCalcType = 0 ; // default to the old (odd) behaviour
  //cache validity flags
  m_A_valid = false;
  m_groundingLineProximity_valid = false;  m_viscousTensor_valid = false;

  zeroFlux* cfptr = new zeroFlux;
  m_basalFluxPtr = cfptr;
  m_offsetTime = 0.0;

}

AmrIce::~AmrIce()
{
  if (s_verbosity > 4)
    {
      pout() << "AmrIce::~AmrIce()" << endl;
    }

  // clean up memory
  for (int lev=0; lev<m_velocity.size(); lev++)
    {
      if (m_velocity[lev] != NULL)
        {
          delete m_velocity[lev];
          m_velocity[lev] = NULL;
        }
    }
 
  // clean up ice fraction
  for (int lev=0; lev<m_iceFrac.size(); lev++)
    {
      if (m_iceFrac[lev] != NULL)
        {
          delete m_iceFrac[lev];
          m_iceFrac[lev] = NULL;
        }
    }
 


  // clean up velocityRHS storage if appropriate
  for (int lev=0; lev<m_velRHS.size(); lev++)
    {
      if (m_velRHS[lev] != NULL)
        {
          delete m_velRHS[lev];
          m_velRHS[lev] = NULL;
        }
    }
 
  // clean up basal C storage if appropriate
  for (int lev=0; lev < m_velBasalC.size(); lev++)
    {
      if (m_velBasalC[lev] != NULL)
        {
          delete m_velBasalC[lev];
          m_velBasalC[lev] = NULL;
        }
    }

  // clean up cell centered mu coefficient storage if appropriate
  for (int lev=0; lev < m_cellMuCoef.size(); lev++)
    {
      if (m_cellMuCoef[lev] != NULL)
        {
          delete m_cellMuCoef[lev];
          m_cellMuCoef[lev] = NULL;
        }
    }

  // clean up face advection velocity storage if appropriate
  for (int lev=0; lev < m_faceVelAdvection.size(); lev++)
    {
      if (m_faceVelAdvection[lev] != NULL)
        {
          delete m_faceVelAdvection[lev];
          m_faceVelAdvection[lev] = NULL;
        }
    }
  
  // clean up face total velocity storage if appropriate
  for (int lev=0; lev < m_faceVelTotal.size(); lev++)
    {
      if (m_faceVelTotal[lev] != NULL)
        {
          delete m_faceVelTotal[lev];
          m_faceVelTotal[lev] = NULL;
        }
    }

  for (int lev=0; lev < m_diffusivity.size(); lev++)
    {
      if (m_diffusivity[lev] != NULL)
	{
	  delete m_diffusivity[lev];
	  m_diffusivity[lev] = NULL;
	}
    }

  // clean up surface thickness storage if appropriate
  for (int lev=0; lev < m_surfaceThicknessSource.size(); lev++)
    {
      if (m_surfaceThicknessSource[lev] != NULL)
	{
	  delete m_surfaceThicknessSource[lev];
	  m_surfaceThicknessSource[lev] = NULL;
	}
      if (m_basalThicknessSource[lev] != NULL)
	{
	  delete m_basalThicknessSource[lev];
	  m_basalThicknessSource[lev] = NULL;
	}
        if (m_volumeThicknessSource[lev] != NULL)
	{
	  delete m_volumeThicknessSource[lev];
	  m_volumeThicknessSource[lev] = NULL;
	}
      
      if (m_calvedIceThickness[lev] != NULL)
	{
	  delete m_calvedIceThickness[lev];
	  m_calvedIceThickness[lev] = NULL;
	}
      if (m_removedIceThickness[lev] != NULL)
	{
	  delete m_removedIceThickness[lev];
	  m_removedIceThickness[lev] = NULL;
	}
      if (m_addedIceThickness[lev] != NULL)
	{
	  delete m_addedIceThickness[lev];
	  m_addedIceThickness[lev] = NULL;
	}
    

    }
  
  for (int lev=0; lev < m_divThicknessFlux.size(); lev++)
    {
      if (m_divThicknessFlux[lev] != NULL)
	{
	  delete m_divThicknessFlux[lev];
	  m_divThicknessFlux[lev] = NULL;
	}
    }

#if BISICLES_Z == BISICLES_LAYERED
  for (int lev=0; lev < m_layerSFaceXYVel.size(); lev++)
    {
      if (m_layerSFaceXYVel[lev] != NULL)
        {
          delete m_layerSFaceXYVel[lev];
          m_layerSFaceXYVel[lev] = NULL;
        }
    }
  
  for (int lev=0; lev < m_layerSFaceSVel.size(); lev++)
    {
      if (m_layerSFaceSVel[lev] != NULL)
        {
          delete m_layerSFaceSVel[lev];
          m_layerSFaceSVel[lev] = NULL;
        }
    }

  for (int lev=0; lev < m_layerXYFaceXYVel.size(); lev++)
    {
      if (m_layerXYFaceXYVel[lev] != NULL)
        {
          delete m_layerXYFaceXYVel[lev];
          m_layerXYFaceXYVel[lev] = NULL;
        }
    }

#endif
  

  // clean up internalEnergy  storage if appropriate
  for (int lev=0; lev < m_internalEnergy.size(); lev++)
    {
      if (m_internalEnergy[lev] != NULL)
        {
          delete m_internalEnergy[lev];
          m_internalEnergy[lev] = NULL;
        }
    }

 // clean up internalEnergy  storage if appropriate
  for (int lev=0; lev < m_temperature.size(); lev++)
    {
      if (m_temperature[lev] != NULL)
        {
          delete m_temperature[lev];
          m_temperature[lev] = NULL;
        }
    }
  
  for (int lev=0; lev < m_A.size(); lev++)
    {
      if (m_A[lev] != NULL)
        {
          delete m_A[lev];
          m_A[lev] = NULL;
        }
    }
#if BISICLES_Z == BISICLES_LAYERED

  for (int lev=0; lev < m_sA.size(); lev++)
    {
      if (m_sA[lev] != NULL)
        {
          delete m_sA[lev];
          m_sA[lev] = NULL;
        }
    }
 for (int lev=0; lev < m_tillWaterDepth.size(); lev++)
    {
      if (m_tillWaterDepth[lev] != NULL)
	{
	  delete m_tillWaterDepth[lev];
	  m_tillWaterDepth[lev] = NULL;
	}
    }
  for (int lev=0; lev < m_sInternalEnergy.size(); lev++)
    {
      if (m_sInternalEnergy[lev] != NULL)
	{
	  delete m_sInternalEnergy[lev];
	  m_sInternalEnergy[lev] = NULL;
	}
    }
  for (int lev=0; lev < m_sTemperature.size(); lev++)
    {
      if (m_sTemperature[lev] != NULL)
        {
          delete m_sTemperature[lev];
          m_sTemperature[lev] = NULL;
        }
    }

  
  for (int lev=0; lev < m_sHeatFlux.size(); lev++)
    {
      if (m_sHeatFlux[lev] != NULL)
	{
	  delete m_sHeatFlux[lev];
	  m_sHeatFlux[lev] = NULL;
	}
    }
  for (int lev=0; lev < m_bInternalEnergy.size(); lev++)
    {
      if (m_bInternalEnergy[lev] != NULL)
        {
          delete m_bInternalEnergy[lev];
          m_bInternalEnergy[lev] = NULL;
        }
    }
  for (int lev=0; lev < m_bTemperature.size(); lev++)
    {
      if (m_bTemperature[lev] != NULL)
        {
          delete m_bTemperature[lev];
          m_bTemperature[lev] = NULL;
        }
    }
  for (int lev=0; lev < m_bHeatFlux.size(); lev++)
    {
      if (m_bHeatFlux[lev] != NULL)
	{
	  delete m_bHeatFlux[lev];
	  m_bHeatFlux[lev] = NULL;
	}
    }
  for (int lev=0; lev < m_bA.size(); lev++)
    {
      if (m_bA[lev] != NULL)
	{
	  delete m_bA[lev];
	  m_bA[lev] = NULL;
	}
    }
#endif

  for (int lev = 0; lev < m_old_thickness.size(); lev++) 
    {
      if (m_old_thickness[lev] != NULL) 
        {
          delete m_old_thickness[lev];
          m_old_thickness[lev] = NULL;
        }
    }

  for (int lev = 0; lev < m_groundingLineProximity.size(); lev++) 
    {
      if (m_groundingLineProximity[lev] != NULL) 
        {
          delete m_groundingLineProximity[lev];
          m_groundingLineProximity[lev] = NULL;
        }
    }

  for (int lev = 0; lev < m_dragCoef.size(); lev++) 
    {
      if (m_dragCoef[lev] != NULL) 
        {
          delete m_dragCoef[lev];
          m_dragCoef[lev] = NULL;
        }
    }
  
  for (int lev = 0; lev < m_viscosityCoefCell.size(); lev++) 
    {
      if (m_viscosityCoefCell[lev] != NULL) 
        {
          delete m_viscosityCoefCell[lev];
          m_viscosityCoefCell[lev] = NULL;
        }
    }

  for (int lev = 0; lev < m_viscousTensorCell.size(); lev++) 
    {
      if (m_viscousTensorCell[lev] != NULL) 
        {
          delete m_viscousTensorCell[lev];
          m_viscousTensorCell[lev] = NULL;
        }
    }
  for (int lev = 0; lev < m_viscousTensorFace.size(); lev++) 
    {
      if (m_viscousTensorFace[lev] != NULL) 
        {
          delete m_viscousTensorFace[lev];
          m_viscousTensorFace[lev] = NULL;
        }
    }

  for (int lev = 0; lev < m_deltaTopography.size(); lev++) 
    {
      if (m_deltaTopography[lev] != NULL) 
        {
          delete m_deltaTopography[lev];
          m_deltaTopography[lev] = NULL;
        }
    }

  if (m_velSolver != NULL)
    {
      // code crashes here. comment out until I figure out the problem...
      delete m_velSolver;
      m_velSolver = NULL;
    }

  if (m_constitutiveRelation != NULL)
    {
      delete m_constitutiveRelation;
      m_constitutiveRelation = NULL;
    }

  if (m_rateFactor != NULL)
    {
      delete m_rateFactor;
      m_rateFactor = NULL;
    }

  if (m_basalFrictionRelation != NULL)
    {
      delete m_basalFrictionRelation;
      m_basalFrictionRelation = NULL;
    }
  
  if (m_basalRateFactor != NULL)
    {
      delete m_basalRateFactor;
      m_basalRateFactor = NULL;
    }

  for (int lev=0; lev<m_thicknessPatchGodVect.size(); lev++)
    {
      if (m_thicknessPatchGodVect[lev] != NULL)
	{
	  delete m_thicknessPatchGodVect[lev];
	  m_thicknessPatchGodVect[lev] = NULL;
	}
    }


  if (m_thicknessPhysPtr != NULL)
    {
      delete m_thicknessPhysPtr;
      m_thicknessPhysPtr = NULL;
    }

  if (m_thicknessIBCPtr != NULL)
    {
      delete m_thicknessIBCPtr;
      m_thicknessIBCPtr = NULL;
    }


  if (m_internalEnergyIBCPtr != NULL)
    {
      delete m_internalEnergyIBCPtr;
      m_internalEnergyIBCPtr = NULL;
    }

  if (m_muCoefficientPtr != NULL)
    {
      delete m_muCoefficientPtr;
      m_muCoefficientPtr = NULL;
    }

  if (m_surfaceFluxPtr != NULL)
    {
      delete m_surfaceFluxPtr;
      m_surfaceFluxPtr = NULL;
    }
  if (m_basalFluxPtr != NULL)
    {
      delete m_basalFluxPtr;
      m_basalFluxPtr = NULL;
    }
  if (m_calvingModelPtr != NULL)
    {
      delete m_calvingModelPtr;
      m_calvingModelPtr = NULL;
    }
  if (m_basalFrictionPtr != NULL)
    {
      delete m_basalFrictionPtr;
      m_basalFrictionPtr = NULL;
    }
  
 if (m_surfaceHeatBoundaryDataPtr != NULL)
    {
      delete m_surfaceHeatBoundaryDataPtr;
      m_surfaceHeatBoundaryDataPtr = NULL;
    }
 
  if (m_basalHeatBoundaryDataPtr != NULL)
    {
      delete m_basalHeatBoundaryDataPtr;
      m_basalHeatBoundaryDataPtr = NULL;
    }
  
#ifdef HAVE_PYTHON
  if (m_tagPythonFunction != NULL)
    {
      Py_DECREF(m_tagPythonFunction);
    }
  if (m_tagPythonModule != NULL)
    {
      Py_DECREF(m_tagPythonModule);
    }
#endif

  // that should be it!

}


void
AmrIce::initialize()
{

  CH_TIME("AmrIce::initialize");
  if (s_verbosity > 3) 
    {
      pout() << "AmrIce::initialize" << endl;
    }

  // set time to be 0 -- do this now in case it needs to be changed later
  m_time = 0.0;
  m_cur_step = 0;

  // first, read in info from parmParse file
  ParmParse ppCon("constants");
  ppCon.query("ice_density",m_iceDensity);
  ppCon.query("sea_water_density",m_seaWaterDensity);
  ppCon.query("gravity",m_gravity);
  ppCon.query("seconds_per_unit_time",m_seconds_per_unit_time);
  IceThermodynamics::setConstants(m_iceDensity,m_seaWaterDensity,m_gravity,m_seconds_per_unit_time);
  
  ParmParse ppAmr("amr");
  Vector<int> ancells(3); 
  // allows for domains with lower indices which are not positive
  Vector<int> domLoIndex(SpaceDim,0); 
  // slc : SpaceDim == 2 implies poor-mans multidim mode, in which we still
  // care about the number of vertical layers. 
  Vector<Real> domsize(SpaceDim);

  // assumption is that domains are not periodic
  bool is_periodic[SpaceDim];
  for (int dir=0; dir<SpaceDim; dir++)
    is_periodic[dir] = false;
  Vector<int> is_periodic_int(SpaceDim, 0);

  ppAmr.get("maxLevel", m_max_level);
  ppAmr.query("tagCap",m_tag_cap);
  
  ppAmr.getarr("num_cells", ancells, 0, ancells.size());
  
  // this one doesn't have a vertical dimension
  ppAmr.queryarr("domainLoIndex", domLoIndex, 0, SpaceDim);




# if 0
  // this is now handled in main and passed in using the
  // setDomainSize function
  if (ppAmr.contains("domain_size"))
    {
      ppAmr.getarr("domain_size", domsize, 0, SpaceDim);
      m_domainSize = RealVect(D_DECL(domsize[0], domsize[1], domsize[2]));
    }
  
#endif
  
 

  ppAmr.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++) 
    {
      is_periodic[dir] = (is_periodic_int[dir] == 1);
    }

  ppAmr.query("velocity_exit", m_velocity_exit);

  ppAmr.query("cfl", m_cfl);

  m_initial_cfl = m_cfl;
  ppAmr.query("initial_cfl", m_initial_cfl);

  ppAmr.query("max_dt_grow_factor", m_max_dt_grow);

  ppAmr.query("time_step_ticks", m_timeStepTicks);
  // max_dt_grow must be at least two in this case - or it will never grow
  if (m_timeStepTicks)
    m_max_dt_grow = std::max(m_max_dt_grow, two);

  ppAmr.query("fixed_dt", m_fixed_dt);
  ppAmr.query("offsetTime", m_offsetTime);


  ppAmr.query("isothermal",m_isothermal);

  setOutputOptions(ppAmr);
  

  ppAmr.query("evolve_thickness", m_evolve_thickness);
  ppAmr.query("evolve_topography_fix_surface", m_evolve_topography_fix_surface);
  ppAmr.query("evolve_velocity", m_evolve_velocity);
  ppAmr.query("velocity_solve_interval", m_velocity_solve_interval);

  ppAmr.query("grounded_ice_stable", m_grounded_ice_stable);
  ppAmr.query("floating_ice_stable", m_floating_ice_stable);
  ppAmr.query("floating_ice_basal_flux_is_dhdt", m_floating_ice_basal_flux_is_dhdt);
  ppAmr.query("floating_ice_basal_flux_is_min_dhdt", m_floating_ice_basal_flux_is_min_dhdt);
  CH_assert( ! (m_floating_ice_basal_flux_is_dhdt && m_floating_ice_basal_flux_is_min_dhdt) ); 
  CH_assert( ! (m_floating_ice_basal_flux_is_dhdt && m_floating_ice_stable) );
  CH_assert( ! (m_floating_ice_stable && m_floating_ice_basal_flux_is_min_dhdt) );

  ppAmr.query("grounded_ice_basal_flux_is_dhdt",m_grounded_ice_basal_flux_is_dhdt);
  ppAmr.query("mask_sources", m_frac_sources);

  // there are several options that make no sense with ice fraction evolution
  bool thickness_evolution_inhibited = (!m_evolve_thickness) ||
    m_grounded_ice_stable ||  m_floating_ice_stable || m_floating_ice_basal_flux_is_dhdt ||
    m_floating_ice_basal_flux_is_min_dhdt || m_grounded_ice_basal_flux_is_dhdt;
  m_evolve_ice_frac = ! thickness_evolution_inhibited;
  //but allow the user to modify that anyway
  ppAmr.query("evolve_ice_frac",m_evolve_ice_frac);
  
  ppAmr.query("grounding_line_proximity_scale",m_groundingLineProximityScale);
  ppAmr.query("grounding_line_proximity_calc_type", m_groundingLineProximityCalcType);
  bool tempBool = m_write_dHDt;
  ppAmr.query("write_dHDt", tempBool);
  m_write_dHDt = (tempBool == 1);
  ppAmr.query("write_mask", m_write_mask);
  ppAmr.query("write_solver_rhs", m_write_solver_rhs);

  if (m_max_level > 0)
    {
      //m_refinement_ratios.resize(m_max_level+1,-1);
      ppAmr.getarr("ref_ratio", m_refinement_ratios, 0, m_max_level);
    }
  else
    {
      m_refinement_ratios.resize(1);
      m_refinement_ratios[0] = -1;
    }

  ppAmr.query("verbosity", s_verbosity);
  ppAmr.get("regrid_interval", m_regrid_interval);
  m_n_regrids = 0;
  ppAmr.query("interpolate_zb", m_interpolate_zb);
  ppAmr.query("regrid_thickness_interpolation_method", m_regrid_thickness_interpolation_method);
  ppAmr.get("blockFactor", m_block_factor);
  ppAmr.get("fill_ratio", m_fill_ratio);
  ppAmr.query("nestingRadius", m_nesting_radius);

#ifdef CH_USE_PETSC
  // petsc solvers require nesting radius >= 3
  if (m_nesting_radius < 3)
    {
      MayDay::Warning("PETSC solvers require nesting radius >= 3 -- resetting to 3");
      m_nesting_radius = 3;
    }
#endif

  bool isThereATaggingCriterion = false;
  ppAmr.query("tag_on_grad_velocity", m_tagOnGradVel);
  isThereATaggingCriterion |= m_tagOnGradVel;

  ppAmr.query("tagging_val", m_tagging_val);

  ppAmr.query("tag_on_laplacian_velocity", m_tagOnLapVel);
  isThereATaggingCriterion |= m_tagOnLapVel;

  ppAmr.query("tag_on_grounded_laplacian_velocity", m_tagOnGroundedLapVel);
  isThereATaggingCriterion |= m_tagOnGroundedLapVel;

  // if we set either of these to be true, require that we also provide the threshold
  if (m_tagOnLapVel | m_tagOnGroundedLapVel)
    {
      ppAmr.get("lap_vel_tagging_val", m_laplacian_tagging_val);
      ppAmr.query("lap_vel_tagging_max_basal_friction_coef", m_laplacian_tagging_max_basal_friction_coef);
    }


  ppAmr.query("tag_on_strain_rate_invariant",m_tagOnEpsSqr);
  isThereATaggingCriterion |= m_tagOnEpsSqr;
  // if we set this to be true, require that we also provide the threshold
  if (m_tagOnEpsSqr)
    {
      ppAmr.get("strain_rate_invariant_tagging_val", m_epsSqr_tagVal);
    }


  ppAmr.query("tag_on_velocity_rhs",m_tagOnVelRHS);
  isThereATaggingCriterion |= m_tagOnVelRHS;

  // if we set this to be true, require that we also provide the threshold
  if (m_tagOnVelRHS)
    {
      ppAmr.get("velocity_rhs_tagging_val", m_velRHS_tagVal);
    }
  
  ppAmr.query("tag_grounding_line", m_tagGroundingLine);
  isThereATaggingCriterion |= m_tagGroundingLine;
  // if we set this to be true, require that we also provide the threshold
  if (m_tagGroundingLine)
    {
      ppAmr.get("grounding_line_tagging_min_vel",m_groundingLineTaggingMinVel);
      ppAmr.query("grounding_line_tagging_max_basal_friction_coef", m_groundingLineTaggingMaxBasalFrictionCoef);
    }
  
  
  ppAmr.query("tag_vel_dx", m_tagVelDx);
  isThereATaggingCriterion |= m_tagVelDx;
  // if we set this to be true, require that we also provide the threshold
  if (m_tagVelDx)
    {
      ppAmr.get("vel_dx_tagging_val",m_velDx_tagVal);
      ppAmr.query("vel_dx_finest_level_grounded",m_velDx_tagVal_finestLevelGrounded);
      m_velDx_tagVal_finestLevelFloating = m_velDx_tagVal_finestLevelGrounded;
      ppAmr.query("vel_dx_finest_level_floating",m_velDx_tagVal_finestLevelFloating);
    }

  ppAmr.query("tag_thin_cavity", m_tag_thin_cavity);
  ppAmr.query("tag_thin_cavity_thickness", m_tag_thin_cavity_thickness);

#ifdef HAVE_PYTHON
  ppAmr.query("tag_python", m_tagPython);
  isThereATaggingCriterion |= m_tagPython;
  if (m_tagPython)
    {
      std::string s;
      ppAmr.get("tag_python_module", s);
      PythonInterface::InitializePythonModule
	(&m_tagPythonModule,  s);
      ppAmr.get("tag_python_function",s);
      PythonInterface::InitializePythonFunction
	(&m_tagPythonFunction, m_tagPythonModule , s);
    }
#endif
  
  ppAmr.query("tag_ice_margin", m_tagMargin);
  isThereATaggingCriterion |= m_tagMargin;
  // if we set this to be true, require finest level to tag
  if (m_tagMargin)
    {
      m_margin_tagVal_finestLevel = m_max_level+1;
      ppAmr.query("margin_finest_level",m_margin_tagVal_finestLevel);
    }

  ppAmr.query("tag_all_ice", m_tagAllIce);
  isThereATaggingCriterion |= m_tagAllIce;

  ppAmr.query("tag_all_ice_on_level_0", m_tagAllIceOnLevel0);
  isThereATaggingCriterion |= m_tagAllIceOnLevel0;  

  ppAmr.query("tag_entire_domain", m_tagEntireDomain);
  isThereATaggingCriterion |= m_tagEntireDomain;


  ppAmr.query("tag_on_div_H_grad_vel",m_tagOndivHgradVel);
  isThereATaggingCriterion |= m_tagOndivHgradVel;

  // if we set this to be true, require that we also provide the threshold
  if (m_tagOndivHgradVel)
    {
      ppAmr.get("div_H_grad_vel_tagging_val", m_divHGradVel_tagVal);
    }


  // here is a good place to set default to grad(vel)
  //if ((!m_tagOnGradVel) && (!m_tagOnLapVel) && (!m_tagOnEpsSqr))
  if (!isThereATaggingCriterion)
    {
      m_tagOnGradVel = true;
    }
  
  ppAmr.query("tags_grow", m_tags_grow);
  {
    Vector<int> tgd(SpaceDim,0);
    ppAmr.queryarr("tags_grow_dir", tgd, 0, tgd.size());
    for (int dir =0; dir < SpaceDim; dir++)
      {
	m_tags_grow_dir[dir] = tgd[dir];
      } 
  }
  ppAmr.query("max_box_size", m_max_box_size);

  if (ppAmr.contains("max_base_grid_size") )
    {
      ppAmr.get("max_base_grid_size", m_max_base_grid_size);
    }
  else 
    {
      m_max_base_grid_size = m_max_box_size;
    }

  ppAmr.query("report_discharge", m_report_discharge);
  ppAmr.query("report_time_interval", m_report_time_interval);
  
  ppAmr.query("eliminate_remote_ice", m_eliminate_remote_ice);
  ppAmr.query("eliminate_remote_ice_max_iter", m_eliminate_remote_ice_max_iter);
  ppAmr.query("eliminate_remote_ice_tol", m_eliminate_remote_ice_tol);
  ppAmr.query("eliminate_remote_ice_after_regrid", m_eliminate_remote_ice_after_regrid);

  // get temporal accuracy
  ppAmr.query("temporal_accuracy", m_temporalAccuracy);

  // number of ghost cells depends on what scheme we're using
  if (m_temporalAccuracy < 3)
    {
      m_num_thickness_ghost = 4;
    }
  else 
    {
      m_num_thickness_ghost = 1;      
    }

  // get solver type
  ppAmr.query("velocity_solver_type", m_solverType);
  ppAmr.query("max_solver_iterations",m_maxSolverIterations);
  ppAmr.query("velocity_solver_tolerance", m_velocity_solver_tolerance);
  ppAmr.query("do_initial_velocity_solve", m_doInitialVelSolve);
  ppAmr.query("do_initial_velocity_guess", m_doInitialVelGuess);
  ppAmr.query("initial_velocity_guess_type", m_initialGuessType);
  ppAmr.query("initial_velocity_guess_const_mu", m_initialGuessConstMu);
  ppAmr.query("initial_velocity_guess_solver_type", m_initialGuessSolverType);

  {
    Vector<Real> t(SpaceDim,0.0);
    ppAmr.queryarr("initial_velocity_guess_const_vel", t, 0, SpaceDim);
    m_initialGuessConstVel = RealVect(D_DECL(t[0], t[1], t[2]));
  }

  ppAmr.query("additional_velocity",m_additionalVelocity);

  //thickness diffusion options
  std::string diffusionTreatment = "none";
  ppAmr.query("diffusion_treatment", diffusionTreatment);
  if (diffusionTreatment == "implicit")
    {
      m_diffusionTreatment = IMPLICIT;
    }
  else if (diffusionTreatment == "explicit")
    m_diffusionTreatment = EXPLICIT;
  ppAmr.query("additional_diffusivity",m_additionalDiffusivity);


  //option to advance thickness/internalEnergy only on coarser levels
  ppAmr.query("finest_timestep_level",m_finest_timestep_level);

  ppAmr.query("reset_floating_friction", m_reset_floating_friction_to_zero);
  ppAmr.query("basal_length_scale", m_basalLengthScale);
 
  // now set up problem domains
  {
    IntVect loVect = IntVect(D_DECL(domLoIndex[0], domLoIndex[1], domLoIndex[3]));
    IntVect hiVect(D_DECL(domLoIndex[0]+ancells[0]-1, 
                          domLoIndex[1]+ancells[1]-1, 
                          domLoIndex[2]+ancells[2]-1));
#if BISICLES_Z == BISICLES_LAYERED
    {
      int nLayers = ancells[2];
      Vector<Real> faceSigma(nLayers+1);
      Real dsigma = 1.0 / Real(nLayers);
      for (unsigned int l = 0; l < faceSigma.size(); ++l)
	faceSigma[l] = dsigma * (Real(l));
      {
	ParmParse ppGeo("geometry");
	ppGeo.queryarr("sigma",faceSigma,0,faceSigma.size());
      }
      ppAmr.queryarr("sigma",faceSigma,0,faceSigma.size());
      setLayers(faceSigma);
    }
#endif
    ProblemDomain baseDomain(loVect, hiVect);
    // now set periodicity
    for (int dir=0; dir<SpaceDim; dir++) 
      baseDomain.setPeriodic(dir, is_periodic[dir]);

    // now set up vector of domains
    m_amrDomains.resize(m_max_level+1);
    m_amrDx.resize(m_max_level+1);

    m_amrDomains[0] = baseDomain;
    m_amrDx[0] = m_domainSize[0]/baseDomain.domainBox().size(0);

    for (int lev=1; lev<= m_max_level; lev++)
      {
        m_amrDomains[lev] = refine(m_amrDomains[lev-1],
                                   m_refinement_ratios[lev-1]);
        m_amrDx[lev] = m_amrDx[lev-1]/m_refinement_ratios[lev-1];
      }
  } // leaving problem domain setup scope
  
  std::string tagSubsetBoxesFile = "";
  m_vectTagSubset.resize(m_max_level);
  
  ppAmr.query("tagSubsetBoxesFile",tagSubsetBoxesFile);
  
  if (tagSubsetBoxesFile != "")
    {
      if (procID() == uniqueProc(SerialTask::compute))
  	{
	 
  	  ifstream is(tagSubsetBoxesFile.c_str(), ios::in);
  	  int lineno = 1;
  	  if (is.fail())
  	    {
  	      pout() << "Can't open " << tagSubsetBoxesFile << std::endl;
  	      MayDay::Error("Cannot open refine boxes file");
  	    }

  	  for (int lev = 0; lev < m_max_level; lev++)
  	    {
              // allowable tokens to identify levels in tag subset file
  	      const char level[6] = "level";
              const char domain[7] = "domain";
  	      char s[6];
  	      is >> s;
  	      if (std::string(level) == std::string(s))
  		{
                  int inlev;
                  is >> inlev;
                  if (inlev != lev)
                    {
                      pout() << "expected ' " << lev << "' at line " << lineno << std::endl;
                      MayDay::Error("bad input file");
                    }
                } 
              else if (std::string(domain) == std::string(s))
                {
                  // basic idea here is that we read in domain box
                  // (domains must be ordered from coarse->fine)
                  // until we get to a domain box which matches ours.
                  // This lets us make a single list of subset regions
                  // which we can use for any coarsening/refining of the domain
                  const Box& levelDomainBox = m_amrDomains[lev].domainBox();
                  bool stillLooking = true;
                  while (stillLooking)
                    {
                      Box domainBox;
                      is >> domainBox;
                      if (domainBox == levelDomainBox)
                        {
                          pout() << "Found a domain matching level " << lev << endl;
                          stillLooking = false;
                        }
                      else // move on until we find our level
                        {
                          // read in info for the level we're ignoring
                          //advance to next line
                          while (is.get() != '\n');
			  lineno++;
                          int nboxes;
                          is >> nboxes;
                          if (nboxes > 0)
                            {
                              for (int i = 0; i < nboxes; ++i)
                                {
                                  Box box;
                                  is >> box;
				  while (is.get() != '\n');
				  lineno++;
                                }
                            } 
                          is >> s;
                          if (std::string(domain) != std::string(s))
                            {
                              pout() << "expected '" << domain
                                     << "' at line " << lineno << ", got " 
                                     << s << std::endl;
                              MayDay::Error("bad input file");
                            }                            
                        }
                    }
                }
              else
                {
  		  pout() << "expected '" << level << "' or '" << domain
                         << "' at line " << lineno << ", got " 
                         << s << std::endl;
  		  MayDay::Error("bad input file");
  		}
              //advance to next line
              while (is.get() != '\n');
	      lineno++;
              int nboxes;
              is >> nboxes;
              if (nboxes > 0)
                {
                  for (int i = 0; i < nboxes; ++i)
                    {
                      Box box;
                      is >> box;
		      while (is.get() != '\n');
		      lineno++;
                      m_vectTagSubset[lev] |= box;
                      pout() << " level " << lev << " refine box : " << box << std::endl;
                    }
                }
              //advance to next line
              while (is.get() != '\n');
	      lineno++;
              
              if (lev > 0)
		{
		  //add lower level's subset to this subset
		  IntVectSet crseSet (m_vectTagSubset[lev-1]);
		  if (!crseSet.isEmpty())
		    {
		      crseSet.refine(m_refinement_ratios[lev-1]);
		      // crseSet.nestingRegion(m_block_factor,m_amrDomains[lev]);
		      if (m_vectTagSubset[lev].isEmpty())
			{
			  m_vectTagSubset[lev] = crseSet;
			} 
		      else
			{
			  m_vectTagSubset[lev] &= crseSet;
			} 
		    }
		 
		}
	      
  	    } // end loop over levels

	} // end if serial compute
      for (int lev = 0; lev < m_max_level; lev++)
	broadcast(m_vectTagSubset[lev], uniqueProc(SerialTask::compute));
    }
  /// PatchGodunov used for thickness advection
  if (m_temporalAccuracy < 3)
    {
      // get PatchGodunov options -- first set reasonable defaults.
      // can over-ride from ParmParse
      ParmParse pph ("amr.thickness_advection");
      
      int normalPredOrder = 2;
      pph.query("normal_pred_order",normalPredOrder);
      bool useFourthOrderSlopes = true;
      pph.query("fourth_order_slopes",useFourthOrderSlopes);
      bool usePrimLimiting = true;
      bool useCharLimiting = false;
      bool useFlattening = false;
      Real artificialViscosity = 0.0;
      pph.query("artificial_viscosity",artificialViscosity);
      bool useArtificialViscosity = artificialViscosity > TINY_NORM;     
      
      // define AdvectionPhysics ptr
      // (does this need to be a special Ice-advection pointer due
      // to H factors?      
      
      m_thicknessPhysPtr = new AdvectPhysics;
      m_thicknessPhysPtr->setPhysIBC(m_thicknessIBCPtr);
    
      

      m_thicknessPatchGodVect.resize(m_max_level+1, NULL);

      for (int lev=0; lev<=m_max_level; lev++)
        {


          m_thicknessPatchGodVect[lev] = new PatchGodunov;
          m_thicknessPatchGodVect[lev]->define(m_amrDomains[lev],
                                               m_amrDx[lev],
                                               m_thicknessPhysPtr,
                                               normalPredOrder,
                                               useFourthOrderSlopes,
                                               usePrimLimiting,
                                               useCharLimiting,
                                               useFlattening,
                                               useArtificialViscosity,
                                               artificialViscosity);
        }
      

      //m_internalEnergyIBCPtr = new IceInternalEnergyIBC;
     

    } // end if temporal accuracy < 3
  
 

  // check to see if we're using predefined grids
  bool usePredefinedGrids = false;
  std::string gridFile;
  if (ppAmr.contains("gridsFile"))
    {
      usePredefinedGrids = true;
      ppAmr.get("gridsFile",gridFile);
    }

  
  ParmParse geomPP("geometry");
  if (geomPP.contains("basalSlope") )
    {
      Vector<Real> basalSlope(SpaceDim, 0.0);
      geomPP.getarr("basalSlope", basalSlope, 0, SpaceDim);
      D_TERM(
             m_basalSlope[0] = basalSlope[0];,
             m_basalSlope[1] = basalSlope[1];,
             m_basalSlope[2] = basalSlope[2];);
    }


  // check to see if we're restarting from a checkpoint file
  if (!ppAmr.contains("restart_file"))
    {
      // if we're not restarting
      
      // now set up data holders
      m_old_thickness.resize(m_max_level+1, NULL);
      m_velocity.resize(m_max_level+1, NULL);
      m_iceFrac.resize(m_max_level+1, NULL);
      m_faceVelAdvection.resize(m_max_level+1, NULL);
      m_faceVelTotal.resize(m_max_level+1, NULL);
      m_diffusivity.resize(m_max_level+1);
      m_velBasalC.resize(m_max_level+1,NULL);
      m_cellMuCoef.resize(m_max_level+1,NULL);
      m_velRHS.resize(m_max_level+1, NULL);
      m_surfaceThicknessSource.resize(m_max_level+1, NULL);
      m_basalThicknessSource.resize(m_max_level+1, NULL);
      m_volumeThicknessSource.resize(m_max_level+1, NULL);
      m_calvedIceThickness.resize(m_max_level+1, NULL);
      m_removedIceThickness.resize(m_max_level+1, NULL);
      m_addedIceThickness.resize(m_max_level+1, NULL);
      m_divThicknessFlux.resize(m_max_level+1, NULL);
      m_internalEnergy.resize(m_max_level+1, NULL);
      m_tillWaterDepth.resize(m_max_level+1, NULL);
      m_deltaTopography.resize(m_max_level+1, NULL);
#if BISICLES_Z == BISICLES_LAYERED
      m_layerXYFaceXYVel.resize(m_max_level+1, NULL);
      m_layerSFaceXYVel.resize(m_max_level+1, NULL);
      m_layerSFaceSVel.resize(m_max_level+1, NULL);
      m_sInternalEnergy.resize(m_max_level+1, NULL);
      m_bInternalEnergy.resize(m_max_level+1, NULL);
      m_sHeatFlux.resize(m_max_level+1, NULL);
      m_bHeatFlux.resize(m_max_level+1, NULL);
#endif
      // allocate storage for m_old_thickness,  m_velocity, etc
      for (int lev=0; lev<m_velocity.size(); lev++)
        {
          m_old_thickness[lev] = new LevelData<FArrayBox>;
          m_velocity[lev] = new LevelData<FArrayBox>;
          m_iceFrac[lev] = new LevelData<FArrayBox>;
	  m_faceVelAdvection[lev] = new LevelData<FluxBox>;
	  m_faceVelTotal[lev] = new LevelData<FluxBox>;
	  m_velBasalC[lev] = new LevelData<FArrayBox>;
	  m_cellMuCoef[lev] = new LevelData<FArrayBox>;
	  m_velRHS[lev] = new LevelData<FArrayBox>;
	  m_diffusivity[lev] = new LevelData<FluxBox>;
	  m_internalEnergy[lev] = new LevelData<FArrayBox>;
	  m_tillWaterDepth[lev] = new LevelData<FArrayBox>;
	  m_deltaTopography[lev] = new LevelData<FArrayBox>;
#if BISICLES_Z == BISICLES_LAYERED
	  m_sInternalEnergy[lev] = new LevelData<FArrayBox>;
	  m_bInternalEnergy[lev] = new LevelData<FArrayBox>;
	  m_sHeatFlux[lev] = new LevelData<FArrayBox>;
	  m_bHeatFlux[lev] = new LevelData<FArrayBox>;
	  m_layerXYFaceXYVel[lev] = new LevelData<FluxBox>;
	  m_layerSFaceXYVel[lev] = new LevelData<FArrayBox>;
	  m_layerSFaceSVel[lev] = new LevelData<FArrayBox>;
#endif

        }

      int finest_level = -1;
      if (usePredefinedGrids)
        {
          setupFixedGrids(gridFile);
        } 
      else
        {
          // now create  grids
          initGrids(finest_level);
        }
      
     
    }
  else
    {
      // we're restarting from a checkpoint file
      string restart_file;
      ppAmr.get("restart_file", restart_file);
      m_do_restart = true;
#ifdef CH_USE_HDF5
      restart(restart_file);
#else
      MayDay::Error("restart_file specified but hdf5 support not built"); 
#endif // hdf5

    }

  // last thing to do is to set this to true from here on out...
  m_doInitialVelSolve = true;

  // set up counter of number of cells
  m_num_cells.resize(m_max_level+1, 0);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LayoutIterator lit = levelGrids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& thisBox = levelGrids.get(lit());
          m_num_cells[lev] += thisBox.numPts();
        }
    }


  // finally, set up covered_level flags
  m_covered_level.resize(m_max_level+1, 0);

  // note that finest level can't be covered.
  for (int lev=m_finest_level-1; lev>=0; lev--)
    {

      // if the next finer level is covered, then this one is too.
      if (m_covered_level[lev+1] == 1)
        {
          m_covered_level[lev] = 1;
        }
      else
        {
          // see if the grids finer than this level completely cover it
          IntVectSet fineUncovered(m_amrDomains[lev+1].domainBox());
          const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];

          LayoutIterator lit = fineGrids.layoutIterator();
          for (lit.begin(); lit.ok(); ++lit)
            {
              const Box& thisBox = fineGrids.get(lit());
              fineUncovered.minus_box(thisBox);
            }

          if (fineUncovered.isEmpty()) 
            {
              m_covered_level[lev] = 1;
            }
        }
    } // end loop over levels to determine covered levels

  m_cf_domain_diagnostic_data.initDiagnostics
    (*this, m_vect_coordSys, m_amrGrids, m_refinement_ratios, m_amrDx[0], time() , m_finest_level);

}  
  
/// set BC for thickness advection
void
AmrIce::setThicknessBC( IceThicknessIBC* a_thicknessIBC)
{
  m_thicknessIBCPtr = a_thicknessIBC->new_thicknessIBC(); 
}

/// set BC for internalEnergy advection
void 
AmrIce::setInternalEnergyBC( IceInternalEnergyIBC* a_internalEnergyIBC)
{
  m_internalEnergyIBCPtr = a_internalEnergyIBC->new_internalEnergyIBC();
}

void 
AmrIce::defineSolver()
{
  if (m_solverType == Picard)
    {
      MayDay::Error("PicardSolver is deprecated (for now)");
      
      
    }
  else  if (m_solverType == JFNK)
    {
      // for now, at least, just delete any existing solvers
      // and rebuild them from scratch
     
      JFNKSolver* jfnkSolver;

      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;
      int numLevels = m_finest_level +1;
      
      // make sure that the IBC has the correct grid hierarchy info
      m_thicknessIBCPtr->setGridHierarchy(m_vect_coordSys, m_amrDomains);

      if (m_velSolver != NULL)
	{
	  // assume that any extant solver is also a JFNKSolver
	  jfnkSolver = dynamic_cast<JFNKSolver*>(m_velSolver);
	  CH_assert(jfnkSolver != NULL);
	}
      else {
	jfnkSolver = new JFNKSolver();
      }
      
      jfnkSolver->define(m_amrDomains[0],
			 m_constitutiveRelation,
			 m_basalFrictionRelation,
			 m_amrGrids,
			 m_refinement_ratios,
			 dxCrse,
			 m_thicknessIBCPtr,
			 numLevels);

      m_velSolver = jfnkSolver;

    }
  else  if (m_solverType == KnownVelocity)
    {
      if (m_velSolver != NULL)
        {
          delete m_velSolver;
          m_velSolver = NULL;
        }
      m_velSolver = new KnownVelocitySolver;
      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;
      int numLevels = m_finest_level +1;
      m_velSolver->define(m_amrDomains[0],
                          m_constitutiveRelation,
			  m_basalFrictionRelation,
                          m_amrGrids,
                          m_refinement_ratios,
                          dxCrse,
                          m_thicknessIBCPtr,
                          numLevels);
    }
#ifdef CH_USE_PETSC
  else if (m_solverType == PetscNLSolver)
    {
      // for now, at least, just delete any existing solvers
      // and rebuild them from scratch
      if (m_velSolver != NULL)
        {
          delete m_velSolver;
          m_velSolver = NULL;
        }

      m_velSolver = new PetscIceSolver;
      
      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;

      int numLevels = m_finest_level +1;

      // make sure that the IBC has the correct grid hierarchy info
      m_thicknessIBCPtr->setGridHierarchy(m_vect_coordSys, m_amrDomains);

      m_velSolver->define(m_amrDomains[0],
                          m_constitutiveRelation,
			  m_basalFrictionRelation,
                          m_amrGrids,
                          m_refinement_ratios,
                          dxCrse,
                          m_thicknessIBCPtr,
                          numLevels);
      m_velSolver->setVerbosity(s_verbosity);

      m_velSolver->setTolerance(m_velocity_solver_tolerance);

      if (m_maxSolverIterations > 0)
        {
          m_velSolver->setMaxIterations(m_maxSolverIterations);
        }
    }
#endif
#ifdef CH_USE_FAS
  else if (m_solverType == FASMGAMR)
    {
      // for now, at least, just delete any existing solvers
      // and rebuild them from scratch
      if (m_velSolver != NULL)
        {
          delete m_velSolver;
          m_velSolver = NULL;
        }
      
      FASIceSolver *solver = new FASIceSolver;
      m_velSolver = solver;
      
      solver->setParameters( "FASSolver" );

      RealVect dxCrse = m_amrDx[0]*RealVect::Unit;

      int numLevels = m_finest_level + 1;

      // make sure that the IBC has the correct grid hierarchy info
      m_thicknessIBCPtr->setGridHierarchy( m_vect_coordSys, m_amrDomains );

      solver->define( m_amrDomains[0],
		      m_constitutiveRelation,
		      m_basalFrictionRelation,
		      m_amrGrids,
		      m_refinement_ratios,
		      dxCrse,
		      m_thicknessIBCPtr,
		      numLevels );

      solver->setTolerance( m_velocity_solver_tolerance );

      if (m_maxSolverIterations > 0)
        {
          solver->setMaxIterations( m_maxSolverIterations );
        }
    }
#endif
#ifdef HAVE_PYTHON
  else if (m_solverType == Python)
    {
      // for now, at least, just delete any existing solvers
      // and rebuild them from scratch
      if (m_velSolver != NULL)
        {
          delete m_velSolver;
          m_velSolver = NULL;
        }
      m_velSolver = new PythonInterface::PythonVelocitySolver;
      m_velSolver->define( m_amrDomains[0],
			   m_constitutiveRelation,
			   m_basalFrictionRelation,
			   m_amrGrids,
			   m_refinement_ratios,
			   m_amrDx[0]*RealVect::Unit,
			   m_thicknessIBCPtr,
			   m_finest_level + 1 );
    }
#endif
  else if (m_solverType == InverseVerticallyIntegrated)
    {
      if (m_velSolver != NULL)
        {
	  //not sure if we are OK with rebuilding solvers?
          delete m_velSolver;
          m_velSolver = NULL;
        }
      InverseVerticallyIntegratedVelocitySolver* ptr 
	= new InverseVerticallyIntegratedVelocitySolver;
      
      ptr->define( *this, 
		   m_amrDomains[0],
		   m_constitutiveRelation,
		   m_basalFrictionRelation,
		   m_amrGrids,
		   m_refinement_ratios,
		   m_amrDx[0]*RealVect::Unit,
		   m_thicknessIBCPtr,
		   m_finest_level + 1 );
      m_velSolver = ptr;
      
    }
  else
    {
      MayDay::Error("unsupported velocity solver type");
    }
 
}

//inline 
//Real remainder(Real a, Real b)
//{
//  Real p = a/b; int i(p);
//  return std::min( p - i, p - 1 - i);
//}

void AmrIce::setToZero(Vector<LevelData<FArrayBox>*>& a_data)
{
  for (int lev=0; lev < std::min(int(a_data.size()),finestLevel()+1); lev++)
    {
      LevelData<FArrayBox>& data = *a_data[lev];
      for (DataIterator dit(data.dataIterator()); dit.ok(); ++dit)
	{
	  data[dit].setVal(0.0);
	}
    }
}


void
AmrIce::run(Real a_max_time, int a_max_step)
{

  CH_TIME("AmrIce::run");
  
  if (s_verbosity > 3) 
    {
      pout() << "AmrIce::run -- max_time= " << a_max_time 
             << ", max_step = " << a_max_step << endl;
    }

  Real dt;
  // only call computeInitialDt if we're not doing restart
  if (!m_do_restart)
    {
      dt = computeInitialDt();
    } 
  else
    {
      dt = computeDt();
    }

  // advance solution until done
  if ( !(m_plot_time_interval > TIME_EPS) || m_plot_time_interval > a_max_time) m_plot_time_interval = a_max_time;
  if ( !(m_report_time_interval > TIME_EPS) || m_report_time_interval > a_max_time) m_report_time_interval = a_max_time;
  
  while ( a_max_time > m_time && (m_cur_step < a_max_step))
    {
      Real next_plot_time = m_plot_time_interval * (1.0 + Real(int((m_time/m_plot_time_interval))));
      if ( !(next_plot_time > m_time))
	{
	  //trap case where machine precision results in (effectively)
          // m_plot_time_interval * (1.0 + Real(int((m_time/m_plot_time_interval)))) == m_time
	  next_plot_time += m_plot_time_interval;
	}

      next_plot_time = std::min(next_plot_time, a_max_time); 

      m_next_report_time = m_time;
      m_next_report_time = std::min(m_next_report_time, a_max_time); 


      while ( (next_plot_time > m_time) && (m_cur_step < a_max_step)
	      && (dt > TIME_EPS))
	{

#ifdef CH_USE_HDF5
	  // dump plotfile before regridding  
	  if ( m_plot_interval > 0 && (m_cur_step%m_plot_interval == 0))
	    {
	      writePlotFile();
	    }
	  // dump checkpoint before regridding 
	  if ((m_check_interval > 0) && (m_cur_step%m_check_interval == 0)
	      && (m_cur_step != m_restart_step))
	    {
	      writeCheckpointFile();
	      if (m_cur_step > 0 && m_check_exit)
		{
		  if (s_verbosity > 2)
		    {
		      pout() << "AmrIce::exit on checkpoint" << endl;
		      return;
		    }
		}
	    }	  
#endif	  
	  setToZero(m_calvedIceThickness); 
	  setToZero(m_removedIceThickness);
	  setToZero(m_addedIceThickness);

	  if ((m_cur_step != 0) && (m_cur_step%m_regrid_interval ==0))
	    {
	      regrid();
	    }
	  
	  if (m_cur_step != 0)
	    {
	      // compute dt after regridding in case number of levels has changed
	      dt = computeDt();           
	    }
	  
	  //Real trueDt = dt; //we will need to restore dt if we change it below
	  if (next_plot_time - m_time + TIME_EPS < dt) 
	    dt =  std::max(2 * TIME_EPS, next_plot_time - m_time);
	  
	  timeStep(dt);

	  //m_dt = trueDt; 
	  // restores the correct timestep in cases where it was chosen just to reach a plot file
	  //update CF data mean
	  if (m_plot_style_cf) accumulateCFData(dt);	  
	  
	} // end of plot_time_interval
#ifdef CH_USE_HDF5
      if (m_plot_interval >= 0 && abs(next_plot_time - a_max_time) > TIME_EPS)
	writePlotFile();
	
#endif
    } // end timestepping loop
 
  // dump out final plotfile, if appropriate
#ifdef CH_USE_HDF5
  
  if (m_plot_interval >= 0)
    {
      writePlotFile();
    }


  // dump out final checkpoint file, if appropriate
  if (m_check_interval >= 0)
    {
      writeCheckpointFile();
    }
#endif    
  
  if (s_verbosity > 2)
    {
      pout() << "AmrIce::run finished" << endl;
    }
}


void
AmrIce::timeStep(Real a_dt)
{

  CH_TIME("AmrIce::timestep");
    
  if (s_verbosity >=2) 
    {
      pout() << "Timestep " << m_cur_step 
             << " Advancing solution from time " 
             << m_time << " ( " << time() << ")" " with dt = " << a_dt << endl;
    }

  m_dt = a_dt;

  
  // first copy thickness into old thickness   
  for (int lev=0; lev <= m_finest_level ; lev++)
    {
      
      LevelData<FArrayBox>& oldThickness = *m_old_thickness[lev];
      LevelData<FArrayBox>& currentThickness = m_vect_coordSys[lev]->getH();

      // this way we avoid communication and maintain ghost cells...
      DataIterator dit = oldThickness.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          oldThickness[dit].copy(currentThickness[dit],0, 0, 1);
        }

    }
        
  // assumption here is that we've already computed the current velocity 
  // field, most likely at initialization or at the end of the last timestep...
  // so, we don't need to recompute the velocity at the start.

  // use PatchGodunov hyperbolic solver
  
#if 0 // this doesn't appear to be used anywhere anymore
  // need a grown velocity field
  IntVect grownVelGhost(2*IntVect::Unit);
  Vector<LevelData<FArrayBox>* > grownVel(m_finest_level+1, NULL);
#endif

  // holder for half-time face velocity
  Vector<LevelData<FluxBox>* > H_half(m_finest_level+1,NULL);
  // thickness fluxes 
  Vector<LevelData<FluxBox>* > vectFluxes(m_finest_level+1, NULL);
  
  
  // allocate storage
  for (int lev = finestTimestepLevel() ; lev>=0 ; lev--)
    {
      
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

      IntVect ghostVect = IntVect::Unit;      
      H_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, 
                                           ghostVect);

      // if we're doing AMR, we'll need to average these fluxes
      // down to coarser levels. As things stand now, 
      // CoarseAverageFace requires that the coarse LevelData<FluxBox>
      // have a ghost cell. 
      vectFluxes[lev] = new LevelData<FluxBox>(m_amrGrids[lev],1, ghostVect);

      LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev];
      
      
      
      // ensure that ghost cells for thickness  are filled in
      if (lev > 0)
        {          
          int nGhost = levelOldThickness.ghostVect()[0];
          PiecewiseLinearFillPatch thicknessFiller(levelGrids, 
                                                   m_amrGrids[lev-1],
                                                   1, 
                                                   m_amrDomains[lev-1],
                                                   m_refinement_ratios[lev-1],
                                                   nGhost);
          
          // since we're not subcycling, don't need to interpolate in time
          Real time_interp_coeff = 0.0;
          thicknessFiller.fillInterp(levelOldThickness,
                                     *m_old_thickness[lev-1],
                                     *m_old_thickness[lev-1],
                                     time_interp_coeff,
                                     0, 0, 1);
          
          
          
        }
      // these are probably unnecessary...
      levelOldThickness.exchange();
      
      
      // do we need to also do a coarseAverage for the vel here?
    }
    // compute face-centered thickness (H) at t + dt/2
  computeH_half(H_half, a_dt);
  
  //  compute face- and layer- centered E*H and H  at t + dt/2 (E is internal energy)
  Vector<LevelData<FluxBox>* > layerEH_half(m_finest_level+1,NULL);
  Vector<LevelData<FluxBox>* > layerH_half(m_finest_level+1,NULL);
  if (!m_isothermal)
    computeInternalEnergyHalf(layerEH_half, layerH_half, m_layerXYFaceXYVel, a_dt, m_time); 

  // compute thickness fluxes
  computeThicknessFluxes(vectFluxes, H_half, m_faceVelAdvection);
  
  if (m_report_discharge && (m_next_report_time - m_time) < (a_dt + TIME_EPS))
    {
      m_cf_domain_diagnostic_data.computeDischarge
	(m_vect_coordSys, vectFluxes, m_amrGrids, m_amrDx, m_refinement_ratios, 
	 m_time, time(), m_cur_step, m_finest_level, s_verbosity);
    }

  // update ice fraction through advection
  //advectIceFrac(m_iceFrac, m_volumeThicknessSource, m_faceVelAdvection,vectFluxes, a_dt);
  if (m_evolve_ice_frac)
    {
      advectIceFrac(m_iceFrac, m_faceVelAdvection, a_dt);
    }
  
  // make a copy of m_vect_coordSys before it is overwritten \todo clean up
  Vector<RefCountedPtr<LevelSigmaCS> > vectCoords_old (m_finest_level+1);
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      IntVect sigmaCSGhost = m_vect_coordSys[lev]->ghostVect();
      RealVect dx = m_amrDx[lev]*RealVect::Unit;
      vectCoords_old[lev] = RefCountedPtr<LevelSigmaCS> 
        (new LevelSigmaCS(m_amrGrids[lev], dx, sigmaCSGhost));
      LevelSigmaCS& levelCoords_old = *vectCoords_old[lev];
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      
      
      levelCoords_old.setIceDensity(levelCoords.iceDensity());
      levelCoords_old.setWaterDensity(levelCoords.waterDensity());
      levelCoords_old.setGravity(levelCoords.gravity());
      // todo replace the copies below with a deepCopy of levelCoords
      for (DataIterator dit( m_amrGrids[lev]); dit.ok(); ++dit)
        {
          FArrayBox& oldH = levelCoords_old.getH()[dit];
          const FArrayBox& H = levelCoords.getH()[dit];
          oldH.copy(H);
        }
      levelCoords_old.setTopography(levelCoords.getTopography());
      {
        LevelSigmaCS* crseCoords = (lev > 0)?&(*vectCoords_old[lev-1]):NULL;
        int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
        levelCoords_old.recomputeGeometry( crseCoords, refRatio);
      }
#if BISICLES_Z == BISICLES_LAYERED
      levelCoords_old.setFaceSigma(levelCoords.getFaceSigma());
#endif
    }

  // compute div(F) and update geometry
  updateGeometry(m_vect_coordSys, vectCoords_old, vectFluxes, a_dt);
  
  
  // update internal energy
  if (!m_isothermal)
    {
      // \todo tify up in here: this also updates m_layerSFaceSVel 
      updateInternalEnergy(layerEH_half, layerH_half, m_layerXYFaceXYVel,
			   m_layerSFaceXYVel, a_dt, m_time,
			   m_vect_coordSys, vectCoords_old, 
			   m_surfaceThicknessSource, m_basalThicknessSource, m_volumeThicknessSource);
    }
  
  notifyObservers(Observer::PostGeometryUpdate);
  
  // clean up temporary storage
  for (int lev=0; lev<=m_finest_level; lev++)
    {
          
      if (H_half[lev] != NULL)
        {
          delete H_half[lev];
          H_half[lev] = NULL;
        }
      
      if (layerEH_half[lev] != NULL)
        {
          delete layerEH_half[lev];
          layerEH_half[lev] = NULL;
        }
      if (layerH_half[lev] != NULL)
        {
          delete layerH_half[lev];
          layerH_half[lev] = NULL;
        }

      if (vectFluxes[lev] != NULL)
        {
          delete vectFluxes[lev];
          vectFluxes[lev] = NULL;
        }      
    }
  
  if (m_temporalAccuracy > 2)
    {
      MayDay::Error("AmrIce::timestep -- un-defined temporal accuracy");
    }

  //update time (velocity is to be computed at the step end)
  m_time += a_dt;
  m_cur_step += 1;
  // compute new ice velocity field
  if (m_evolve_velocity )
    {
      if (s_verbosity > 3) 
	{
	  pout() << "AmrIce::timeStep solveVelocityField() (step end) " << endl;
	}
      solveVelocityField();
    }
  
  if ((m_next_report_time - m_time) < (a_dt + TIME_EPS))
    {
       m_cf_domain_diagnostic_data.endTimestepDiagnostics
	 (m_vect_coordSys, m_old_thickness, m_divThicknessFlux, m_basalThicknessSource, m_surfaceThicknessSource, m_volumeThicknessSource,
	  m_calvedIceThickness, m_addedIceThickness, m_removedIceThickness,
	  m_amrGrids, m_refinement_ratios, m_amrDx[0], time(), m_time, m_dt,
	  m_cur_step, m_finest_level, s_verbosity);
    }

  if (s_verbosity > 0) 
    {
      pout () << "AmrIce::timestep " << m_cur_step
              << " --     end time = " 
	      << setiosflags(ios::fixed) << setprecision(6) << setw(12)
              << m_time  << " ( " << time() << " )"
        //<< " (" << m_time/secondsperyear << " yr)"
              << ", dt = " 
        //<< setiosflags(ios::fixed) << setprecision(6) << setw(12)
              << a_dt
        //<< " ( " << a_dt/secondsperyear << " yr )"
	      << resetiosflags(ios::fixed)
              << endl;
    }

  
  int totalCellsAdvanced = 0;
  for (int lev=0; lev<m_num_cells.size(); lev++) 
    {
      totalCellsAdvanced += m_num_cells[lev];
    }
     
  if (s_verbosity > 0) 
    {
      pout() << "Time = " << m_time  
             << " cells advanced = " 
             << totalCellsAdvanced << endl;

      for (int lev=0; lev<m_num_cells.size(); lev++) 
        {
          pout () << "Time = " << m_time 
                  << "  level " << lev << " cells advanced = " 
                  << m_num_cells[lev] << endl;
        }
    }
}


// compute half-time face-centered thickness using unsplit PPM
void
AmrIce::computeH_half(Vector<LevelData<FluxBox>* >& a_H_half, Real a_dt)
{
  CH_TIME("AmrIce::computeH_half");

  setToZero(m_volumeThicknessSource); // \todo : think about making this user settable
  
  for (int lev=0; lev<= finestTimestepLevel();  lev++)
    {
      
      // get AdvectPhysics object from PatchGodunov object
      PatchGodunov* patchGod = m_thicknessPatchGodVect[lev];
      AdvectPhysics* advectPhysPtr = dynamic_cast<AdvectPhysics*>(patchGod->getGodunovPhysicsPtr());
      if (advectPhysPtr == NULL)
        {
          MayDay::Error("AmrIce::timestep -- unable to upcast GodunovPhysics to AdvectPhysics");
        }
      
      patchGod->setCurrentTime(m_time);
      
      // loop over grids on this level and compute H-Half
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFaceVel = *m_faceVelAdvection[lev];
      LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev];
      LevelData<FluxBox>& levelHhalf = *a_H_half[lev];
      
      LevelData<FArrayBox>& levelSTS = *m_surfaceThicknessSource[lev];
      LevelData<FArrayBox>& levelBTS = *m_basalThicknessSource[lev];
      LevelData<FArrayBox>& levelVTS = *m_volumeThicknessSource[lev];
      CH_assert(m_surfaceFluxPtr != NULL);
      
      // set surface thickness source
      m_surfaceFluxPtr->surfaceThicknessFlux(levelSTS, *this, lev, a_dt);
      
      // set basal thickness source
      m_basalFluxPtr->surfaceThicknessFlux(levelBTS, *this, lev, a_dt);
 
     
      
      
      LevelData<FArrayBox> levelCCVel(levelGrids, SpaceDim, IntVect::Unit);
      EdgeToCell( levelFaceVel, levelCCVel);
      
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          patchGod->setCurrentBox(levelGrids[dit]);
          advectPhysPtr->setVelocities(&(levelCCVel[dit]), 
                                       &(levelFaceVel[dit]));
          
          FArrayBox advectiveSource(levelSTS[dit].box(),1);
          advectiveSource.copy(levelSTS[dit]);
          advectiveSource.plus(levelBTS[dit]);
	  advectiveSource.plus(levelVTS[dit]); 
          // add a diffusive source term div(D grad H)) to  advectiveSource
          if (m_diffusionTreatment == IMPLICIT)
            {
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  Box faceBox = levelGrids[dit].surroundingNodes(dir);
                  FArrayBox flux(faceBox,1);
                  FORT_FACEDERIV(CHF_FRA1(flux,0),
                                 CHF_CONST_FRA1(levelOldThickness[dit],0),
                                 CHF_BOX(faceBox),
                                 CHF_CONST_REAL(dx(lev)[dir]),
                                 CHF_INT(dir),
                                 CHF_INT(dir));
                  CH_assert(flux.norm(0) < HUGE_NORM);
                  flux *= (*m_diffusivity[lev])[dit][dir];
                  CH_assert(flux.norm(0) < HUGE_NORM);
                  FORT_DIVERGENCE(CHF_CONST_FRA(flux),
                                  CHF_FRA(advectiveSource),
                                  CHF_BOX(levelGrids[dit]),
                                  CHF_CONST_REAL(dx(lev)[dir]),
                                  CHF_INT(dir));
                  
                }
            }

	  // for this part, we need the *effective* thickness
	  FArrayBox he(levelOldThickness[dit].box(), 1);
	  FArrayBox f(levelOldThickness[dit].box(), 1);
	  f.copy( (*m_iceFrac[lev])[dit]);
	  f += TINY_NORM;
	  he.copy(levelOldThickness[dit]);
	  he /= f;
          
          patchGod->computeWHalf(levelHhalf[dit],
                                  levelOldThickness[dit],
                                 advectiveSource,
                                 a_dt,
                                 levelGrids[dit]);
          
          
        } //end loop over grids
      
      
    } // end loop over levels for computing Whalf
  
  // coarse average new H-Half to covered regions
  for (int lev= finestTimestepLevel(); lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],
                                     1, m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*a_H_half[lev-1], *a_H_half[lev]);
    }
  

}




void 
AmrIce::computeThicknessFluxes(Vector<LevelData<FluxBox>* >& a_vectFluxes,
                               const Vector<LevelData<FluxBox>* >& a_H_half,
                               const Vector<LevelData<FluxBox>* >& a_faceVelAdvection)
{

  CH_TIME("AmrIce::computeThicknessFluxes");
  
  for (int lev=0; lev<=finestTimestepLevel(); lev++)
    {
      LevelData<FluxBox>& levelFaceVel = *a_faceVelAdvection[lev];
      LevelData<FluxBox>& levelFaceH = *a_H_half[lev];
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& faceVel = levelFaceVel[dit];
          FluxBox& faceH = levelFaceH[dit];
          FluxBox& flux = (*a_vectFluxes[lev])[dit];
	  
          const Box& gridBox = levelGrids[dit];
          
          for (int dir=0; dir<SpaceDim; dir++)
            {
              Box faceBox(gridBox);
              faceBox.surroundingNodes(dir);
              flux[dir].copy(faceH[dir], faceBox);
              flux[dir].mult(faceVel[dir], faceBox, 0, 0, 1);
            }
        }
    } // end loop over levels
  
  // average fine fluxes down to coarse levels
  for (int lev=finestTimestepLevel(); lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],
                                     1, m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*a_vectFluxes[lev-1], *a_vectFluxes[lev]);
    }
  
}


void AmrIce::setStableSources(FArrayBox& a_smb,
			      FArrayBox& a_bmb,
			      FArrayBox& a_vmb,
			      const FArrayBox& a_divuh,
			      const BaseFab<int>& a_mask,
			      const Box& a_box) const
{

  bool floating_ice_stable = m_floating_ice_stable || (!m_evolve_thickness);
  bool grounded_ice_stable = m_grounded_ice_stable || (!m_evolve_thickness);
  
  //modify thickness sources if any part of the ice sheet is to be held stable, or dh/dt is set  
  for (BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      
      if (a_mask(iv) == FLOATINGMASKVAL)
	{
	  //keep floating ice stable if required, assuming that any balance is due to basal freezing/melting
	  if ( floating_ice_stable )
	    {
	      a_bmb(iv) = a_divuh(iv) - a_smb(iv) ;
	      a_vmb(iv) = 0.0;
	    }
	  else if (m_floating_ice_basal_flux_is_dhdt)
	    {
	      a_bmb(iv) += ( a_divuh(iv) - a_smb(iv) );
	      a_vmb(iv) = 0.0;
	    }
	  else if (m_floating_ice_basal_flux_is_min_dhdt)
	    {
	      a_bmb(iv) -= std::max(0.0,-a_divuh(iv)+a_smb(iv));
	      a_vmb(iv) = 0.0;
	    }
	} // end if (a_mask(iv) == FLOATINGMASKVAL)
      else if (a_mask(iv) == GROUNDEDMASKVAL)
	{
	  //keep grounded ice stable if required
	  if (grounded_ice_stable )
	    {
	      if (a_divuh(iv) < 0)
	      	{
	      	  // need to thin the ice. attempt to balance between isothermal compression and flux out of the base 
		  //a_bmb(iv) = a_divuh(iv);
		  //a_vmb(iv) = - a_smb(iv);
		  a_vmb(iv) = a_divuh(iv) - a_smb(iv) - a_bmb(iv);
	      	}
	      else
	      	{
	      	  //need to thicken the ice. distrubuted new energy evenly through the thickness, i.e isothermal expansion
	      	  a_vmb(iv) = a_divuh(iv) - a_smb(iv) - a_bmb(iv);
	      	}
	    }
	  else if (m_grounded_ice_basal_flux_is_dhdt)
	    {
	      a_vmb(iv) = a_divuh(iv) - a_smb(iv);
	    }	  
	} // end else if (a_mask(iv) == GROUNDEDMASKVAL)
      else
	{
	  // ice free regions.
	  if (!m_evolve_thickness)
	    {
	      a_vmb(iv) = a_divuh(iv) - a_smb(iv) - a_bmb(iv);
	    }
	}
    } //  end for (BoxIterator bit(a_box); bit.ok(); ++bit)
}

  

// update  ice thickness *and* bedrock elevation
void
AmrIce::updateGeometry(Vector<RefCountedPtr<LevelSigmaCS> >& a_vect_coordSys_new, 
		       Vector<RefCountedPtr<LevelSigmaCS> >& a_vect_coordSys_old, 
		       const Vector<LevelData<FluxBox>* >& a_vectFluxes, 
		       Real a_dt)
{

  CH_TIME("AmrIce::updateGeometry");
  for (int lev=0; lev <= finestTimestepLevel() ; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFlux = *a_vectFluxes[lev];
      LevelSigmaCS& levelCoords = *(a_vect_coordSys_new[lev]);
      LevelData<FArrayBox>& levelNewH = levelCoords.getH();
      LevelData<FArrayBox>& levelOldH = (*a_vect_coordSys_old[lev]).getH();
      LevelData<FArrayBox>& levelDivThckFlux = *m_divThicknessFlux[lev];
      const RealVect& dx = levelCoords.dx();              
      
      DataIterator dit = levelGrids.dataIterator();          
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          FArrayBox& newH = levelNewH[dit];
          FArrayBox& oldH = levelOldH[dit];//(*m_old_thickness[lev])[dit];
          FluxBox& thisFlux = levelFlux[dit];
          newH.setVal(0.0);
          
          // loop over directions and increment with div(F)
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // use the divergence from 
              // Chombo/example/fourthOrderMappedGrids/util/DivergenceF.ChF
              FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                              CHF_FRA(newH),
                              CHF_BOX(gridBox),
                              CHF_CONST_REAL(dx[dir]),
                              CHF_INT(dir));
              
              
            }
          
	  levelDivThckFlux[dit].copy(newH);

	  // // modify sources as required to balance div(uh)
	  setStableSources((*m_surfaceThicknessSource[lev])[dit],
	  		   (*m_basalThicknessSource[lev])[dit],
	  		   (*m_volumeThicknessSource[lev])[dit],
	  		   levelDivThckFlux[dit],
	  		   levelCoords.getFloatingMask()[dit],
	  		   gridBox);
	  
          // add in thickness source
          // if there are still diffusive fluxes to deal
          // with, the source term will be included then
 
	  //if (m_diffusionTreatment != IMPLICIT)
            {
              if (m_frac_sources)
                {
                  // scale surface fluxes by mask values
                  const FArrayBox& thisFrac = (*m_iceFrac[lev])[dit];
                  FArrayBox sources(gridBox,1);
                  sources.setVal(0.0);
                  sources.plus((*m_surfaceThicknessSource[lev])[dit], gridBox,
                               0, 0, 1);
                  sources.plus((*m_basalThicknessSource[lev])[dit], gridBox, 
                               0, 0, 1);
		  sources.plus((*m_volumeThicknessSource[lev])[dit], gridBox, 
                               0, 0, 1);
                  
                  sources.mult(thisFrac, gridBox, 0, 0, 1);
                  newH.minus(sources, gridBox, 0, 0, 1);
                  
                }
              else 
                {

                  // just add in sources directly
                  newH.minus((*m_surfaceThicknessSource[lev])[dit], gridBox,0,0,1);
                  newH.minus((*m_basalThicknessSource[lev])[dit], gridBox,0,0,1);
		  newH.minus((*m_volumeThicknessSource[lev])[dit], gridBox,0,0,1);
                }
            }
          
	          
	  newH *= -1*a_dt;
          newH.plus(oldH, 0, 0, 1);

	  for (BoxIterator bit(gridBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();

	      Real H=newH(iv);
	      // Remove negative thickness by limiting low bmb and/or smb
	      // Calculate the effective basal and surface mass fluxes
	      if (H < 0.0)
		{
		  Real excessSource=0.0;
		  if (a_dt > 0.0)
		    {
		      excessSource=H/a_dt;
		    }
		  Real bts = (*m_basalThicknessSource[lev])[dit](iv);
		  Real sts = (*m_surfaceThicknessSource[lev])[dit](iv);
		  Real oldBTS=bts;
		  if (bts < 0.0)
		    {
		      bts = std::min(bts-excessSource,0.0);
		    }
		  sts = sts+oldBTS - excessSource - bts;
		  (*m_basalThicknessSource[lev])[dit](iv) = bts;
		  (*m_surfaceThicknessSource[lev])[dit](iv) = sts;
		  
                  newH(iv)=0.0;
		}
	    }

          
        } // end loop over grids 
    } // end loop over levels
  
  
  
  //include any diffusive fluxes
  if (m_evolve_thickness && m_diffusionTreatment == IMPLICIT)
    {
      if (m_grounded_ice_stable || m_floating_ice_stable)
	{
	  CH_assert( !(m_grounded_ice_stable || m_floating_ice_stable));
	  MayDay::Error("implicit diffusion not implemented with grounded_ice_stable or floating_ice_stable ");
	}
      //MayDay::Error("m_diffusionTreatment == IMPLICIT no yet implemented");
      //implicit thickness correction
      if (m_frac_sources)
        {
          MayDay::Error("scaling sources by ice fraction values not implemented yet");
        }
      implicitThicknessCorrection(a_dt, m_surfaceThicknessSource,  m_basalThicknessSource, m_volumeThicknessSource);
    }

  //update the topography (fixed surface case)
  if (m_evolve_topography_fix_surface)
    {
      // update the bedrock so that the surface remains constant on grounded ice
      for (int lev=0; lev <= finestTimestepLevel() ; lev++)
	{
	  for (DataIterator dit(m_amrGrids[lev]);dit.ok();++dit)
	    {
	      FArrayBox& newH = a_vect_coordSys_new[lev]->getH()[dit];
	      FArrayBox& oldH = a_vect_coordSys_old[lev]->getH()[dit];
	      FArrayBox& topg = a_vect_coordSys_new[lev]->getTopography()[dit];
	      
	      const BaseFab<int>& mask = a_vect_coordSys_old[lev]->getFloatingMask()[dit];
	      FORT_EVOLVEGROUNDEDBED(CHF_FRA1(newH,0), CHF_FRA1(oldH,0), 
				     CHF_FRA1(topg,0), CHF_CONST_FIA1(mask,0), 
				     CHF_BOX(topg.box()));
	      FArrayBox& deltaTopg = (*m_deltaTopography[lev])[dit];
	      deltaTopg += topg;
	      deltaTopg -= a_vect_coordSys_old[lev]->getTopography()[dit];	      
	    }
	}
    }


  //update the topography (gia)
  if (m_topographyFluxPtr != NULL)
    {
      for (int lev=0; lev <= finestTimestepLevel() ; lev++)
	{
	  LevelSigmaCS& levelCoords = *(a_vect_coordSys_new[lev]);
	  DisjointBoxLayout& levelGrids = m_amrGrids[lev];
	  LevelData<FArrayBox>& levelTopg = levelCoords.getTopography();
	  LevelData<FArrayBox>& levelDeltaTopg = *m_deltaTopography[lev];
	  LevelData<FArrayBox> levelSrc(levelGrids,1,IntVect::Zero);
	  m_topographyFluxPtr->surfaceThicknessFlux(levelSrc, *this, lev, a_dt);
	  for (DataIterator dit(levelGrids);dit.ok();++dit)
	    {
	      FArrayBox& src = levelSrc[dit];
	      src *= a_dt;
	      levelTopg[dit] += src;
	      levelDeltaTopg[dit] += src;
	    }
	}
    }
  
  


  // average down thickness and topography to coarser levels and fill in ghost cells
  // before calling recomputeGeometry. 
  int Hghost = 2;
  Vector<LevelData<FArrayBox>* > vectH(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectB(m_finest_level+1, NULL);
  for (int lev=0; lev<vectH.size(); lev++)
    {
      IntVect HghostVect = Hghost*IntVect::Unit;
      LevelSigmaCS& levelCoords = *(a_vect_coordSys_new[lev]);
      vectH[lev] = &levelCoords.getH();
      vectB[lev] = &levelCoords.getTopography();
    }
  
  //average from the finest level down
  for (int lev =  finestTimestepLevel() ; lev > 0 ; lev--)
    {
      CoarseAverage averager(m_amrGrids[lev],
                             1, m_refinement_ratios[lev-1]);
      averager.averageToCoarse(*vectH[lev-1],
                               *vectH[lev]);
      averager.averageToCoarse(*vectB[lev-1],
                               *vectB[lev]);
    }

  // now pass back over and do PiecewiseLinearFillPatch
  for (int lev=1; lev<vectH.size(); lev++)
    {
      
      PiecewiseLinearFillPatch filler(m_amrGrids[lev],
                                      m_amrGrids[lev-1],
                                      1, 
                                      m_amrDomains[lev-1],
                                      m_refinement_ratios[lev-1],
                                      Hghost);
      
      Real interp_coef = 0;
      filler.fillInterp(*vectH[lev],
                        *vectH[lev-1],
                        *vectH[lev-1],
                        interp_coef,
                        0, 0, 1);
      filler.fillInterp(*vectB[lev],
                        *vectB[lev-1],
                        *vectB[lev-1],
                        interp_coef,
                        0, 0, 1);
    }
  
  
  //re-fill ghost cells ouside the domain
  for (int lev=0; lev <= finestTimestepLevel()  ; ++lev)
    {
      RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
      m_thicknessIBCPtr->setGeometryBCs(*a_vect_coordSys_new[lev],
                                        m_amrDomains[lev],levelDx, m_time, m_dt);
    }
  
  //allow calving model to modify geometry and velocity
  applyCalvingCriterion(CalvingModel::PostThicknessAdvection);
  
  //dont allow thickness or iceFrac to be negative
  for (int lev=0; lev<= m_finest_level; lev++)
    {

      
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelSigmaCS& levelCoords = *(a_vect_coordSys_new[lev]);
      LevelData<FArrayBox>& levelH = levelCoords.getH();

      if (m_evolve_ice_frac)
	{
	  updateIceFrac(levelH, lev);
	}
      else
	{
	  setIceFrac(levelH, lev);
	}
      
      DataIterator dit = levelGrids.dataIterator();

      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
        {
          Real lim = 0.0;
          FORT_MAXFAB1(CHF_FRA(levelH[dit]), 
                       CHF_CONST_REAL(lim), 
                       CHF_BOX(levelH[dit].box()));
        }
    }

  
  
  //average from the finest level down
  for (int lev =  finestTimestepLevel() ; lev > 0 ; lev--)
    {
      CoarseAverage averager(m_amrGrids[lev],
                             1, m_refinement_ratios[lev-1]);
      averager.averageToCoarse(*vectH[lev-1],
                               *vectH[lev]);      
    }
  
  for (int lev=1; lev<vectH.size(); lev++)
    {      
      PiecewiseLinearFillPatch filler(m_amrGrids[lev],
                                      m_amrGrids[lev-1],
                                      1, 
                                      m_amrDomains[lev-1],
                                      m_refinement_ratios[lev-1],
                                      Hghost);
      
      Real interp_coef = 0;
      filler.fillInterp(*vectH[lev],
                        *vectH[lev-1],
                        *vectH[lev-1],
                        interp_coef,
                        0, 0, 1);
    }
  
  //interpolate levelSigmaCS to any levels above finestTimestepLevel()
  for (int lev = finestTimestepLevel()+1 ; lev<= m_finest_level; lev++)
    {
      m_vect_coordSys[lev]->interpFromCoarse(*m_vect_coordSys[lev-1],
                                             m_refinement_ratios[lev-1],
                                             false , true);
    }
  
  
  // recompute thickness-derived data in SigmaCS
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
      LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
      int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
      levelCoords.recomputeGeometry(crseCoords, refRatio);            
    }

  //the grounding line may have moved.
  m_groundingLineProximity_valid = false;
}


void
AmrIce::levelSetup(int a_level, const DisjointBoxLayout& a_grids)
{
  IntVect ghostVect = IntVect::Unit;
  // 4 ghost cells needed for advection. Could later go back and
  // make this a temporary if the additional storage becomes an issue...
  IntVect thicknessGhostVect = m_num_thickness_ghost*IntVect::Unit;
  m_old_thickness[a_level]->define(a_grids, 1,
                                   thicknessGhostVect);

  if (a_level == 0 || m_velocity[a_level] == NULL)
    {
      m_velocity[a_level] = new LevelData<FArrayBox>(a_grids, SpaceDim,
                                                     ghostVect);
    }
  else
    {
      // do velocity a bit differently in order to use previously 
      // computed velocity field as an initial guess
      {
        LevelData<FArrayBox>* newVelPtr = new LevelData<FArrayBox>(a_grids,
                                                                   SpaceDim,
                                                                   ghostVect);
        
        // first do interp from coarser level
        FineInterp velInterp(a_grids, SpaceDim, 
                             m_refinement_ratios[a_level-1],
                             m_amrDomains[a_level]);
        
        velInterp.interpToFine(*newVelPtr, *m_velocity[a_level-1]);
        
        // can only copy from existing level if we're not on the
        // newly created level
        //if (a_level != new_finest_level)
        if (m_velocity[a_level]->isDefined())
          {
            m_velocity[a_level]->copyTo(*newVelPtr);
          }
        
        // finally, do an exchange (this may wind up being unnecessary)
        newVelPtr->exchange();
        
        delete (m_velocity[a_level]);
        m_velocity[a_level] = newVelPtr;
      }
    } // end interpolate/copy new velocity

  levelAllocate(&m_faceVelAdvection[a_level] ,a_grids,1,IntVect::Unit);
  levelAllocate(&m_faceVelTotal[a_level],a_grids,1,IntVect::Unit);
  levelAllocate(&m_diffusivity[a_level],a_grids, 1, IntVect::Zero);
  levelAllocate(&m_iceFrac[a_level],a_grids, 1, IntVect::Unit);

#if BISICLES_Z == BISICLES_LAYERED
  levelAllocate(&m_layerXYFaceXYVel[a_level] ,a_grids, m_nLayers, IntVect::Unit);
  levelAllocate(&m_layerSFaceXYVel[a_level], a_grids, SpaceDim*(m_nLayers + 1), IntVect::Unit);
  levelAllocate(&m_layerSFaceSVel[a_level], a_grids, m_nLayers + 1, IntVect::Unit);
#endif

  levelAllocate(&m_velBasalC[a_level],a_grids, 1, ghostVect);
  levelAllocate(&m_cellMuCoef[a_level],a_grids, 1, ghostVect);
  levelAllocate(&m_velRHS[a_level],a_grids, SpaceDim,  IntVect::Zero);
  levelAllocate(&m_surfaceThicknessSource[a_level], a_grids,   1, IntVect::Unit) ;
  levelAllocate(&m_basalThicknessSource[a_level], a_grids,  1, IntVect::Unit) ;
  levelAllocate(&m_volumeThicknessSource[a_level], a_grids,  1, IntVect::Unit) ;
  levelAllocate(&m_divThicknessFlux[a_level],a_grids,   1, IntVect::Zero) ;
  levelAllocate(&m_calvedIceThickness[a_level],a_grids, 1, IntVect::Unit);
  levelAllocate(&m_removedIceThickness[a_level],a_grids, 1, IntVect::Unit);
  levelAllocate(&m_addedIceThickness[a_level],a_grids, 1, IntVect::Unit);
  levelAllocate(&m_deltaTopography[a_level],a_grids, 1, IntVect::Zero);
  // probably eventually want to do this differently
  RealVect dx = m_amrDx[a_level]*RealVect::Unit;

  //IntVect sigmaCSGhost = IntVect::Unit;
  IntVect sigmaCSGhost = thicknessGhostVect;
  m_vect_coordSys.resize(m_max_level+1);
  m_vect_coordSys[a_level] = RefCountedPtr<LevelSigmaCS >(new LevelSigmaCS(a_grids, 
									   dx,
									   sigmaCSGhost));
  m_vect_coordSys[a_level]->setIceDensity(m_iceDensity);
  m_vect_coordSys[a_level]->setWaterDensity(m_seaWaterDensity);
  m_vect_coordSys[a_level]->setGravity(m_gravity);

#if BISICLES_Z == BISICLES_LAYERED
  //in poor-man's multidim mode, use one FArrayBox component per layer
  //to hold the 3D internalEnergy field
  levelAllocate(&m_internalEnergy[a_level],a_grids, m_nLayers,thicknessGhostVect);
  levelAllocate(&m_tillWaterDepth[a_level], a_grids, 1, thicknessGhostVect);
  levelAllocate(&m_sInternalEnergy[a_level], a_grids, 1, thicknessGhostVect);
  levelAllocate(&m_bInternalEnergy[a_level], a_grids, 1, thicknessGhostVect);
  levelAllocate(&m_sHeatFlux[a_level], a_grids, 1, thicknessGhostVect);
  levelAllocate(&m_bHeatFlux[a_level], a_grids, 1, thicknessGhostVect);
  m_vect_coordSys[a_level]->setFaceSigma(getFaceSigma());

#elif BISICLES_Z == BISICLES_FULLZ
  levelAllocate(&m_internalEnergy[a_level],a_grids, 1, thicknessGhostVect);
#endif




}

void
AmrIce::initData(Vector<RefCountedPtr<LevelSigmaCS> >& a_vectCoordSys,
                 Vector<LevelData<FArrayBox>* >& a_velocity)
{

  CH_TIME("AmrIce::initData");
  //  for (int i=0; i < 
  //m_diagnostic_values[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::initData" << endl;
    }

  m_groundingLineProximity_valid = false;
  m_A_valid = false;

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
      m_thicknessIBCPtr->define(m_amrDomains[lev],levelDx[0]);
      LevelSigmaCS* crsePtr = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
      int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:0;
      
      m_thicknessIBCPtr->initializeIceGeometry(*a_vectCoordSys[lev],
					       levelDx,
					       m_domainSize,
					       m_time,
					       crsePtr,
					       refRatio);
      

	


      a_vectCoordSys[lev]->recomputeGeometry(crsePtr, refRatio);



      const LevelData<FArrayBox>& levelThickness = m_vect_coordSys[lev]->getH();
      setIceFrac(levelThickness, lev);
      a_vectCoordSys[lev]->recomputeGeometry(crsePtr, refRatio);

      // initialize oldH to be the current value
      LevelData<FArrayBox>& currentH = a_vectCoordSys[lev]->getH();
      currentH.copyTo(*m_old_thickness[lev]);

#if BISICLES_Z == BISICLES_LAYERED
      m_internalEnergyIBCPtr->initializeIceInternalEnergy
	(*m_internalEnergy[lev], *m_tillWaterDepth[lev], *m_sInternalEnergy[lev], *m_bInternalEnergy[lev], *this, lev, 0.0);
#elif BISICLES_Z == BISICLES_FULLZ
      m_internalEnergyIBCPtr->initializeIceInternalEnergy(*m_temperature[lev],*this, lev, 0.0);
#endif
     
    }
  // tempearture depends on internal energy
  updateTemperature();

  // this is a good time to check for remote ice
  // (don't bother if we're doing it as a matter of course, since we'd
  // wind up doing it 2x)
  if ((m_eliminate_remote_ice_after_regrid) && !(m_eliminate_remote_ice))
    eliminateRemoteIce();
  
  setToZero(m_deltaTopography);

  applyCalvingCriterion(CalvingModel::Initialization);
  

  
  // now call velocity solver to initialize velocity field, force a solve no matter what the time step
  solveVelocityField(true);

  // may be necessary to average down here
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverage avgDown(m_amrGrids[lev],
                            SpaceDim, m_refinement_ratios[lev-1]);
      avgDown.averageToCoarse(*m_velocity[lev-1], *m_velocity[lev]);
    }


  //#define writePlotsImmediately
#ifdef  writePlotsImmediately
  if (m_plot_interval >= 0)
    {
#ifdef CH_USE_HDF5
      writePlotFile();
#endif
    }
#endif


}

/// solve for velocity field
void
AmrIce::solveVelocityField(bool a_forceSolve, Real a_convergenceMetric)
{

  CH_TIME("AmrIce::solveVelocityField");

  notifyObservers(Observer::PreVelocitySolve);

  if (m_eliminate_remote_ice)
    eliminateRemoteIce();

  //ensure A is up to date
#if BISICLES_Z == BISICLES_LAYERED
  if (!m_A_valid)
    {
      computeA(m_A,m_sA,m_bA,m_internalEnergy,m_sInternalEnergy,m_bInternalEnergy,
	       m_vect_coordSys);
      m_A_valid = true;
    }
#else
  MayDay::Error("AmrIce::SolveVelocityField full z calculation of A not done"); 
#endif
  //certainly the viscous tensr field will need re-computing
  m_viscousTensor_valid = false;
  

  
  // define basal friction
  Vector<LevelData<FArrayBox>* > vectC(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectC0(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectRhs(m_finest_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      vectRhs[lev] = m_velRHS[lev];
      vectC[lev] = m_velBasalC[lev];
      vectC0[lev] = new LevelData<FArrayBox>; vectC0[lev]->define(*vectC[lev]);
    }

  //
  setMuCoefficient(m_cellMuCoef);

  // set basal friction coeffs C,C0. C = 0 for floating ice. C0 != 0 at walls
  setBasalFriction(vectC, vectC0);

  m_basalFrictionRelation->updateThermodynamics(&m_tillWaterDepth);
  
  // right hand side of the stress-balance
  defineVelRHS(vectRhs);




  
  // write out sumRhs if appropriate
  if (s_verbosity > 3)
    {
      Real sumRhs = computeSum(vectRhs,
                               m_refinement_ratios,
                               m_amrDx[0],
                               Interval(0,0),
                               0);

      pout() << "Sum(rhs) for velocity solve = " << sumRhs << endl;

    }
  
  // put this in place to catch runs where plotfile writing is
  // going to hang _before_ I waste a few hours waiting for the 
  // velocity solve
  //#define writeTestPlots
#ifdef  writeTestPlots
  if (m_plot_interval >= 0 && m_cur_step == 0)
    {
      writePlotFile();
    }
#endif

  if (m_doInitialVelSolve) 
    {      
      if (m_finest_level == 0 && m_doInitialVelGuess)
        {
          // only really want or need to do this once
          m_doInitialVelGuess = false;
	  if (m_initialGuessType == SlidingLaw)
	    {
	      pout() << "computing an initial guess via a sliding law u = rhs/C "  << endl;
	      // compute initial guess as rhs/beta
	      LevelData<FArrayBox>& vel = *m_velocity[0];
	      LevelData<FArrayBox>& C = *m_velBasalC[0];
	      LevelData<FArrayBox>& rhs = *m_velRHS[0];
	      const DisjointBoxLayout& levelGrids = m_amrGrids[0];
	      
	      DataIterator dit = vel.dataIterator();
	      for (dit.begin(); dit.ok(); ++dit)
		{
		  FORT_VELINITIALGUESS(CHF_FRA(vel[dit]),
				       CHF_FRA(rhs[dit]),
				       CHF_FRA1(C[dit],0),
				       CHF_BOX(levelGrids[dit]));
		}
	    }
	  else if (m_initialGuessType == ConstMu)
	    {
	      if (s_verbosity > 3) 
		{
		  pout() << "computing an initial guess by solving the velocity equations "
			 <<" with constant mu = " 
			 << m_initialGuessConstMu   
			 << " and constant initial velocity = " << m_initialGuessConstVel
			 << endl;
		}

	      // compute initial guess by solving a linear problem with a
	      // modest constant viscosity
	      constMuRelation* newPtr = new constMuRelation;
	      newPtr->setConstVal(m_initialGuessConstMu);
	      for (int lev=0; lev < m_finest_level + 1; lev++)
		{
		  for (DataIterator dit(m_amrGrids[lev]);dit.ok();++dit)
		    {
		      for (int dir = 0; dir < SpaceDim; dir++)
			{
			  (*m_velocity[lev])[dit].setVal(m_initialGuessConstVel[dir],dir);
			}
		    }

		}

              // do this by saving the exisiting velSolver and 
              // constitutiveRelation, re-calling defineSolver, then 
              // doing solve.
              IceVelocitySolver* velSolverSave = m_velSolver;
              ConstitutiveRelation* constRelSave = m_constitutiveRelation;
              int solverTypeSave = m_solverType;

              // new values prior to calling defineSolver
             
              m_constitutiveRelation = static_cast<ConstitutiveRelation*>(newPtr);

	      Real finalNorm = 0.0, initialNorm = 0.0, convergenceMetric = -1.0;
	      //Vector<LevelData<FArrayBox>* > muCoef(m_finest_level + 1,NULL);
	      int rc = -1;
	      if (m_initialGuessSolverType == JFNK)
		{
		  //JFNK can be instructed to assume a linear solve
		  m_solverType = JFNK;
                  // (DFM 2/4/14) this is not a memory leak -- velSolver is 
                  // saved in velSolverSave and will be swapped back after 
                  // the initial guess solve
		  m_velSolver = NULL;
		  defineSolver();
		  JFNKSolver* jfnkSolver = dynamic_cast<JFNKSolver*>(m_velSolver);
		  CH_assert(jfnkSolver != NULL);
		  const bool linear = true;
		  rc = jfnkSolver->solve( m_velocity, 
					  m_calvedIceThickness, m_addedIceThickness, m_removedIceThickness,
					  initialNorm,finalNorm,convergenceMetric,
					  linear , m_velRHS, m_velBasalC, vectC0, m_A, m_cellMuCoef,
					  m_vect_coordSys, m_time, 0, m_finest_level);
		}
	      else if (m_initialGuessSolverType == Picard)
		{
		  ParmParse pp("picardSolver");
		  Real tol = 1.e-4; int nits = 1;
		  // since the constant-viscosity solve is a linear solve,
		  // Picard is the best option.
		  m_solverType = Picard;
		  m_velSolver = NULL;
		  defineSolver();

		  pp.query("linearsolver_tolerance", tol );
		  pp.query("max_picard_iterations", nits );
		  m_velSolver->setTolerance(nits);		  
		  m_velSolver->setMaxIterations(tol);

		  rc = m_velSolver->solve(m_velocity, 
					  m_calvedIceThickness, m_addedIceThickness, m_removedIceThickness,
					  initialNorm,finalNorm,convergenceMetric,
					  m_velRHS, m_velBasalC, vectC0, m_A, m_cellMuCoef,
					  m_vect_coordSys, m_time, 0, m_finest_level);
		}
	      else
		{
		  MayDay::Error("unknown initial guess solver type");
		}


	      if (rc != 0)
		{
		  MayDay::Warning("constant mu solve failed");
		}

              // now put everything back the way it was...
	      delete m_constitutiveRelation;
              delete m_velSolver;
              m_velSolver = velSolverSave;
              m_constitutiveRelation = constRelSave;
              m_solverType = solverTypeSave;

#if 0	      
	      //put the solver back how it was
	      m_velSolver->define(m_amrDomains[0],
				  m_constitutiveRelation,
				  m_basalFrictionRelation,
				  m_amrGrids,
				  m_refinement_ratios,
				  dxCrse,
				  m_thicknessIBCPtr,
				  numLevels);
#endif
      
	    }
	  else if (m_initialGuessType == Function)
	    {
	      ParmParse pp("amr");
	      std::string functionType = "constant";
	      pp.query("initial_velocity_function_type", functionType );
	      if (functionType == "flowline")
		{
		  Real dx; 
		  std::string file, set;
		  pp.get("initial_velocity_function_flowline_dx", dx);
		  pp.get("initial_velocity_function_flowline_file", file);
		  pp.get("initial_velocity_function_flowline_set", set);
		  ExtrudedPieceWiseLinearFlowline f(file,set,dx);
		  for (int lev = 0; lev <  m_finest_level + 1; lev++)
		    {
		      LevelData<FArrayBox>& levelVel = *m_velocity[lev];
		      for (DataIterator dit(levelVel.disjointBoxLayout());
			   dit.ok(); ++dit)
			{
			  FArrayBox& vel = levelVel[dit];
			  const Box& box = vel.box(); 
			  for (BoxIterator bit(box); bit.ok(); ++bit)
			    {
			      const IntVect& iv = bit();
			      RealVect x = RealVect(iv) * m_amrDx[lev] 
				+ 0.5 * m_amrDx[lev];
			      vel(iv,0) = f(x);
			    }
			}
		    }
		}
	      
	    }
	  else
	    {
	      MayDay::Error("AmrIce::SolveVelocityField unknown initial guess type");
	    }
	}
#ifdef CH_USE_HDF5
      if (m_write_presolve_plotfiles)
        {
          string save_prefix = m_plot_prefix;
          m_plot_prefix.append("preSolve.");
	  bool t_write_fluxVel = m_write_fluxVel;
	  m_write_fluxVel = false; // turning this off in preSolve files for now
          writePlotFile();
          m_write_fluxVel = t_write_fluxVel;
	  m_plot_prefix = save_prefix;
        }
#endif

      int solverRetVal; 
    
      //set u = 0 in ice free cells
      for (int lev=0; lev <= m_finest_level ; ++lev)
	{
	  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
	  LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	      FArrayBox& vel = (*m_velocity[lev])[dit];
	      for (BoxIterator bit(levelGrids[dit]); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		  if (mask(iv) == OPENSEAMASKVAL || 
		      mask(iv) == OPENLANDMASKVAL )
		    {
		      D_TERM(vel(iv,0) = 0.0;, vel(iv,1) = 0.0;,)
		    } 
		}
	    }
	}

	if (a_forceSolve || ((m_cur_step+1)%m_velocity_solve_interval == 0))
	  {

	    //any thickness change within m_velSolver->solve is assumed to be calving.
	    notifyObservers(Observer::PreCalving);
	    solverRetVal = m_velSolver->solve(m_velocity, 
					      m_calvedIceThickness, 
					      m_addedIceThickness,
					      m_removedIceThickness,
					      m_velocitySolveInitialResidualNorm, 
					      m_velocitySolveFinalResidualNorm,
					      a_convergenceMetric,
					      m_velRHS, m_velBasalC, vectC0,
					      m_A, m_cellMuCoef,
					      m_vect_coordSys,
					      m_time,
					      0, m_finest_level);

	    notifyObservers(Observer::PostCalving);
	    
	    if (solverRetVal != 0)
	      {
		pout() << " solver return value = "
		       << solverRetVal << std::endl;
		//MayDay::Warning("solver return value != 0"); 
	      }
	    for (int lev = 0; lev <= m_finest_level; lev++)
	      {
		m_thicknessIBCPtr->velocityGhostBC
		  (*m_velocity[lev],*m_vect_coordSys[lev],
		   m_amrDomains[lev], m_time);
	      }
	    
	    //special case for inverse problems : read back C and muCoef if they are ready
	    InverseIceVelocitySolver* invPtr = dynamic_cast<InverseIceVelocitySolver*>(m_velSolver);
	    if (invPtr)
	      {
		BasalFriction* bfptr = invPtr->basalFriction();
		if (bfptr)
		  {
		    if (m_basalFrictionPtr) delete m_basalFrictionPtr; 
		    m_basalFrictionPtr = bfptr;
		  }

		MuCoefficient* mcptr = invPtr->muCoefficient();
		if (mcptr)
		  {
		    if (m_muCoefficientPtr) delete m_muCoefficientPtr;
		    m_muCoefficientPtr = mcptr;
		  } 
	      } // end special case for inverse problems

	    
	  } // end if (a_forceSolve || ((m_cur_step+1)%m_velocity_solve_interval == 0))
  
	// update the viscous tensor  
	updateViscousTensor();
	//allow calving model to modify geometry 
	applyCalvingCriterion(CalvingModel::PostVelocitySolve);
	
	for (int lev=0; lev<=m_finest_level; lev++)
	  {
	    if (vectC0[lev] != NULL)
	      {
		delete vectC0[lev]; vectC0[lev] = NULL;
	      }
	  }
  
	/// This is probably the most useful notification, as a velocity 
	/// solve is carried out at the end of every major stage
	notifyObservers(Observer::PostVelocitySolve);
	
    } // end if (m_doInitialSolve)

  // update the viscous tensor  
  updateViscousTensor();
  //calculate the face centred (flux) velocity and diffusion coefficients
  computeFaceVelocity(m_faceVelAdvection,m_faceVelTotal,m_diffusivity,m_layerXYFaceXYVel, m_layerSFaceXYVel);

}

	  

    


void AmrIce::defineVelRHS(Vector<LevelData<FArrayBox>* >& a_vectRhs)
{

  Vector<RealVect> dx;
  Vector<LevelData<FArrayBox>*> rhs;
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      dx.push_back(m_vect_coordSys[lev]->dx());
      rhs.push_back(a_vectRhs[lev]);
    }

  IceUtility::defineRHS(rhs, m_vect_coordSys,  m_amrGrids, dx);

  //\todo : move this into IceUtility::defineRHS
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      // finally, modify RHS in problem-dependent ways,
      m_thicknessIBCPtr->modifyVelocityRHS(*a_vectRhs[lev],  *m_vect_coordSys[lev],
                                           m_amrDomains[lev],m_time, m_dt);
    }

}


/// set mu coefficient (phi) prior to velocity solve
void
AmrIce::setMuCoefficient(Vector<LevelData<FArrayBox>* >& a_cellMuCoef)
{
  CH_assert(m_muCoefficientPtr != NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      m_muCoefficientPtr->setMuCoefficient(*a_cellMuCoef[lev],
					   *m_vect_coordSys[lev],
                                           this->time(),
                                           m_dt);
      if (lev > 0)
	{
	  PiecewiseLinearFillPatch ghostFiller
	    (m_amrGrids[lev],m_amrGrids[lev-1],1,m_amrDomains[lev-1],
	     m_refinement_ratios[lev-1],1);
	  
	  ghostFiller.fillInterp(*a_cellMuCoef[lev],
				 *a_cellMuCoef[lev-1],
				 *a_cellMuCoef[lev-1],
				 1.0,0,0,1);
	}
      a_cellMuCoef[lev]->exchange();
    }
}


/// set basal friction coefficients C,C0 prior to velocity solve
void
AmrIce::setBasalFriction(Vector<LevelData<FArrayBox>* >& a_vectC,Vector<LevelData<FArrayBox>* >& a_vectC0)
{

  // first, compute C and C0 as though there was no floating ice
  CH_assert(m_basalFrictionPtr != NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      m_basalFrictionPtr->setBasalFriction(*a_vectC[lev], *m_vect_coordSys[lev],
                                           this->time(),m_dt); 
      if (m_basalRateFactor != NULL)
	{
	  //basal temperature dependence
	  LevelData<FArrayBox>& C = *a_vectC[lev];
	  Vector<Real> bSigma(1,1.0);
	  LevelData<FArrayBox> A(C.disjointBoxLayout(),1,C.ghostVect());
	  IceUtility::computeA(A, bSigma,*m_vect_coordSys[lev],  
			       m_basalRateFactor, *m_bInternalEnergy[lev]);
	  for (DataIterator dit = C.dataIterator(); dit.ok(); ++dit)
	    {
	      C[dit] /= A[dit];
	    }
	}
     
      a_vectC[lev]->exchange();
    }

  // compute C0
  // C0 include a term (wall drag) that depends on C,
  // so needs to be computed before setting C = 0 in floating regions
  IceUtility::computeC0(a_vectC0, a_vectC, m_amrGrids, m_vect_coordSys, m_amrDx, m_finest_level);

  if ( m_reset_floating_friction_to_zero )
    {
      //set C = 0 in floating region, possibly employing a thickness-above-flotation interpolation 
      for (int lev=0; lev<=m_finest_level; lev++)
	{
	  IceUtility::setFloatingBasalFriction(*a_vectC[lev], *m_vect_coordSys[lev], m_amrGrids[lev]);
	}
    }

}




/// given the current cell centred velocity field, compute a face centred velocity field
void 
AmrIce::computeFaceVelocity(Vector<LevelData<FluxBox>* >& a_faceVelAdvection, 
			    Vector<LevelData<FluxBox>* >& a_faceVelTotal,
			    Vector<LevelData<FluxBox>* >& a_diffusivity,
			    Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
			    Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel) 
{
  CH_assert(m_constitutiveRelation != NULL);

  LevelData<FArrayBox>* cellDiffusivity = NULL;

  for (int lev = 0; lev <= m_finest_level; lev++)
    {
      LevelData<FArrayBox>* crseVelPtr = (lev > 0)?m_velocity[lev-1]:NULL;
      int nRefCrse = (lev > 0)?m_refinement_ratios[lev-1]:1;

      

      LevelData<FArrayBox>* crseCellDiffusivityPtr = 
	(lev > 0)?cellDiffusivity:NULL;

      cellDiffusivity = new LevelData<FArrayBox>(m_amrGrids[lev],1,IntVect::Unit);
      
      CH_assert(cellDiffusivity != NULL);
      CH_assert(a_faceVelAdvection[lev] != NULL);
      CH_assert(a_faceVelTotal[lev] != NULL);
      CH_assert(a_diffusivity[lev] != NULL);
      CH_assert(a_layerXYFaceXYVel[lev] != NULL);
      CH_assert(a_layerSFaceXYVel[lev] != NULL);
      CH_assert(m_velocity[lev] != NULL);
      CH_assert(m_vect_coordSys[lev] != NULL);
      CH_assert(m_A[lev] != NULL);
      CH_assert(m_sA[lev] != NULL);
      CH_assert(m_bA[lev] != NULL);
      
      IceUtility::computeFaceVelocity
       	(*a_faceVelAdvection[lev], *a_faceVelTotal[lev], *a_diffusivity[lev],
	 *cellDiffusivity,*a_layerXYFaceXYVel[lev], *a_layerSFaceXYVel[lev],
	 *m_velocity[lev],*m_vect_coordSys[lev], m_thicknessIBCPtr, 
	 *m_A[lev], *m_sA[lev], *m_bA[lev], 
	 crseVelPtr,crseCellDiffusivityPtr, nRefCrse, 
	 m_constitutiveRelation, m_additionalVelocity, m_diffusionTreatment == IMPLICIT);

      if (crseCellDiffusivityPtr != NULL)
	delete crseCellDiffusivityPtr;

    }

  if (cellDiffusivity != NULL)
    delete cellDiffusivity;
}

// /// compute div(vel*H) at a given time
// void
// AmrIce::computeDivThicknessFlux(Vector<LevelData<FArrayBox>* >& a_divFlux,
//                                 Vector<LevelData<FluxBox>* >& a_flux,
//                                 Vector<LevelData<FArrayBox>* >& a_thickness,
//                                 Real a_time, Real a_dt)
// {

//   //Vector<LevelData<FluxBox>* > faceVel(m_finest_level+1, NULL);
  
//   for (int lev=0; lev<= m_finest_level; lev++)
//     {
//       const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
//       // construct face-centered velocity field
//       LevelData<FluxBox> faceVel(levelGrids, 1, IntVect::Unit);
                                 
//       if (lev > 0)
//         {
//           int nVelGhost = m_velocity[lev]->ghostVect()[0];
          
//           PiecewiseLinearFillPatch velFiller(levelGrids, 
//                                              m_amrGrids[lev-1],
//                                              m_velocity[0]->nComp(), 
//                                              m_amrDomains[lev-1],
//                                              m_refinement_ratios[lev-1],
//                                              nVelGhost);
          
//           // since we're not subcycling, don't need to interpolate in time
//           Real time_interp_coeff = 0.0;
//           velFiller.fillInterp(*m_velocity[lev],
//                                *m_velocity[lev-1],
//                                *m_velocity[lev-1],
//                                time_interp_coeff,
//                                0, 0, m_velocity[0]->nComp());
          
//         } // end if lev > 0
//       m_velocity[lev]->exchange();

//       // average velocities to faces
//       CellToEdge(*m_velocity[lev],faceVel);
//       faceVel.exchange();
      
//       // flux = faceVel*faceH      
//       LevelData<FluxBox>& levelFlux = *a_flux[lev];
//       LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
//       LevelData<FluxBox>& faceH = levelCoords.getFaceH();

//       DataIterator dit = levelGrids.dataIterator();
//       for (dit.begin(); dit.ok(); ++dit)
//         {
//           FluxBox& thisflux = levelFlux[dit];
//           FluxBox& thisVel = faceVel[dit];
//           FluxBox& thisH = faceH[dit];

//           for (int dir=0; dir<SpaceDim; dir++)
//             {
//               thisflux[dir].copy(thisVel[dir]);
//               thisflux[dir].mult(thisH[dir]);
//             }
//         }

//       // average fluxes to coarser levels, if needed
//       if (lev>0)
//         {
//           CoarseAverageFace faceAverager(m_amrGrids[lev],
//                                          1, m_refinement_ratios[lev-1]);
//           faceAverager.averageToCoarse(*a_flux[lev-1], *a_flux[lev]);
          
//         }
      
//     } // end loop over levels

//   // now compute div(flux)
  
//   // compute div(F) and add source term
//   setToZero(m_volumeThicknessFlux);
//   for (int lev=0; lev<= m_finest_level; lev++)
//     {
//       DisjointBoxLayout& levelGrids = m_amrGrids[lev];
//       LevelData<FluxBox>& levelFlux = *a_flux[lev];
//       LevelData<FArrayBox>& levelDiv = *a_divFlux[lev];
//       LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);

//       LevelData<FArrayBox>& surfaceThicknessSource = *m_surfaceThicknessSource[lev];
//       m_surfaceFluxPtr->surfaceThicknessFlux(surfaceThicknessSource, *this, lev, a_dt);
      
//       LevelData<FArrayBox>& basalThicknessSource = *m_basalThicknessSource[lev];
//       m_basalFluxPtr->surfaceThicknessFlux(basalThicknessSource, *this, lev, a_dt);

//       const RealVect& dx = levelCoords.dx();          

//       DataIterator dit = levelGrids.dataIterator();
      
//       for (dit.begin(); dit.ok(); ++dit)
//         {
//           const Box& gridBox = levelGrids[dit];
//           FArrayBox& thisDiv = levelDiv[dit];
          
//           FluxBox& thisFlux = levelFlux[dit];
//           thisDiv.setVal(0.0);
          
//           // loop over directions and increment with div(F)
//           for (int dir=0; dir<SpaceDim; dir++)
//             {
//               // use the divergence from 
//               // Chombo/example/fourthOrderMappedGrids/util/DivergenceF.ChF
//               FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
//                               CHF_FRA(thisDiv),
//                               CHF_BOX(gridBox),
//                               CHF_CONST_REAL(dx[dir]),
//                               CHF_INT(dir));
//             }

//           // add in thickness source here
//           thisDiv.minus(surfaceThicknessSource[dit], gridBox,0,0,1);
// 	  thisDiv.minus(basalThicknessSource[dit], gridBox,0,0,1);
// 	  thisDiv.minus(volumeThicknessSource[dit], gridBox,0,0,1);
//           //thisDiv *= -1*a_dt;
//         } // end loop over grids
//     } // end loop over levels
  
// }

// increment phi := phi + dt*dphi
void
AmrIce::incrementWithDivFlux(Vector<LevelData<FArrayBox>* >& a_phi,
                             const Vector<LevelData<FArrayBox>* >& a_dphi,
                             Real a_dt)
{
  for (int lev=0; lev<a_phi.size(); lev++)
    {
      LevelData<FArrayBox>& levelPhi = *a_phi[lev];
      const LevelData<FArrayBox>& level_dPhi = *a_dphi[lev];

      DataIterator dit = levelPhi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          levelPhi[dit].plus(level_dPhi[dit], a_dt);
        }
    }
}
 

// increment coordSys with new thickness
void
AmrIce::updateCoordSysWithNewThickness(const Vector<LevelData<FArrayBox>* >& a_thickness)
{
  CH_assert(a_thickness.size() >= m_finest_level);
  
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      const LevelData<FArrayBox>& levelH = *a_thickness[lev];
      LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      LevelData<FArrayBox>& levelCS_H = levelCS.getH();
      DataIterator dit = levelH.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisH = levelCS_H[dit];
          thisH.copy(levelH[dit]);
        }
      {
	LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
	int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
	levelCS.recomputeGeometry(crseCoords, refRatio);
      }
    } // end loop over levels      
}

void
AmrIce::setIceFrac(const LevelData<FArrayBox>& a_thickness, int a_level)
{
  // initialize fraction to 1 if H>0, 0 o/w...
  DataIterator dit = m_iceFrac[a_level]->dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisFrac = (*m_iceFrac[a_level])[dit];
      thisFrac.setVal(0.0);
      const FArrayBox& thisH = a_thickness[dit];
      BoxIterator bit(thisFrac.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          if (thisH(iv,0) > 0) thisFrac(iv,0) = 1.0;
        }
    }
}

void
AmrIce::updateIceFrac(LevelData<FArrayBox>& a_thickness, int a_level)
{
  // set ice fraction to 0 if no ice in cell...

  // "zero" thickness value
  // Check that this rountine doesn't interfer with diagnostics (endTimestepDiagnostics).
  Real ice_eps = 1.0e-6;
  DataIterator dit = m_iceFrac[a_level]->dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisFrac = (*m_iceFrac[a_level])[dit];
      FArrayBox& thisH = a_thickness[dit];
      BoxIterator bit(thisFrac.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
	  if (thisH(iv,0) < ice_eps || thisFrac(iv,0) < ice_eps) 
            {
              thisFrac(iv,0) = 0.0;
              thisH(iv,0) = 0.0;
            }          
        }
    }
}

/// update covered cells, ghost cells etc, of a_iceFrac
void
AmrIce::updateInvalidIceFrac(Vector<LevelData<FArrayBox> * > a_iceFrac)
{
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      LevelData<FArrayBox>& levelFrac = *a_iceFrac[lev];
      const DisjointBoxLayout& fracGrids = levelFrac.getBoxes();
      Real levelDx = m_amrDx[lev];
      
      if (lev < m_finest_level)
	{
	  CoarseAverage av(m_amrGrids[lev+1], 1, m_refinement_ratios[lev]);
	  av.averageToCoarse(levelFrac,*a_iceFrac[lev+1]);	    }
      
      if (lev > 0)
	{
	  int nGhost = levelFrac.ghostVect()[0];
	  PiecewiseLinearFillPatch pwl(m_amrGrids[lev], m_amrGrids[lev-1],
				       1, m_amrDomains[lev-1], m_refinement_ratios[lev-1],
				       nGhost);
	  pwl.fillInterp(levelFrac,*m_iceFrac[lev-1],*m_iceFrac[lev-1],0.0, 0, 0, 1);  
	}
      m_thicknessIBCPtr->setIceFractionBCs(levelFrac,m_amrDomains[lev],dx(lev), time(), m_dt);
      levelFrac.exchange();
    }
}

#define ICE_FRAC_EPS 1.0e-6

// /// update real-valued ice fraction and thickness source according to the calving rate
// void
// AmrIce::advectIceFracUpstream
// (Vector<LevelData<FArrayBox>* >& a_iceFrac,
//  Vector<LevelData<FArrayBox>* >& a_thicknessSource,
//  const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
//  const Vector<LevelData<FArrayBox>* >& a_calvingRate,
//  Real a_subDt, Real a_dt)
// {

//   // for (int lev=0; lev<= m_finest_level; lev++)
//   //   {
//   //     LevelData<FArrayBox>& levelFrac = *a_iceFrac[lev];
//   //     const DisjointBoxLayout& fracGrids = levelFrac.getBoxes();
//   //     Real levelDx = m_amrDx[lev];
//   //     //LevelData<FluxBox> faceRate( m_amrGrids[lev], 1, IntVect::Unit);
//   //     //CellToEdge ( (*a_calvingRate[lev]), faceRate);

      
//   //     //Compute a source term from the calving rate: only applies at front cells
//   //     //where iceFrac > epsilon and iceFrac(a neighbour) < epsilon,
//   //     //and add to the volume sourc
//   //     //modify iceFrac accordingly
//   //     for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
//   // 	{
//   // 	  const Box& box = m_amrGrids[lev][dit];
//   // 	  FArrayBox& src = (*a_thicknessSource[lev])[dit];
//   // 	  //src.setVal(0.0);
//   // 	  FArrayBox& frac = levelFrac[dit];
//   // 	  const FArrayBox& thk = (*m_old_thickness[lev])[dit];
//   // 	  FArrayBox dh(frac.box(), 1); dh.copy(frac);
//   // 	  const FArrayBox& oldfrac = dh; //might as well
//   // 	  for (int dir = 0; dir < SpaceDim; dir++)
//   // 	    {
//   // 	      //const FArrayBox& rate = faceRate[dit][dir];
//   // 	      const FArrayBox& rate = (*a_calvingRate[lev])[dit];
//   // 	      //FArrayBox erate(rate.box(),1.0);
//   // 	      //erate.setVal(1.0);
//   // 	      Real eps = 0.5 * m_cfl;
	      
//   // 	      FORT_ABLATEFRAC(CHF_FRA1(frac,0),CHF_CONST_FRA1(oldfrac,0),
//   // 			      CHF_CONST_FRA1(rate,0),
//   // 			      CHF_REAL(levelDx),
//   // 			      CHF_REAL(a_subDt),
//   // 			      CHF_REAL(eps),
//   // 			      CHF_BOX(box),
//   // 			      CHF_INT(dir));
//   // 	      Real et = a_subDt * 1.0e-3;
//   // 	      FArrayBox t(frac.box(),1);
//   // 	      t.copy(frac);
//   // 	      FORT_ABLATEFRAC(CHF_FRA1(frac,0),CHF_CONST_FRA1(t,0),
//   // 			      CHF_CONST_FRA1(rate,0),
//   // 			      CHF_REAL(levelDx),
//   // 			      CHF_REAL(et),
//   // 			      CHF_REAL(eps),
//   // 			      CHF_BOX(box),
//   // 			      CHF_INT(dir));


	      
//   // 	    } // loop over directions
//   // 	  dh -= frac; //ice fraction lost due to ablation
//   // 	  dh *= thk; // ice lost to ablation
//   // 	  dh /= a_dt; // rate of ice lost to ablation
//   // 	  src -= dh; // include rate of ice loss in source term
//   // 	} // loop over boxes
//   //   } // end loop over levels
// }

// /// update real-valued ice fraction through advection from neighboring cells
// void
// AmrIce::advectIceFracDownstream(Vector<LevelData<FArrayBox>* >& a_iceFrac,
// 				Vector<LevelData<FArrayBox>* >& a_thicknessSource,
// 				const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
// 				const Vector<LevelData<FluxBox>* >& a_thicknessFlux,
// 				const Vector<LevelData<FArrayBox>* >& a_calvingRate,
// 				Real a_subDt, Real a_dt)
// {


  

  
//   for (int lev=0; lev<= m_finest_level; lev++)
//     {
//       LevelData<FArrayBox>& levelFrac = *a_iceFrac[lev];
//       const LevelData<FluxBox>& levelFaceVel = *a_faceVelAdvection[lev];
//       const DisjointBoxLayout& fracGrids = levelFrac.getBoxes();
//       Real levelDx = m_amrDx[lev];
      
//       for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
// 	{
// 	  // only update valid cells
// 	  const Box& gridBox = fracGrids[dit];
// 	  FArrayBox& thisFrac = levelFrac[dit];
// 	  FArrayBox oldfrac(thisFrac.box(),1); oldfrac.copy(thisFrac);
// 	  FArrayBox dh(gridBox,1); dh.setVal(0.0);
// 	  const FluxBox& thisFaceVel = levelFaceVel[dit];
// 	  const FluxBox& faceHFlux = levelFaceVel[dit];
// 	  int dbg = 0;
// 	  for (int dir=0; dir<SpaceDim; dir++)
// 	    {
// 	      const FArrayBox& u_ice_face = thisFaceVel[dir];
// 	      const FArrayBox& hflux = (*a_thicknessFlux[lev])[dit][dir];
// 	      const FArrayBox& h = (*m_old_thickness[lev])[dit];
// 	      Box faceBox = gridBox;
// 	      faceBox.surroundingNodes(dir);
// 	      FArrayBox u_calv_face(faceBox, 1);
// 	      u_calv_face.setVal(0.0);
	      
// 	      Real eps = 0.0005*m_cfl;
// 	      eps = 1.0e-8;
// 	      const FArrayBox& rate = (*a_calvingRate[lev])[dit];
	      
	      
// 	      // FORT_ABLATERATE(CHF_FRA1(u_calv_face,0),
// 	      // 		      CHF_CONST_FRA1(oldfrac,0),
// 	      // 		      CHF_CONST_FRA1(rate,0),
// 	      // 		      CHF_REAL(eps),
// 	      // 		      CHF_BOX(faceBox),
// 	      // 		      CHF_INT(dir));
	      
// 	      FORT_ADVECTFRAC(CHF_FRA1(thisFrac,0),
// 			      CHF_FRA1(dh,0),
// 			      CHF_CONST_FRA1(oldfrac,0),
// 			      CHF_CONST_FRA1(u_ice_face,0),
// 			      CHF_CONST_FRA1(u_calv_face,0),
// 			      CHF_CONST_FRA1(hflux,0),
// 			      CHF_CONST_FRA1(h,0),
// 			      CHF_REAL(levelDx),
// 			      CHF_REAL(a_subDt),
// 			      CHF_REAL(eps),
// 			      CHF_BOX(gridBox),
// 			      CHF_INT(dir));
// 	    }  // end loop over directions
// 	  // thickness sink
// 	  // on return from FORT_ADVECTFRAC, dh = - df/dt * dt due to void growth,
// 	  // we want dh/dt = - h / f df/dt (h/f is constant)
// 	  // dh *=  (*m_old_thickness[lev])[dit];
// 	  // oldfrac += TINY_NORM;
// 	  // dh /= oldfrac;
// 	  if (m_frac_sources)
// 	    {
// 	      // source will be multiplied by f, so
// 	      //thisFrac += TINY_NORM;
// 	      dh /= thisFrac;
// 	    }
// 	  dh /= a_dt;
// 	  //(*a_thicknessSource[lev])[dit] += dh;
	  
// 	  dbg++;
// 	} // end loop over boxes
//     } // end loop over levels 
// }

// /// update real-valued ice fraction through advection from neighboring cells
// void
// AmrIce::advectIceFrac(Vector<LevelData<FArrayBox>* >& a_iceFrac,
// 		      Vector<LevelData<FArrayBox>* >& a_thicknessSource,
//                       const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
// 		      const Vector<LevelData<FluxBox>* >& a_thicknessFlux,
//                       Real a_dt)
// {

//   ///compute and store the cell-centered calving rate
//   Vector<LevelData<FArrayBox> *> calvingRate(m_finest_level + 1, NULL);
//   for (int lev=0; lev<= m_finest_level; lev++)
//     {
//       calvingRate[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Unit);
//       (*m_calvingModelPtr).getCalvingRate(*calvingRate[lev], *this, lev);
      
//       if (lev > 0)
// 	{
// 	  PiecewiseLinearFillPatch pwl(m_amrGrids[lev], m_amrGrids[lev-1],
// 				       1, m_amrDomains[lev-1], m_refinement_ratios[lev-1],
// 				       1);
// 	  pwl.fillInterp(*calvingRate[lev],*calvingRate[lev-1],*calvingRate[lev-1],0.0, 0, 0, 1);  
// 	}
     
//       calvingRate[lev]->exchange();
      
//       //while we are here, initialize thickness source
//       for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
// 	    {
// 	      FArrayBox& src = (*a_thicknessSource[lev])[dit];
// 	      src.setVal(0.0);
// 	    }	  
//     }
//   Real max_calving_rate = computeNorm(calvingRate, refRatios() , m_amrDx[0], Interval(0,0), 0);
//   Real dt_max = m_cfl * m_amrDx[finestLevel()] / max_calving_rate;
  
//   //advection step. 
//   int n_sub_step = int (a_dt / dt_max);
//   for (int sub_step = 0; sub_step < n_sub_step; sub_step++)
//     {
//       Real sub_step_dt = a_dt / Real(n_sub_step);
//       updateInvalidIceFrac(a_iceFrac);
//       applyCalvingCriterion(CalvingModel::PostThicknessAdvection);
//       advectIceFracDownstream(a_iceFrac, a_thicknessSource,
// 			      a_faceVelAdvection,  a_thicknessFlux, calvingRate, sub_step_dt, a_dt);
//     }
    
 
//   // clean up
//   updateInvalidIceFrac(a_iceFrac);
//   for (int lev=0; lev<= m_finest_level; lev++)
//     {    
//       if ( calvingRate[lev] ) delete calvingRate[lev]; 
//     }

// }



/// update real-valued ice fraction through advection from neighboring cells
void
AmrIce::advectIceFrac(Vector<LevelData<FArrayBox>* >& a_iceFrac,
                      const Vector<LevelData<FluxBox>* >& a_faceVelAdvection,
                      Real a_dt)
{
  // for now, set fill threshold to be (1-cfl) on the theory 
  // that we want to declare a cell full before it actually over-fills
  Real fillThreshold = (1.0 - m_cfl);
  Real epsFrac = 1.0e-8;
  Real epsVel = 1.0e-8;
  int calvingNormalToFront = 1;

  for (int lev=0; lev<= m_finest_level; lev++)
    {
      LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      LevelData<FArrayBox>& levelFrac = *a_iceFrac[lev];
      levelFrac.exchange();
      const LevelData<FluxBox>& levelFaceVel = *a_faceVelAdvection[lev];
      const LevelData<FArrayBox>& levelCCVel = *m_velocity[lev];
      //LevelData<FArrayBox>& levelFrontVel = *a_frontVel[lev];
      const DisjointBoxLayout& fracGrids = levelFrac.getBoxes();
      const LevelData<FArrayBox>& levelTopg = levelCoords.getTopography();
      const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
      const LevelData<FArrayBox>& levelThck =  levelCoords.getH();
      Real levelDx = m_amrDx[lev];
      
      LevelData<FArrayBox> levelOldFrac(fracGrids, 1, levelFrac.ghostVect());
      LevelData<FArrayBox> levelCopyFrac(fracGrids, 1, levelFrac.ghostVect());
      LevelData<FArrayBox> levelTmpVel(fracGrids, 2, levelFrac.ghostVect());
      LevelData<FArrayBox> levelFrontVel(fracGrids, 2, levelFrac.ghostVect());
      
      for (DataIterator dit(fracGrids); dit.ok(); ++dit)
        {
	  //FArrayBox& oldFrac= levelFrac[dit];
	  levelOldFrac[dit].copy(levelFrac[dit]);
	  levelCopyFrac[dit].setVal(0.0);
	  levelTmpVel[dit].setVal(0.0);
	}
      levelOldFrac.exchange();
      LevelData<FArrayBox> levelCalvingRate(fracGrids, 1, levelFrac.ghostVect());
      (*m_calvingModelPtr).getCalvingRate(levelCalvingRate, *this, lev);

      DataIterator dit = levelFrac.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          // only update valid cells
          const Box& gridBox = fracGrids[dit];
	  const FArrayBox& oldFrac = levelOldFrac[dit];
          const FluxBox& thisFaceVel = levelFaceVel[dit];
          const FArrayBox& thisCCVel = levelCCVel[dit];

	  FArrayBox& calvingRate = levelCalvingRate[dit];

	  bool matt_option = false;
	  ParmParse pp("amr");
	  pp.query("matt_option",  matt_option);
	  if (matt_option)
	    {
	      const FArrayBox& thck = levelThck[dit];
	      const FArrayBox& frac = levelFrac[dit];
	      FArrayBox t(calvingRate.box(), 1); // allocate some memory for the calculation
	      t.copy(thck); // t = h
	      t += 1.0e-10; // t = h + epsilon
	      calvingRate *= frac; // Uc = Uc * f
	      calvingRate /= t; // Uc = Uc * f/(h+ epsilon)
	    }
	  
	  
	  const FArrayBox& topg = levelTopg[dit];
	  const BaseFab<int>& mask = levelMask[dit];
          FArrayBox& thisFrontVel = levelFrontVel[dit];
          FArrayBox& thisTmpVel = levelTmpVel[dit];

          // (DFM -- 1/4/21) do it this way so that this will still work in 1D
          const FArrayBox& thisFaceVelx = thisFaceVel[0];
          // this points at the x-faces in 1D, at the y-faces in 2D
          const FArrayBox& thisFaceVely = thisFaceVel[SpaceDim-1];
          
	  FORT_FINDFRONTVEL(CHF_CONST_FRA1(oldFrac,0),
			      CHF_CONST_FRA1(thisFaceVelx,0),
			      CHF_CONST_FRA1(thisFaceVely,0),
			      CHF_CONST_FRA1(calvingRate,0),
			      CHF_CONST_FRA1(topg,0),
			      CHF_CONST_FIA1(mask,0),
			      CHF_CONST_FRA(thisCCVel),
			      CHF_FRA(thisTmpVel),
			      CHF_FRA(thisFrontVel),
			      CHF_CONST_INT(calvingNormalToFront),
			      CHF_CONST_REAL(epsFrac),
			      CHF_CONST_REAL(epsVel),
			      CHF_CONST_REAL(levelDx),
			      CHF_CONST_INT(m_cur_step),
			      CHF_BOX(gridBox));
	  

	} // end loop over boxes

      levelFrontVel.exchange();

      m_thicknessIBCPtr->velocityGhostBC( levelFrontVel , *m_vect_coordSys[lev],m_amrDomains[lev],m_time);

      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = fracGrids[dit];
          FArrayBox& thisFrac = levelFrac[dit];
	  const FArrayBox& oldFrac = levelOldFrac[dit];
          const FArrayBox& thisFrontVel = levelFrontVel[dit];

          for (int dir=0; dir<SpaceDim; dir++)
            {

	    FORT_ADVECTFRAC(CHF_FRA1(thisFrac,0),
			      CHF_CONST_FRA1(oldFrac,0),
                              CHF_CONST_FRA1(thisFrontVel,dir),
                              CHF_CONST_REAL(levelDx),
                              CHF_CONST_REAL(a_dt),
			      CHF_CONST_REAL(epsFrac),
			      CHF_CONST_REAL(epsVel),
                              CHF_BOX(gridBox),
                              CHF_INT(dir));

            } // end loop over directions
	} // end loop over boxes
      levelFrac.exchange();
      m_thicknessIBCPtr->setIceFractionBCs(levelFrac,m_amrDomains[lev],dx(lev), time(), m_dt);

      int sumOutOfRange=0;
      int MaxIt=m_max_box_size*(lev+1);
      for (int ii=0; ii<MaxIt; ii++)
	{
	  if (sumOutOfRange > 0 || ii==0)
	    {
	      sumOutOfRange=0;
	      for (dit.begin(); dit.ok(); ++dit)
		{
		  const Box& gridBox = fracGrids[dit];
		  FArrayBox& thisFrac = levelFrac[dit];
		  const FArrayBox& oldFrac = levelOldFrac[dit];
		  const FArrayBox& thisFrontVel = levelFrontVel[dit];

		  FArrayBox& copyFrac=levelCopyFrac[dit];
		  copyFrac.copy(thisFrac);
		  int outOfRange=0;

		  FORT_GETFRAC(CHF_FRA1(thisFrac,0),
			       CHF_CONST_FRA1(oldFrac,0),
			       CHF_CONST_FRA1(copyFrac,0),
			       CHF_CONST_FRA(thisFrontVel),
			       CHF_CONST_REAL(epsFrac),
			       CHF_CONST_REAL(epsVel),
			       CHF_CONST_INT(m_cur_step),
			       CHF_CONST_INT(ii),
			       CHF_INT(outOfRange),
			       CHF_CONST_INT(lev),
                               CHF_CONST_INT(MaxIt),
			       CHF_BOX(gridBox));
		  sumOutOfRange += outOfRange;
		} // end loop over boxes
	    }

	  m_thicknessIBCPtr->setIceFractionBCs(levelFrac,m_amrDomains[lev],dx(lev), time(), m_dt);
	} // end loop over ii
      levelFrac.exchange();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  const Box& gridBox = fracGrids[dit];
	  FArrayBox& thisFrac = levelFrac[dit];
	  const FArrayBox& oldFrac = levelOldFrac[dit];
	  const FArrayBox& thisFrontVel = levelFrontVel[dit];
	  //FArrayBox& thisFrontVel = levelFrontVel[dit];
	  FArrayBox& thisTmpVel = levelTmpVel[dit];

	  FORT_ISOLATEDFRAC(CHF_FRA1(thisFrac,0),
			    CHF_CONST_FRA1(oldFrac,0),
			    CHF_CONST_FRA(thisFrontVel),
			    CHF_CONST_INT(m_cur_step),
			    CHF_CONST_INT(lev),
			    CHF_BOX(gridBox));
	  //thisFrontVel.copy(thisTmpVel);
	} // end loop over boxes
      levelFrac.exchange();
    } // end loop over levels

  //coarse average from finer levels & exchange
  for (int lev = m_finest_level; lev >= 0 ; --lev)
    {
      if (lev > 0)
	{
	  CoarseAverage averager(m_amrGrids[lev],
				 a_iceFrac[lev]->nComp(),
				 m_refinement_ratios[lev-1]);
	  
	  averager.averageToCoarse(*a_iceFrac[lev-1], *a_iceFrac[lev]);
	}
    }
  
}


// compute timestep
Real 
AmrIce::computeDt()
{
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::computeDt" << endl;
    }

  if (m_fixed_dt > TIME_EPS)
    return m_fixed_dt;

  Real dt = 1.0e50;
  Real maxVelAll = 0.0;
  for (int lev=0; lev<= finestTimestepLevel(); lev++)
    {

      Real dtLev = dt;
      Real maxVelLev = maxVelAll;
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      const LevelData<FluxBox>& levelVel = *m_faceVelAdvection[lev]; 
      DataIterator levelDit = levelVel.dataIterator();
      for (levelDit.reset(); levelDit.ok(); ++levelDit)
	{
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      int p = 0;
	      Box faceBox = levelGrids[levelDit];
	      faceBox.surroundingNodes(dir);
	      Real maxVel = 1.0e-10 + levelVel[levelDit][dir].norm(faceBox,p, 0, 1);
	      Real localDt = m_amrDx[lev]/maxVel;
	      dtLev = min(dtLev, localDt);
	      maxVelLev = std::max(maxVelLev, maxVel);
	    }
	}
      
      if (m_diffusionTreatment == EXPLICIT){
	MayDay::Error("diffusion_treatment == explicit not supported now : use none");
      }
      dt = min(dt, dtLev);
      maxVelAll = std::max(maxVelAll, maxVelLev);
    }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL,
			     MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("communication error on norm");
    }
  dt = tmp;

  result = MPI_Allreduce(&maxVelAll, &tmp, 1, MPI_CH_REAL,
			 MPI_MAX, Chombo_MPI::comm);
  
if (result != MPI_SUCCESS)
    {
      MayDay::Error("communication error on norm");
    }
  maxVelAll = tmp;
#endif

  if (m_cur_step == 0)
    {
      dt *= m_initial_cfl;
    } 
  else 
    {
      dt *= m_cfl;
    }

  // also check to see if max grow rate applies
  // (m_dt > 0 test screens out initial time, when we set m_dt to a negative 
  // number by default)
  // Use the value stored in m_stable_dt in case dt was altered to hit a plot interval
  // m_max_dt_grow < 0 implies that we don't enforce this.
  if ((m_max_dt_grow > 0) && (dt > m_max_dt_grow*m_stable_dt) && (m_stable_dt > 0) )
    dt = m_max_dt_grow*m_stable_dt;
  
  if (m_timeStepTicks){
    // reduce time step to integer power of two
    dt = std::pow(2.0, std::floor(std::log(dt)/std::log(two)));
    
  }
  
  m_stable_dt = dt;
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::computeDt dt = " << dt << ", max vel = " << maxVelAll <<  endl;
    }

  if (maxVelAll > m_velocity_exit || dt < TIME_EPS)
    {
      std::stringstream ss;
      ss << " max(vel) = " << maxVelAll
	 << " > amr.velocity_exit =  " << m_velocity_exit 
	 << " || dt = " << dt
	 << " < TIME_EPS = " << TIME_EPS;
        
      m_plot_prefix = m_plot_prefix + std::string("_error_.");
      writePlotFile();
      pout() << " AmrIce::computeDt exit because " << ss.str() << endl;
      MayDay::Error(ss.str().c_str());
    }

  return dt;

}

Real 
AmrIce::computeInitialDt()
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::computeInitialDt" << endl;
    }


  // for now, just call computeDt;
  Real dt = computeDt();
  return dt;
}



//determine the grouding line proximity 
/**

   Solves the elliptic problem 
   a * phi - b* grad^2 phi = 0;
   with natural boundary conditions.

   for grounded ice, a = 10^5 and b = 1
   for floating ice, s = 0 and b = 1
*/
void AmrIce::updateGroundingLineProximity() const
{

  CH_TIME("AmrIce::updateGroundingLineProximity");

  if (m_groundingLineProximity_valid)
    return;

  if (m_groundingLineProximity.size() < m_finest_level + 1)
    {
      m_groundingLineProximity.resize(m_finest_level + 1, NULL);
    }

  if (s_verbosity > 0)
    {
      pout() << "AmrIce::updateGroundingLineProximity() max level = " << m_finest_level << " " << endl; 
    }

  //Natural boundary conditions
  BCHolder bc(ConstDiriNeumBC(IntVect::Zero, RealVect::Zero,
  			      IntVect::Zero, RealVect::Zero));

  //BCHolder bc(ConstDiriNeumBC(IntVect(0,0), RealVect(-1.0,-1.0),
  //			      IntVect(0,0), RealVect(1.0,1.0)));

  Vector<RefCountedPtr<LevelData<FArrayBox> > > a(m_finest_level + 1);
  Vector<RefCountedPtr<LevelData<FluxBox> > > b(m_finest_level + 1);
  Vector<LevelData<FArrayBox>* > rhs(m_finest_level+ 1,NULL);
  Vector<DisjointBoxLayout> grids(finestTimestepLevel() + 1);
  Vector<ProblemDomain> domains(finestTimestepLevel() + 1);
  Vector<RealVect> dx(finestTimestepLevel() + 1);

  for (int lev=0; lev <= m_finest_level; ++lev)
    {
      dx[lev] = m_amrDx[lev]*RealVect::Unit;
      domains[lev] = m_amrDomains[lev];

      const LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
      const LevelData<BaseFab<int> >& levelMask = levelCS.getFloatingMask();
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      a[lev] = RefCountedPtr<LevelData<FArrayBox> >
 	(new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero));
      b[lev] = RefCountedPtr<LevelData<FluxBox> >
 	(new LevelData<FluxBox>(levelGrids, 1, IntVect::Zero));
      rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
      
      grids[lev] = levelGrids;

     

      if (m_groundingLineProximity[lev] != NULL)
	{
	  delete m_groundingLineProximity[lev];
	  m_groundingLineProximity[lev] = NULL;
	}
      m_groundingLineProximity[lev] =  new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
      
      LevelData<FArrayBox>& levelPhi = *m_groundingLineProximity[lev];

      const Real& crseDx = m_amrDx[0];
      Real crseDxSq = crseDx*crseDx;

      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
 	{
	  FluxBox& B = (*b[lev])[dit];
	  
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      B[dir].setVal(crseDxSq);
	    }

 	  FArrayBox& r =  (*rhs[lev])[dit];
 	  r.setVal(0.0);
	  FArrayBox& A =  (*a[lev])[dit];
	  A.setVal(0.0);
	  FArrayBox& phi = levelPhi[dit];
	  phi.setVal(0.0);

 	  const BaseFab<int>& mask = levelMask[dit];
 	  const Box& gridBox = levelGrids[dit];
	  //	  const FArrayBox& u = (*m_velocity[lev])[dit];

	  Real AcoefF = crseDx / m_groundingLineProximityScale;
	  Real AcoefG = 1.0 ;
	  if (m_groundingLineProximityCalcType > 0)
	    {
	      AcoefF = crseDx / m_groundingLineProximityScale;
	      AcoefF *= AcoefF;
	      
	    }

 	  for (BoxIterator bit(gridBox);bit.ok();++bit)
 	    {
 	      const IntVect& iv = bit();
 	      if (mask(iv) == GROUNDEDMASKVAL )
 		{
 		  A(iv) = AcoefG;
		  r(iv) = AcoefG;
		  
 		} 
	      else
		{
		  
		  A(iv) = AcoefF;
		  r(iv) = 0.0;
		}
	      
 	    }
	  phi.copy(r);
 	}

      rhs[lev]->exchange();
      levelPhi.exchange();
      m_groundingLineProximity[lev]->exchange();
      a[lev]->exchange();
      b[lev]->exchange();
    }


  VCAMRPoissonOp2Factory* poissonOpFactory = new VCAMRPoissonOp2Factory;
  poissonOpFactory->define(domains[0], grids , m_refinement_ratios,
 			   m_amrDx[0], bc, 1.0, a,  1.0 , b);
  RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > 
    opFactoryPtr(poissonOpFactory);

  MultilevelLinearOp<FArrayBox> poissonOp;
  poissonOp.define(grids, m_refinement_ratios, domains, dx, opFactoryPtr, 0);
    
  RelaxSolver<Vector<LevelData<FArrayBox>* > >* relaxSolver
    = new RelaxSolver<Vector<LevelData<FArrayBox>* > >();

  relaxSolver->define(&poissonOp,false);
  relaxSolver->m_verbosity = s_verbosity;
  relaxSolver->m_normType = 0;
  relaxSolver->m_eps = 1.0e-8;
  relaxSolver->m_imax = 12;
  relaxSolver->m_hang = 0.05;
  relaxSolver->solve(m_groundingLineProximity,rhs);

  delete(relaxSolver);

#ifdef DUMP_PROXIMITY
  std::string file("proximity.2d.hdf5");
  Real dt = 0.0; 
  Real time = 0.0;
  Vector<std::string> names(1,"proximity");
  WriteAMRHierarchyHDF5(file ,grids, m_groundingLineProximity ,names, m_amrDomains[0].domainBox(),
  			m_amrDx[0], dt, m_time, m_refinement_ratios, m_groundingLineProximity.size());
#endif
  
  for (int lev=0; lev <= m_finest_level ; ++lev)
    {
      if (rhs[lev] != NULL)
 	{
 	  delete rhs[lev];
	  rhs[lev] = NULL;
 	}
    }

  m_groundingLineProximity_valid = true;
}

//access the viscous tensor (cell-centered)
const LevelData<FArrayBox>* AmrIce::viscousTensor(int a_level) const
{
  //updateViscousTensor();
  CH_assert(m_viscousTensor_valid);

  if (!(m_viscousTensorCell.size() > a_level))
    {
      std::string msg("AmrIce::viscousTensor !(m_viscousTensorCell.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_viscousTensorCell.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_viscousTensorCell[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::viscousTensor m_viscousTensorCell[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;

}

//access the viscous tensor (cell-centered)
const LevelData<FArrayBox>* AmrIce::viscosityCoefficient(int a_level) const
{
  CH_assert(m_viscousTensor_valid);
  //updateViscousTensor();
  if (!(m_viscosityCoefCell.size() > a_level))
    {
      std::string msg("AmrIce::viscosityCoef !(m_viscosityCoefCell.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_viscosityCoefCell.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_viscosityCoefCell[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::viscosityCoef m_viscosityCoefCell[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;

}

const LevelData<FArrayBox>* AmrIce::surfaceThicknessSource(int a_level) const
{
  if (!(m_surfaceThicknessSource.size() > a_level))
    {
      std::string msg("AmrIce::surfaceThicknessSource !(m_surfaceThicknessSource.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_surfaceThicknessSource.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_surfaceThicknessSource[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::surfaceThicknessSource m_surfaceThicknessSource[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;
}

const LevelData<FArrayBox>* AmrIce::basalThicknessSource(int a_level) const
{
  if (!(m_basalThicknessSource.size() > a_level))
    {
      std::string msg("AmrIce::basalThicknessSource !(m_basalThicknessSource.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_basalThicknessSource.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_basalThicknessSource[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::basalThicknessSource m_basalThicknessSource[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;
}


//access the drag coefficient (cell-centered)
const LevelData<FArrayBox>* AmrIce::dragCoefficient(int a_level) const
{
  //updateViscousTensor();
  CH_assert(m_viscousTensor_valid);
  if (!(m_dragCoef.size() > a_level))
    {
      std::string msg("AmrIce::dragCoef !(m_dragCoef.size() > a_level))");
      pout() << msg << endl;
      CH_assert((m_dragCoef.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_dragCoef[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::dragCoef m_dragCoefCell[a_level] == NULL ");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;
}



//update the viscous tensor components
void AmrIce::updateViscousTensor() const
{
  CH_TIME("AmrIce::updateViscousTensor");

  if (m_viscousTensor_valid)
    return;
  
  if (m_viscousTensorCell.size() < m_finest_level + 1)
    {
      m_viscousTensorCell.resize(m_finest_level + 1, NULL);
    }
  if (m_viscosityCoefCell.size() < m_finest_level + 1)
    {
      m_viscosityCoefCell.resize(m_finest_level + 1, NULL);
    }
  if (m_dragCoef.size() < m_finest_level + 1)
    {
      m_dragCoef.resize(m_finest_level + 1, NULL);
    }
  // if (m_magnitudeBasalDrag.size() < m_finest_level + 1)
  //   {
  //     m_magnitudeBasalDrag..resize(m_finest_level + 1, NULL);
  //   }
  
  if (m_viscousTensorFace.size() < m_finest_level + 1)
    {
      m_viscousTensorFace.resize(m_finest_level + 1, NULL);
    }

 
  Vector<LevelData<FluxBox>*> faceA(m_finest_level + 1,NULL);
  Vector<RefCountedPtr<LevelData<FluxBox> > > viscosityCoef;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > dragCoef;
  Vector<LevelData<FArrayBox>* > C0(m_finest_level + 1,  NULL);

  Vector<RealVect> vdx(m_finest_level + 1);
  for (int lev =0; lev <= m_finest_level; lev++)
    {
      faceA[lev] = new LevelData<FluxBox>(m_amrGrids[lev],m_A[lev]->nComp(),IntVect::Unit);
      CellToEdge(*m_A[lev],*faceA[lev]);

      if (m_viscousTensorFace[lev] != NULL)
	{
	  delete m_viscousTensorFace[lev];m_viscousTensorFace[lev]=NULL;
	}
      m_viscousTensorFace[lev] = new LevelData<FluxBox>(m_amrGrids[lev],SpaceDim,IntVect::Unit);

      if (m_viscousTensorCell[lev] != NULL)
	{
	  delete m_viscousTensorCell[lev];m_viscousTensorCell[lev]=NULL;
	}
      m_viscousTensorCell[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],SpaceDim*SpaceDim,IntVect::Unit);
      
      if (m_dragCoef[lev] != NULL)
	{
	  delete m_dragCoef[lev]; m_dragCoef[lev] = NULL;
	}
      m_dragCoef[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],SpaceDim,IntVect::Zero);

      // if (m_magnitudeBasalDrag[lev] != NULL)
      // 	{
      // 	  delete m_[lev]; m_dragCoef[lev] = NULL;
      // 	}
      // m_dragCoef[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],SpaceDim,IntVect::Zero);

      
      if (m_viscosityCoefCell[lev] != NULL)
	{
	  delete m_viscosityCoefCell[lev]; m_viscosityCoefCell[lev] = NULL;
	}
      m_viscosityCoefCell[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],SpaceDim,IntVect::Zero);

      if (C0[lev] != NULL)
	{
	  delete C0[lev];C0[lev] = NULL;
	}
      C0[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1,m_velBasalC[0]->ghostVect());
      DataIterator dit = m_amrGrids[lev].dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          (*C0[lev])[dit].setVal(0.0);
        }
      vdx[lev] = RealVect::Unit*m_amrDx[lev];
    }

  //these parameters don't matter because we don't solve anything here. 
  Real vtopSafety = 1.0;
  int vtopRelaxMinIter = 4;
  Real vtopRelaxTol = 1.0;
  Real muMin = 0.0; 
  Real muMax = 1.23456789e+300;

  int numLevels = m_finest_level + 1;
  IceNonlinearViscousTensor state(m_amrGrids, m_refinement_ratios, m_amrDomains, vdx, m_vect_coordSys, 
				  m_velocity, m_velBasalC, C0, numLevels-1, 
				  *m_constitutiveRelation,  *m_basalFrictionRelation, *m_thicknessIBCPtr,  
				  m_A, faceA, m_time, vtopSafety, vtopRelaxMinIter, vtopRelaxTol, 
				  muMin, muMax);
  state.setState(m_velocity);
  viscosityCoef = state.mu();
  dragCoef = state.alpha();
  state.computeViscousTensorFace(m_viscousTensorFace);
  
  for (int lev =0; lev < numLevels; lev++)
    {
      
      //If a cell is adjacent to a calving front, we  set the (vertically integrated)
      //viscous tensor components at the intervening face to zero. That works well enough for the velocity solves,
      //but causes pain here because the cell average (in EdgeToCell) will end up half the value at the other face.
      for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
      	{

      	  const FArrayBox& thck = m_vect_coordSys[lev]->getH()[dit];
      	  //const FArrayBox& dsdx = m_vect_coordSys[lev]->getGradSurface()[dit];
	  const FArrayBox& usrf = m_vect_coordSys[lev]->getSurfaceHeight()[dit];
      	  const BaseFab<int>& mask = m_vect_coordSys[lev]->getFloatingMask()[dit];
      	  const Real& rhoi = m_vect_coordSys[lev]->iceDensity();
      	  //const Real& rhoo = m_vect_coordSys[lev]->waterDensity();
      	  const Real& gravity = m_vect_coordSys[lev]->gravity();
      	  //const Real rgr = rhoi * gravity * (1.0-rhoi/rhoo);
      	  //const RealVect& dx = m_vect_coordSys[lev]->dx();

      	  for (int dir = 0; dir < SpaceDim; dir++)
      	    {
      	      FArrayBox& facevt = (*m_viscousTensorFace[lev])[dit][dir];
      	      Real factor = rhoi * gravity;
      	      FORT_SETFRONTFACEVT(CHF_FRA1(facevt,dir),
      				  CHF_CONST_FRA1(thck,0),
      				  CHF_CONST_FRA1(usrf,0),
      				  CHF_CONST_FIA1(mask,0),
      				  CHF_CONST_INT(dir),
      				  CHF_CONST_REAL(factor),
      				  CHF_BOX(m_amrGrids[lev][dit]));
      	    }
      	}


      EdgeToCell(*m_viscousTensorFace[lev],*m_viscousTensorCell[lev]);
      if (lev > 0)
	{
	  PiecewiseLinearFillPatch ghostFiller
	    (m_amrGrids[lev],
	     m_amrGrids[lev-1],
	     m_viscousTensorCell[lev-1]->nComp(),
	     m_amrDomains[lev-1],
	     m_refinement_ratios[lev-1],
	     m_viscousTensorCell[lev-1]->ghostVect()[0]);
	  
	  ghostFiller.fillInterp(*m_viscousTensorCell[lev], 
				 *m_viscousTensorCell[lev-1], 
				 *m_viscousTensorCell[lev-1],1.0,0,0,
				 m_viscousTensorCell[lev-1]->nComp());

	}
      m_viscousTensorCell[lev]->exchange();

      EdgeToCell(*viscosityCoef[lev],*m_viscosityCoefCell[lev]);

      for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
      	{
	  const BaseFab<int>& mask = m_vect_coordSys[lev]->getFloatingMask()[dit];
	  FArrayBox& cellvt = (*m_viscousTensorCell[lev])[dit];
	  const Real z = 0.0;
	  for (int comp = 0; comp < SpaceDim * SpaceDim; comp++)
	    {
	      FORT_SETICEFREEVAL(CHF_FRA1(cellvt,comp), 
				 CHF_CONST_FIA1(mask,0),
				 CHF_CONST_REAL(z),
				 CHF_BOX(m_amrGrids[lev][dit]));
	    }
	}

      dragCoef[lev]->copyTo(Interval(0,0),*m_dragCoef[lev],Interval(0,0));

      if (faceA[lev] != NULL)
	{
	  delete faceA[lev]; faceA[lev] = NULL;
	}
      if (C0[lev] != NULL)
	{
	  delete C0[lev]; C0[lev] = NULL;
	}
    }

  m_viscousTensor_valid = true;

}

//access the grounding line proximity
const LevelData<FArrayBox>* AmrIce::groundingLineProximity(int a_level) const
{

  updateGroundingLineProximity();
  
  if (!(m_groundingLineProximity.size() > a_level))
    {
      std::string msg("AmrIce::groundingLineProximity !(m_groundingLineProximity.size() > a_level)");
      pout() << msg << endl;
      CH_assert((m_groundingLineProximity.size() > a_level));
      MayDay::Error(msg.c_str());
    }


  LevelData<FArrayBox>* ptr = m_groundingLineProximity[a_level];
  if (ptr == NULL)
    {
      std::string msg("AmrIce::groundingLineProximity m_groundingLineProximity[a_level] == NULL)");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }

  return ptr;
}


void AmrIce::applyCalvingCriterion(CalvingModel::Stage a_stage)
{

  CH_TIME("AmrIce::applyCalvingCriterion");
  CH_assert(m_calvingModelPtr != NULL);
  
  // observers (e.g AMRMelange) may care about the calved ice
  if (a_stage != CalvingModel::Initialization)
    notifyObservers(Observer::PreCalving);

  //allow calving model to modify geometry 
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      LevelData<FArrayBox>& thck = m_vect_coordSys[lev]->getH();
      LevelData<FArrayBox>& frac = *m_iceFrac[lev];
      LevelData<FArrayBox>& calvedIce = *m_calvedIceThickness[lev];
      LevelData<FArrayBox>& addedIce = *m_addedIceThickness[lev];
      LevelData<FArrayBox>& removedIce = *m_removedIceThickness[lev];
      m_calvingModelPtr->applyCriterion(thck, calvedIce, addedIce, removedIce, frac, *this, lev, a_stage);	  
 
   }
  
  // observers (e.g AMRMelange) may care about the calved ice
  if (a_stage != CalvingModel::Initialization)
  notifyObservers(Observer::PostCalving);
  
  // usually a good time to eliminate remote ice
  if (m_eliminate_remote_ice) eliminateRemoteIce();
  
}


///Identify regions of floating ice that are remote
///from grounded ice and eliminate them.
void AmrIce::eliminateRemoteIce()
{
  
  //any thickness change in eliminateRemoteIce is assumed to be calving: observers may care
  notifyObservers(Observer::PreCalving);
  
  IceUtility::eliminateRemoteIce(m_vect_coordSys, m_velocity, 
				 m_calvedIceThickness, m_addedIceThickness,
				 m_removedIceThickness,
				 m_amrGrids, m_amrDomains, 
				 m_refinement_ratios, m_amrDx[0], 
				 m_finest_level, m_eliminate_remote_ice_max_iter,
				 m_eliminate_remote_ice_tol,s_verbosity);

  //any thickness change in eliminateRemoteIce is assumed to be calving: observers may care
  notifyObservers(Observer::PostCalving);
  
}




void 
AmrIce::implicitThicknessCorrection(Real a_dt,
				    const Vector<LevelData<FArrayBox>* >& a_sts,
				    const Vector<LevelData<FArrayBox>* >& a_bts,
				    const Vector<LevelData<FArrayBox>* >& a_vts
				    )
{

  CH_TIME("AmrIce::implicitThicknessCorrection");
  if (s_verbosity > 3)
    {
      pout() << "AmrIce::implicitThicknessCorrection" << std::endl;
    }

  if  (m_temporalAccuracy == 1)
    {  
      //implicit Euler : solve (I - dt P) H = H_pred + dt * S
      
      //slc: at the moment, I'm setting eveything up every time-step,
      //pretending that diffusion is constant in time, and using the multi-grid
      //solver only. All these things are to be improved 

      //Natural boundary conditions - OK for now, but ought to get 
      //moved into subclasses of IceThicknessIBC
      BCHolder bc(ConstDiriNeumBC(IntVect::Zero, RealVect::Zero,
      				  IntVect::Zero, RealVect::Zero));

      Vector<RefCountedPtr<LevelData<FArrayBox> > > I(finestTimestepLevel() + 1);
      Vector<RefCountedPtr<LevelData<FluxBox> > > D(finestTimestepLevel() + 1);
      Vector<LevelData<FArrayBox>* > H(finestTimestepLevel() + 1);
      Vector<LevelData<FArrayBox>* > rhs(finestTimestepLevel()+ 1);
      Vector<DisjointBoxLayout> grids(finestTimestepLevel() + 1);

      for (int lev=0; lev <= finestTimestepLevel(); ++lev)
	{
	  LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
	  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

	  I[lev] = RefCountedPtr<LevelData<FArrayBox> >
	    (new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit));

	  H[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);

	  D[lev] = RefCountedPtr<LevelData<FluxBox> >(m_diffusivity[lev]);
	  D[lev].neverDelete();

	  rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);

	  grids[lev] = levelGrids;

	  const LevelData<FArrayBox>& levelSTS = *a_sts[lev];
	  const LevelData<FArrayBox>& levelBTS = *a_bts[lev];
	  const LevelData<FArrayBox>& levelVTS = *a_vts[lev];

	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      
	      (*I[lev])[dit].setVal(one);
	      (*H[lev])[dit].copy(levelCoords.getH()[dit] , 0 , 0, 1);
	      (*rhs[lev])[dit].copy( (*H[lev])[dit] , 0 , 0, 1);
	      //(*rhs[lev])[dit].plus(levelSTS[dit],a_dt);
	      //(*rhs[lev])[dit].plus(levelBTS[dit],a_dt);
	      //(*rhs[lev])[dit].plus(levelVTS[dit],a_dt);
	      (*D[lev])[dit][0].plus(m_additionalDiffusivity);
	      (*D[lev])[dit][1].plus(m_additionalDiffusivity);
	      
	    }
	  rhs[lev]->exchange();
	  H[lev]->exchange();
	  m_diffusivity[lev]->exchange();
	  I[lev]->exchange();
	}

      VCAMRPoissonOp2Factory poissonOpFactory;//= new VCAMRPoissonOp2Factory;
      poissonOpFactory.define(m_amrDomains[0], grids , m_refinement_ratios,
			      m_amrDx[0], bc, 1.0, I,  a_dt, D);
    
      //Plain MG
      BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
      AMRMultiGrid<LevelData<FArrayBox> > mgSolver;
      mgSolver.define(m_amrDomains[0], poissonOpFactory , &bottomSolver, finestTimestepLevel()+1);
      //parse these
      mgSolver.m_eps = 1.0e-10;
      mgSolver.m_normThresh = 1.0e-10;
    
      int numMGSmooth = 4;
      mgSolver.m_pre = numMGSmooth;
      mgSolver.m_post = numMGSmooth;
      mgSolver.m_bottom = numMGSmooth;
      
      mgSolver.solve(H, rhs, finestTimestepLevel(), 0,  false);
   
      for (int lev=0; lev <= finestTimestepLevel()  ; ++lev)
	{
	  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
	  LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
          LevelData<FArrayBox>& levelCoord_H = levelCoords.getH();
	  
	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      CH_assert( (*H[lev])[dit].norm(0,0,1) < HUGE_THICKNESS);
	      levelCoord_H[dit].copy( (*H[lev])[dit], 0, 0, 1);

	      //put sensible values into the corners.
	      FArrayBox &thisH = levelCoord_H[dit];
	      Box sbox = thisH.box();
	      sbox.grow(-levelCoord_H.ghostVect()[0]);
	      FORT_EXTRAPCORNER2D(CHF_FRA(thisH),
      				  CHF_BOX(sbox));

	    }

	  if (rhs[lev] != NULL)
	    {
	      delete rhs[lev];
	      rhs[lev] = NULL;
	    }
	  if (H[lev] != NULL)
	    {
	      delete H[lev];
	      H[lev] = NULL;
	    }
	}
    }
  else 
    {    
      MayDay::Error("AmrIce::implicitThicknessCorrection, invalid temporal accuracy");
    }


  

}




void AmrIce::helmholtzSolve
(Vector<LevelData<FArrayBox>* >& a_phi,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 Real a_alpha, Real a_beta) const
{

  // AMRPoissonOp supports only one component of phi
  // if m_finest_level > 0 (its LevelFluxRegisters are
  // defined with only one component, so we will do
  // one component at a time. should try to avoid some
  // of the rhs copies...
  for (int icomp = 0; icomp < a_phi[0]->nComp(); ++icomp)
    {
      
      //make a copy of a_phi with one ghost cell
      Vector<LevelData<FArrayBox>* > phi(m_finest_level + 1, NULL);
      Vector<LevelData<FArrayBox>* > rhs(m_finest_level + 1, NULL);
      Vector<DisjointBoxLayout> grids(m_finest_level + 1);
      for (int lev=0; lev < m_finest_level + 1; ++lev)
	{
	  grids[lev] = m_amrGrids[lev];

	  const LevelData<FArrayBox>& levelPhi = *a_phi[lev]; 
	  phi[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 
					      1, IntVect::Unit);
	  levelPhi.copyTo(Interval(icomp,icomp),*phi[lev], Interval(0,0));
	  phi[lev]->exchange();
	  

	  const LevelData<FArrayBox>& levelRhs = *a_rhs[lev];
	  rhs[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 
					      1, IntVect::Zero);
	  levelRhs.copyTo(Interval(icomp,icomp),*rhs[lev], Interval(0,0));
	  rhs[lev]->exchange();
      
	}


      //Natural boundary conditions
      BCHolder bc(ConstDiriNeumBC(IntVect::Zero, RealVect::Zero,
				  IntVect::Zero, RealVect::Zero));
      
      
      AMRPoissonOpFactory opf;
      opf.define(m_amrDomains[0],  grids , m_refinement_ratios,
		 m_amrDx[0], bc, a_alpha, -a_beta );
      
      AMRMultiGrid<LevelData<FArrayBox> > mgSolver;
      BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
      mgSolver.define(m_amrDomains[0], opf, &bottomSolver, m_finest_level+1);
      mgSolver.m_eps = TINY_NORM;
      mgSolver.m_normThresh = TINY_NORM;
      mgSolver.m_iterMax = 1;
      int numMGSmooth = 1;
      mgSolver.m_pre = numMGSmooth;
      mgSolver.m_post = numMGSmooth;
      mgSolver.m_bottom = numMGSmooth;
      mgSolver.m_verbosity = s_verbosity - 1;
      
      mgSolver.solve(phi, rhs, m_finest_level, 0,  false);
      
      for (int lev=0; lev < m_finest_level + 1; ++lev)
	{
	  LevelData<FArrayBox>& levelPhi = *a_phi[lev];
	  phi[lev]->copyTo(Interval(0,0), levelPhi, Interval(icomp, icomp));
	  
	  if (phi[lev] != NULL){
	    delete phi[lev];
	    phi[lev] = NULL;
	  }
	}
    }
}


void AmrIce::helmholtzSolve
(Vector<LevelData<FArrayBox>* >& a_phi, Real a_alpha, Real a_beta) const
{
  
  Vector<LevelData<FArrayBox>* > rhs(m_finest_level + 1, NULL);
 
  for (int lev=0; lev < m_finest_level + 1; ++lev)
    {
      const LevelData<FArrayBox>& levelPhi = *a_phi[lev]; 
      rhs[lev] = new LevelData<FArrayBox>
	(m_amrGrids[lev], levelPhi.nComp(), IntVect::Zero);
      levelPhi.copyTo(*rhs[lev]);
    }
 
  helmholtzSolve(a_phi, rhs, a_alpha, a_beta);

  for (int lev=0; lev < m_finest_level + 1; ++lev)
    {
      if (rhs[lev] != NULL){
	delete rhs[lev];
	rhs[lev] = NULL;
      }
    }

}


#if BISICLES_Z == BISICLES_LAYERED

/// update the flow law coefficient A
void AmrIce::computeA(Vector<LevelData<FArrayBox>* >& a_A, 
		      Vector<LevelData<FArrayBox>* >& a_sA,
		      Vector<LevelData<FArrayBox>* >& a_bA,
		      const Vector<LevelData<FArrayBox>* >& a_internalEnergy, 
		      const Vector<LevelData<FArrayBox>* >& a_sInternalEnergy,
		      const Vector<LevelData<FArrayBox>* >& a_bInternalEnergy,
		      const Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys) const
		      
{
  if (s_verbosity > 0)
    {
      pout() <<  "AmrIce::computeA" <<  std::endl;
    }

  //for now, throw a_A etc away and recompute
  for (int lev = 0; lev < a_A.size(); ++lev)
    {
      if (a_A[lev] != NULL)
	{
	  delete a_A[lev]; a_A[lev] = NULL;
	}
      
      if (a_sA[lev] != NULL)
	{
	  delete a_sA[lev]; a_sA[lev] = NULL;
	}
      if (a_bA[lev] != NULL)
	{
	  delete a_bA[lev]; a_bA[lev] = NULL;
	}
      
    }
  a_A.resize(m_finest_level+1,NULL);
  a_sA.resize(m_finest_level+1,NULL);
  a_bA.resize(m_finest_level+1,NULL);
	
  for (int lev = 0; lev <= m_finest_level; ++lev)
    {
      const LevelSigmaCS& levelCoords = *a_coordSys[lev];
      
      const Vector<Real>& sigma = levelCoords.getSigma();
      a_A[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], m_nLayers, IntVect::Unit);
      IceUtility::computeA(*a_A[lev], sigma, levelCoords,  m_rateFactor, *a_internalEnergy[lev] );
      
      Vector<Real> sSigma(1,0.0);
      a_sA[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1, IntVect::Unit);
      IceUtility::computeA(*a_sA[lev], sSigma, levelCoords,  
			   m_rateFactor, *a_sInternalEnergy[lev]);
      Vector<Real> bSigma(1,1.0);
      a_bA[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],1, IntVect::Unit);
      IceUtility::computeA(*a_bA[lev], bSigma, levelCoords,  
			   m_rateFactor, *a_bInternalEnergy[lev]);
   
    }//end loop over AMR levels

  if (s_verbosity > 0)
    {
      Real Amin = computeMin(a_A,  m_refinement_ratios, Interval(0,a_A[0]->nComp()-1));
      Real Amax = computeMax(a_A,  m_refinement_ratios, Interval(0,a_A[0]->nComp()-1));
      pout() << Amin << " <= A(x,y,sigma) <= " << Amax << std::endl;

    }
}

#endif

/// diagnostic function -- integrates thickness over domain
Real
AmrIce::computeTotalIce() const
{
  Vector<LevelData<FArrayBox>* > thickness(m_finest_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      // need a const_cast to make things all line up right
      // (but still essentially const)
      thickness[lev] = const_cast<LevelData<FArrayBox>* >(&levelCoords.getH());
    }

  Interval thicknessInt(0,0);
  Real totalIce = computeSum(thickness, m_refinement_ratios,
                             m_amrDx[0], thicknessInt, 0);


  return totalIce;

}

void AmrIce::computeAreaFraction(LevelData<FArrayBox>& a_area,
				 int a_maskVal,
				 int a_level) const
{
  CH_assert(a_level <= m_finest_level)
    {
      const LevelData<BaseFab<int> >& levelMask = m_vect_coordSys[a_level]->getFloatingMask();
      // set a_area = 1.0 where we have grounded ice
      DataIterator dit=a_area.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  const BaseFab<int>& thisMask = levelMask[dit];
	  FArrayBox& thisIce = a_area[dit];
	  thisIce.setVal(0.0);
	  BoxIterator bit(thisIce.box());
	  for (bit.begin(); bit.ok(); ++bit)
	    {
	      IntVect iv = bit();
	      if (thisMask(iv,0) == a_maskVal)
		{
		  thisIce(iv,0) = 1.0;
		}
	    }
	}
    }
}

void AmrIce::setLayers(const Vector<Real>& a_sigma)
{
   
  if (m_sigmaSet)
    MayDay::Error("AmrIce::SetLayers should only be called once");
  m_faceSigma = a_sigma;
  m_nLayers = m_faceSigma.size()-1;
  CH_assert(m_nLayers > 0 && m_faceSigma[0] < TINY_NORM && abs(m_faceSigma[m_nLayers] - 1.0) < TINY_NORM);
  m_sigmaSet = true;

}


void AmrIce::incrementIceThickness
(Vector<LevelData<FArrayBox>*> a_delta_thk)
{

  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {
      LevelSigmaCS& coords = *m_vect_coordSys[lev];
      LevelData<FArrayBox>& thk = coords.getH();
      for (DataIterator dit( m_amrGrids[lev]); dit.ok(); ++dit)
	{
	  thk[dit] += ( (*a_delta_thk[lev])[dit] );
	}
    }

  // average down thickness and topography to coarser levels, fill in ghost cells
  int Hghost = 2;
  Vector<LevelData<FArrayBox>* > vectH(m_finest_level+1, NULL);
 
  for (int lev=0; lev<vectH.size(); lev++)
    {
      IntVect HghostVect = Hghost*IntVect::Unit;
      LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
      vectH[lev] = &levelCoords.getH();
    }
  
  //average from the finest level down
  for (int lev =  finestTimestepLevel() ; lev > 0 ; lev--)
    {
      CoarseAverage averager(m_amrGrids[lev],
                             1, m_refinement_ratios[lev-1]);
      averager.averageToCoarse(*vectH[lev-1],
                               *vectH[lev]);
    }

  // now pass back over and do PiecewiseLinearFillPatch
  for (int lev=1; lev<vectH.size(); lev++)
    {
      
      PiecewiseLinearFillPatch filler(m_amrGrids[lev],
                                      m_amrGrids[lev-1],
                                      1, 
                                      m_amrDomains[lev-1],
                                      m_refinement_ratios[lev-1],
                                      Hghost);
      
      Real interp_coef = 0;
      filler.fillInterp(*vectH[lev],
                        *vectH[lev-1],
                        *vectH[lev-1],
                        interp_coef,
                        0, 0, 1);
    }
 
  //re-fill ghost cells ouside the domain
  for (int lev=0; lev <= finestTimestepLevel()  ; ++lev)
    {
      RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
      m_thicknessIBCPtr->setGeometryBCs(*m_vect_coordSys[lev],
                                        m_amrDomains[lev],levelDx, m_time, m_dt);
    }
  
  //allow calving model to modify geometry and velocity
  //applyCalvingCriterion(CalvingModel::PostThicknessAdvection);
  
  //dont allow thickness to be negative
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
      LevelData<FArrayBox>& levelH = levelCoords.getH();
      DataIterator dit = levelGrids.dataIterator();          
      
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
        {
          Real lim = 0.0;
          FORT_MAXFAB1(CHF_FRA(levelH[dit]), 
                       CHF_CONST_REAL(lim), 
                       CHF_BOX(levelH[dit].box()));
        }
    }

  //average from the finest level down
  for (int lev =  finestTimestepLevel() ; lev > 0 ; lev--)
    {
      CoarseAverage averager(m_amrGrids[lev],
                             1, m_refinement_ratios[lev-1]);
      averager.averageToCoarse(*vectH[lev-1],
                               *vectH[lev]);      
    }
  
  for (int lev=1; lev<vectH.size(); lev++)
    {      
      PiecewiseLinearFillPatch filler(m_amrGrids[lev],
                                      m_amrGrids[lev-1],
                                      1, 
                                      m_amrDomains[lev-1],
                                      m_refinement_ratios[lev-1],
                                      Hghost);
      
      Real interp_coef = 0;
      filler.fillInterp(*vectH[lev],
                        *vectH[lev-1],
                        *vectH[lev-1],
                        interp_coef,
                        0, 0, 1);
    }
  
  // recompute thickness-derived data in SigmaCS
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      LevelSigmaCS& levelCoords = *(m_vect_coordSys[lev]);
      LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
      int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
      levelCoords.recomputeGeometry(crseCoords, refRatio);            
    }
      
}



#include "NamespaceFooter.H"
