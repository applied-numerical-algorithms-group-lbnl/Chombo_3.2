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
#include <functional>
#ifdef HAVE_UDUNITS2
#include <udunits2.h>
#endif

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
#include "FASIceSolverI.H"
#endif
#include "KnownVelocitySolver.H"
#include "VCAMRPoissonOp2.H"
#include "AMRPoissonOpF_F.H"
#include "CH_HDF5.H"
#include "IceUtility.H"
#include "LevelMappedDerivatives.H"
#include "FillFromReference.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif

#include "NamespaceHeader.H"


void AmrIce::setOutputOptions(ParmParse& a_pp)
{

  m_plot_prefix = "plot.";
  m_plot_style_cf = false;
  m_plot_style_amr = true;
  m_output_file_numbering = time_step;
  m_plot_interval = 10000000;
  m_plot_time_interval = 1.0e+12;
  m_reduced_plot = false;
  m_write_presolve_plotfiles = false;
  m_write_solver_rhs = false;
  m_write_dHDt = true;
  m_write_fluxVel = true;
  m_write_viscousTensor = false;
  m_write_baseVel = true;
  m_write_internal_energy = false;
  m_write_map_file = false;
  m_write_thickness_sources = false;
  m_write_ismip6 = false;
  m_write_layer_velocities = false;
  m_write_mask = false;
  m_check_prefix = "chk.";
  m_check_interval = -1;
  m_check_overwrite = true;
  m_check_exit = false;

  a_pp.query("check_interval", m_check_interval);
  a_pp.query("check_prefix", m_check_prefix);
  a_pp.query("check_overwrite", m_check_overwrite);
  a_pp.query("check_exit", m_check_exit);
  
  a_pp.query("plot_interval", m_plot_interval);
  a_pp.query("plot_time_interval", m_plot_time_interval);
  CH_assert(m_plot_time_interval > 0);
  a_pp.query("write_preSolve_plotfiles", m_write_presolve_plotfiles);
  a_pp.query("plot_prefix", m_plot_prefix);
  a_pp.query("plot_style_cf", m_plot_style_cf);
  a_pp.query("plot_style_amr", m_plot_style_amr);
  a_pp.query("write_flux_velocities", m_write_fluxVel);
  a_pp.query("write_viscous_tensor", m_write_viscousTensor);
  a_pp.query("write_base_velocities", m_write_baseVel);
  a_pp.query("write_internal_energy", m_write_internal_energy);
  a_pp.query("write_thickness_sources", m_write_thickness_sources);
  a_pp.query("write_ismip6", m_write_ismip6);
  a_pp.query("write_layer_velocities", m_write_layer_velocities);
  a_pp.query("write_map_file", m_write_map_file);
  a_pp.query("write_mask",m_write_mask);
  a_pp.query("reduced_plot", m_reduced_plot);

  
  std::string pn = "time_step";
  a_pp.query("output_file_numbering",pn);
  if (pn == "time_step")
    {
      m_output_file_numbering = time_step;
    }
  else if  (pn == "time_seconds")
    {
      m_output_file_numbering = time_seconds;
    }
  else if  (pn == "time_years")
    {
      m_output_file_numbering= time_years;
    }
  else if  (pn == "time_yyyymmdd_360")
    {
      m_output_file_numbering = time_yyyymmdd_360 ;
    }
  else
    {
      MayDay::Error("Unknown output_file_numbering");
    }
}

///
long AmrIce::outputNumbering()
{

  long rc = 0;
  Real second = time() * m_seconds_per_unit_time;
  Real year = second/SECONDS_PER_TROPICAL_YEAR;
  int yyyy = int(year);
  Real month = (year - Real(yyyy)) * 12.0;
  int mm = int(month);
  Real day = (year - Real(yyyy) - Real(mm)/12.0) * 30.0;
  int dd = int(day);

  /// nobody likes new year on 00/00/yyyy
  mm++;dd++;
  
  switch (m_output_file_numbering)
    {
    case(time_seconds):
      rc = long(second);
      break;
    case(time_years):
      rc = long(year);
      break;
    case(time_yyyymmdd_360):
      rc = long(dd + 100*mm + 10000*yyyy);
      break;
    default:
      rc = long(m_cur_step);
    }

  return rc;
}

void 
AmrIce::writePlotFile()
{
  if (m_plot_style_cf && m_cur_step > 0 ) flushCFData();
  if (m_plot_style_amr) writeAMRPlotFile();
}

/// write hdf5 plotfile to the standard location
void 
AmrIce::writeAMRPlotFile() 
{
  CH_TIME("AmrIce::writeAMRPlotFile");

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::writeAMRPlotFile" << endl;
    }
  
  // plot comps: thickness + horizontal velocity + z_bottom + z_surface + z_base 
  int numPlotComps = 4 + SpaceDim;

  if (m_write_ismip6)
    {
      CH_assert(!m_reduced_plot);
      m_write_baseVel=true;
      m_write_mask=true;
      // include surface and basal temperatures.
      numPlotComps += 2;
      // include drag coefficient.
      if (!m_write_viscousTensor && !m_write_solver_rhs)
	{
	  numPlotComps += 1;
	}
    }

  if (m_reduced_plot)
    {
      // plot comps: thickness + horizontal velocity + zb + zs
      numPlotComps = 3 + SpaceDim;
    }

  // may need a zvel for visit to do "3d" streamlines correctly
  bool writeZvel = !m_reduced_plot;
  if (writeZvel) numPlotComps+=1;

  if (m_write_fluxVel)
    {
      numPlotComps += SpaceDim;
      if (writeZvel) numPlotComps+=1;
    }
  
 
  if (m_write_baseVel)
    {
      numPlotComps += SpaceDim;
      if (writeZvel) numPlotComps+=1;
    }
  // write both integer and real-valued masks
  if (m_write_mask) numPlotComps += 2;
  if (m_write_dHDt) numPlotComps += 1;
  if (m_write_solver_rhs) 
    {
      numPlotComps += SpaceDim ;
      // include basal_friction and C0 iff !m_reduced_plot 
      if (!m_reduced_plot)
        {
          numPlotComps += 2;
        }
    }
  
  

  if (m_write_internal_energy) 
    numPlotComps += m_internalEnergy[0]->nComp() + 1; // layer energies plus till water depth
#if BISICLES_Z == BISICLES_LAYERED
  if (m_write_internal_energy)
    numPlotComps += 4;// surface and basal internalEnergys and heat fluxes

  //layer velocities
  if (m_write_layer_velocities)
    {
      numPlotComps += SpaceDim * (m_nLayers+1);
      if (writeZvel) numPlotComps += (m_nLayers+1);
    }

  if (m_write_viscousTensor)
    {
      // effective drag and viscosity coefficients
      numPlotComps += 2;
      if (!m_reduced_plot)
	{
	  numPlotComps += SpaceDim * SpaceDim; // viscous tensor components                        
	}
    }
  
  if (m_write_thickness_sources)
    {
      numPlotComps += 2;  // surface and basal sources
      if (!m_reduced_plot)
	{
	  numPlotComps += 4; // divThicknessFlux, calving flux and accumulated calving
	  numPlotComps += 2; // calving rate and water depth
	}
    }


#endif
  // generate data names

  string thicknessName("thickness");
  string xVelName("xVel");
  string yVelName("yVel");
  string zVelName("zVel");
  string zsName("Z_surface");
  string zbName("Z_base");
  string zbottomName("Z_bottom");
  string dthicknessName("dThickness/dt");
  string betaName("basal_friction");
  string solverRhsxName("xRhs");
  string solverRhsyName("yRhs");
  string C0Name("C0");
  string maskName("mask");
  string fracName("iceFrac");
  string xfVelName("xfVel");
  string yfVelName("yfVel");
  string zfVelName("zfVel");
  string xbVelName("xbVel");
  string ybVelName("ybVel");
  string zbVelName("zbVel");

  string internalEnergyName("internalEnergy");
  string heatFluxName("heatflux");
  string tillWaterDepthName("tillWaterDepth");

  string sTempName("sTemp");
  string bTempName("bTemp");

#if BISICLES_Z == BISICLES_LAYERED
   string xlayerVelName("xlayerVel");
  string ylayerVelName("ylayerVel");
  string zlayerVelName("zlayerVel");
#endif

  string xxVTname("xxViscousTensor");
  string xyVTname("xyViscousTensor");
  string xzVTname("xzViscousTensor");
  string yxVTname("yxViscousTensor");
  string yyVTname("yyViscousTensor");
  string yzVTname("yzViscousTensor");
  string zxVTname("zxViscousTensor");
  string zyVTname("zyViscousTensor");
  string zzVTname("zzViscousTensor");
  string viscosityCoefName("viscosityCoef");
  //string yViscosityCoefName("yViscosityCoef");
  //string zViscosityCoefName("zViscosityCoef");
  string dragCoefName("dragCoef");

  string activeBasalThicknessSourceName("activeBasalThicknessSource");
  string activeSurfaceThicknessSourceName("activeSurfaceThicknessSource");
  string basalThicknessSourceName("basalThicknessSource");
  string surfaceThicknessSourceName("surfaceThicknessSource");
  string divergenceThicknessFluxName("divergenceThicknessFlux");
  string calvedIceThicknessName("calvingFlux");
  string calvedThicknessSourceName("calvedThicknessSource");
  string calvingRateName("calvingRate");
  string waterDepthName("waterDepth");

  Vector<string> vectName(numPlotComps);
  //int dThicknessComp;

  vectName[0] = thicknessName;
  vectName[1] = xVelName;
  if (SpaceDim > 1)
    vectName[2] = yVelName;
  int comp = SpaceDim+1;
  if (writeZvel) 
    {
      vectName[comp] = zVelName;
      comp++;
    }

  vectName[comp] = zsName;
  comp++;

  if (!m_reduced_plot)
    {
      vectName[comp] = zbottomName;
      comp++;
    }

  vectName[comp] = zbName;
  comp++;

  if (m_write_solver_rhs)
    {
      if (!m_reduced_plot)
	{
	  vectName[comp] = betaName;
	  comp++;
	  
	  vectName[comp] = C0Name;
	  comp++;
	}

      if (SpaceDim == 1)
        {
          vectName[comp] = solverRhsxName;
          comp++;
        }
      else if (SpaceDim == 2)
        {
          vectName[comp] = solverRhsxName;
          comp++;
          vectName[comp] = solverRhsyName;
          comp++;
        }
      else
        {
          MayDay::Error("writeSolverRHS undefined for this dimensionality");
        }
    }

  if (m_write_dHDt)
    {
      vectName[comp] = dthicknessName;      
      comp++;
    } 

  if (m_write_mask)
    {
      vectName[comp] = maskName;      
      comp++;
      vectName[comp] = fracName;
      comp++;
    } 


  if (m_write_fluxVel)
    {
      vectName[comp] = xfVelName;
      comp++;
      if (SpaceDim > 1)
        {
          vectName[comp] = yfVelName;
          comp++;
        }
      
      if (writeZvel) 
        {
          vectName[comp] = zfVelName;
          comp++;
        }
    }

 

  if (m_write_baseVel)
    {
      vectName[comp] = xbVelName;
      comp++;
      if (SpaceDim > 1)
        {
          vectName[comp] = ybVelName;
          comp++;
        }
      
      if (writeZvel) 
        {
          vectName[comp] = zbVelName;
          comp++;
        }
    }

  if (m_write_internal_energy)
    {
#if BISICLES_Z == BISICLES_LAYERED
      vectName[comp] = internalEnergyName + string("Surface");
      comp++;
#endif    
      for (int l = 0; l < m_internalEnergy[0]->nComp(); ++l)
	{
	  char idx[12]; sprintf(idx, "%04d", l);
	  vectName[comp] = internalEnergyName + string(idx);
	  comp++;
	}
      vectName[comp] = tillWaterDepthName; 
      comp++;
#if BISICLES_Z == BISICLES_LAYERED
      vectName[comp] = internalEnergyName + string("Base");
      comp++;
      vectName[comp] = heatFluxName + string("Surface");
      comp++;
      vectName[comp] = heatFluxName + string("Base");
      comp++;
#endif 
    }

#if BISICLES_Z == BISICLES_LAYERED
  if (m_write_layer_velocities){
    for (int l = 0; l < m_nLayers + 1; ++l)
      {
	char idx[12]; sprintf(idx, "%04d", l);
	vectName[comp] = xlayerVelName + string(idx);
	comp++;
	vectName[comp] = ylayerVelName + string(idx);
	comp++;
	if (writeZvel) 
	  {
	    vectName[comp] = zlayerVelName + string(idx);
	    comp++;
	  }
      }
  }
#endif

  if (m_write_viscousTensor)
    {
      vectName[comp] = dragCoefName; comp++;
      vectName[comp] = viscosityCoefName; comp++;
      if (!m_reduced_plot)
	{
	  vectName[comp] = xxVTname;comp++;
	  if (SpaceDim > 1)
	    {
	      vectName[comp] = yxVTname;comp++;
	      if (SpaceDim > 2)
		{
		  vectName[comp] = zzVTname;comp++;
		}
	      vectName[comp] = xyVTname;comp++;
	      vectName[comp] = yyVTname;comp++;
	      
	      if (SpaceDim > 2)
		{
		  vectName[comp] = zyVTname;comp++;
		  vectName[comp] = xzVTname;comp++;
		  vectName[comp] = yzVTname;comp++;
		  vectName[comp] = zzVTname;comp++;
		}
	    }
	}
    }

  if (m_write_thickness_sources)
    {
      vectName[comp] = activeBasalThicknessSourceName; comp++;
      vectName[comp] = activeSurfaceThicknessSourceName; comp++;
      
      if (!m_reduced_plot)
	{
	  vectName[comp] = divergenceThicknessFluxName; comp++;	
	  vectName[comp] = basalThicknessSourceName; comp++;
	  vectName[comp] = surfaceThicknessSourceName; comp++;
	  vectName[comp] = calvedIceThicknessName; comp++;
	  vectName[comp] = calvingRateName; comp++;
	  vectName[comp] = waterDepthName; comp++;
	}
    }

  if (m_write_ismip6)
    {      
      vectName[comp] = sTempName; comp++;
      vectName[comp] = bTempName; comp++;
      if (!m_write_viscousTensor && !m_write_solver_rhs)
	{
	  vectName[comp] = dragCoefName; comp++;
	}
    }

  // allow observers to add variables to the plot file
  for (int i = 0; i < m_observers.size(); i++)
    m_observers[i]->addPlotVars(vectName);
  numPlotComps = vectName.size();


  Box domain = m_amrDomains[0].domainBox();
  int numLevels = m_finest_level +1;
  // compute plot data
  Vector<LevelData<FArrayBox>* > plotData(m_velocity.size(), NULL);

  // temp storage for C0
  Vector<LevelData<FArrayBox>* > vectC0(m_velocity.size(), NULL);

  // ghost vect makes things simpler
  IntVect ghostVect(IntVect::Unit);
  
  for (int lev=0; lev<numLevels; lev++)
    {
      // first allocate storage
      plotData[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                               numPlotComps,
                                               ghostVect);

      vectC0[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                             1,
                                             m_velBasalC[0]->ghostVect());
      DataIterator dit = m_amrGrids[lev].dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          (*vectC0[lev])[dit].setVal(0.0);
        }
    }

  if (m_write_solver_rhs)
    {
      setBasalFriction(m_velBasalC, vectC0);
      defineVelRHS(m_velRHS);
      
    }


  for (int lev=0; lev<numLevels; lev++)
    {
      // now copy new-time solution into plotData
      Interval thicknessComps(0,0);
      Interval velocityComps(1,SpaceDim);

      LevelData<FArrayBox>& plotDataLev = *plotData[lev];

      const LevelSigmaCS& levelCS = (*m_vect_coordSys[lev]);
      const LevelData<FArrayBox>& levelH = levelCS.getH();
      const LevelData<FArrayBox>& levelZbase = levelCS.getTopography();
      LevelData<FArrayBox> levelZsurf(m_amrGrids[lev], 1, ghostVect);
      levelCS.getSurfaceHeight(levelZsurf);
      LevelData<FArrayBox> levelSTS (m_amrGrids[lev], 1, ghostVect);
      LevelData<FArrayBox> levelBTS (m_amrGrids[lev], 1, ghostVect);
      LevelData<FArrayBox> levelCalvingRate(m_amrGrids[lev], 1, ghostVect);
      LevelData<FArrayBox> levelWaterDepth(m_amrGrids[lev], 1, ghostVect);
      
      if (m_write_thickness_sources)
	{
	  m_surfaceFluxPtr->surfaceThicknessFlux(levelSTS, *this, lev, m_dt);
	  m_basalFluxPtr->surfaceThicknessFlux(levelBTS, *this, lev, m_dt); 
	  (*m_calvingModelPtr).getCalvingRate(levelCalvingRate, *this, lev);
	  (*m_calvingModelPtr).getWaterDepth(levelWaterDepth, *this, lev);
	}


      // set these up here because calls to dragCoefficient may lead
      // to calls to things like CoarseAverage which contain dataIterator loops
      LevelData<FArrayBox>* dragCoefficientPtr = NULL;
      LevelData<FArrayBox>* viscosityCoefficientPtr = NULL;
      LevelData<FArrayBox>* viscousTensorPtr = NULL;
      
      if (m_write_viscousTensor)      
        {
          dragCoefficientPtr = const_cast<LevelData<FArrayBox>* >(dragCoefficient(lev));
          viscosityCoefficientPtr = const_cast<LevelData<FArrayBox>*>(viscosityCoefficient(lev));
          viscousTensorPtr = const_cast<LevelData<FArrayBox>*>(viscousTensor(lev));
        }
      
      DataIterator dit = m_amrGrids[lev].dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = m_amrGrids[lev][dit];
          FArrayBox& thisPlotData = plotDataLev[dit];
          comp = 0;
          const FArrayBox& thisH = levelH[dit];
          
          thisPlotData.copy(thisH, 0, comp, 1);

          comp++;
          const FArrayBox& thisVel = (*m_velocity[lev])[dit];
          thisPlotData.copy(thisVel, 0, comp, SpaceDim);
          
          comp += SpaceDim;
	 
          if (writeZvel) 
            {
              // use zVel = zero for the moment
              Real zVel = 0.0;
              thisPlotData.setVal(zVel, comp);
              ++comp;
            }

          const FArrayBox& zBase = levelZbase[dit];
          
          // account for background slope of base 
          FArrayBox backgroundBase(thisPlotData.box(), 1);
          BoxIterator bit(thisPlotData.box());
	  const RealVect& basalSlope = levelCS.getBackgroundSlope();
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc(iv);
              loc += 0.5*RealVect::Unit;
              loc *= m_amrDx[lev];

              backgroundBase(iv,0) = D_TERM(loc[0]*basalSlope[0],
                                            +loc[1]*basalSlope[1],
                                            +loc[2]*basalSlope[2]);
            }
          // zsurface
          FArrayBox& zSurf = levelZsurf[dit];
          thisPlotData.copy(zSurf, 0, comp, 1);
          thisPlotData.plus(backgroundBase, 0, comp, 1);
          ++comp;

	  if (!m_reduced_plot)
	    {
	      // zbottom (bottom of ice
	      thisPlotData.copy(zSurf, 0, comp, 1);
	      thisPlotData.minus(thisH, 0, comp, 1);
	      thisPlotData.plus(backgroundBase, 0, comp, 1);
	      ++comp;
	    }

          // zbase 
          thisPlotData.copy(zBase, 0, comp, 1);
          thisPlotData.plus(backgroundBase, 0, comp, 1);
          ++comp;

          if (m_write_solver_rhs)
            {
	      if (!m_reduced_plot)
		{
		  thisPlotData.copy((*m_velBasalC[lev])[dit],0,comp,1);
		  comp++;
		  thisPlotData.copy((*vectC0[lev])[dit],0,comp,1);
		  comp++;
		}
	      thisPlotData.copy((*m_velRHS[lev])[dit],0,comp,SpaceDim);
              comp += SpaceDim;
            }

          // now copy for dthickness/dt 
          if (m_write_dHDt)
            {
              const FArrayBox& thisOldH = (*m_old_thickness[lev])[dit];
              thisPlotData.copy(thisH, 0, comp, 1);
              thisPlotData.minus(thisOldH, 0, comp, 1);
              if (m_dt > 0)
                {
                  thisPlotData.divide(m_dt, comp, 1);

                }              
              ++comp;

	    } // end if we are computing dHDt
      
	  if (m_write_mask)
            {
	      const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	      FArrayBox tmp(mask.box(),1);
	      for (BoxIterator bit(mask.box());bit.ok();++bit)
		{
		  tmp(bit()) = Real( mask(bit()) ) ;
		}
	      thisPlotData.copy(tmp,0,comp,1);
	      comp++;
              // now copy real-valued ice fraction
              const FArrayBox& iceFracFab = (*m_iceFrac[lev])[dit];
              thisPlotData.copy(iceFracFab,0,comp,1);
	      comp++;

	    }

	  // const FArrayBox& thisSurfaceVel = (*m_velocity[lev])[dit];
          // thisPlotData.copy(thisSurfaceVel, 0, comp, 2);
          
          // comp += 2;

	  

          if (m_write_fluxVel)
            {
              for (int dir = 0; dir < SpaceDim; ++dir)
                {
                  
                  const FArrayBox& thisVel = (*m_faceVelTotal[lev])[dit][dir];
                  for (BoxIterator bit(gridBox); bit.ok(); ++bit)
                    {
                      const IntVect& iv = bit();
                      const IntVect ivp = iv + BASISV(dir);
                      thisPlotData(iv,comp) = half*(thisVel(iv) + thisVel(ivp));
                    }
                  comp++;
                }            

              if (writeZvel) 
                {
                  // use zVel = zero for the moment
                  Real zVel = 0.0;
                  thisPlotData.setVal(zVel, comp);
                  ++comp;
                }
            }



          if (m_write_baseVel)
            {
              const FArrayBox& thisBaseVel = (*m_velocity[lev])[dit];
              thisPlotData.copy(thisBaseVel, 0, comp, SpaceDim);
              
              comp += SpaceDim;
              
              if (writeZvel) 
                {
                  // use zVel = zero for the moment
                  Real zVel = 0.0;
                  thisPlotData.setVal(zVel, comp);
                  ++comp;
                }	  
            }
	  if (m_write_internal_energy)
	    {
#if BISICLES_Z == BISICLES_LAYERED
	      {
		const FArrayBox& thisTemp = (*m_sInternalEnergy[lev])[dit];
		thisPlotData.copy(thisTemp, 0, comp, thisTemp.nComp());
		comp++;
	      }
#endif
	      {
		const FArrayBox& thisTemp = (*m_internalEnergy[lev])[dit];
		thisPlotData.copy(thisTemp, 0, comp, thisTemp.nComp());
  		comp += thisTemp.nComp();
	      }

	      {
		const FArrayBox& thisTemp = (*m_tillWaterDepth[lev])[dit];
		thisPlotData.copy(thisTemp, 0, comp, 1);
		comp ++;
	       }
	      
#if BISICLES_Z == BISICLES_LAYERED
	      {
		const FArrayBox& thisTemp = (*m_bInternalEnergy[lev])[dit];
		thisPlotData.copy(thisTemp, 0, comp, thisTemp.nComp());
		comp++;
		thisPlotData.copy((*m_sHeatFlux[lev])[dit], 0, comp, 1);
		comp++;
		thisPlotData.copy((*m_bHeatFlux[lev])[dit], 0, comp, 1);
		comp++;
	      }
#endif
	    }
#if BISICLES_Z == BISICLES_LAYERED
	  if (m_write_layer_velocities)
	    {
	      const FArrayBox& thisVel = (*m_layerSFaceXYVel[lev])[dit];
	      const FArrayBox& thisSVel = (*m_layerSFaceSVel[lev])[dit];
	      for (int j = 0; j < m_nLayers + 1; ++j)
		{
		  thisPlotData.copy(thisVel, j*SpaceDim, comp, SpaceDim);
		  
		  comp+= SpaceDim;
		  // end loop over components
		  if (writeZvel) 
		    {
		      thisPlotData.copy(thisSVel, j, comp++, 1);
		      //++comp;
		    } 
	      
		} // end loop over layers
	    }
#endif
	  if (m_write_viscousTensor)
	    {
	      thisPlotData.copy( (*dragCoefficientPtr)[dit],0,comp);
	      comp++;
	      thisPlotData.copy( (*viscosityCoefficientPtr)[dit],0,comp);
	      comp++;
	      if (!m_reduced_plot)
		{
		  thisPlotData.copy( (*viscousTensorPtr)[dit],0,comp, SpaceDim*SpaceDim);
		  comp += SpaceDim * SpaceDim;
		}
	    }
	 
	  if (m_write_thickness_sources)
	    {
	      thisPlotData.copy((*m_basalThicknessSource[lev])[dit], 0, comp, 1);
	      if (m_frac_sources)
		{
		  thisPlotData.mult( (*m_iceFrac[lev])[dit],0,comp,1);
		}
	      comp++;
		  
	 
	      thisPlotData.copy((*m_surfaceThicknessSource[lev])[dit], 0, comp, 1);
	      if (m_frac_sources)
		{
		  // scale by ice fraction
		  thisPlotData.mult( (*m_iceFrac[lev])[dit],0,comp,1);
		}
	      comp++;

	      if (!m_reduced_plot)
		{
		  thisPlotData.copy((*m_divThicknessFlux[lev])[dit], 0, comp, 1);
		  comp++;

		  thisPlotData.copy(levelBTS[dit], 0, comp, 1);
		  if (m_frac_sources)
		    {
		      thisPlotData.mult( (*m_iceFrac[lev])[dit],0,comp,1);
		    }
		  comp++;

		  thisPlotData.copy(levelSTS[dit], 0, comp, 1);
		  if (m_frac_sources)
		    {
		      // scale by ice fraction
		      thisPlotData.mult( (*m_iceFrac[lev])[dit],0,comp,1);
		    }
		  comp++;
	      
		  thisPlotData.copy((*m_calvedIceThickness[lev])[dit], 0, comp, 1);
		  if (m_dt > 0)
		    {
		      thisPlotData.divide(m_dt, comp, 1);
		    }              
		  comp++;

		  thisPlotData.copy(levelCalvingRate[dit], 0, comp, 1);
		  comp++;

		  thisPlotData.copy(levelWaterDepth[dit], 0, comp, 1);
		  comp++;
		  
		}
	    }

	  if (m_write_ismip6)
	    {      
	      thisPlotData.copy( (*m_sTemperature[lev])[dit],0,comp);
	      comp++;
	      thisPlotData.copy( (*m_bTemperature[lev])[dit],0,comp);
	      comp++;

	      if (!m_write_viscousTensor && !m_write_solver_rhs)
		{
		  thisPlotData.copy( (*dragCoefficient(lev))[dit],0,comp);
		  comp++;
		}
	    }
		
	} // end loop over boxes on this level

      
      //allow observers to write data. 
      for (int i = 0; i < m_observers.size(); i++)
	{
	  Vector<std::string> vars;
	  m_observers[i]->addPlotVars(vars);
	  if (vars.size() > 0)
	    {
	      Interval interval(comp, comp + vars.size() - 1);
	      LevelData<FArrayBox> obsPlotData;
	      aliasLevelData( obsPlotData, plotData[lev], interval);
	      m_observers[i]->writePlotData(obsPlotData, lev);
	    }
	}


      // this is just so that visit surface plots look right
      // fill coarse-fine ghost-cell values with interpolated data
      if (lev > 0)
        {
          PiecewiseLinearFillPatch interpolator(m_amrGrids[lev],
                                                m_amrGrids[lev-1],
                                                numPlotComps,
                                                m_amrDomains[lev-1],
                                                m_refinement_ratios[lev-1],
                                                ghostVect[0]);
          
          // no interpolation in time
          Real time_interp_coeff = 0.0;
          interpolator.fillInterp(*plotData[lev],
                                  *plotData[lev-1],
                                  *plotData[lev-1],
                                  time_interp_coeff,
                                  0, 0,  numPlotComps);
        }
      // just in case...
      plotData[lev]->exchange();
    } // end loop over levels for computing plot data
  
  // generate plotfile name
  std::string fs;
  if (m_output_file_numbering ==  time_yyyymmdd_360)
    {
      fs.assign("%s%08ld.");
    }
  else
    {
      fs.assign("%s%06ld.");
    }
  char* iter_str = new char[m_plot_prefix.size() + fs.size() + 16];
  sprintf(iter_str, fs.c_str(), m_plot_prefix.c_str(), outputNumbering() );
  string filename(iter_str);

  delete[] iter_str;
 
  // need to pull out SigmaCS pointers:
  Vector<const LevelSigmaCS* > vectCS(m_vect_coordSys.size(), NULL);
  for (int lev=0; lev<numLevels; lev++)
    {
      vectCS[lev] = dynamic_cast<const LevelSigmaCS* >(&(*m_vect_coordSys[lev]));
    }
  if (m_write_map_file)
    {
      WriteSigmaMappedAMRHierarchyHDF5(filename, m_amrGrids, plotData, vectName, 
                                       vectCS, domain, m_dt, m_time,
                                       m_refinement_ratios,
                                       numLevels);
    }
  else
    {
      if (SpaceDim == 1)
        {
          filename.append("1d.hdf5");
        } 
      else if (SpaceDim == 2)
        {
          filename.append("2d.hdf5");
        } 
      else if (SpaceDim == 3)
        {
          filename.append("3d.hdf5");
        }

      this->writeAMRHierarchyHDF5(filename, m_amrGrids, plotData, vectName, 
				  domain, m_amrDx[0], m_dt, time(), m_refinement_ratios, 
				  numLevels);

    }

  // need to delete plotData
  for (int lev=0; lev<numLevels; lev++)
    {
      if (plotData[lev] != NULL)
        {
          delete plotData[lev];
          plotData[lev] = NULL;
        }      

      if (vectC0[lev] != NULL)
        {
          delete vectC0[lev];
          vectC0[lev] = NULL;
        } 

      // if (faceA[lev] != NULL)
      // 	{
      // 	  delete faceA[lev];
      // 	  faceA[lev] = NULL;
      // 	}
    
      // if (viscousTensor[lev] != NULL)
      // 	{
      // 	  delete viscousTensor[lev];
      // 	  viscousTensor[lev] = NULL;
      // 	}
    }
}

/// write checkpoint file out for later restarting
void 
AmrIce::writeCheckpointFile() 
{
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::writeCheckpointfile" << endl;
    }

  CH_TIME("AmrIce::writeCheckpointFile");

  // generate checkpointfile name
  char* iter_str;
  if (m_check_overwrite)
    {
      // overwrite the same checkpoint file, rather than re-writing them
      std::string fs("%s.%dd.hdf5");
      iter_str = new char[m_check_prefix.size() + fs.size() + 16];
      sprintf(iter_str, "%s.%dd.hdf5", m_check_prefix.c_str(), SpaceDim);
      
    }
  else 
    {
      // or hang on to them, if you are a bit sentimental. It's better than keeping
      // every core dump you generate.
      std::string fs;
      if (m_output_file_numbering ==  time_yyyymmdd_360)
	{
	  fs.assign("%s%08ld.%dd.hdf5");
	}
      else
	{
	  fs.assign("%s%06ld.%dd.hdf5");
	}

      iter_str = new char[m_check_prefix.size() + fs.size() + 16];
      sprintf(iter_str, fs.c_str(), m_check_prefix.c_str(), outputNumbering() , SpaceDim );
    }

  CH_assert(iter_str != NULL);

  if (s_verbosity > 3) 
    {
      pout() << "checkpoint file name = " << iter_str << endl;
    }

  writeCheckpointFile(std::string(iter_str));
  delete[] iter_str;
}

/// write checkpoint file out for later restarting
void 
AmrIce::writeCheckpointFile(const string& a_file) 
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::writeCheckpointfile((const string& a_file)" << endl;
      pout() << "checkpoint file name a_file = " << a_file << endl;
    }

#ifdef CH_USE_HDF5

  string thicknessName("thickness");
  Vector<string> vectName(1);
  for (int comp=0; comp<1; comp++)
    {
      char idx[12]; sprintf(idx, "%d", comp);
      vectName[comp] = thicknessName+string(idx);
    } 
  Box domain = m_amrDomains[0].domainBox();
  //int numLevels = m_finest_level +1;      

  
  HDF5Handle handle(a_file.c_str(), HDF5Handle::CREATE);

  // write amr data -- only dump out things which are essential
  // to restarting the computation (i.e. max_level, finest_level, 
  // time, refinement ratios, etc.).  Other paramters (regrid 
  // intervals, block-factor, etc can be changed by the inputs
  // file of the new run.
  // At the moment, the maximum level is not allowed to change,
  // although in principle, there is no real reason why it couldn't
  // 
  HDF5HeaderData&  header = m_headerData;
  header.m_int["max_level"] = m_max_level;
  header.m_int["finest_level"] = m_finest_level;
  header.m_int["current_step"] = m_cur_step;
  header.m_real["time"] = m_time;
  header.m_real["dt"] = m_dt;
  header.m_int["num_comps"] = 2 +  m_velocity[0]->nComp() 
    + m_internalEnergy[0]->nComp() + m_deltaTopography[0]->nComp();
#if BISICLES_Z == BISICLES_LAYERED



  header.m_int["num_comps"] +=2; // surface and base internalEnergys
#endif
  // at the moment, save cfl, but it can be changed by the inputs
  // file if desired.
  header.m_real["cfl"] = m_cfl;

  // periodicity info
  D_TERM(
         if (m_amrDomains[0].isPeriodic(0))
	   header.m_int["is_periodic_0"] = 1;
         else
	   header.m_int["is_periodic_0"] = 0; ,

         if (m_amrDomains[0].isPeriodic(1))
	   header.m_int["is_periodic_1"] = 1;
         else
	   header.m_int["is_periodic_1"] = 0; ,

         if (m_amrDomains[0].isPeriodic(2))
	   header.m_int["is_periodic_2"] = 1;
         else
	   header.m_int["is_periodic_2"] = 0; 
         );
         

  // set up component names
  char compStr[30];
  //string thicknessName("thickness");
  string compName;
  int nComp = 0;
  for (int comp=0; comp < 1; comp++)
    {
      // first generate component name
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = thicknessName + string(idx);
      sprintf(compStr, "component_%04d", comp);
      header.m_string[compStr] = compName;
     
    }
  nComp++;

  string baseHeightName("bedHeight");
  for (int comp=0; comp < 1; comp++)
    {
      // first generate component name
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = baseHeightName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
      
    }
  nComp++;

  string baseDeltaName("deltaBedHeight");
  for (int comp=0; comp < 1; comp++)
    {
      // first generate component name
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = baseDeltaName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
      
    }
  nComp++;


  string iceFracName("iceFrac");
  for (int comp=0; comp < 1; comp++)
    {
      // first generate component name
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = iceFracName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
      
    }
  nComp++;

  string velocityName("velocity");
  for (int comp=0; comp < m_velocity[0]->nComp() ; comp++) 
    {
      // first generate component name
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = velocityName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
    }
  nComp += m_velocity[0]->nComp() ;

  string basalFrictionName("basalFriction");
  for (int comp=0; comp < 1; comp++)
    {
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = basalFrictionName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
    }
  nComp++;

  string muCoefName("muCoef");
  for (int comp=0; comp < 1; comp++)
    {
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = muCoefName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
    }
  nComp++;


  string internalEnergyName("internalEnergy");
  for (int comp=0; comp < m_internalEnergy[0]->nComp() ; comp++) 
    {
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = internalEnergyName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
    }
  
  nComp += m_internalEnergy[0]->nComp() ;

  string tillWaterDepthName("tillWaterDepth");
  for (int comp=0; comp < 1  ; comp++) 
    {
      char idx[12]; sprintf(idx, "%04d", comp);
      compName = tillWaterDepthName + string(idx);
      sprintf(compStr, "component_%04d", comp + nComp);
      header.m_string[compStr] = compName;
    }
  
  nComp += m_tillWaterDepth[0]->nComp() ;

  
#if BISICLES_Z == BISICLES_LAYERED
  {
    sprintf(compStr, "component_%04d", nComp);
    compName = "sInternalEnergy";
    header.m_string[compStr] = compName;
    nComp += 1;
    sprintf(compStr, "component_%04d", nComp);
    compName = "bInternalEnergy";
    header.m_string[compStr] = compName;
    nComp += 1;
    //layer data
    const Vector<Real>& sigma = getFaceSigma();
    string s("sigma");
    for (int l =0; l < sigma.size(); ++l)
      {
	char idx[12]; sprintf(idx, "%04d", l);
	header.m_real[s + string(idx)] = sigma[l];  
      }
  }
#endif

  //allow observers to add checkpoint variables
  for (int i = 0; i < m_observers.size(); i++)
    {
      Vector<std::string> vars;
      m_observers[i]->addCheckVars(vars);
      for (int j = 0; j < vars.size(); j++)
	{
	  nComp++;
	  sprintf(compStr, "component_%04d", nComp);
	  header.m_string[compStr] = vars[j].c_str();
	}
    }

  header.writeToFile(handle);

  // now loop over levels and write out each level's data
  // note that we loop over all allowed levels, even if they
  // are not defined at the moment.
  for (int lev=0; lev<= m_max_level; lev++)
    {
      // set up the level string
      char levelStr[20];
      sprintf(levelStr, "%d", lev);
      const std::string label = std::string("level_") + levelStr;
      
      handle.setGroup(label);
      
      // set up the header info
      HDF5HeaderData levelHeader;
      if (lev < m_max_level)
        {
          levelHeader.m_int["ref_ratio"] = m_refinement_ratios[lev];
        }
      levelHeader.m_real["dx"] = m_amrDx[lev];
      levelHeader.m_box["prob_domain"] = m_amrDomains[lev].domainBox();
      
      levelHeader.writeToFile(handle);
      
      // now write the data for this level
      // only try to write data if level is defined.
      if (lev <= m_finest_level)
        {
          write(handle, m_amrGrids[lev]);

	  const IntVect ghost = IntVect::Unit*2;
          const LevelSigmaCS& levelCS = *m_vect_coordSys[lev];

	  write(handle, levelCS.getH() , "thicknessData", levelCS.getH().ghostVect());

	  write(handle, levelCS.getTopography() , "bedHeightData",
		levelCS.getTopography().ghostVect()  );

	  write(handle, *m_deltaTopography[lev] , "deltaBedHeightData",
		m_deltaTopography[lev]->ghostVect()  );

	  write(handle, *m_iceFrac[lev] , "iceFracData",
		m_iceFrac[lev]->ghostVect()  );

	  write(handle, *m_velocity[lev], "velocityData", 
		m_velocity[lev]->ghostVect());

	  write(handle, *m_velBasalC[lev], "basalFrictionData", 
		m_velBasalC[lev]->ghostVect());

	  write(handle, *m_cellMuCoef[lev], "muCoefData", 
		m_cellMuCoef[lev]->ghostVect());
	  
	  write(handle, *m_internalEnergy[lev], "internalEnergyData", 
		m_internalEnergy[lev]->ghostVect());

	  write(handle, *m_tillWaterDepth[lev], "tillWaterDepthData", 
		m_tillWaterDepth[lev]->ghostVect());

#if BISICLES_Z == BISICLES_LAYERED
	  write(handle, *m_sInternalEnergy[lev], "sInternalEnergyData", 
		m_sInternalEnergy[lev]->ghostVect());
	  write(handle, *m_bInternalEnergy[lev], "bInternalEnergyData", 
		m_bInternalEnergy[lev]->ghostVect());
#endif
	  //allow observers to write to the checkpoint
	  for (int i = 0; i < m_observers.size(); i++)
	    {
	      m_observers[i]->writeCheckData(handle, lev);
	    }
        }
    }// end loop over levels
  
  handle.close();
#endif
}


/// read checkpoint file for restart 
void 
AmrIce::readCheckpointFile(HDF5Handle& a_handle)
{

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::readCheckpointFile" << endl;
    }

#ifndef CH_USE_HDF5
  MayDay::Error("code must be compiled with HDF5 to read checkpoint files");
#endif

#ifdef CH_USE_HDF5

 
  HDF5HeaderData& header = m_headerData;
  header.readFromFile(a_handle);
  ParmParse ppAmr("amr");
  
  //check for various components. Maybe rethink this when HDF5::SetGroup
  //is fixed...
  bool containsDeltaBedHeight(false);
  bool containsInternalEnergy(false);
  bool containsTillWaterDepth(false);
  bool containsTemperature(false);
  bool containsIceFrac(false);
  bool containsBasalFriction(false);
  bool containsMuCoef(false);

  map<std::string, std::string>::const_iterator i;
  for (i = header.m_string.begin(); i!= header.m_string.end(); ++i)
    {
      if (i->second == "deltaBedHeight0000")
	{
	  containsDeltaBedHeight = true;
	}
      if (i->second == "temperature0000")
	{
	  containsTemperature = true;
	}
      if (i->second == "internalEnergy0000")
	{
	  containsInternalEnergy = true;
	}
      if (i->second == "tillWaterDepth0000")
	{
	  containsTillWaterDepth = true;
	}
      if (i->second == "iceFrac0000")
	{
	  containsIceFrac = true;
	  ppAmr.query("restart_read_ice_frac", containsIceFrac); // backward compatibility option.
	}
      if (i->second == "basalFriction0000")
	{
	  containsBasalFriction = true;
	}
      if (i->second == "muCoef0000")
	{
	  containsMuCoef = true;
	}
    }

  if (s_verbosity >= 3)
    {
      pout() << "hdf5 header data: " << endl;
      pout() << header << endl;
    }

  // read max level
  if (header.m_int.find("max_level") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain max_level");
    }
  // we can change max level upon restart
  int max_level_check = header.m_int["max_level"];
  if (max_level_check != m_max_level)
    {
      if (s_verbosity > 0)
        {
          pout() << "Restart file has a different max level than inputs file"
                 << endl;
          pout() << "     max level from inputs file = " 
                 << m_max_level << endl;
          pout() << "     max level in checkpoint file = " 
                 << max_level_check << endl;                 
          pout() << "Using max level from inputs file" << endl;
        }
    }
  // read finest level
  if (header.m_int.find("finest_level") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain finest_level");
    }

  m_finest_level = header.m_int["finest_level"];
  if (m_finest_level > m_max_level)
    {
      MayDay::Error("finest level in restart file > max allowable level!");
    }

  // read current step
  if (header.m_int.find("current_step") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain current_step");
    }

  m_cur_step = header.m_int["current_step"];
  m_restart_step = m_cur_step;

  // optionally, over-ride the step number in the restart checkpoint file with one specified in the inputs
 
  if (ppAmr.contains("restart_step") )
    {
      int restart_step;
      ppAmr.get("restart_step", restart_step);
      m_cur_step = restart_step;
      m_restart_step = restart_step;
    }

  
  // read time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain time");
    }

  m_time = header.m_real["time"];

  //optionally, over-ride the time in the restart checkpoint file with one specified in the inputs      
  if (ppAmr.contains("restart_time") )
    {
      bool set_time = true;
      ppAmr.query("restart_set_time",set_time); // set amr.restart_set_time = false to prevent time reset
      if (set_time){
	Real restart_time;
	ppAmr.get("restart_time", restart_time);
	m_time = restart_time;
      }
    }
  

  m_dt = header.m_real["dt"];

  // read num comps
  if (header.m_int.find("num_comps") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain num_comps");
    }
  
  //int numComps = header.m_int["num_comps"];

  // read cfl
  if (header.m_real.find("cfl") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain cfl");
    }

  Real check_cfl = header.m_real["cfl"];
  ParmParse ppCheck("amr");

  if (ppCheck.contains("cfl"))
    { 
      // check for consistency and warn if different
      if (check_cfl != m_cfl)
	{
	  if (s_verbosity > 0)
	    {
	      pout() << "CFL in checkpoint file different from inputs file" 
		     << endl;
	      pout() << "     cfl in inputs file = " << m_cfl << endl;
	      pout() << "     cfl in checkpoint file = " << check_cfl 
		     << endl;
	      pout() << "Using cfl from inputs file" << endl;                
	    }
	}  // end if cfl numbers differ
    } // end if cfl present in inputs file
  else
    {
      m_cfl = check_cfl;
    }          

  // read periodicity info
  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility 
  // bool isPeriodic[SpaceDim];
  // D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
  //          isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
  //        else
  //          isPeriodic[0] = false; ,

  //        if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
  //          isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
  //        else
  //          isPeriodic[1] = false; ,

  //        if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
  //          isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
  //        else
  //          isPeriodic[2] = false;);

#if BISICLES_Z == BISICLES_LAYERED
  {
    //retrieve and *check* sigma data
    Vector<Real> sigma;
    int l = 0;
    string s("sigma");
    bool found = false;
    do {
      char idx[12]; sprintf(idx, "%04d", l);
      string ss = s + string(idx);
      map<std::string, Real>::const_iterator it = header.m_real.find(ss);
      found = (it != header.m_real.end());
      if (found)
	{
	  sigma.push_back(it->second);
	  CH_assert( Abs( sigma[l] - m_faceSigma[l] ) < TINY_NORM); 
	  ++l;
	}
      
    } while (found);
    CH_assert(m_nLayers == sigma.size() - 1); // layer number has not changed
  }
#endif

  // now resize stuff 
  m_amrDomains.resize(m_max_level+1);
  m_amrGrids.resize(m_max_level+1);
  m_amrDx.resize(m_max_level+1);
  m_old_thickness.resize(m_max_level+1, NULL);
  m_iceFrac.resize(m_max_level+1, NULL);
  m_velocity.resize(m_max_level+1, NULL);
  m_diffusivity.resize(m_max_level+1);
  m_vect_coordSys.resize(m_max_level+1);
  m_velRHS.resize(m_max_level+1);
  m_surfaceThicknessSource.resize(m_max_level+1,NULL);
  m_volumeThicknessSource.resize(m_max_level+1,NULL);
  m_basalThicknessSource.resize(m_max_level+1,NULL);
  m_calvedIceThickness.resize(m_max_level+1, NULL);
  m_removedIceThickness.resize(m_max_level+1, NULL);
  m_addedIceThickness.resize(m_max_level+1, NULL);
  m_deltaTopography.resize(m_max_level+1, NULL);
  m_divThicknessFlux.resize(m_max_level+1,NULL);
  m_velBasalC.resize(m_max_level+1,NULL);
  m_cellMuCoef.resize(m_max_level+1,NULL);
  m_faceVelAdvection.resize(m_max_level+1,NULL);
  m_faceVelTotal.resize(m_max_level+1,NULL);
  m_internalEnergy.resize(m_max_level+1,NULL);
  m_tillWaterDepth.resize(m_max_level+1,NULL);
#if BISICLES_Z == BISICLES_LAYERED
  m_sInternalEnergy.resize(m_max_level+1,NULL);
  m_bInternalEnergy.resize(m_max_level+1,NULL);
  m_sHeatFlux.resize(m_max_level+1,NULL);
  m_bHeatFlux.resize(m_max_level+1,NULL);
  m_layerSFaceXYVel.resize(m_max_level+1,NULL);
  m_layerSFaceSVel.resize(m_max_level+1,NULL);
  m_layerXYFaceXYVel.resize(m_max_level+1,NULL);
#endif
  IntVect sigmaCSGhost = m_num_thickness_ghost*IntVect::Unit;
	 

  // now read in level-by-level data
  for (int lev=0; lev<= m_max_level; lev++)
    {
      // set up the level string
      char levelStr[20];
      sprintf(levelStr, "%d", lev);
      const std::string label = std::string("level_") + levelStr;
      
      a_handle.setGroup(label);

      // read header info
      HDF5HeaderData levheader;
      levheader.readFromFile(a_handle);
      
      if (s_verbosity >= 3)
        {
          pout() << "level " << lev << " header data" << endl;
          pout() << levheader << endl;
        }

      // Get the refinement ratio
      if (lev < max_level_check)
        {
          int checkRefRatio;
          if (levheader.m_int.find("ref_ratio") == levheader.m_int.end())
            {
              MayDay::Error("checkpoint file does not contain ref_ratio");
            }
          checkRefRatio = levheader.m_int["ref_ratio"];

          // check for consistency
          if (checkRefRatio != m_refinement_ratios[lev])
            {
	      
	      MayDay::Error("inputs file and checkpoint file ref ratios inconsistent");
            }
        }
      
      // read dx
      if (levheader.m_real.find("dx") == levheader.m_real.end())
        {
          MayDay::Error("checkpoint file does not contain dx");
        }
      
      if (lev <= max_level_check)
	{
	  if ( Abs(m_amrDx[lev] - levheader.m_real["dx"]) > TINY_NORM )
	    {
	      MayDay::Error("restart file dx != input file dx");
	    }
	  
	  // read problem domain box
	  if (levheader.m_box.find("prob_domain") == levheader.m_box.end())
	    {
	      MayDay::Error("checkpoint file does not contain prob_domain");
	    }
	  Box domainBox = levheader.m_box["prob_domain"];

	  if (m_amrDomains[lev].domainBox() != domainBox)
	    { 
	      MayDay::Error("restart file domain != input file domain");
	    }
	}

      // the rest is only applicable if this level is defined
      if (lev <=  m_finest_level && lev <= max_level_check)
        {
          // read grids          
          Vector<Box> grids;
          const int grid_status = read(a_handle, grids);
          if (grid_status != 0) 
            {
              MayDay::Error("checkpoint file does not contain a Vector<Box>");
            }
          // do load balancing
          int numGrids = grids.size();
          Vector<int> procIDs(numGrids);
          LoadBalance(procIDs, grids);
          DisjointBoxLayout levelDBL(grids, procIDs, m_amrDomains[lev]);
          m_amrGrids[lev] = levelDBL;

          // allocate this level's storage
	  // 4 ghost cells needed for advection.
          m_old_thickness[lev] = new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
#if BISICLES_Z == BISICLES_LAYERED
	  m_internalEnergy[lev] =  new LevelData<FArrayBox>
	    (levelDBL, m_nLayers, m_num_thickness_ghost*IntVect::Unit);
	  m_tillWaterDepth[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
	  m_sInternalEnergy[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
	  m_bInternalEnergy[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
	  m_sHeatFlux[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
	  m_bHeatFlux[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
#elif BISICLES_Z == BISICLES_FULLZ
	  m_internalEnergy[lev] =  new LevelData<FArrayBox>
	    (levelDBL, 1, m_num_thickness_ghost*IntVect::Unit);
#endif
	  // other quantities need only one;
	  IntVect ghostVect(IntVect::Unit);
          m_velocity[lev] = new LevelData<FArrayBox>(levelDBL, SpaceDim, 
                                                     ghostVect);

          m_iceFrac[lev] = new LevelData<FArrayBox>(levelDBL, 1, ghostVect);

	  m_faceVelAdvection[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Unit);
	  m_faceVelTotal[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Unit);
#if BISICLES_Z == BISICLES_LAYERED
	  m_layerXYFaceXYVel[lev] = new LevelData<FluxBox>
	    (m_amrGrids[lev], m_nLayers, IntVect::Unit);
	  m_layerSFaceXYVel[lev] = new LevelData<FArrayBox>
	    (m_amrGrids[lev], SpaceDim*(m_nLayers + 1), IntVect::Unit);
	  m_layerSFaceSVel[lev] = new LevelData<FArrayBox>
	    (m_amrGrids[lev], m_nLayers + 1, IntVect::Unit);
	  
#endif

	  m_velBasalC[lev] = new LevelData<FArrayBox>(levelDBL, 1, ghostVect);
	  m_cellMuCoef[lev] = new LevelData<FArrayBox>(levelDBL, 1, ghostVect);
	  m_velRHS[lev] = new LevelData<FArrayBox>(levelDBL, SpaceDim, IntVect::Zero);
	  m_surfaceThicknessSource[lev] =  new LevelData<FArrayBox>(levelDBL,   1, IntVect::Unit) ;
	  m_volumeThicknessSource[lev] = new LevelData<FArrayBox>(levelDBL, 1, IntVect::Unit);
	  m_basalThicknessSource[lev] = new LevelData<FArrayBox>(levelDBL,   1, IntVect::Unit) ;
	  m_calvedIceThickness[lev] =  new LevelData<FArrayBox>(levelDBL,   1, IntVect::Unit) ;
	  m_removedIceThickness[lev] =  new LevelData<FArrayBox>(levelDBL,   1, IntVect::Unit) ;
	  m_addedIceThickness[lev] =  new LevelData<FArrayBox>(levelDBL,   1, IntVect::Unit) ;
	  m_deltaTopography[lev] =  new LevelData<FArrayBox>(levelDBL,   1, IntVect::Zero) ;
	  m_divThicknessFlux[lev] =  new LevelData<FArrayBox>(levelDBL,   1, IntVect::Zero) ;
	  m_diffusivity[lev] = new LevelData<FluxBox>(levelDBL, 1, IntVect::Zero);

          // read this level's data
          LevelData<FArrayBox>& old_thickness = *m_old_thickness[lev];  
	  LevelData<FArrayBox> tmpThickness;
	  tmpThickness.define(old_thickness);

          int dataStatus = read<FArrayBox>(a_handle, tmpThickness, "thicknessData", levelDBL);
	  for (DataIterator dit(levelDBL);dit.ok();++dit)
	    {
	      old_thickness[dit].copy(tmpThickness[dit]);
	    }

          if (dataStatus != 0)
            {
              MayDay::Error("checkpoint file does not contain thickness data");
            }

	  LevelData<FArrayBox> bedHeight;
	  bedHeight.define(old_thickness);
	  dataStatus = read<FArrayBox>(a_handle, bedHeight, "bedHeightData", levelDBL);

	  if (dataStatus != 0)
            {
              MayDay::Error("checkpoint file does not contain bed height data");
            }
	  
	  LevelData<FArrayBox> deltaBedHeight;
	  deltaBedHeight.define(old_thickness);
	  
	  if (containsDeltaBedHeight)
	    {
	      dataStatus = read<FArrayBox>(a_handle,deltaBedHeight,"deltaBedHeightData",levelDBL);
	      if (dataStatus != 0)
		{
		  MayDay::Error("checkpoint file does not contain delta bed height data");
		}
	    }
	  else
	    {
	      for (DataIterator dit = deltaBedHeight.disjointBoxLayout(); dit.ok();++dit)
		{
		  deltaBedHeight[dit].setVal(0.0);
		}
	    }

	  //having read thickness and base data, we can define
          //the co-ordinate system 
	  RealVect dx = m_amrDx[lev]*RealVect::Unit;
          m_vect_coordSys[lev] = RefCountedPtr<LevelSigmaCS >
            (new LevelSigmaCS(m_amrGrids[lev], dx, sigmaCSGhost));
	  m_vect_coordSys[lev]->setIceDensity(m_iceDensity);
	  m_vect_coordSys[lev]->setWaterDensity(m_seaWaterDensity);
	  m_vect_coordSys[lev]->setGravity(m_gravity);
#if BISICLES_Z == BISICLES_LAYERED
	  m_vect_coordSys[lev]->setFaceSigma(getFaceSigma());
#endif
          LevelSigmaCS& levelCS = *m_vect_coordSys[lev];
          LevelData<FArrayBox>& levelH = levelCS.getH();

          DataIterator dit = levelH.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              levelH[dit].copy((*m_old_thickness[lev])[dit]);
            }
          levelCS.setTopography(bedHeight);
	  if (deltaBedHeight.isDefined())
	    {
	      for (dit.begin(); dit.ok(); ++dit)
		{
		  (*m_deltaTopography[lev])[dit].copy(deltaBedHeight[dit]);
		} 
	    }

	  {
	    LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
	    int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
	    levelCS.recomputeGeometry(crseCoords, refRatio);
	  }
          LevelData<FArrayBox>& velData = *m_velocity[lev];
	  dataStatus = read<FArrayBox>(a_handle,
                                       velData,
                                       "velocityData",
				       levelDBL);
	  m_velocitySolveInitialResidualNorm = 1.0e+6; //\todo fix this

	  //this check doesn't work, because HDF5::SetGroup attempts
	  //to create a group if it doesn't exist. And since the file has been
	  // opened readonly, the previous call generates an exception
 			       
          if (dataStatus != 0)
            {
	      MayDay::Error("checkpoint file does not contain velocity data");
	      
            }


	  if (containsIceFrac)
	    {
	      LevelData<FArrayBox>& iceFracData = *m_iceFrac[lev];
	      dataStatus = read<FArrayBox>(a_handle,
					   iceFracData,
					   "iceFracData",
				       levelDBL);
          
	      /// note that although this check appears to work, it makes a mess of a_handle and the next lot of data are not read...
	      if (dataStatus != 0)
		{
		  MayDay::Warning("checkpoint file does not contain ice fraction data -- initializing based on current ice thicknesses"); 
		  const LevelData<FArrayBox>& levelThickness = m_vect_coordSys[lev]->getH();
		  setIceFrac(levelThickness, lev);
		} // end if no ice fraction in data
	      else
		{
		  // ensure that ice fraction is set to zero where there's no ice
		  // or not, since this should have been done before writing
		  // updateIceFrac(m_vect_coordSys[lev]->getH(), lev);
		}
	    } 
	  else
	    {
	      MayDay::Warning("checkpoint file does not contain ice fraction data -- initializing based on current ice thicknesses"); 
	      const LevelData<FArrayBox>& levelThickness = m_vect_coordSys[lev]->getH();
	      setIceFrac(levelThickness, lev);
	    }

	  if (containsBasalFriction)
	    {
	      dataStatus = read<FArrayBox>(a_handle, *m_velBasalC[lev],
					   "basalFrictionData", levelDBL);
	    }
	  else
	    {
	      MayDay::Warning("checkpoint file does not basal friction coefficient data"); 
	    }

	  if (containsMuCoef)
	    {
	      dataStatus = read<FArrayBox>(a_handle, *m_cellMuCoef[lev],
					   "muCoefData", levelDBL);
	    }
	  else
	    {
	      MayDay::Warning("checkpoint file does not mu coefficient data"); 
	    }

	  {
	    // read internal energy , or read temperature and convert to internal energy
 	    std::string dataName, sDataName, bDataName;

	    if (containsInternalEnergy)
	      {
		dataName = "internalEnergyData";
		sDataName = "sInternalEnergyData";
		bDataName = "bInternalEnergyData";
	      }
	    else if (containsTemperature)
	      {
		dataName = "temperatureData";
		sDataName = "sTemperatureData";
		bDataName = "bTemperatureData";	
	      }
	    else
	      {
		MayDay::Error("checkpoint file does not contain internal energy or temperature data"); 
	      }
	    

	    LevelData<FArrayBox>& internalEnergyData = *m_internalEnergy[lev];
	    dataStatus = read<FArrayBox>(a_handle, internalEnergyData, dataName,levelDBL);
	    if (dataStatus != 0)
	      {
		MayDay::Error("checkpoint file does not contain internal energy data"); 
	      }
	    

#if BISICLES_Z == BISICLES_LAYERED	  
	    LevelData<FArrayBox>& sInternalEnergyData = *m_sInternalEnergy[lev];
	    dataStatus = read<FArrayBox>(a_handle,  sInternalEnergyData, sDataName,levelDBL);	  
	    if (dataStatus != 0)
	      {
		MayDay::Error("checkpoint file does not contain surface internal energy data"); 
	      }
	    
	    LevelData<FArrayBox>& bInternalEnergyData = *m_bInternalEnergy[lev];
	    dataStatus = read<FArrayBox>(a_handle, bInternalEnergyData, bDataName, levelDBL);	  
	    if (dataStatus != 0)
	      {
		MayDay::Error("checkpoint file does not contain basal internal energy data"); 
	      }
#endif
	  
	    if (containsTemperature)
	      {
		// need to covert tempearture to internal energy
		for (DataIterator dit(levelDBL);dit.ok();++dit)
		  {
		    FArrayBox& E = (*m_internalEnergy[lev])[dit];
		    FArrayBox T(E.box(),E.nComp()); T.copy(E);
		    IceThermodynamics::composeInternalEnergy(E,T,E.box() );
#if BISICLES_Z == BISICLES_LAYERED
		    FArrayBox& sE = (*m_sInternalEnergy[lev])[dit];
		    FArrayBox sT(sE.box(),sE.nComp()); sT.copy(sE);
		    IceThermodynamics::composeInternalEnergy(sE,sT,sE.box(),false );
		    // it is possible for sT to be meaningless, so avoid test 
		    FArrayBox& bE = (*m_bInternalEnergy[lev])[dit];
		    FArrayBox bT(bE.box(),bE.nComp()); bT.copy(bE);
		    IceThermodynamics::composeInternalEnergy(bE,bT,bE.box(),false); 
		    // it is possible for bT to be meaningless, so avoid test 
#endif	    
		  }
	      }
	    CH_assert(m_internalEnergy[lev]->nComp() == m_nLayers);

	    LevelData<FArrayBox>& tillWaterDepthData = *m_tillWaterDepth[lev];
	    if (containsTillWaterDepth)
	      {
		dataStatus = read<FArrayBox>(a_handle, tillWaterDepthData, "tillWaterDepthData" ,levelDBL);
		if (dataStatus != 0)
		  {
		    MayDay::Error("checkpoint file does not contain till water depth data"); 
		  }
	      }
	    else
	      {
		for (DataIterator dit = tillWaterDepthData.disjointBoxLayout(); dit.ok();++dit)
		  {
		    tillWaterDepthData[dit].setVal(0.0);
		  }
	      }
	  }

	  //allow observers to read from the checkpoint
	  for (int i = 0; i < m_observers.size(); i++)
	    {
	      m_observers[i]->readCheckData(a_handle, header,  lev, levelDBL);
	    }


	} // end if this level is defined
    } // end loop over levels                                    
          
  // do we need to close the handle?
  
  //update the temperature since it depends on the internal energy
  updateTemperature();
  
  
  //this is just to make sure the diffusivity is computed
  //(so I should improve that)
  defineSolver();
  m_doInitialVelSolve = false; // since we have just read the velocity field
  m_doInitialVelGuess = false; // ditto

  
  if (dynamic_cast<InverseIceVelocitySolver*>(m_velSolver))
    {
      //special case for inverse problems : in most cases the basal friction and
      //mu coeffcient data read from the checkpoint file will be overwritten
      //in solveVelocityField(). Here we want them to be input data 
      if (containsBasalFriction)
	{	
	  if (m_basalFrictionPtr) delete m_basalFrictionPtr; 
	  m_basalFrictionPtr = InverseIceVelocitySolver::basalFriction(m_velBasalC, refRatios(), dx(0));
	}
      
      if (containsMuCoef)
	{ 	
	  if (m_muCoefficientPtr) delete m_muCoefficientPtr;
	  m_muCoefficientPtr = InverseIceVelocitySolver::muCoefficient(m_cellMuCoef, refRatios(), dx(0));
	}
    }

  if (s_verbosity > 3) 
    {
      pout() << "AmrIce::readCheckPointFile solveVelocityField() " << endl;
    }
  solveVelocityField();
  m_doInitialVelSolve = true;

#endif
  
}

#ifdef CH_USE_HDF5

void AmrIce::writeMetaDataHDF5(HDF5Handle& a_handle) const
{
  //Additional data (BISICLES specific)
  HDF5HeaderData headerData;
  headerData.m_int["max_level"] = m_max_level;
  headerData.m_int["finest_level"] = m_finest_level;
  headerData.m_int["current_step"] = m_cur_step; 
  headerData.m_real["time"] = time();
  headerData.m_real["dt"] = m_dt;
  headerData.m_string["svn_version"] = SVN_REV;
  headerData.m_string["svn_repository"] = SVN_REP;
  headerData.m_string["svn_url"] = SVN_URL;
  headerData.m_int["bisicles_version_major"] = BISICLES_VERSION_MAJOR;
  headerData.m_int["bisicles_version_minor"] = BISICLES_VERSION_MINOR;
  headerData.m_int["bisicles_patch_number"] = BISICLES_PATCH_NUMBER;
  headerData.m_int["chombo_version_major"] = CHOMBO_VERSION_MAJOR;
  headerData.m_int["chombo_version_minor"] = CHOMBO_VERSION_MINOR;
#ifdef CHOMBO_TRUNK
  headerData.m_int["chombo_patch_number"] = -1;
#else
  headerData.m_int["chombo_patch_number"] = CHOMBO_PATCH_NUMBER;
#endif
  headerData.m_real["density_of_ice"] = m_iceDensity;
  headerData.m_real["seconds_per_unit_time"] = m_seconds_per_unit_time;

  // Additional data that doesn't mean much to bisicles but might to other programs,
  // coordinate systsem data
  {
    ParmParse pp("CRS");
    int epsg = 0; // EPSG is the format I like, becuuse it has a single integer value \todo WKT? PROJ.4? 
    pp.query("EPSG", epsg);
    if (epsg != 0)
      {
	headerData.m_int["crs_EPSG"] = epsg;
      }
    // x-,y-cordinate of BISICLES origin in the CRS
    Real x0 = 0.0; 
    pp.query("origin_x", x0);
    headerData.m_real["crs_origin_x"] = x0;
    Real y0 = 0.0; 
    pp.query("origin_y", y0);
    headerData.m_real["crs_origin_y"] = y0;
  }
  
  headerData.writeToFile(a_handle);
}


void AmrIce::writeAMRHierarchyHDF5(HDF5Handle& a_handle, 
				   const Vector<DisjointBoxLayout>& a_grids,
				   const Vector<LevelData<FArrayBox>* > & a_data,
				   const Vector<string>& a_name,
				   const Box& a_domain,
				   const Real& a_dx,
				   const Real& a_dt,
				   const Real& a_time,
				   const Vector<int>& a_ratio,
				   const int& a_numLevels) const
{
 
     
  //Chombo AMR data (VisIt compatible)
  WriteAMRHierarchyHDF5(a_handle, a_grids, a_data, a_name, 
			a_domain, a_dx, a_dt, a_time, a_ratio, 
			a_numLevels);

  writeMetaDataHDF5(a_handle);
  
}

void AmrIce::writeAMRHierarchyHDF5(const string& filename, 
				   const Vector<DisjointBoxLayout>& a_grids,
				   const Vector<LevelData<FArrayBox>* > & a_data,
				   const Vector<string>& a_name,
				   const Box& a_domain,
				   const Real& a_dx,
				   const Real& a_dt,
				   const Real& a_time,
				   const Vector<int>& a_ratio,
				   const int& a_numLevels) const
{
 
  HDF5Handle handle(filename.c_str(), HDF5Handle::CREATE);
  this->writeAMRHierarchyHDF5(handle, a_grids, a_data, a_name, 
			      a_domain, a_dx, a_dt, a_time, a_ratio, 
			      a_numLevels);  
  handle.close();
}


/// set up for restart
void 
AmrIce::restart(const string& a_restart_file)
{
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::restart" << endl;
    }

  HDF5Handle handle(a_restart_file, HDF5Handle::OPEN_RDONLY);
  // first read in data from checkpoint file
  readCheckpointFile(handle);
  handle.close();
  // don't think I need to do anything else, do I?
  
  //special case for inverse problems
  InverseIceVelocitySolver* invPtr = dynamic_cast<InverseIceVelocitySolver*>(m_velSolver);
  if (invPtr)
    {
      invPtr->setPreviousTime(m_time);
    }
}


#define CFIO_MAX_BOX_SIZE 32
#define CFIO_BLOCK_FACTOR 2
#define CFIO_PARM_PARSE_PREFIX "CFIO"
#define CFIO_PARM_PARSE_FIELD_PREFIX "CFIO"
/*
  see 
  http://cfconventions.org/Data/cf-standard-names/52/build/cf-standard-name-table.html
  for cf land ice names
*/
#define CFIO_FIELD_LAND_ICE_THICKNESS_SHORT_NAME "lithk"
#define CFIO_FIELD_SURFACE_ALTITUDE_SHORT_NAME "orog"
#define CFIO_FIELD_BEDROCK_ALTITUDE_SHORT_NAME "topg"
#define CFIO_FIELD_LAND_ICE_BASAL_VELOCITY_SHORT_NAME "velbase"
#define CFIO_FIELD_LAND_ICE_SURFACE_VELOCITY_SHORT_NAME "velsurf"
#define CFIO_FIELD_LAND_ICE_BASAL_X_VELOCITY_SHORT_NAME "xvelbase"
#define CFIO_FIELD_LAND_ICE_BASAL_Y_VELOCITY_SHORT_NAME "yvelbase"
#define CFIO_FIELD_LAND_ICE_SURFACE_X_VELOCITY_SHORT_NAME "xvelsurf"
#define CFIO_FIELD_LAND_ICE_SURFACE_Y_VELOCITY_SHORT_NAME "yvelsurf"
#define CFIO_FIELD_LAND_ICE_AREA_FRACTION_SHORT_NAME "sftgif"
#define CFIO_FIELD_GROUNDED_ICE_AREA_FRACTION_SHORT_NAME "sftgrf"
#define CFIO_FIELD_FLOATING_ICE_AREA_FRACTION_SHORT_NAME "sftflf"
#define CFIO_FIELD_TEMPERATURE_ICE_TOP_SHORT_NAME "litemptop"
#define CFIO_FIELD_TEMPERATURE_ICE_BASE_SHORT_NAME "litempbot"
#define CFIO_FIELD_TEMPERATURE_GROUNDED_ICE_BASE_SHORT_NAME "litempbotgr"
#define CFIO_FIELD_TEMPERATURE_FLOATING_ICE_BASE_SHORT_NAME "litempbotfl"
#define CFIO_FIELD_MAGNITUDE_BASAL_DRAG_SHORT_NAME "strbasemag"
 
#define CFIO_FIELD_LAND_ICE_THICKNESS_IMBALANCE_SHORT_NAME "dlithkdt"
#define CFIO_FIELD_LAND_ICE_CALVING_FLUX_SHORT_NAME "licalvf"
#define CFIO_FIELD_SURFACE_MASS_BALANCE_FLUX_SHORT_NAME "acabf"
#define CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_SHORT_NAME "libmassbf"
#define CFIO_FIELD_GROUNDED_ICE_BASAL_MASS_BALANCE_FLUX_SHORT_NAME "libmassbfgr"
#define CFIO_FIELD_FLOATING_ICE_BASAL_MASS_BALANCE_FLUX_SHORT_NAME "libmassbffl"
#define CFIO_FIELD_ICE_FRONT_CALVING_AND_MELT_SHORT_NAME "lifmassbf"
#define CFIO_FIELD_GEOTHERMAL_FLUX_SHORT_NAME "hfgeoubed"

#define CFIO_FIELD_LAND_ICE_THICKNESS_CF_NAME "land_ice_thickness"
#define CFIO_FIELD_LAND_ICE_BASAL_VELOCITY_CF_NAME "land_ice_basal_velocity"
#define CFIO_FIELD_LAND_ICE_SURFACE_VELOCITY_CF_NAME "land_ice_surface_velocity"
#define CFIO_FIELD_SURFACE_ALTITUDE_CF_NAME "surface_altitude"
#define CFIO_FIELD_BEDROCK_ALTITUDE_CF_NAME "bedrock_altitude"
#define CFIO_FIELD_LAND_ICE_BASAL_X_VELOCITY_CF_NAME "land_ice_basal_x_velocity"
#define CFIO_FIELD_LAND_ICE_BASAL_Y_VELOCITY_CF_NAME "land_ice_basal_y_velocity"
#define CFIO_FIELD_LAND_ICE_SURFACE_X_VELOCITY_CF_NAME "land_ice_surface_x_velocity"
#define CFIO_FIELD_LAND_ICE_SURFACE_Y_VELOCITY_CF_NAME "land_ice_surface_y_velocity"
#define CFIO_FIELD_LAND_ICE_AREA_FRACTION_CF_NAME "land_ice_area_fraction"
#define CFIO_FIELD_GROUNDED_ICE_AREA_FRACTION_CF_NAME "grounded_ice_sheet_area_fraction"
#define CFIO_FIELD_FLOATING_ICE_AREA_FRACTION_CF_NAME "floating_ice_sheet_area_fraction"
#define CFIO_FIELD_TEMPERATURE_ICE_TOP_CF_NAME "temperature_at_top_of_ice_sheet_model"
#define CFIO_FIELD_TEMPERATURE_ICE_BASE_CF_NAME "temperature_at_base_of_ice_sheet_model"
#define CFIO_FIELD_TEMPERATURE_GROUNDED_ICE_BASE_CF_NAME "temperature_at_base_of_ice_sheet_model"
#define CFIO_FIELD_TEMPERATURE_FLOATING_ICE_BASE_CF_NAME "temperature_at_base_of_ice_sheet_model" 
#define CFIO_FIELD_MAGNITUDE_BASAL_DRAG_CF_NAME "land_ice_basal_drag" 

#define CFIO_FIELD_LAND_ICE_CALVING_FLUX_CF_NAME "land_ice_specific_mass_flux_due_to_calving"
#define CFIO_FIELD_SURFACE_MASS_BALANCE_FLUX_CF_NAME "land_ice_surface_specific_mass_balance_flux"
#define CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_CF_NAME "land_ice_basal_specific_mass_balance_flux"
#define CFIO_FIELD_GROUNDED_ICE_BASAL_MASS_BALANCE_FLUX_CF_NAME "land_ice_basal_specific_mass_balance_flux_ice"
#define CFIO_FIELD_FLOATING_ICE_BASAL_MASS_BALANCE_FLUX_CF_NAME "land_ice_basal_specific_mass_balance_flux_ice"

#define CFIO_FIELD_LAND_ICE_THICKNESS_IMBALANCE_CF_NAME "tendency_of_land_ice_thickness"
#define CFIO_FIELD_ICE_FRONT_CALVING_AND_MELT_CF_NAME "land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting"
#define CFIO_FIELD_GEOTHERMAL_FLUX_CF_NAME "upward_geothermal_heat_flux_at_ground_level_in_land_ice"

#define CFIO_FIELD_LAND_ICE_THICKNESS_LONG_NAME "Ice-sheet thickness"
#define CFIO_FIELD_SURFACE_ALTITUDE_LONG_NAME "Ice-sheet altitude"
#define CFIO_FIELD_BEDROCK_ALTITUDE_LONG_NAME "Bedrock altitude"
#define CFIO_FIELD_LAND_ICE_BASAL_X_VELOCITY_LONG_NAME "X-component of land ice basal velocity"
#define CFIO_FIELD_LAND_ICE_BASAL_Y_VELOCITY_LONG_NAME "Y-component of land ice basal velocity"
#define CFIO_FIELD_LAND_ICE_SURFACE_X_VELOCITY_LONG_NAME "X-component of land ice surface velocity"
#define CFIO_FIELD_LAND_ICE_SURFACE_Y_VELOCITY_LONG_NAME "Y-component of land ice surface velocity"
#define CFIO_FIELD_LAND_ICE_AREA_FRACTION_LONG_NAME "Land ice area fraction"
#define CFIO_FIELD_GROUNDED_ICE_AREA_FRACTION_LONG_NAME "Grounded ice area fraction"
#define CFIO_FIELD_FLOATING_ICE_AREA_FRACTION_LONG_NAME "Floating ice-sheet area fraction"
#define CFIO_FIELD_TEMPERATURE_ICE_TOP_LONG_NAME "Temperature at top of ice-sheet model"
#define CFIO_FIELD_TEMPERATURE_ICE_BASE_LONG_NAME "Basal temperature of ice sheet"
#define CFIO_FIELD_TEMPERATURE_GROUNDED_ICE_BASE_LONG_NAME "Basal temperature of grounded ice sheet"
#define CFIO_FIELD_TEMPERATURE_FLOATING_ICE_BASE_LONG_NAME "Basal temperature of floating ice shelf"
#define CFIO_FIELD_MAGNITUDE_BASAL_DRAG_LONG_NAME "Land ice basal drag" 

#define CFIO_FIELD_LAND_ICE_THICKNESS_IMBALANCE_LONG_NAME "Ice thickness imbalance"
#define CFIO_FIELD_LAND_ICE_CALVING_FLUX_LONG_NAME "Land ice calving flux"
#define CFIO_FIELD_SURFACE_MASS_BALANCE_FLUX_LONG_NAME "Surface mass balance flux"
#define CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_LONG_NAME "Basal mass balance of ice sheet"
#define CFIO_FIELD_GROUNDED_ICE_BASAL_MASS_BALANCE_FLUX_LONG_NAME "Basal mass balance of grounded ice sheet"
#define CFIO_FIELD_FLOATING_ICE_BASAL_MASS_BALANCE_FLUX_LONG_NAME "Basal mass balance of floating ice shelf"
#define CFIO_FIELD_ICE_FRONT_CALVING_AND_MELT_LONG_NAME "Ice frontcalving and melt flux"
#define CFIO_FIELD_GEOTHERMAL_FLUX_LONG_NAME "Bedrock geothermal heat flux"

void AmrIce::initCFData()
{

  CH_TIME("AmrIce::initCFData");


  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::initCFData" << endl;
    }

  
  ParmParse pp (CFIO_PARM_PARSE_PREFIX);
  m_cf_level = 0;
  pp.query("level",m_cf_level);
  CH_assert(m_cf_level <= m_finest_level);
  
  // need to chose a dbl if cf_level > 0, which will usually be the case
  ProblemDomain pd = grids(m_cf_level).physDomain();
  Vector<int> proc;
  Vector<Box> boxes;
  domainSplit(pd, boxes, CFIO_MAX_BOX_SIZE , CFIO_BLOCK_FACTOR);
  LoadBalance(proc, boxes);
  DisjointBoxLayout dbl(boxes, proc, pd);
  
  //need to work out the component count; set up functions that evaluate the CF fields at the same time
  ParmParse ppf(CFIO_PARM_PARSE_FIELD_PREFIX);
  auto test = [](const char* name, const ParmParse& a_pp, bool a_default)
    {
      a_pp.query(name, a_default);
      pout() << CFIO_PARM_PARSE_FIELD_PREFIX << "." << name << " =  " << a_default << endl;
      return a_default;
    };

									      

/* 
   Fields from Table A3 from Nowicki et al, Geosci. Model Dev., 9, 4521–4545, 2016 
   -----------------------------------------------------------------------------
   
   Tier 1 (Mandatory)
   ------------------
   surface_altitude (m), land_ice_thickness(m), bedrock_altitude(m), 
   
   land_ice_surface_x_velocity ( m/a), land_ice_surface_y_velocity ( m/a),
   land_ice_basal_x_velocity ( m/a), land_ice_basal_y_velocity ( m/a),
   
   temperature_at_top_of_ice_sheet_model (K), temperature_at_base_of_ice_sheet_model (K)
   
   land_ice_area_fraction,  grounded_ice_sheet_area_fraction, floating_ice_sheet_area_fraction
   
   also, 'surface mass balance and its components' (kg m^-2 s^-1). 

*/


  // ice thickness
  if (test(CFIO_FIELD_LAND_ICE_THICKNESS_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF thickness" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_THICKNESS_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_THICKNESS_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_THICKNESS_LONG_NAME);
      m_uniform_cf_units.push_back("m");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return &(geometry(a_lev)->getH());});
      m_cf_field_interval.push_back( Interval(0,0));
    }

  // ice surface elevation
  if (test(CFIO_FIELD_SURFACE_ALTITUDE_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF surface" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_SURFACE_ALTITUDE_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_SURFACE_ALTITUDE_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_SURFACE_ALTITUDE_LONG_NAME);
      m_uniform_cf_units.push_back("m");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return &geometry(a_lev)->getSurfaceHeight();});
      m_cf_field_interval.push_back( Interval(0,0));
    }

  // bedrock elevation
  if (test(CFIO_FIELD_BEDROCK_ALTITUDE_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF bedrock" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_BEDROCK_ALTITUDE_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_BEDROCK_ALTITUDE_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_BEDROCK_ALTITUDE_LONG_NAME);
      m_uniform_cf_units.push_back("m");
      m_cf_field_function.push_back( [this](int a_lev, LevelData<FArrayBox>& a_buf){return &(geometry(a_lev)->getTopography());} );
      m_cf_field_interval.push_back( Interval(0,0));
    }
  
  // x and y components of basal velocity. Ask for x, get y too
  if (test(CFIO_FIELD_LAND_ICE_BASAL_VELOCITY_SHORT_NAME,ppf,false))
    { 
      pout() << "AmrIceIO:: CF basal veolcity" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_BASAL_X_VELOCITY_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_BASAL_X_VELOCITY_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_BASAL_X_VELOCITY_LONG_NAME);
      m_uniform_cf_units.push_back("m yr^-1");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return &(*m_layerSFaceXYVel[a_lev]);} );
      m_cf_field_interval.push_back( Interval(m_nLayers,m_nLayers+1)); // when x is requested, supply x and y
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_BASAL_Y_VELOCITY_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_BASAL_Y_VELOCITY_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_BASAL_Y_VELOCITY_LONG_NAME);
      m_uniform_cf_units.push_back("m yr^-1");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return  &(*m_layerSFaceXYVel[a_lev]) ;} );
      m_cf_field_interval.push_back( Interval(m_nLayers,m_nLayers+1)); // when y is requested, supply x and y. Shouldn't happen
    } 

  // x and y components of surface velocity. Ask for x, get y too
  if (test(CFIO_FIELD_LAND_ICE_SURFACE_VELOCITY_SHORT_NAME,ppf,false))
    { 
      pout() << "AmrIceIO:: CF surface velocity" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_SURFACE_X_VELOCITY_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_SURFACE_X_VELOCITY_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_SURFACE_X_VELOCITY_LONG_NAME);
      m_uniform_cf_units.push_back("m yr^-1");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return  &(*m_layerSFaceXYVel[a_lev]);} );
      m_cf_field_interval.push_back( Interval(0,1)); // when x is requested, supply x and y
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_SURFACE_Y_VELOCITY_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_SURFACE_Y_VELOCITY_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_SURFACE_Y_VELOCITY_LONG_NAME);
      m_uniform_cf_units.push_back("m yr^-1");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return  &(*m_layerSFaceXYVel[a_lev]) ;} );
      m_cf_field_interval.push_back( Interval(0,1)); // when y is requested, supply x and y. Shouldn't happen
    }

  // temperature at the top of the ice sheet model
  if (test(CFIO_FIELD_TEMPERATURE_ICE_TOP_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF top temperature" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_TOP_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_TOP_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_TOP_LONG_NAME);
      m_uniform_cf_units.push_back("K");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return  &(*m_sTemperature[a_lev]);});
      m_cf_field_interval.push_back( Interval(0,0));
    }
  
  // temperature at the base of the ice sheet model
  if (test(CFIO_FIELD_TEMPERATURE_ICE_BASE_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF basal temperature" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_BASE_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_BASE_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_BASE_LONG_NAME);
      m_uniform_cf_units.push_back("K");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf){return  &(*m_bTemperature[a_lev]);});
      m_cf_field_interval.push_back( Interval(0,0));
    }

  // temperature at the base of the ice sheet model - grounded ice
  if (test(CFIO_FIELD_TEMPERATURE_GROUNDED_ICE_BASE_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF basal temperature for grounded ice" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_TEMPERATURE_GROUNDED_ICE_BASE_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_BASE_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_TEMPERATURE_GROUNDED_ICE_BASE_LONG_NAME);
      m_uniform_cf_units.push_back("K");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf)
				    {
				      computeAreaFraction(a_buf, GROUNDEDMASKVAL, a_lev);
				      const LevelData<FArrayBox>& t = (*m_bTemperature[a_lev]);
				      for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
					{
					  a_buf[dit] *= t[dit];
					}
				      return &a_buf;
				    } );
      m_cf_field_interval.push_back( Interval(0,0));

    }

  // temperature at the base of the ice sheet model - floating ice shelf
  if (test(CFIO_FIELD_TEMPERATURE_FLOATING_ICE_BASE_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF basal temperature for ice shelves" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_TEMPERATURE_FLOATING_ICE_BASE_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_TEMPERATURE_ICE_BASE_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_TEMPERATURE_FLOATING_ICE_BASE_LONG_NAME);
      m_uniform_cf_units.push_back("K");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf)
				    {
				      computeAreaFraction(a_buf, FLOATINGMASKVAL, a_lev);
				      const LevelData<FArrayBox>& t = (*m_bTemperature[a_lev]);
				      for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
					{
					  a_buf[dit] *= t[dit];
					}
				      return &a_buf;
				    } );
      m_cf_field_interval.push_back( Interval(0,0));
    }

  
  // // ice area fraction
  if (test(CFIO_FIELD_LAND_ICE_AREA_FRACTION_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF ice fraction" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_AREA_FRACTION_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_AREA_FRACTION_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_AREA_FRACTION_LONG_NAME);
      m_uniform_cf_units.push_back("1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   const DisjointBoxLayout& levelGrids = m_amrGrids[a_lev];
	   LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	       a_buf[dit].setVal(0.0);
	       Real a = 1.0;
	       int b = GROUNDEDMASKVAL;
	       FORT_SETONMASK(CHF_FRA1(a_buf[dit],0),
			      CHF_CONST_FIA1(mask,0),
			      CHF_CONST_INT(b),CHF_CONST_REAL(a),
			      CHF_BOX(levelGrids[dit]));
	       b = FLOATINGMASKVAL;
	       FORT_SETONMASK(CHF_FRA1(a_buf[dit],0),
			      CHF_CONST_FIA1(mask,0),
			      CHF_CONST_INT(b),CHF_CONST_REAL(a),
			      CHF_BOX(levelGrids[dit]));
	      
	     }
	   return &a_buf;
	 } );
      m_cf_field_interval.push_back( Interval(0,0));
    }

  //grounded ice area fraction
  if (test(CFIO_FIELD_GROUNDED_ICE_AREA_FRACTION_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF grounded ice fraction" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_GROUNDED_ICE_AREA_FRACTION_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_GROUNDED_ICE_AREA_FRACTION_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_GROUNDED_ICE_AREA_FRACTION_LONG_NAME);
      m_uniform_cf_units.push_back("1");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf)
				    {
				      computeAreaFraction(a_buf, GROUNDEDMASKVAL, a_lev);
				      return &a_buf;
				    } );
      m_cf_field_interval.push_back( Interval(0,0));
    }
  
  //floating ice area fraction
  if (test(CFIO_FIELD_FLOATING_ICE_AREA_FRACTION_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF floating ice fraction" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_FLOATING_ICE_AREA_FRACTION_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_FLOATING_ICE_AREA_FRACTION_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_FLOATING_ICE_AREA_FRACTION_LONG_NAME);
      m_uniform_cf_units.push_back("1");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf)
       				    {computeAreaFraction(a_buf, FLOATINGMASKVAL, a_lev);
				      return &a_buf;
				    } );
      m_cf_field_interval.push_back( Interval(0,0));
    }

 // Ice thickness imbalance, dH/dt
  if (test(CFIO_FIELD_LAND_ICE_THICKNESS_IMBALANCE_SHORT_NAME,ppf,false) || m_write_ismip6)
    {
      pout() << "AmrIceIO:: CF ice thickness imbalance" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_THICKNESS_IMBALANCE_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_THICKNESS_IMBALANCE_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_THICKNESS_IMBALANCE_LONG_NAME);
      m_uniform_cf_units.push_back("m yr^-1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   const LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];
	   const LevelData<FArrayBox>& levelH = levelCS.getH();
	   const LevelData<FArrayBox>& levelOldH = (*m_old_thickness[a_lev]);

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       a_buf[dit].copy(levelH[dit]);
	       a_buf[dit].minus(levelOldH[dit]);
	       if (m_dt > 0)
		 {
		   a_buf[dit] /= m_dt;
		 }
	     }
	   return &a_buf;
	 } );
      m_cf_field_interval.push_back( Interval(0,0));
    }
  
 // surface mass balance
  if (test(CFIO_FIELD_SURFACE_MASS_BALANCE_FLUX_SHORT_NAME,ppf,false) || m_write_ismip6)
    {
      pout() << "AmrIceIO:: CF smb" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_SURFACE_MASS_BALANCE_FLUX_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_SURFACE_MASS_BALANCE_FLUX_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_SURFACE_MASS_BALANCE_FLUX_LONG_NAME);
      m_uniform_cf_units.push_back("kg m^-2 yr^-1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];
	   Real rhoi = levelCS.iceDensity();

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       a_buf[dit].copy((*m_surfaceThicknessSource[a_lev])[dit]);
	       a_buf[dit] *= rhoi;
	     }
	   return &a_buf;
	} );
      m_cf_field_interval.push_back( Interval(0,0));
    }
 

/* 
   Fields from Table A3 from Nowicki et al, Geosci. Model Dev., 9, 4521–4545, 2016 
   -----------------------------------------------------------------------------

   Tier 2 (Desirable)
   ------------------

   land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting (kg m^-2 s_2)
   land_ice_basal_specific_mass_balance_flux (kg m^-2 s_2)

   land_ice_surface_upward_velocity ( m/a),land_ice_basal_upward_velocity (m/a)
   
   land_ice_vertical_mean_x_velocity (m/a), land_ice_vertical_mean_y_velocity (m/a)
   
   surface_snow_area_fraction

*/

 // basal mass balance
  if (test(CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF bmb" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_LONG_NAME);
      m_uniform_cf_units.push_back("kg m^-2 yr^-1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];
	   Real rhoi = levelCS.iceDensity();

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       a_buf[dit].copy((*m_basalThicknessSource[a_lev])[dit]);
	       a_buf[dit] *= rhoi;
	     }
	   return &a_buf;
	} );

      m_cf_field_interval.push_back( Interval(0,0));
    }

 // basal mass balance beneath grounded ice
  if (test(CFIO_FIELD_GROUNDED_ICE_BASAL_MASS_BALANCE_FLUX_SHORT_NAME,ppf,false) || m_write_ismip6)
    {
      pout() << "AmrIceIO:: CF bmb for grounded ice" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_GROUNDED_ICE_BASAL_MASS_BALANCE_FLUX_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_GROUNDED_ICE_BASAL_MASS_BALANCE_FLUX_LONG_NAME);
      m_uniform_cf_units.push_back("kg m^-2 yr^-1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   computeAreaFraction(a_buf, GROUNDEDMASKVAL, a_lev);
	   const LevelData<FArrayBox>& b = (*m_basalThicknessSource[a_lev]);
	   LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];
	   Real rhoi = levelCS.iceDensity();

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       a_buf[dit] *= b[dit];
	       a_buf[dit] *= rhoi;
	     }
	   return &a_buf;
	 } );
      m_cf_field_interval.push_back( Interval(0,0));
    }

 // basal mass balance beneath floating ice
  if (test(CFIO_FIELD_FLOATING_ICE_BASAL_MASS_BALANCE_FLUX_SHORT_NAME,ppf,false) || m_write_ismip6)
    {
      pout() << "AmrIceIO:: CF bmb for floating ice" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_FLOATING_ICE_BASAL_MASS_BALANCE_FLUX_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_BASAL_MASS_BALANCE_FLUX_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_FLOATING_ICE_BASAL_MASS_BALANCE_FLUX_LONG_NAME);
      m_uniform_cf_units.push_back("kg m^-2 yr^-1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   computeAreaFraction(a_buf, FLOATINGMASKVAL, a_lev);
	   const LevelData<FArrayBox>& b = (*m_basalThicknessSource[a_lev]);
	   LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];
	   Real rhoi = levelCS.iceDensity();

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       a_buf[dit] *= b[dit];
	       a_buf[dit] *= rhoi;
	     }
	   return &a_buf;
	 } );
      m_cf_field_interval.push_back( Interval(0,0));
    }

 // vertical front mass balance
  if (test(CFIO_FIELD_ICE_FRONT_CALVING_AND_MELT_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF vertical front mb" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_ICE_FRONT_CALVING_AND_MELT_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_ICE_FRONT_CALVING_AND_MELT_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_ICE_FRONT_CALVING_AND_MELT_LONG_NAME);
      m_uniform_cf_units.push_back("kg m^-2 yr^-1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];
	   Real rhoi = levelCS.iceDensity();

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       a_buf[dit].copy((*m_volumeThicknessSource[a_lev])[dit]);
	       a_buf[dit].plus((*m_calvedIceThickness[a_lev])[dit]);
	       a_buf[dit] *= rhoi;
	     }
	   return &a_buf;
	} );

      m_cf_field_interval.push_back( Interval(0,0));
    }

 /* 
    Fields from Table A3 from Nowicki et al, Geosci. Model Dev., 9, 4521–4545, 2016 
    ----------------------------------------------------------------------------- 

    Tier 3 (Experimental)
    ---------------------
    
    magnitude_of_basal_drag_at_land_ice_base (Pa)
    land_ice_specific_mass_flux_due_to_calving (kg m^-2 s_2)
    upward_geothermal_heat_flux_at_ground_level_in_land_ice (W m^−2)
    
 */
  
 // ice calving
  if (test(CFIO_FIELD_LAND_ICE_CALVING_FLUX_SHORT_NAME,ppf,false) || m_write_ismip6)
    {
      pout() << "AmrIceIO:: CF calving flux" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_LAND_ICE_CALVING_FLUX_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_LAND_ICE_CALVING_FLUX_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_LAND_ICE_CALVING_FLUX_LONG_NAME);
      m_uniform_cf_units.push_back("kg m^-1 yr^-1");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   LevelSigmaCS& levelCS = *m_vect_coordSys[a_lev];
	   Real rhoi = levelCS.iceDensity();

	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {
	       a_buf[dit].copy((*m_calvedIceThickness[a_lev])[dit]);
	       a_buf[dit] *= rhoi;
	     }
	   return &a_buf;
	} );

      m_cf_field_interval.push_back( Interval(0,0));
    }


  if (test(CFIO_FIELD_MAGNITUDE_BASAL_DRAG_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF basal drag" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_MAGNITUDE_BASAL_DRAG_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_MAGNITUDE_BASAL_DRAG_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_MAGNITUDE_BASAL_DRAG_LONG_NAME);
      m_uniform_cf_units.push_back("Pa");
      m_cf_field_function.push_back
	([this](int a_lev, LevelData<FArrayBox>& a_buf)
	 {
	   for (DataIterator dit=a_buf.dataIterator(); dit.ok(); ++dit)
	     {

	       const FArrayBox& dragC = (*dragCoefficient(a_lev))[dit];
	       const FArrayBox& vel = (*m_velocity[a_lev])[dit];

	       for (BoxIterator bit(dragC.box()); bit.ok(); ++bit)
		 {
		   const IntVect& iv = bit();
		   const Real& u = vel(iv,0);
		   const Real& v = vel(iv,1);
		   a_buf[dit](iv) = dragC(iv) * std::sqrt(u*u + v*v);
		 }
	     }
	   return &a_buf;
	 } );

      m_cf_field_interval.push_back( Interval(0,0));
    }

  if (test(CFIO_FIELD_GEOTHERMAL_FLUX_SHORT_NAME,ppf,false))
    {
      pout() << "AmrIceIO:: CF geothermal flux" << endl;
      m_uniform_cf_data_name.push_back(CFIO_FIELD_GEOTHERMAL_FLUX_SHORT_NAME);
      m_uniform_cf_standard_name.push_back(CFIO_FIELD_GEOTHERMAL_FLUX_CF_NAME);
      m_uniform_cf_long_name.push_back(CFIO_FIELD_GEOTHERMAL_FLUX_LONG_NAME);
      m_uniform_cf_units.push_back("W m^-2 ?");
      m_cf_field_function.push_back([this](int a_lev, LevelData<FArrayBox>& a_buf)
				    {
				      (*this).basalHeatBoundaryData().evaluate(a_buf, *this, a_lev, m_dt);
				      return &a_buf;
				    } );
      m_cf_field_interval.push_back( Interval(0,0));
    }

  int n_comp =  m_uniform_cf_data_name.size();
  

  if (n_comp > 0)
    {
      m_uniform_cf_data.define(dbl, n_comp, IntVect::Zero);
    }
  else
    {
      pout() << " no " << CFIO_PARM_PARSE_FIELD_PREFIX << ".<field> entries specified or write_ismip6 " << endl;
      CH_assert(m_uniform_cf_data.nComp() > 0);
    }
}
  


void AmrIce::accumulateCFData(Real a_dt, bool a_reset)
{
  CH_TIME("AmrIce::accumulateCFData");

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::accumulateCFData" << endl;
    }

  bool reset = a_reset || !m_uniform_cf_data.isDefined();

  if (! m_uniform_cf_data.isDefined())
    {
      initCFData();
    }

  if (reset)
    {
      m_cf_start_time = time();
      m_cf_dt = 0.0;
      for (DataIterator dit(m_uniform_cf_data.dataIterator()); dit.ok(); ++dit)
	{
	  m_uniform_cf_data[dit].setVal(0.0);
	}
    }

  if (a_reset)
    {
      m_cf_domain_diagnostic_data.reset();
    }
  else
    {
      //record diagnostic data
      m_cf_domain_diagnostic_data.record(*this);
    }
 
  // put all the fields we want to average into a LevelData;
  LevelData<FArrayBox> current(m_uniform_cf_data.disjointBoxLayout(),
			       m_uniform_cf_data.nComp(),
			       IntVect::Unit);

  // I'm so excited: since Chombo now seems to require C++11, I can use lambda expressions.
  // fun below is not a very nice one, admittedly, since it has side effects.   

  // flatten one or components of the data returned by the function a_ld to
  // to the LevelData<FArrayBox> current, starting at component a_comp.
  // return the index of the next available component of current.
  //auto fun = [&current, this](int a_comp,
  //			      std::function<const LevelData<FArrayBox>*(int)> a_ld,
  //			      Interval a_interval = Interval(0,0))
    
  auto fun = [&current, this](int a_comp,
  			      std::function<const LevelData<FArrayBox>*(int,LevelData<FArrayBox>&) > a_ld,
  			      Interval a_interval = Interval(0,0))


    {
      Vector<LevelData<FArrayBox> *> tmp(m_finest_level+1);
      Vector<RealVect> dx_src(m_finest_level + 1);
      for (int lev = 0; lev <= m_finest_level; lev++)
	{
	  dx_src[lev] = dx(lev);
	  // creating new (empty) LevelData<FArrayBox> and aliasing is safer than the obvious methods
	  tmp[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],a_interval.size(),IntVect::Unit);
	  //only minor evil here, since there will be no writes to ld
	  LevelData<FArrayBox>* ld = const_cast<LevelData<FArrayBox>*>( a_ld(lev, *tmp[lev] ) );
	  if (tmp[lev] != ld) aliasLevelData(*tmp[lev],  ld , a_interval);
	}

      int dcomp = a_interval.size() - 1;
      LevelData<FArrayBox> alias;
      aliasLevelData(alias, &current, Interval(a_comp,a_comp + dcomp));
      flattenCellData(alias, dx(m_cf_level), tmp, dx_src, false);
      
      for (int lev = 0; lev <= m_finest_level; lev++)
	delete tmp[lev];

      return a_comp + a_interval.size();
    };

  // evaluate all the individual fields
  int comp = 0;
  while (comp < m_uniform_cf_data_name.size())
    {
     comp =  fun(comp, m_cf_field_function[comp], m_cf_field_interval[comp]);
    }

  // update the integral with respect to time
  Real w = std::max(a_dt, 1.0e-10);
  m_cf_dt += w;
  for (DataIterator dit(m_uniform_cf_data.dataIterator()); dit.ok(); ++dit)
    {
      current[dit] *= w;
      m_uniform_cf_data[dit] += current[dit];
    }

}

void AmrIce::flushCFData()
{

  CH_TIME("AmrIce::flushCFData");

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::flushCFData" << endl;
    }

  // do we end up in this state after reading a checkpoint?
  if (! m_uniform_cf_data.isDefined())
    {
      if (s_verbosity > 3) 
	{ 
	  pout() << "AmrIce::flushCFData no data to flush" << endl;
	}
      return;
      //initCFData();
    }

  // time integral -> time mean.
  Real cf_end_time = time(); 
  for (DataIterator dit(m_uniform_cf_data.dataIterator()); dit.ok(); ++dit)
    {
      m_uniform_cf_data[dit] /= m_cf_dt;
    }
 

  std::string fs;
  if (m_output_file_numbering ==  time_yyyymmdd_360)
    {
      fs.assign("%s%08ld.CF.%dd.hdf5");
    }
  else
    {
      fs.assign("%s%06ld.CF.%dd.hdf5");
    }
  char* iter_str = new char[m_plot_prefix.size() + 64 + fs.size()];
  sprintf(iter_str, fs.c_str(), m_plot_prefix.c_str(), outputNumbering() , SpaceDim );
  string filename(iter_str);
  delete iter_str;

  
  const DisjointBoxLayout& dbl = m_uniform_cf_data.disjointBoxLayout();

  HDF5Handle handle(filename.c_str(), HDF5Handle::CREATE);
  Vector<DisjointBoxLayout> grids(1,dbl);
  Vector<LevelData<FArrayBox>* > cf_data(1,&m_uniform_cf_data);
  Real dx = m_amrDx[m_cf_level];
  Vector<int> ref_ratio(1,2);

  this->writeAMRHierarchyHDF5(handle,
			      grids,cf_data,
			      m_uniform_cf_data_name,
			      dbl.physDomain().domainBox(),
			      dx,m_dt,
			      cf_end_time,
			      ref_ratio,1);

  std:: string group = handle.getGroup();
  for (int i = 0; i < m_uniform_cf_data_name.size(); ++i)
    { 
      HDF5HeaderData attributeInfo;

      attributeInfo.m_string["Short name"] = m_uniform_cf_data_name[i];
      attributeInfo.m_string["Long name"] = m_uniform_cf_long_name[i];
      attributeInfo.m_string["Units"] = m_uniform_cf_units[i];
      attributeInfo.m_string["Standard name"] = m_uniform_cf_standard_name[i];
      handle.setGroup(group + "/" + m_uniform_cf_data_name[i] + "_attribute");
      attributeInfo.writeToFile(handle);
      handle.setGroup(group);
    }

  // //   HDF5Handle handle(filename, HDF5Handle::mode::OPEN_RDWR);
  m_cf_domain_diagnostic_data.write(handle);

  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::flushCFData close handle" << endl;
    }
  
  handle.close();
  
  //reset CF data following a write
  if (s_verbosity > 3) 
    { 
      pout() << "AmrIce::flushCFData reset" << endl;
    }
  accumulateCFData(0.0, true);
}

#endif


#include "NamespaceFooter.H"
