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
#include "FASIceSolverI.H"
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


#if BISICLES_Z == BISICLES_LAYERED
void AmrIce::updateInternalEnergy(Vector<LevelData<FluxBox>* >& a_layerEH_half, 
				  Vector<LevelData<FluxBox>* >& a_layerH_half,
				  const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel,
				  const Vector<LevelData<FArrayBox>* >& a_layerSFaceXYVel, 
				  const Real a_dt, const Real a_time,
				  Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysNew,
				  Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSysOld,
				  const Vector<LevelData<FArrayBox>*>& a_surfaceThicknessSource,
				  const Vector<LevelData<FArrayBox>*>& a_basalThicknessSource,
				  const Vector<LevelData<FArrayBox>*>& a_volumeThicknessSource)
{

CH_TIME("AmrIce::updateInternalEnergy");

  
  //update the internalEnergy fields, 2D case
  Vector<LevelData<FluxBox>* > vectLayerFluxes(m_finest_level+1, NULL);
  Vector<LevelData<FluxBox>* > vectLayerThicknessFluxes(m_finest_level+1, NULL);
  //Vector<LevelData<FArrayBox>* > vectUSigma(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > vectDivUHxy(m_finest_level+1, NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      LevelData<FluxBox>& levelXYFaceXYVel = *a_layerXYFaceXYVel[lev];
      LevelData<FluxBox>& levelFaceEH = *a_layerEH_half[lev];
      LevelData<FluxBox>& levelFaceH = *a_layerH_half[lev];
      IntVect ghostVect = IntVect::Unit;//CoarseAverageFace requires a ghost cell

      //vectUSigma[lev] = new LevelData<FArrayBox>
      //	(m_amrGrids[lev], m_nLayers + 1 , IntVect::Zero);
       
      vectDivUHxy[lev] = new LevelData<FArrayBox>
	(m_amrGrids[lev], m_nLayers + 1 , ghostVect);
     
      vectLayerFluxes[lev] = new LevelData<FluxBox>
	(m_amrGrids[lev], levelXYFaceXYVel.nComp() , ghostVect);

      vectLayerThicknessFluxes[lev] = new LevelData<FluxBox>
	(m_amrGrids[lev], levelXYFaceXYVel.nComp() , ghostVect);

      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  FluxBox& faceVel = levelXYFaceXYVel[dit];
	  FluxBox& faceEH = levelFaceEH[dit];
	  FluxBox& faceH = levelFaceH[dit];
	  FluxBox& flux = (*vectLayerFluxes[lev])[dit];
	  FluxBox& thicknessFlux = (*vectLayerThicknessFluxes[lev])[dit];

	  const Box& gridBox = levelGrids[dit];
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      Box faceBox(gridBox);
	      faceBox.surroundingNodes(dir);
	      flux[dir].copy(faceEH[dir], faceBox);
	      flux[dir].mult(faceVel[dir], faceBox, 0, 0, faceVel[dir].nComp());

	      

	      thicknessFlux[dir].copy(faceH[dir],faceBox);
	      thicknessFlux[dir].mult(faceVel[dir], faceBox, 0, 0, faceVel[dir].nComp());
	        
	      // CH_assert(flux[dir].norm(faceBox,0,0,flux[dir].nComp()) < HUGE_NORM);
	      // CH_assert(thicknessFlux[dir].norm(faceBox,0,0,thicknessFlux[dir].nComp()) < HUGE_NORM);
		
	    }
	}
    }
  // average fine fluxes down to coarse levels
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],
				     vectLayerFluxes[lev]->nComp(), m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*vectLayerFluxes[lev-1], *vectLayerFluxes[lev]);
      faceAverager.averageToCoarse(*vectLayerThicknessFluxes[lev-1], *vectLayerThicknessFluxes[lev]);
    }
 

     
  //cross-layer velocity u^sigma
  for (int lev=0; lev <= m_finest_level; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
     
      LevelSigmaCS& levelCoordsNew = *(a_coordSysNew[lev]);
      LevelSigmaCS& levelCoordsOld = *(a_coordSysOld[lev]);
      LevelData<FArrayBox> dHdt(levelGrids,1,IntVect::Zero);
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  const FArrayBox& oldH = levelCoordsOld.getH()[dit];
	  const FArrayBox& newH = levelCoordsNew.getH()[dit];
	  dHdt[dit].copy(newH);
	  dHdt[dit] -= oldH;
	  dHdt[dit] *= 1.0/a_dt;
	}
      
      IceUtility::computeSigmaVelocity(*m_layerSFaceSVel[lev],
				       *vectLayerThicknessFluxes[lev],
				       *a_layerSFaceXYVel[lev],
				       dHdt,
				       m_amrGrids[lev],
				       *a_surfaceThicknessSource[lev],
				       *a_basalThicknessSource[lev],
				       levelCoordsNew.getDSigma(),
				       levelCoordsNew.dx(),
				       a_dt);
    }


  // compute rhs =a_dt *(H*dissipation - div(u H T)) and update solution
  for (int lev=0; lev <= m_finest_level; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFlux = *vectLayerFluxes[lev];
      LevelSigmaCS& levelCoordsOld = *(a_coordSysOld[lev]);
      LevelSigmaCS& levelCoordsNew = *(a_coordSysNew[lev]);
      const Vector<Real>& dSigma = levelCoordsNew.getDSigma();
      //caculate dissipation due to internal stresses
      LevelData<FArrayBox> dissipation(levelGrids,m_nLayers,IntVect::Zero);
      {
	LevelData<FArrayBox>* crseVelPtr = NULL;
	int nRefCrse = -1;
	if (lev > 0)
	  {
	    crseVelPtr = m_velocity[lev-1];
	    nRefCrse = m_refinement_ratios[lev-1];
	  }

	m_velocity[lev]->exchange();

	m_constitutiveRelation->computeDissipation
	  (dissipation,*m_velocity[lev],  crseVelPtr,
	   nRefCrse, *m_A[lev],
	   levelCoordsOld , m_amrDomains[lev],  IntVect::Zero);
      }
      
      LevelData<FArrayBox>& surfaceHeatFlux = *m_sHeatFlux[lev];
      if (surfaceHeatBoundaryDirichlett())
	{
	  surfaceHeatBoundaryData().evaluate(*m_sInternalEnergy[lev], *this, lev, a_dt);
	  if (surfaceHeatBoundaryTemperature())
	    {
	      //convert surface temperature to internal energy
	      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
		{
		  FArrayBox& E = (*m_sInternalEnergy[lev])[dit];
		  FArrayBox T(E.box(),1);
		  T.copy(E);
		  FArrayBox W(E.box(),1);
		  W.setVal(0.0);
		  IceThermodynamics::composeInternalEnergy(E, T, W, E.box());
		}
	    }
	}
      else
	{
	  m_surfaceHeatBoundaryDataPtr->evaluate(surfaceHeatFlux, *this, lev, a_dt);
	}
      
      LevelData<FArrayBox>& basalHeatFlux = *m_bHeatFlux[lev];
      basalHeatBoundaryData().evaluate(basalHeatFlux, *this, lev, a_dt);


      /// compute (spatially variable) till water drainage rate
      LevelData<FArrayBox> tillWaterDrainFactor(levelGrids, 1, IntVect::Unit);
      {
	// on the fly creation seems lazy, but we don't need this anywhere else
	SurfaceFlux* ptr = SurfaceFlux::parse("tillWaterDrainFactor");
	if (ptr == NULL)
	  {
	    // default: use the constant rate set in IceThermodynamics
	    ptr = new constantFlux(IceThermodynamics::m_till_water_drain_factor);
	  }
	ptr->evaluate(tillWaterDrainFactor, *this, lev, a_dt);
	delete ptr;
      }
	  
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  const Box& box = levelGrids[dit];
	  const FArrayBox& oldH = levelCoordsOld.getH()[dit];
	  const FArrayBox& newH = levelCoordsNew.getH()[dit];
	  FArrayBox& E = (*m_internalEnergy[lev])[dit];
	  FArrayBox& Wt = (*m_tillWaterDepth[lev])[dit];
	  FArrayBox& sT = (*m_sInternalEnergy[lev])[dit];	
	  FArrayBox& bT = (*m_bInternalEnergy[lev])[dit];
	  

	  // first, do the ordinary fluxes : if we just had
	  // horizontal advection and grad(H) = grad(S) = 0., 
	  // this would be the lot
	       
	  FArrayBox rhs(box, m_nLayers);
	  rhs.setVal(0.0);
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      Real dx = levelCoordsOld.dx()[dir];              
	   
	      FORT_DIVERGENCE(CHF_CONST_FRA(levelFlux[dit][dir]),
			      CHF_FRA(rhs),
			      CHF_BOX(box),
			      CHF_CONST_REAL(dx),
			      CHF_INT(dir));
	     

	    }
	  for (int layer = 0; layer < dissipation.nComp(); ++layer)
	    {
	      dissipation[dit].mult(newH,0,layer,1);
	    }
	   dissipation[dit] /= levelCoordsNew.iceDensity();
	  //add a layer heat source proportional to the volume thickness source
	  for (int layer = 0; layer < dissipation.nComp(); ++layer)
	    {
	      FArrayBox dE(box,1);
	      dE.copy(E, layer, 0, 1);
	      dE *= (*a_volumeThicknessSource[lev])[dit];
	      //dE *= dSigma[layer];
	      dissipation[dit].plus(dE,0,layer,1); 
	    }
	      

	  
	  dissipation[dit] /= levelCoordsOld.iceDensity();
	  rhs -= dissipation[dit]; 
	  rhs *= -a_dt;

	  //compute heat flux across base due to basal dissipation
	  FArrayBox basalDissipation(rhs.box(),1);
	  m_basalFrictionRelation->computeDissipation
	    (basalDissipation , (*m_velocity[lev])[dit] , (*m_velBasalC[lev])[dit], 1.0, 
	     levelCoordsOld , dit ,lev, rhs.box());
	  

	  //add to user set (e.g geothermal) heat flux
	  basalHeatFlux[dit] += basalDissipation;

	  //zero heat flux outside grounded ice
	  for (BoxIterator bit(rhs.box());bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (levelCoordsOld.getFloatingMask()[dit](iv) != GROUNDEDMASKVAL)
		{
		  basalHeatFlux[dit](iv) = 0.0;
		}
	    }
	  
	  //basalHeatFlux[dit] /= levelCoordsNew.iceDensity()); // scale conversion
	  FArrayBox scaledBasalHeatFlux(basalHeatFlux[dit].box(),basalHeatFlux[dit].nComp());
	  scaledBasalHeatFlux.copy(basalHeatFlux[dit]);
	  scaledBasalHeatFlux /=  levelCoordsNew.iceDensity();

	  //surfaceHeatFlux[dit] /= (levelCoordsNew.iceDensity()); // scale conversion
	  FArrayBox scaledSurfaceHeatFlux(surfaceHeatFlux[dit].box(),surfaceHeatFlux[dit].nComp());
	  scaledSurfaceHeatFlux.copy(surfaceHeatFlux[dit]);
	  scaledSurfaceHeatFlux /= levelCoordsNew.iceDensity();

	  //solve H(t+dt)E(t+dt) + vertical transport terms = H(t)E(t) - rhs(t+/dt)
	  //with either a Dirichlett or flux boundary condition at the upper surface and a flux condition at base
	 
	  Real halftime = time() + 0.5*a_dt;
	  int nLayers = m_nLayers;
	  const Real& rhoi = levelCoordsNew.iceDensity();
	  const Real& rhoo = levelCoordsNew.waterDensity();
	  const Real& gravity = levelCoordsNew.gravity();
	 
	  int surfaceTempDirichlett = surfaceHeatBoundaryDirichlett()?1:0;
	  const FArrayBox& uSigma = (*m_layerSFaceSVel[lev])[dit];

	  IceThermodynamics::timestep(E, Wt, sT, bT,
				      scaledSurfaceHeatFlux,
				      scaledBasalHeatFlux,
				      tillWaterDrainFactor[dit],
				      levelCoordsOld.getFloatingMask()[dit],
				      levelCoordsNew.getFloatingMask()[dit],
				      oldH,
				      newH,
				      uSigma,
				      rhs,
				      levelCoordsOld.getFaceSigma(),
				      dSigma,
				      halftime,
				      a_dt,
				      nLayers,
				      surfaceHeatBoundaryDirichlett(),
				      box);
	  
	  scaledBasalHeatFlux *= (levelCoordsNew.iceDensity());
	  basalHeatFlux[dit].copy(scaledBasalHeatFlux);
	  scaledSurfaceHeatFlux *= (levelCoordsNew.iceDensity());
	  surfaceHeatFlux[dit].copy(scaledSurfaceHeatFlux);	    
	} // end update internal energy loop over grids
    } // end update internal energy loop over levels


  // horizontal conduction / smoothing. Occasional...
  
  ParmParse pp("Thermodynamics");

  int skip = -1;
  pp.query("smooth_interval",skip); 


  Real lambda = std::sqrt( IceThermodynamics::IceConductivity()/IceThermodynamics::IceHeatCapacity() );
  if (skip > 0 && int(m_time) % skip == 0)
    {
      pp.query("length_scale",  lambda);
      helmholtzSolve(m_internalEnergy, 1.0, lambda*lambda * double(skip));      
    }
  
  pp.query("till_water_smooth_interval",skip); 
  if (skip > 0 && int(m_time) % skip == 0)
    {
      pp.query("till_water_length_scale",  lambda);
      helmholtzSolve(m_tillWaterDepth, 1.0, lambda * lambda * double(skip));      
    }
  
  
  //coarse average from finer levels & exchange
  for (int lev = m_finest_level; lev >= 0 ; --lev)
    {
      if (lev > 0)
	{
	  CoarseAverage avN(m_amrGrids[lev],
			    m_amrGrids[lev-1],
			    m_internalEnergy[lev]->nComp(),
			    m_refinement_ratios[lev-1], 
			    IntVect::Zero);
	  
	  
	  
	  avN.averageToCoarse(*m_internalEnergy[lev-1], *m_internalEnergy[lev]);
	
	  
	  CoarseAverage avOne(m_amrGrids[lev],m_amrGrids[lev-1],
			      1,m_refinement_ratios[lev-1], IntVect::Zero);
	  
	  avOne.averageToCoarse(*m_tillWaterDepth[lev-1], *m_tillWaterDepth[lev]);
	  avOne.averageToCoarse(*m_sInternalEnergy[lev-1], *m_sInternalEnergy[lev]);
	  avOne.averageToCoarse(*m_bInternalEnergy[lev-1], *m_bInternalEnergy[lev]);
	  avOne.averageToCoarse(*m_sHeatFlux[lev-1], *m_sHeatFlux[lev]);
	  avOne.averageToCoarse(*m_bHeatFlux[lev-1], *m_bHeatFlux[lev]);
	}
      
      m_internalEnergy[lev]->exchange();
      m_tillWaterDepth[lev]->exchange();
      m_sInternalEnergy[lev]->exchange();
      m_bInternalEnergy[lev]->exchange();
      m_sHeatFlux[lev]->exchange();
      m_bHeatFlux[lev]->exchange();
    }
  
  for (int lev = 0; lev < vectLayerFluxes.size(); ++lev)
    {
      // if (vectUSigma[lev] != NULL)
      // 	{
      // 	  delete vectUSigma[lev]; vectUSigma[lev] = NULL;
      // 	}

      if (vectDivUHxy[lev] != NULL)
	{
	  delete  vectDivUHxy[lev]; vectDivUHxy[lev] = NULL;
	}

      if (vectLayerFluxes[lev] != NULL)
	{
	  delete vectLayerFluxes[lev];vectLayerFluxes[lev] = NULL;
	}

      if (vectLayerThicknessFluxes[lev] != NULL)
	{
	  delete vectLayerThicknessFluxes[lev];vectLayerThicknessFluxes[lev] = NULL;
	}

    }
  
  //update the temperature since it depends on the internal energy
  updateTemperature();
  
  //finally, A is no longer valid 
  m_A_valid = false;
  //#endif

}
#endif

#if BISICLES_Z == BISICLES_LAYERED
void AmrIce::updateTemperature()
{
  
  //update the temperature (a derived field);
  for (int lev=0; lev < m_temperature.size() ; lev++)
    {
      if (m_temperature[lev])
	{delete m_temperature[lev];m_temperature[lev]=NULL;}
      if (m_bTemperature[lev])
	{delete m_bTemperature[lev];m_bTemperature[lev]=NULL;}
      if (m_sTemperature[lev])
	{ delete m_sTemperature[lev];m_sTemperature[lev]=NULL;}
    }
  
  (m_temperature.resize(m_internalEnergy.size()));
  (m_sTemperature.resize(m_sInternalEnergy.size()));
  (m_bTemperature.resize(m_sInternalEnergy.size()));
  
  for (int lev=0; lev<=m_finest_level; lev++)
    {  
      m_temperature[lev] = new LevelData<FArrayBox>
	(m_internalEnergy[lev]->disjointBoxLayout(),m_internalEnergy[lev]->nComp(), IntVect::Unit );
      m_sTemperature[lev] = new LevelData<FArrayBox>
	(m_sInternalEnergy[lev]->disjointBoxLayout(),m_sInternalEnergy[lev]->nComp(), IntVect::Unit );
      m_bTemperature[lev] = new LevelData<FArrayBox>
	(m_bInternalEnergy[lev]->disjointBoxLayout(),m_sInternalEnergy[lev]->nComp(), IntVect::Unit );

      const LevelSigmaCS& coordSys = *m_vect_coordSys[lev];
      const Vector<Real>& sigma = coordSys.getSigma();
      
      for (DataIterator dit(m_amrGrids[lev]); dit.ok(); ++dit)
	{
	  FArrayBox& T = (*m_temperature[lev]) [dit];
	  FArrayBox& sT = (*m_sTemperature[lev]) [dit];
	  FArrayBox& bT = (*m_bTemperature[lev]) [dit];
	  const FArrayBox& E = (*m_internalEnergy[lev]) [dit];
	  const FArrayBox& sE = (*m_sInternalEnergy[lev]) [dit];
	  const FArrayBox& bE = (*m_bInternalEnergy[lev]) [dit];
	  const Box& box = bT.box();
	  FArrayBox pressure(box,1);
	  FArrayBox w(box,1);
	  
	  //Surface temperature
	  pressure.setVal(0.0);
	  IceThermodynamics::decomposeInternalEnergy(sT, w, sE , pressure, box);
	  
	  //Base temperature  
	  pressure.copy(coordSys.getH()[dit]);
	  pressure *= coordSys.iceDensity() * coordSys.gravity();
	  IceThermodynamics::decomposeInternalEnergy(bT, w, bE , pressure, box);
	  
	  // Bulk temperature
	  for (int layer = 0; layer < sigma.size(); ++layer)
	    {
	      pressure.copy(coordSys.getH()[dit]);
	      pressure *= coordSys.iceDensity() * coordSys.gravity() * sigma[layer];
	      FArrayBox layerE(box,1);
	      layerE.copy(E,layer,0,1); // an alias would be better?
	      FArrayBox layerT(box,1); 
	      IceThermodynamics::decomposeInternalEnergy(layerT, w, layerE , pressure, box);
	      T.copy(layerT,0,layer,1);
	    }
	}
    }
 
}
#endif

#if BISICLES_Z == BISICLES_LAYERED
//compute the face- and layer- centered internal energy (a_layerEH_half)
//and thickness (a_layerH_half) at time a_time + 1/2 * a_dt
void AmrIce::computeInternalEnergyHalf(Vector<LevelData<FluxBox>* >& a_layerEH_half,
				       Vector<LevelData<FluxBox>* >& a_layerH_half,
				       const Vector<LevelData<FluxBox>* >& a_layerXYFaceXYVel, 
				       const Real a_dt, const Real a_time)
{

  CH_TIME("AmrIce::computeInternalEnergyHalf");
  
  //delete and re-create storage for a_layerEH_half and a_layerH_half.
  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {
      
      if (a_layerEH_half[lev] != NULL)
	delete(a_layerEH_half[lev]);

      a_layerEH_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 
						   m_internalEnergy[lev]->nComp(), 
						   IntVect::Unit);
      if (a_layerH_half[lev] != NULL)
	delete(a_layerH_half[lev]);
      
      a_layerH_half[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 
						  m_internalEnergy[lev]->nComp(), 
						  IntVect::Unit);
    }


  //assume the ghost regions of m_internalEnergy are not correct
  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {
      if (lev > 0)
	{
	  PiecewiseLinearFillPatch pwl(m_amrGrids[lev],
				       m_amrGrids[lev-1],
				       m_internalEnergy[lev]->nComp(),
				       m_amrDomains[lev-1],
				       m_refinement_ratios[lev-1],
				       m_internalEnergy[lev]->ghostVect()[0]);
	  pwl.fillInterp(*m_internalEnergy[lev],*m_internalEnergy[lev-1],
			 *m_internalEnergy[lev-1],1.0,0,0,m_internalEnergy[lev]->nComp());
	}
      m_internalEnergy[lev]->exchange();
    }
   
  for (int lev = 0 ; lev <= m_finest_level; lev++)
    {
      //in the 2D case (with poor man's multidim) this
      //is a little pained using AdvectPhysics, but for the time being
      //we need to construct a single component thisLayerEH_Half for each layer, 
      //given a internalEnergy and horizontal velocity and then copy it into 
      //the multicomponent EH_half[lev]
     
      // PatchGodunov object for layer thickness/energy advection
      PatchGodunov patchGodunov;
      {
	int normalPredOrder = 0;
	bool useFourthOrderSlopes = false;
	bool usePrimLimiting = false;
	bool useCharLimiting = false;
	bool useFlattening = false;
	bool useArtificialViscosity = false;
	Real artificialViscosity = 0.0;
	AdvectPhysics advectPhys;
	advectPhys.setPhysIBC(m_internalEnergyIBCPtr);

	patchGodunov.define(m_amrDomains[lev], m_amrDx[lev],
			    &advectPhys, normalPredOrder,
			    useFourthOrderSlopes,usePrimLimiting,
			    useCharLimiting,useFlattening,
			    useArtificialViscosity,artificialViscosity);
	patchGodunov.setCurrentTime(m_time);
      }
      
      AdvectPhysics* advectPhysPtr = dynamic_cast<AdvectPhysics*>(patchGodunov.getGodunovPhysicsPtr());
      if (advectPhysPtr == NULL)
	{
	  MayDay::Error("AmrIce::computeInternalEnergyHalf -- unable to upcast GodunovPhysics to AdvectPhysics");
	}

      const LevelData<FArrayBox>& levelInternalEnergy = *m_internalEnergy[lev]; 
      const LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev]; 
      const LevelData<FluxBox>& levelLayerXYFaceXYVel = *a_layerXYFaceXYVel[lev]; 
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

      for (int layer = 0; layer < m_nLayers; ++layer)
	{
	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      const Box& box = levelInternalEnergy[dit].box(); // grid box plus ghost cells
	      
	      FluxBox layerXYFaceXYVel(box,1);
	      layerXYFaceXYVel.setVal(0.0);
	      for (int dir = 0; dir < SpaceDim; ++dir){
		layerXYFaceXYVel[dir].copy(levelLayerXYFaceXYVel[dit][dir],layer,0,1);
		Box faceBox = levelGrids[dit].surroundingNodes(dir);
		CH_assert(layerXYFaceXYVel[dir].norm(faceBox,0) < HUGE_NORM);
	      }

	      FArrayBox layerCellXYVel(box,SpaceDim);
	      EdgeToCell(layerXYFaceXYVel,layerCellXYVel);

	      //\todo compute bulk heat sources
	      FArrayBox heatSource(levelGrids[dit], 1);
	      heatSource.setVal(0.0);

	      patchGodunov.setCurrentBox(levelGrids[dit]);
	      advectPhysPtr->setVelocities(&layerCellXYVel,&layerXYFaceXYVel);

	      FArrayBox WGdnv(box,1);

	      //HE at half time and cell faces
	      WGdnv.copy(levelInternalEnergy[dit],layer,0,1);
	      WGdnv *= levelOldThickness[dit];
	      Box grownBox = levelGrids[dit];
	      grownBox.grow(1);
	      FluxBox HEhalf(grownBox,1);
	      HEhalf.setVal(0.0);
	      patchGodunov.computeWHalf(HEhalf,
					WGdnv,
					heatSource,
					a_dt,
					levelGrids[dit]);
	      for (int dir = 0; dir < SpaceDim; ++dir)
		{
		  Box faceBox(levelGrids[dit]);
		  faceBox.surroundingNodes(dir);
		  (*a_layerEH_half[lev])[dit][dir].copy(HEhalf[dir],0,layer,1);

		  Real hemax = HEhalf[dir].norm(faceBox,0) + 0.0;
		  if (hemax > HUGE_NORM)
		    {
		      pout() << "AmrIce::computeInternalEnergyHalf max(HE[face]) = "
			     << hemax << std::endl;
		    }
		}
	      
	      //H at half time and cell faces
	      WGdnv.copy(levelOldThickness[dit]);
	      FluxBox Hhalf(grownBox,1);
	      //\todo compute layer thickness sources
	      FArrayBox HSource(levelGrids[dit], 1);
	      HSource.setVal(0.0);
	      patchGodunov.computeWHalf(Hhalf,
					WGdnv,
					HSource,
					a_dt,
					levelGrids[dit]);
	      for (int dir = 0; dir < SpaceDim; ++dir)
		{
		  Box faceBox(levelGrids[dit]);
		  faceBox.surroundingNodes(dir);
		  (*a_layerH_half[lev])[dit][dir].copy(Hhalf[dir],0,layer,1);
		}
	      
	    }
	  
	}
    }
      
  // coarse average new EH-Half to covered regions
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverageFace faceAverager(m_amrGrids[lev],a_layerEH_half[lev]->nComp(), m_refinement_ratios[lev-1]);
      faceAverager.averageToCoarse(*a_layerEH_half[lev-1], *a_layerEH_half[lev]);
      faceAverager.averageToCoarse(*a_layerH_half[lev-1], *a_layerH_half[lev]);
    }

}




#endif

#include "NamespaceFooter.H"
