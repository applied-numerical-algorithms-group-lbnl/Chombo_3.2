#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//===========================================================================
// diagnostics.cpp
// read in a bisicles plotfile (which must) 
// and an optional mask, and write out a bunch of diagnostics about the ice sheet. 
// These are
// 0. time
// 1. Volume of ice
// 2. Volume of ice above flotation
// dh/dt
// SMB
// BMB
// Volume of calved ice
// Calving flux  
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LevelSigmaCS.H"
#include "IceConstants.H"
#include "computeSum.H"
#include "CellToEdge.H"
#include "amrIceF_F.H"
#include "AdvectPhysics.H"
#include "PatchGodunov.H"
#include "SigmaCSF_F.H"
#include "PiecewiseLinearFillPatch.H"
#include "FillFromReference.H"
#include "IceThicknessIBC.H"
#include "LevelDataIBC.H"

void createDEM( Vector<LevelData<FArrayBox>* >& topography,  
		Vector<LevelData<FArrayBox>* >& thickness, 
		Vector<std::string>& name, 
		Vector<LevelData<FArrayBox>* >& data,
		Vector<int>& ratio,
		Vector<Real>& dx,
		Real mcrseDx)
{
  CH_TIME("createDEM");
  int numLevels = data.size();

  for (int lev = 0; lev < numLevels; lev++)
    {

      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  (*topography[lev])[dit].setVal(0.0);
	  (*thickness[lev])[dit].setVal(0.0);
	}

      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "Z_base")
	    {
	      data[lev]->copyTo(Interval(j,j),*topography[lev],Interval(0,0));
	    }
	  else if (name[j] == "thickness")
	    {
	      data[lev]->copyTo(Interval(j,j),*thickness[lev],Interval(0,0));
	    }	  
	}

      

      if (lev > 0)
	{
	  
	  const DisjointBoxLayout& crseGrids = topography[lev-1]->disjointBoxLayout();
	  PiecewiseLinearFillPatch filler(grids , crseGrids, 1, 
					  crseGrids.physDomain(), ratio[lev-1], 2);
	  Real time_interp_coeff = 0.0;
	  filler.fillInterp(*topography[lev],*topography[lev-1] ,*topography[lev-1],
			    time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*thickness[lev],*thickness[lev-1] ,*thickness[lev-1],
			    time_interp_coeff,0, 0, 1);

	}
      thickness[lev] -> exchange();
      topography[lev] -> exchange();
    }
}

void createSigmaCS(Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		   Vector<LevelData<FArrayBox>* >& topography, 
		   Vector<LevelData<FArrayBox>* >& thickness, 
		   Vector<Real>& dx, Vector<int>& ratio,
		   Real iceDensity, Real waterDensity, Real gravity)
{

  int numLevels = topography.size();

  IntVect sigmaCSGhost(2*IntVect::Unit);
       
  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& levelGrids = topography[lev]->disjointBoxLayout();
	   
      coords[lev] = RefCountedPtr<LevelSigmaCS> 
	(new LevelSigmaCS(levelGrids, RealVect::Unit*dx[lev], sigmaCSGhost));
      coords[lev]->setIceDensity(iceDensity);
      coords[lev]->setWaterDensity(waterDensity);
      coords[lev]->setGravity(gravity);
      topography[lev]->copyTo(Interval(0,0),coords[lev]->getTopography(),Interval(0,0));   
      thickness[lev]->copyTo(Interval(0,0),coords[lev]->getH(),Interval(0,0));
      coords[lev]->recomputeGeometry(NULL,0);
	   
    }
       
}

void getThicknessSource(Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
			Vector<LevelData<FArrayBox>* >& basalThicknessSource,
			Vector<LevelData<FArrayBox>* >& deltaThickness,
			Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
			Vector<LevelData<FArrayBox>* >& calvingFlux,
			Vector<LevelData<FArrayBox>* >& melangeThickness,
			Vector<LevelData<FArrayBox>* >& topography,
			Vector<Real>& dx, Vector<int>& ratio, 
			Vector<std::string>& name, 
			Vector<LevelData<FArrayBox>* >& data)

{
  int numLevels = data.size();

  for (int lev = 0; lev < numLevels; lev++)
    {

      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
          (*deltaThickness[lev])[dit].setVal(0.0);
          (*surfaceThicknessSource[lev])[dit].setVal(0.0);
          (*basalThicknessSource[lev])[dit].setVal(0.0);
          (*divergenceThicknessFlux[lev])[dit].setVal(0.0);
          (*melangeThickness[lev])[dit].setVal(0.0);
          (*calvingFlux[lev])[dit].setVal(0.0);
	}

      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "dThickness/dt")
	    {
	      data[lev]->copyTo(Interval(j,j),*deltaThickness[lev],Interval(0,0));
	    }
	  else if (name[j] == "activeSurfaceThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*surfaceThicknessSource[lev],Interval(0,0));
	    }	  
	  else if (name[j] == "activeBasalThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*basalThicknessSource[lev],Interval(0,0));
	    }
	  else if (name[j] == "divergenceThicknessFlux")
	    {
	      data[lev]->copyTo(Interval(j,j),*divergenceThicknessFlux[lev],Interval(0,0));
	    }
	  else if (name[j] == "calvedThicknessSource" || name[j] == "calvingFlux")
	    {
	      data[lev]->copyTo(Interval(j,j),*calvingFlux[lev],Interval(0,0));
	    }
	  else if (name[j] == "calvedIceThickness" || name[j] == "melangeThickness")
	    {
	      data[lev]->copyTo(Interval(j,j),*melangeThickness[lev],Interval(0,0));
	    }
	  
	}
     

      if (lev > 0)
	{
	  
	  const DisjointBoxLayout& crseGrids = topography[lev-1]->disjointBoxLayout();
	  PiecewiseLinearFillPatch filler(grids , crseGrids, 1, 
					  crseGrids.physDomain(), ratio[lev-1], 1);
	  Real time_interp_coeff = 0.0;

	  filler.fillInterp(*deltaThickness[lev],*deltaThickness[lev-1],
			    *deltaThickness[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*surfaceThicknessSource[lev],*surfaceThicknessSource[lev-1],
			    *surfaceThicknessSource[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*basalThicknessSource[lev],*basalThicknessSource[lev-1],
			    *basalThicknessSource[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*divergenceThicknessFlux[lev],*divergenceThicknessFlux[lev-1],
			    *divergenceThicknessFlux[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*calvingFlux[lev],*calvingFlux[lev-1],
			    *calvingFlux[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*melangeThickness[lev],*melangeThickness[lev-1],
			    *melangeThickness[lev-1],time_interp_coeff,0, 0, 1);
	}
      deltaThickness[lev]->exchange();
      surfaceThicknessSource[lev]->exchange();
      basalThicknessSource[lev]->exchange();
      divergenceThicknessFlux[lev]->exchange();
      calvingFlux[lev]->exchange();
      melangeThickness[lev]->exchange();
    }
} 

void computeFlux(Vector<LevelData<FluxBox>* >& fluxOfIce, 
		 const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		 Vector<LevelData<FArrayBox>* >& topography,
		 Vector<LevelData<FArrayBox>* >& thickness,
		 Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
		 Vector<LevelData<FArrayBox>* >& basalThicknessSource, 
		 Vector<Real>& dx, Vector<int>& ratio, 
		 Vector<std::string>& name, 
		 Vector<LevelData<FArrayBox>* >& data)
{

  
  int numLevels = data.size();
  Vector<LevelData<FArrayBox>* > ccVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcThck(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      int comp = 0;
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccVel[lev] = new LevelData<FArrayBox>(grids,SpaceDim,IntVect::Unit);
      fcVel[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      fcThck[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      const LevelSigmaCS& levelCS = *coords[lev];

      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  (*ccVel[lev])[dit].setVal(0.0);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      (*fcVel[lev])[dit][dir].setVal(0.0);
	      (*fluxOfIce[lev])[dit][dir].setVal(0.0);
	    }
	}

      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "xVel")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccVel[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "yVel")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccVel[lev],Interval(1,1));
	      comp++;
	    }
	}
      CH_assert(comp == 2);

      ccVel[lev]->exchange();
      if (lev > 0)
	{
	  const DisjointBoxLayout& crseGrids = topography[lev-1]->disjointBoxLayout();
	  PiecewiseLinearFillPatch velFiller(grids , crseGrids, ccVel[lev]->nComp(), 
					     crseGrids.physDomain(), ratio[lev-1], 1);
	  Real time_interp_coeff = 0.0;
	  velFiller.fillInterp(*ccVel[lev],*ccVel[lev-1] ,*ccVel[lev-1],
			       time_interp_coeff,0, 0, ccVel[lev]->nComp());
	}

      CellToEdge(*ccVel[lev], *fcVel[lev]);

      //modification to fluxes at the margins, that is where mask changes to open sea or land.
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{

	  for (int dir = 0; dir < SpaceDim; ++dir)
	      {
		Box faceBox = grids[dit];
		faceBox.surroundingNodes(dir);
		FArrayBox& faceVel = (*fcVel[lev])[dit][dir];
		Box grownFaceBox = faceBox;
		CH_assert(faceVel.box().contains(grownFaceBox));
		FArrayBox vface(faceBox,1);
		FArrayBox faceVelCopy(faceVel.box(), 1); faceVelCopy.copy(faceVel);
		const FArrayBox& cellVel = (*ccVel[lev])[dit];
		const FArrayBox& usrf = levelCS.getSurfaceHeight()[dit];
		const FArrayBox& thck = (*thickness[lev])[dit];
		const FArrayBox& topg = (*topography[lev])[dit];

		FORT_EXTRAPTOMARGIN(CHF_FRA1(faceVel,0),
                                    CHF_FRA1(vface,0),
				    CHF_CONST_FRA1(faceVelCopy,0),
				    CHF_CONST_FRA1(cellVel,dir),
				    CHF_CONST_FRA1(usrf,0),
				    CHF_CONST_FRA1(topg,0),
				    CHF_CONST_FRA1(thck,0),
				    CHF_CONST_INT(dir),
				    CHF_BOX(faceBox));
	      }
	}

      // face centered thickness from PPM
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  AdvectPhysics advectPhys;
	  RealVect levelDx = RealVect::Unit * dx[lev];
	  RefCountedPtr<LevelData<FArrayBox> > levelThck(thickness[0]); levelThck.neverDelete();
	  RefCountedPtr<LevelData<FArrayBox> > levelTopg(topography[0]); levelTopg.neverDelete();
	  LevelDataIBC thicknessIBC(levelThck,levelTopg,levelDx);

	 
	  FluxBox& fcvel = (*fcVel[lev])[dit];
	  FArrayBox& ccvel = (*ccVel[lev])[dit];
	  advectPhys.setPhysIBC(&thicknessIBC);
	  
	  
	  FluxBox& faceh = (*fcThck[lev])[dit];
	  FArrayBox& cch = (*thickness[lev])[dit];
	  
	  FArrayBox src(grids[dit],1); 
	  src.copy( (*surfaceThicknessSource[lev])[dit]) ;
	  src.plus( (*basalThicknessSource[lev])[dit]) ;
	  
	  Real dt = 1.0e-3;
	  int normalPredOrder = 2;
	  bool useFourthOrderSlopes = true;
	  bool usePrimLimiting = true;
	  bool useCharLimiting = false;
	  bool useFlattening = false;
	  bool useArtificialViscosity = false;
	  Real artificialViscosity = 0.0;
	  PatchGodunov patchGod;
	  patchGod.define(grids.physDomain(), dx[lev], &advectPhys,
			  normalPredOrder,
			  useFourthOrderSlopes,
			  usePrimLimiting,
			  useCharLimiting,
			  useFlattening,
			  useArtificialViscosity,
			  artificialViscosity);
	  AdvectPhysics* advectPhysPtr = dynamic_cast<AdvectPhysics*>(patchGod.getGodunovPhysicsPtr());
	  CH_assert(advectPhysPtr != NULL);
	  
	  advectPhysPtr->setVelocities(&ccvel,&fcvel);

	  patchGod.setCurrentTime(0.0);
	  patchGod.setCurrentBox(grids[dit]);
	  patchGod.computeWHalf(faceh, cch, src, dt, grids[dit]);

	}

      //work out uh
      for (DataIterator dit(grids);dit.ok();++dit)
    	{	  
	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      FArrayBox& flux = (*fluxOfIce[lev])[dit][dir];
	      const FArrayBox& vel = (*fcVel[lev])[dit][dir];
   
	      CH_assert(vel.norm(0) < 1.0e+12);

	      flux.copy((*fcVel[lev])[dit][dir]);
	      flux.mult((*fcThck[lev])[dit][dir]);
	     
    	    } // end loop over direction
    	} // end loop over grids

      fluxOfIce[lev] -> exchange();

    } //end loop over levels


  for (int lev = 0; lev < numLevels; lev++)
    {
      if (ccVel[lev] != NULL)
	{
	  delete ccVel[lev];ccVel[lev]= NULL;
	}
      if (fcVel[lev] != NULL)
	{
	  delete fcVel[lev];fcVel[lev]= NULL;
	}
      if (fcThck[lev] != NULL)
	{
	  delete fcThck[lev];fcThck[lev]= NULL;
	}
      
    }


}

void computeIceStats(Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		     Vector<LevelData<FArrayBox>* >& topography, 
		     Vector<LevelData<FArrayBox>* >& thickness, 
		     Vector<LevelData<FArrayBox>* >& melangeThickness, 
		     Vector<LevelData<FArrayBox>* >& sectorMask,
		     Vector<Real>& dx, Vector<int>& ratio,
		     int maskNo)
{ 

  CH_TIMERS("computeIceStats");
  CH_TIMER("createSigmaCS",t1);
  CH_TIMER("integrateH",t2);
  CH_TIMER("integrateHab",t3);
  CH_TIMER("integrateGA",t4);
  CH_TIMER("integrateTA",t5);
  
  int numLevels = topography.size();


  CH_START(t2);
   
  Real iceVolumeAll = 0.0;
  Real melangeVolume = 0.0;

  Vector<LevelData<FArrayBox>* > tmp(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > tmp2(numLevels, NULL);

  for (int lev=0; lev< numLevels; lev++)
    {
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      tmp[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      tmp2[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  FArrayBox& sectorThck = (*tmp[lev])[dit];
	  sectorThck.setVal(0.0);
	  FArrayBox& sectorMelange = (*tmp2[lev])[dit];
	  sectorMelange.setVal(0.0);

	  const FArrayBox& thck = (*thickness[lev])[dit];
	  const FArrayBox& melange = (*melangeThickness[lev])[dit];

	  const Box& b = grids[dit];

	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		{
		  sectorThck(iv)=thck(iv);
		  sectorMelange(iv)=melange(iv);
		}
	    }
	}    
    }
  iceVolumeAll = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
  melangeVolume = computeSum(tmp2, ratio, dx[0], Interval(0,0), 0);

  CH_STOP(t2);
  Real iceVolumeAbove = 0.0;
  Real groundedArea = 0.0;
  Real groundedPlusOpenLandArea = 0.0;
  Real floatingArea = 0.0;
  Real totalArea = 0.0;


  if (iceVolumeAll > 1.0e-10)
    {
      CH_START(t1);
             
      Vector<LevelData<FArrayBox>* > tmp3(numLevels, NULL);
      CH_STOP(t1);
       
      CH_START(t3);
       
      {
	//Compute the total thickness above flotation;
	for (int lev=0; lev< numLevels; lev++)
	  {
	    const DisjointBoxLayout& grids = coords[lev]->grids();

	    for (DataIterator dit(grids);dit.ok();++dit)
	      {
		const FArrayBox& Hab =  coords[lev]->getThicknessOverFlotation()[dit];
		const Box& b = grids[dit];
		FArrayBox& sectorHab = (*tmp[lev])[dit];
		sectorHab.setVal(0.0);
		for (BoxIterator bit(b);bit.ok();++bit)
		  {
		    const IntVect& iv = bit();
		    if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		      {
			sectorHab(iv) = Hab(iv);
		      }
		
		  }
	      }
	  }
	iceVolumeAbove = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
	 
      }
      CH_STOP(t3);
      CH_START(t4);
       
      {
	//grounded and floating area
	 
	for (int lev=0; lev< numLevels; lev++)
	  {	     
	    const DisjointBoxLayout& grids = coords[lev]->grids();
	    tmp3[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
	    for (DataIterator dit(grids);dit.ok();++dit)
	      {
		const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
		const Box& b = grids[dit];
		FArrayBox& a = (*tmp[lev])[dit];
		FArrayBox& a2 = (*tmp2[lev])[dit];
		FArrayBox& a3 = (*tmp3[lev])[dit];
		a.setVal(0.0);
		a2.setVal(0.0);
                a3.setVal(0.0);
		for (BoxIterator bit(b);bit.ok();++bit)
		  {
		    const IntVect& iv = bit();

		    if ( std::abs ( (*sectorMask[lev])[dit](iv) - maskNo) < 1.0e-6 || maskNo == -1)
		      {
			if (mask(iv) == GROUNDEDMASKVAL)
			  {
			    a(iv) = 1.0;
			  }

			if (mask(iv) == FLOATINGMASKVAL)
			  {
			    a2(iv) = 1.0;
			  }

			// grounded plus open land
			if ((mask(iv) == GROUNDEDMASKVAL) || (mask(iv) == OPENLANDMASKVAL))
			  {
			    a3(iv) = 1.0;
			  }
		      }
		
		  }
	      }
	  }
	groundedArea = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
	floatingArea = computeSum(tmp2, ratio, dx[0], Interval(0,0), 0);
        groundedPlusOpenLandArea = computeSum(tmp3, ratio, dx[0], Interval(0,0), 0);

    
      }
      CH_STOP(t4);
      CH_START(t5);
      {
	//total area
	for (int lev=0; lev< numLevels; lev++)
	  {
	     
	    const DisjointBoxLayout& grids = coords[lev]->grids();
	    for (DataIterator dit(grids);dit.ok();++dit)
	      {
		const BaseFab<int>& mask =  coords[lev]->getFloatingMask()[dit];
		const Box& b = grids[dit];
		FArrayBox& a = (*tmp[lev])[dit];
		a.setVal(0.0);
		for (BoxIterator bit(b);bit.ok();++bit)
		  {
		    const IntVect& iv = bit();
		    if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		      {
			if ((mask(iv) == GROUNDEDMASKVAL) || (mask(iv) == FLOATINGMASKVAL))
			  {
			    a(iv) = 1.0;
			  }
		      }
		
		  }
	      }
	  }
	totalArea = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
      }
      CH_STOP(t5);

      for (int lev=0; lev< numLevels; lev++)
	{
	  if (tmp[lev] != NULL)
	    {
	      delete tmp[lev]; tmp[lev];
	    } 

	  if (tmp2[lev] != NULL)
	    {
	      delete tmp2[lev]; tmp2[lev];
	    } 

	  if (tmp3[lev] != NULL)
	    {
	      delete tmp3[lev]; tmp3[lev];
	    } 
	}
    }
    
  pout() << " iceVolumeAll = " << iceVolumeAll << " ";
  pout() << " iceVolumeAbove = " << iceVolumeAbove << " ";
  pout() << " melangeVolume = " << melangeVolume << " ";
  pout() << " groundedArea = " << groundedArea << " ";
  pout() << " floatingArea = " << floatingArea << " ";
  pout() << " totalAreaOfIce = " << totalArea << " ";
  pout() << " groundedPlusOpenLandArea = " << groundedPlusOpenLandArea << " ";
  
}

void computeVolCons(const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		    Vector<LevelData<FluxBox>* >& fluxOfIce,
		    Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
		    Vector<LevelData<FArrayBox>* >& basalThicknessSource,
		    Vector<LevelData<FArrayBox>* >& deltaThickness,
		    Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
		    Vector<LevelData<FArrayBox>* >& calvingFlux,
		    Vector<LevelData<FArrayBox>* >& topography,
		    Vector<LevelData<FArrayBox>* >& thickness, 
		    Vector<LevelData<FArrayBox>* >& sectorMask,
		    Vector<Real>& dx, Vector<int>& ratio, 
		    int maskNo)

{

  
  int numLevels = topography.size();
  Vector<LevelData<FArrayBox>* > ccDH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccFlxDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccSMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccBMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccCMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccDH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccFlxDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccSMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccBMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccCMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

      //work out (dh/dt)_c = smb - bmb + div(uh) and discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  FArrayBox& smb = (*ccSMB[lev])[dit];
    	  FArrayBox& bmb = (*ccBMB[lev])[dit];
    	  FArrayBox& cmb = (*ccCMB[lev])[dit];
    	  FArrayBox& dhdt = (*ccDH[lev])[dit];
	  FArrayBox& flxDiv = (*ccFlxDiv[lev])[dit];
	  FArrayBox& div = (*ccDiv[lev])[dit];
	  smb.setVal(0.0);
	  bmb.setVal(0.0);
	  cmb.setVal(0.0);
	  dhdt.setVal(0.0);
	  flxDiv.setVal(0.0);
	  div.setVal(0.0);

    	  const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		{
		  smb(iv) = (*surfaceThicknessSource[lev])[dit](iv);
		  bmb(iv) = (*basalThicknessSource[lev])[dit](iv);
		  cmb(iv) = (*calvingFlux[lev])[dit](iv);
		  dhdt(iv) = (*deltaThickness[lev])[dit](iv); 
		  div(iv) = (*divergenceThicknessFlux[lev])[dit](iv);
		}
	    }

	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& flux = (*fluxOfIce[lev])[dit][dir];
   
	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		    {
		      Real epsThck = 1.0e-10;
		      
		      flxDiv(iv) += (flux(iv + BASISV(dir)) - flux(iv))/dx[lev]; // flux divergence
		      
		      if ((thck(iv) < epsThck) || (mask(iv) != GROUNDEDMASKVAL && mask(iv) != FLOATINGMASKVAL))
			{

			  if (thck(iv + BASISV(dir)) > epsThck && (mask(iv + BASISV(dir)) == GROUNDEDMASKVAL || mask(iv + BASISV(dir)) == FLOATINGMASKVAL) ) 
			    {
			      discharge(iv) -= flux(iv + BASISV(dir)) / dx[lev];
			    }
			  if (thck(iv - BASISV(dir)) > epsThck && (mask(iv - BASISV(dir)) == GROUNDEDMASKVAL || mask(iv - BASISV(dir)) == FLOATINGMASKVAL) )
			    {
			      discharge(iv) += flux(iv) / dx[lev];
			    }
			} //end if (thck(iv) < epsThck) 

		    } 
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels

  Real sumDischarge = computeSum(ccDischarge, ratio, dx[0], Interval(0,0), 0);
  pout() << " discharge = " << sumDischarge << " ";
  
  Real sumFlxDiv = computeSum(ccFlxDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceReconstr = " << sumFlxDiv << " ";

  Real sumDiv = computeSum(ccDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceFromFile = " << sumDiv << " ";

  Real sumDH = computeSum(ccDH, ratio, dx[0], Interval(0,0), 0);
  pout() << " dhdt = " << sumDH << " ";

  Real sumSMB = computeSum(ccSMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " smb = " << sumSMB << " ";

  Real sumBMB = computeSum(ccBMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " bmb = " << sumBMB << " ";

  Real sumCMB = computeSum(ccCMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " calving = " << sumCMB << " ";


  Real err = sumSMB + sumBMB - sumCMB - sumDH;
  pout() << " (smb+bmb-calving-dhdt) = " << err << " ";



  for (int lev = 0; lev < numLevels; lev++)
    {
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDH[lev] != NULL)
	{
	  delete ccDH[lev];ccDH[lev]= NULL;
	}
      if (ccFlxDiv[lev] != NULL)
	{
	  delete ccFlxDiv[lev];ccFlxDiv[lev]= NULL;
	}
      if (ccDiv[lev] != NULL)
	{
	  delete ccDiv[lev];ccDiv[lev]= NULL;
	}
       if (ccSMB[lev] != NULL)
	{
	  delete ccSMB[lev];ccSMB[lev]= NULL;
	}
      if (ccBMB[lev] != NULL)
	{
	  delete ccBMB[lev];ccBMB[lev]= NULL;
	}
      if (ccCMB[lev] != NULL)
	{
	  delete ccCMB[lev];ccCMB[lev]= NULL;
	}
      
    }

}

void computeIceSheet(const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		    Vector<LevelData<FluxBox>* >& fluxOfIce,
		    Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
		    Vector<LevelData<FArrayBox>* >& basalThicknessSource,
		    Vector<LevelData<FArrayBox>* >& deltaThickness,
		    Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
		    Vector<LevelData<FArrayBox>* >& calvingFlux,
		    Vector<LevelData<FArrayBox>* >& topography,
		    Vector<LevelData<FArrayBox>* >& thickness, 
		    Vector<LevelData<FArrayBox>* >& sectorMask,
		    Vector<Real>& dx, Vector<int>& ratio, 
		    int maskNo)

{

  
  int numLevels = topography.size();
  Vector<LevelData<FArrayBox>* > ccDH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccFlxDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccSMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccBMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccCMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccDH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccFlxDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccSMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccBMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccCMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

      //work out (dh/dt)_c = smb - bmb + div(uh) and discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  FArrayBox& smb = (*ccSMB[lev])[dit];
    	  FArrayBox& bmb = (*ccBMB[lev])[dit];
    	  FArrayBox& cmb = (*ccCMB[lev])[dit];
    	  FArrayBox& dhdt = (*ccDH[lev])[dit];
	  FArrayBox& flxDiv = (*ccFlxDiv[lev])[dit];
	  FArrayBox& div = (*ccDiv[lev])[dit];
	  smb.setVal(0.0);
	  bmb.setVal(0.0);
	  cmb.setVal(0.0);
	  dhdt.setVal(0.0);
	  flxDiv.setVal(0.0);
	  div.setVal(0.0);

    	  const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  Real epsThck = 1.0e-10;

	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		{
		  if (thck(iv) > epsThck)
		    {
		      smb(iv) = (*surfaceThicknessSource[lev])[dit](iv);
		      bmb(iv) = (*basalThicknessSource[lev])[dit](iv);
		      cmb(iv) = (*calvingFlux[lev])[dit](iv);
		      dhdt(iv) = (*deltaThickness[lev])[dit](iv); 
		      div(iv) = (*divergenceThicknessFlux[lev])[dit](iv); 
		    }
		}
	    }

	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& flux = (*fluxOfIce[lev])[dit][dir];
   
	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		    {
		      
		      if (thck(iv) > epsThck)
			{
			  flxDiv(iv) += (flux(iv + BASISV(dir)) - flux(iv))/dx[lev]; // flux divergence
			}
		      else
			{

			  if (thck(iv + BASISV(dir)) > epsThck)
			    {
			      discharge(iv) -= flux(iv + BASISV(dir)) / dx[lev];
			    }
			  if (thck(iv - BASISV(dir)) > epsThck)
			    {
			      discharge(iv) += flux(iv) / dx[lev];
			    }
			} //end if (thck(iv) <= epsThck) 

		    } 
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels

  Real sumDischarge = computeSum(ccDischarge, ratio, dx[0], Interval(0,0), 0);
  pout() << " discharge = " << sumDischarge << " ";
  
  Real sumDH = computeSum(ccDH, ratio, dx[0], Interval(0,0), 0);
  pout() << " dhdt = " << sumDH << " ";

  Real sumSMB = computeSum(ccSMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " smb = " << sumSMB << " ";

  Real sumBMB = computeSum(ccBMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " bmb = " << sumBMB << " ";

  Real sumCMB = computeSum(ccCMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " calving = " << sumCMB << " ";

  Real sumFlxDiv = computeSum(ccFlxDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceReconstr = " << sumFlxDiv << " ";

  Real sumDiv = computeSum(ccDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceFromFile = " << sumDiv << " ";

  Real err = sumSMB + sumBMB - sumFlxDiv - sumDH;
  pout() << " (smb+bmb-flxDivReconstr-dhdt) = " << err << " ";

  err = sumSMB + sumBMB - sumDiv - sumDH;
  pout() << " (smb+bmb-flxDivFromFile-dhdt) = " << err << " ";

  for (int lev = 0; lev < numLevels; lev++)
    {
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDH[lev] != NULL)
	{
	  delete ccDH[lev];ccDH[lev]= NULL;
	}
      if (ccFlxDiv[lev] != NULL)
	{
	  delete ccFlxDiv[lev];ccFlxDiv[lev]= NULL;
	}
      if (ccDiv[lev] != NULL)
	{
	  delete ccDiv[lev];ccDiv[lev]= NULL;
	}
       if (ccSMB[lev] != NULL)
	{
	  delete ccSMB[lev];ccSMB[lev]= NULL;
	}
      if (ccBMB[lev] != NULL)
	{
	  delete ccBMB[lev];ccBMB[lev]= NULL;
	}
      if (ccCMB[lev] != NULL)
	{
	  delete ccCMB[lev];ccCMB[lev]= NULL;
	}
      
    }

}

void computeOutsideIce(const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		    Vector<LevelData<FluxBox>* >& fluxOfIce,
		    Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
		    Vector<LevelData<FArrayBox>* >& basalThicknessSource,
		    Vector<LevelData<FArrayBox>* >& deltaThickness,
		    Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
		    Vector<LevelData<FArrayBox>* >& calvingFlux,
		    Vector<LevelData<FArrayBox>* >& topography,
		    Vector<LevelData<FArrayBox>* >& thickness, 
		    Vector<LevelData<FArrayBox>* >& sectorMask,
		    Vector<Real>& dx, Vector<int>& ratio, 
		    int maskNo)

{

  
  int numLevels = topography.size();
  Vector<LevelData<FArrayBox>* > ccDH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccFlxDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccSMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccBMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccCMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccDH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccFlxDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccSMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccBMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccCMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

      //work out (dh/dt)_c = smb - bmb + div(uh) and discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  FArrayBox& smb = (*ccSMB[lev])[dit];
    	  FArrayBox& bmb = (*ccBMB[lev])[dit];
    	  FArrayBox& cmb = (*ccCMB[lev])[dit];
    	  FArrayBox& dhdt = (*ccDH[lev])[dit];
	  FArrayBox& flxDiv = (*ccFlxDiv[lev])[dit];
	  FArrayBox& div = (*ccDiv[lev])[dit];
	  smb.setVal(0.0);
	  bmb.setVal(0.0);
	  cmb.setVal(0.0);
	  dhdt.setVal(0.0);
	  flxDiv.setVal(0.0);
	  div.setVal(0.0);

    	  const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  Real epsThck = 1.0e-10;

	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		{
		  if (thck(iv) < epsThck)
		    {
		      smb(iv) = (*surfaceThicknessSource[lev])[dit](iv);
		      bmb(iv) = (*basalThicknessSource[lev])[dit](iv);
		      cmb(iv) = (*calvingFlux[lev])[dit](iv);
		      dhdt(iv) = (*deltaThickness[lev])[dit](iv); 
		      div(iv) = (*divergenceThicknessFlux[lev])[dit](iv); 
		    }
		}
	    }

	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& flux = (*fluxOfIce[lev])[dit][dir];
   
	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		    {
		      
		      if (thck(iv) < epsThck)
			{
			  flxDiv(iv) += (flux(iv + BASISV(dir)) - flux(iv))/dx[lev]; // flux divergence
			}
		      else
			{

			  if (thck(iv + BASISV(dir)) < epsThck)
			    {
			      discharge(iv) -= flux(iv + BASISV(dir)) / dx[lev];
			    }
			  if (thck(iv - BASISV(dir)) < epsThck)
			    {
			      discharge(iv) += flux(iv) / dx[lev];
			    }
			} //end if (thck(iv) <= epsThck) 

		    } 
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels

  Real sumDischarge = computeSum(ccDischarge, ratio, dx[0], Interval(0,0), 0);
  pout() << " discharge = " << sumDischarge << " ";
  
  Real sumDH = computeSum(ccDH, ratio, dx[0], Interval(0,0), 0);
  pout() << " dhdt = " << sumDH << " ";

  Real sumSMB = computeSum(ccSMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " smb = " << sumSMB << " ";

  Real sumBMB = computeSum(ccBMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " bmb = " << sumBMB << " ";

  Real sumCMB = computeSum(ccCMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " calving = " << sumCMB << " ";

  Real sumFlxDiv = computeSum(ccFlxDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceReconstr = " << sumFlxDiv << " ";

  Real sumDiv = computeSum(ccDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceFromFile = " << sumDiv << " ";

  // err is not zero if the calving model maintains the front position - check runtime diagnostics (endTimestepDiagnotics in AmrIce.cpp)  
  Real err = sumSMB + sumBMB - sumFlxDiv - sumCMB - sumDH;
  pout() << " (smb+bmb-flxDivReconstr-calving-dhdt) = " << err << " ";

  err = sumSMB + sumBMB - sumDiv - sumCMB - sumDH;
  pout() << " (smb+bmb-flxDivFromFile-calving-dhdt) = " << err << " ";

  for (int lev = 0; lev < numLevels; lev++)
    {
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDH[lev] != NULL)
	{
	  delete ccDH[lev];ccDH[lev]= NULL;
	}
      if (ccFlxDiv[lev] != NULL)
	{
	  delete ccFlxDiv[lev];ccFlxDiv[lev]= NULL;
	}
      if (ccDiv[lev] != NULL)
	{
	  delete ccDiv[lev];ccDiv[lev]= NULL;
	}
       if (ccSMB[lev] != NULL)
	{
	  delete ccSMB[lev];ccSMB[lev]= NULL;
	}
      if (ccBMB[lev] != NULL)
	{
	  delete ccBMB[lev];ccBMB[lev]= NULL;
	}
      if (ccCMB[lev] != NULL)
	{
	  delete ccCMB[lev];ccCMB[lev]= NULL;
	}
      
    }

}


void computeGrounded(const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		     Vector<LevelData<FluxBox>* >& fluxOfIce,
		     Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
		     Vector<LevelData<FArrayBox>* >& basalThicknessSource,
		     Vector<LevelData<FArrayBox>* >& deltaThickness,
		     Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux,
		     Vector<LevelData<FArrayBox>* >& topography,
		     Vector<LevelData<FArrayBox>* >& thickness, 
		     Vector<LevelData<FArrayBox>* >& sectorMask,
		     const Real Hmin,
		     Vector<Real>& dx, Vector<int>& ratio, 
		     int maskNo)

{

  
  int numLevels = topography.size();
  Vector<LevelData<FArrayBox>* > ccDH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccFlxDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccSMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccBMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > tmp(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccDH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccFlxDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccSMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccBMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      tmp[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

      //work out (dh/dt) = smb - bmb + div(uh) and discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  FArrayBox& smb = (*ccSMB[lev])[dit];
    	  FArrayBox& bmb = (*ccBMB[lev])[dit];
    	  FArrayBox& dhdt = (*ccDH[lev])[dit];
	  FArrayBox& flxDiv = (*ccFlxDiv[lev])[dit];
	  FArrayBox& div = (*ccDiv[lev])[dit];
	  smb.setVal(0.0);
	  bmb.setVal(0.0);
	  dhdt.setVal(0.0);
	  flxDiv.setVal(0.0);
	  div.setVal(0.0);
	  FArrayBox& a = (*tmp[lev])[dit];
	  a.setVal(0.0);

    	  const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		{
		  if (thck(iv) > Hmin && mask(iv) == GROUNDEDMASKVAL)
		    {
		      smb(iv) = (*surfaceThicknessSource[lev])[dit](iv);
		      bmb(iv) = (*basalThicknessSource[lev])[dit](iv);
		      dhdt(iv) = (*deltaThickness[lev])[dit](iv);
 		      div(iv) = (*divergenceThicknessFlux[lev])[dit](iv); 
		      a(iv) = 1.0;
		    }
		}
	    }

	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& flux = (*fluxOfIce[lev])[dit][dir];
   
	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		    {
		      if (thck(iv) > Hmin && mask(iv) == GROUNDEDMASKVAL)
			{
			  flxDiv(iv) += (flux(iv + BASISV(dir)) - flux(iv))/dx[lev]; // flux divergence
			}
		      else
			{

			  if (thck(iv + BASISV(dir)) > Hmin && (mask(iv + BASISV(dir)) == GROUNDEDMASKVAL))
			    {
			      discharge(iv) -= flux(iv + BASISV(dir)) / dx[lev];
			    }
			  if (thck(iv - BASISV(dir)) > Hmin && (mask(iv - BASISV(dir)) == GROUNDEDMASKVAL))
			    {
			      discharge(iv) += flux(iv) / dx[lev];
			    }
			} //end if 
		    } // sector loops
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels


  pout() << " minimumThickness = " << Hmin << " ";

  Real groundedArea = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
  pout() << " groundedArea = " << groundedArea << " ";

  Real sumDischarge = computeSum(ccDischarge, ratio, dx[0], Interval(0,0), 0);
  pout() << " discharge = " << sumDischarge << " ";
  
  Real sumDH = computeSum(ccDH, ratio, dx[0], Interval(0,0), 0);
  pout() << " dhdt = " << sumDH << " ";

  Real sumSMB = computeSum(ccSMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " smb = " << sumSMB << " ";

  Real sumBMB = computeSum(ccBMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " bmb = " << sumBMB << " ";

  Real sumFlxDiv = computeSum(ccFlxDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceReconstr = " << sumFlxDiv << " ";

 Real sumDiv = computeSum(ccDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceFromFile = " << sumDiv << " ";

  Real err = sumSMB + sumBMB - sumFlxDiv - sumDH;
  pout() << " (smb+bmb-flxDivReconstr-dhdt) = " << err << " ";

  err = sumSMB + sumBMB - sumDiv - sumDH;
  pout() << " (smb+bmb-flxDivFromFile-dhdt) = " << err << " ";


  for (int lev = 0; lev < numLevels; lev++)
    {
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDH[lev] != NULL)
	{
	  delete ccDH[lev];ccDH[lev]= NULL;
	}
      if (ccFlxDiv[lev] != NULL)
	{
	  delete ccFlxDiv[lev];ccFlxDiv[lev]= NULL;
	}
      if (ccDiv[lev] != NULL)
	{
	  delete ccDiv[lev];ccDiv[lev]= NULL;
	}
       if (ccSMB[lev] != NULL)
	{
	  delete ccSMB[lev];ccSMB[lev]= NULL;
	}
      if (ccBMB[lev] != NULL)
	{
	  delete ccBMB[lev];ccBMB[lev]= NULL;
	}
      if (tmp[lev] != NULL)
	{
	  delete tmp[lev]; tmp[lev];
	} 
      
    }

}

void computeFloating(const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
		     Vector<LevelData<FluxBox>* >& fluxOfIce,
		     Vector<LevelData<FArrayBox>* >& surfaceThicknessSource, 
		     Vector<LevelData<FArrayBox>* >& basalThicknessSource,
		     Vector<LevelData<FArrayBox>* >& deltaThickness,
		     Vector<LevelData<FArrayBox>* >& divergenceThicknessFlux, 
		     Vector<LevelData<FArrayBox>* >& calvingFlux,
		     Vector<LevelData<FArrayBox>* >& topography,
		     Vector<LevelData<FArrayBox>* >& thickness, 
		     Vector<LevelData<FArrayBox>* >& sectorMask,
		     Vector<Real>& dx, Vector<int>& ratio, 
		     int maskNo)

{

  
  int numLevels = topography.size();
  Vector<LevelData<FArrayBox>* > ccDH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccFlxDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDiv(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccSMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccBMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccCMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccDH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccFlxDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDiv[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccSMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccBMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccCMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

      //work out (dh/dt) = smb - bmb + div(uh) and discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  FArrayBox& smb = (*ccSMB[lev])[dit];
    	  FArrayBox& bmb = (*ccBMB[lev])[dit];
    	  FArrayBox& cmb = (*ccCMB[lev])[dit];
    	  FArrayBox& dhdt = (*ccDH[lev])[dit];
	  FArrayBox& flxDiv = (*ccFlxDiv[lev])[dit];
	  FArrayBox& div = (*ccDiv[lev])[dit];
	  smb.setVal(0.0);
	  bmb.setVal(0.0);
	  cmb.setVal(0.0);
	  dhdt.setVal(0.0);
	  flxDiv.setVal(0.0);
	  div.setVal(0.0);

    	  const BaseFab<int>& mask = coords[lev]->getFloatingMask()[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		{
		  if (mask(iv) == FLOATINGMASKVAL)
		    {
		      smb(iv) = (*surfaceThicknessSource[lev])[dit](iv);
		      bmb(iv) = (*basalThicknessSource[lev])[dit](iv);
		      dhdt(iv) = (*deltaThickness[lev])[dit](iv);
		      div(iv) = (*divergenceThicknessFlux[lev])[dit](iv);
		    }

		}
	    }

	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& flux = (*fluxOfIce[lev])[dit][dir];
   
	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  if ( std::abs ( (*sectorMask[lev])[dit](iv) - static_cast<Real>(maskNo)) < 1.0e-6 || maskNo == -1)
		    {
		      Real epsThck = 1.0e-10;
		      
		      if (mask(iv)==FLOATINGMASKVAL)
			{
			  flxDiv(iv) += (flux(iv + BASISV(dir)) - flux(iv))/dx[lev]; // flux divergence
			}
		  		     
		      else if (thck(iv) < epsThck)
			{
			  cmb(iv) = (*calvingFlux[lev])[dit](iv);
			  if (thck(iv + BASISV(dir)) > epsThck && (mask(iv + BASISV(dir)) == FLOATINGMASKVAL) ) 
			    {
			      discharge(iv) -= flux(iv + BASISV(dir)) / dx[lev];
			    }
			  if (thck(iv - BASISV(dir)) > epsThck && (mask(iv - BASISV(dir)) == FLOATINGMASKVAL) )
			    {
			      discharge(iv) += flux(iv) / dx[lev];
			    }
			} //end if (thck(iv) < epsThck) 

		    } 
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels

  Real sumDischarge = computeSum(ccDischarge, ratio, dx[0], Interval(0,0), 0);
  pout() << " dischargeFloating = " << sumDischarge << " ";
  
  Real sumCMB = computeSum(ccCMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " calving = " << sumCMB << " ";

  Real sumDH = computeSum(ccDH, ratio, dx[0], Interval(0,0), 0);
  pout() << " dhdt = " << sumDH << " ";

  Real sumSMB = computeSum(ccSMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " smb = " << sumSMB << " ";

  Real sumBMB = computeSum(ccBMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " bmb = " << sumBMB << " ";

  Real sumFlxDiv = computeSum(ccFlxDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceReconstr = " << sumFlxDiv << " ";

  Real sumDiv = computeSum(ccDiv, ratio, dx[0], Interval(0,0), 0);
  pout() << " fluxDivergenceFromFile = " << sumDiv << " ";

  Real err = sumSMB + sumBMB - sumFlxDiv - sumDH;
  pout() << " (smb+bmb-flxDivReconstr-dhdt) = " << err << " ";

  err = sumSMB + sumBMB - sumDiv - sumDH;
  pout() << " (smb+bmb-flxDivFromFile-dhdt) = " << err << " ";


  for (int lev = 0; lev < numLevels; lev++)
    {
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDH[lev] != NULL)
	{
	  delete ccDH[lev];ccDH[lev]= NULL;
	}
      if (ccFlxDiv[lev] != NULL)
	{
	  delete ccFlxDiv[lev];ccFlxDiv[lev]= NULL;
	}
      if (ccDiv[lev] != NULL)
	{
	  delete ccDiv[lev];ccDiv[lev]= NULL;
	}
       if (ccSMB[lev] != NULL)
	{
	  delete ccSMB[lev];ccSMB[lev]= NULL;
	}
      if (ccBMB[lev] != NULL)
	{
	  delete ccBMB[lev];ccBMB[lev]= NULL;
	}
      if (ccCMB[lev] != NULL)
	{
	  delete ccCMB[lev];ccCMB[lev]= NULL;
	}
      
    }

}


int main(int argc, char* argv[]) {

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif
   
    CH_TIMERS("stats");
    CH_TIMER("loadplot",tp);
    if(argc < 5) 
      { 
	std::cerr << " usage: " << argv[0] << " <plot file> <ice_density> <water_density> <gravity> [mask_file] [mask_no_start = 0] [mask_no_end = mask_no_start] " << std::endl; 
	exit(0); 
      }
    char* plotFile = argv[1];
    Real iceDensity = atof(argv[2]);
    Real waterDensity = atof(argv[3]);
    Real gravity = atof(argv[4]);
    char* maskFile = (argc > 5)?argv[5]:NULL;
    int maskNoStart = 0;
    if (maskFile && argc > 6)
      {
    	maskNoStart = atoi(argv[6]);
      }
    int maskNoEnd = maskNoStart;
    if (maskFile && argc > 7)
      {
    	maskNoEnd = atoi(argv[7]);
      }
    Real Hmin = 50.0;

    Box domainBox;
    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    pout() << " Reading plot file " << endl  ;

    Real dt ,crseDx, time;
    CH_START(tp);
    ReadAMRHierarchyHDF5(std::string(plotFile), grids, data, name , 
			 domainBox, crseDx, dt, time, ratio, numLevels);

    Vector<ProblemDomain> domain(numLevels,domainBox);
    Vector<RealVect> vdx(numLevels,RealVect::Unit*crseDx);
    Vector<Real> dx(numLevels,crseDx);
    for (int lev=1;lev<numLevels;++lev)
      {
	dx[lev] = dx[lev-1] / Real(ratio[lev-1]);
	vdx[lev] = vdx[lev-1] / Real(ratio[lev-1]);
	domain[lev] = domain[lev-1];
	domain[lev].refine(ratio[lev-1]);
      }

    //load the sector mask, if it exists
    //Vector<RefCountedPtr<LevelData<FArrayBox> > > sectorMask;
    
    Box mdomainBox;
    Vector<std::string> mname;
    Vector<LevelData<FArrayBox>* > mdata;
    Vector<DisjointBoxLayout > mgrids;
    Vector<int > mratio;
    int mnumLevels;
    Real mdt ,mcrseDx, mtime;

    if (maskFile)
      {
	pout() << " Reading mask " << endl  ;    	
	ReadAMRHierarchyHDF5(std::string(maskFile), mgrids, mdata, mname , 
			     mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);
      }
    
    CH_STOP(tp);

    Vector<LevelData<FArrayBox>* > sectorMask(numLevels,NULL);
    for (int lev=0;lev<numLevels;++lev)
      {
	sectorMask[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);

	for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
	  {
	    (*sectorMask[lev])[dit].setVal(0.0);
	  }

	if (maskFile)
	  {
	    FillFromReference(*sectorMask[lev], *mdata[0], RealVect::Unit*dx[lev], RealVect::Unit*mcrseDx,true);
	  }

      }

    // Test mask
    for (int lev=0; lev< numLevels; lev++)
      {
	for (DataIterator dit(grids[lev]);dit.ok();++dit)
	  {
	    const Box& b = grids[lev][dit];

	    for (BoxIterator bit(b);bit.ok();++bit)
	      {
		const IntVect& iv = bit();
		Real msk=(*sectorMask[lev])[dit](iv);
		bool maskValExists=false;

		for (int maskNo = maskNoStart; maskNo<= maskNoEnd; maskNo++)
		  {
		    if ( std::abs ( msk - static_cast<Real>(maskNo)) < 1.0e-6)
		      {
			maskValExists=true;
		      }
		  }
		if (!maskValExists)
		  {
		    if (std::abs (msk) > 1.0e-6)
		      {
			// Assume non-masked areas are set to zero.
			Real tmp=std::round(msk);
			if (tmp < static_cast<Real>(maskNoStart))
			  {
			    tmp=std::max(static_cast<Real>(maskNoStart),tmp);
			  }
			(*sectorMask[lev])[dit](iv)=tmp;
			//pout() << " Mask value " << msk << " adjusted value " << tmp << "  lev " << lev << endl;
		      }
		  }
	      }    
	  }
      }

    Vector<LevelData<FArrayBox>* > thickness(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > topography(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > deltaThickness(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > surfaceThicknessSource(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > basalThicknessSource(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > divergenceThicknessFlux(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > calvingFlux(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > melangeThickness(numLevels,NULL);
    Vector<LevelData<FluxBox>* > fluxOfIce(numLevels, NULL);


    for (int lev=0;lev<numLevels;++lev)
      {
	// Extra ghost cells are needed to calculate advection using 4th order Godunov. 
	thickness[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
	topography[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
	deltaThickness[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
	surfaceThicknessSource[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
	basalThicknessSource[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
	divergenceThicknessFlux[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
	calvingFlux[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
	melangeThickness[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
	fluxOfIce[lev] = new LevelData<FluxBox>(grids[lev],1,IntVect::Unit);

      }

    
    pout().setf(ios_base::scientific,ios_base::floatfield); 
    pout().precision(12);

    pout() << " Create DEM " << endl  ;

    createDEM(topography, thickness, name, data, ratio, dx, mcrseDx);

    pout() << " Create Sigma CS " << endl  ;

    Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
    createSigmaCS(coords,topography, thickness, 
		  dx, ratio, iceDensity, waterDensity,gravity);

    pout() << " Get thickness equation components " << endl  ;

    getThicknessSource(surfaceThicknessSource, basalThicknessSource, deltaThickness, 
		       divergenceThicknessFlux, calvingFlux, melangeThickness,
		       topography, dx, ratio, name, data);

    pout() << " Calculate flux " << endl  ;

    computeFlux(fluxOfIce,coords,topography,thickness,surfaceThicknessSource,basalThicknessSource,
		dx,ratio,name,data);

    if (maskFile)
      {
	pout() << " Do sector calculations " << endl;
      }
    pout() << endl;

    for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
      {
	pout() << " Diagnostics for computational domain";
	if (maskFile)
	  {
	    pout() << ", sector = " << maskNo;
	  }
	pout() << endl;
	pout() << " time = " << time  ;

	computeIceStats(coords,topography, thickness, melangeThickness, sectorMask, dx, ratio, maskNo);
	pout() << endl;
	computeVolCons(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
		       deltaThickness, divergenceThicknessFlux, calvingFlux,
		       topography, thickness, sectorMask, dx, ratio, maskNo);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for grounded ice, excluding thin ice";
	if (maskFile)
	  {
	    pout() << ", sector = " << maskNo;
	  }
	pout() << endl;
	computeGrounded(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, topography, thickness, 
			sectorMask, Hmin, dx, ratio, maskNo);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for floating ice";
	if (maskFile)
	  {
	    pout() << ", sector = " << maskNo;
	  }
	pout() << endl;
	computeFloating(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, calvingFlux,
			topography, thickness, sectorMask, dx, ratio, maskNo);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for ice sheet";
	if (maskFile)
	  {
	    pout() << ", sector = " << maskNo;
	  }
	pout() << endl;
	computeIceSheet(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, calvingFlux,
			topography, thickness, sectorMask, dx, ratio, maskNo);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for outside ice sheet";
	if (maskFile)
	  {
	    pout() << ", sector = " << maskNo;
	  }
	pout() << endl;
	computeOutsideIce(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, calvingFlux,
			topography, thickness, sectorMask, dx, ratio, maskNo);
	pout() << endl;
	pout() << endl;

      }

    // now compute stats for all sectors
    if (maskNoStart != maskNoEnd)
      {
	pout() << " Diagnostics for computational domain";
	if (maskFile)
	  {
	    pout() << ", all sectors " << endl;
	  }
	pout() << " time = " << time  ;

	computeIceStats(coords,topography, thickness, melangeThickness, sectorMask, dx, ratio, -1);
	pout() << endl;
	computeVolCons(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
		       deltaThickness, divergenceThicknessFlux, calvingFlux,
		       topography, thickness, sectorMask, dx, ratio, -1);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for grounded ice, excluding thin ice";
	if (maskFile)
	  {
	    pout() << ", all sectors ";
	  }
	pout() << endl;
	computeGrounded(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, topography, thickness, 
			sectorMask, Hmin, dx, ratio, -1);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for floating ice";
	if (maskFile)
	  {
	    pout() << ", all sectors ";
	  }
	pout() << endl;
	computeFloating(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, calvingFlux,
			topography, thickness, sectorMask, dx, ratio, -1);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for ice sheet";
	if (maskFile)
	  {
	    pout() << ", all sectors ";
	  }
	pout() << endl;
	computeIceSheet(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, calvingFlux,
			topography, thickness, sectorMask, dx, ratio, -1);
	pout() << endl;
	pout() << endl;

	pout() << " Diagnostics for outside ice sheet";
	if (maskFile)
	  {
	    pout() << ", all sectors ";
	  }
	pout() << endl;
	computeOutsideIce(coords, fluxOfIce, surfaceThicknessSource, basalThicknessSource, 
			deltaThickness, divergenceThicknessFlux, calvingFlux,
			topography, thickness, sectorMask, dx, ratio, -1);
	pout() << endl;
	pout() << endl;

      }
    
    for (int lev=0;lev<numLevels;++lev)
      {
	if (sectorMask[lev] != NULL) delete sectorMask[lev];
	if (thickness[lev] != NULL) delete thickness[lev];
	if (topography[lev] != NULL) delete topography[lev];
	if (deltaThickness[lev] != NULL) delete deltaThickness[lev];
	if (surfaceThicknessSource[lev] != NULL) delete surfaceThicknessSource[lev];
	if (basalThicknessSource[lev] != NULL) delete basalThicknessSource[lev];
	if (divergenceThicknessFlux[lev] != NULL) delete divergenceThicknessFlux[lev];
	if (calvingFlux[lev] != NULL) delete calvingFlux[lev];
	if (melangeThickness[lev] != NULL) delete melangeThickness[lev];
	if (fluxOfIce[lev] != NULL) delete fluxOfIce[lev];
      }

		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
