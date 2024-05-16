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
// gia-stats.cpp
// version of stats.cpp which is specialized for instances where GIA is active
// 
// read in a bisicles plotfile (which must include thickness and Z_base) 
// and an optional mask, and write out a bunch of stats about the ice sheet. 
// These are
// 0. time
// 1. Volume of ice
// 2. Volume of ice above flotation
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

void createMaskedDEM( Vector<LevelData<FArrayBox>* >& topography,  
		      Vector<LevelData<FArrayBox>* >& thickness, 
		      Vector<LevelData<FArrayBox>* >& basalThicknessSrc, 
		      Vector<LevelData<FArrayBox>* >& mdata, 
		      Vector<std::string>& name, 
		      Vector<LevelData<FArrayBox>* >& data,
		      Vector<int>& ratio,
		      Vector<Real>& dx,
		      Real mcrseDx,
		      int maskStart,
                      int maskEnd)
{
  CH_TIME("createMaskedDEM");
  int numLevels = data.size();
  bool has_src = false;
  if (basalThicknessSrc[0] != NULL) has_src = true;

  for (int lev = 0; lev < numLevels; lev++)
    {

      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  (*topography[lev])[dit].setVal(0.0);
	  (*thickness[lev])[dit].setVal(0.0);
          if (has_src) (*basalThicknessSrc[lev])[dit].setVal(0.0);
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
	  else if (name[j] == "basalThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*basalThicknessSrc[lev],Interval(0,0));
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

          if (has_src)
            {
              filler.fillInterp(*basalThicknessSrc[lev],
                                *basalThicknessSrc[lev-1],
                                *basalThicknessSrc[lev-1],
                                time_interp_coeff,0, 0, 1);
            }
	}
      thickness[lev] -> exchange();
      topography[lev] -> exchange();
      if (has_src) basalThicknessSrc[lev]->exchange();
    }
  if (mdata.size() != 0)
    {
      // there is mask data, so set thickness to zero outside the selected region
      for (int lev =0; lev < numLevels; lev++)
	{
	  const DisjointBoxLayout& levelGrids = topography[lev]->disjointBoxLayout();
	  LevelData<FArrayBox> levelMask(levelGrids,1,IntVect::Zero);
	  FillFromReference(levelMask, *mdata[0], RealVect::Unit*dx[lev], RealVect::Unit*mcrseDx,true);
	  LevelData<FArrayBox>& levelData = *thickness[lev];
	  for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	    {
	      const Box& b = levelGrids[dit];
	      FArrayBox coef(b,1); coef.setVal(0.0);
	      for (BoxIterator bit(b); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		    
                  for (int maskNo = maskStart; maskNo<= maskEnd; maskNo++)
                    {
                      if ( std::abs ( levelMask[dit](iv) - maskNo) < 1.0e-6)
                        {
                          coef(iv) = 1.0;
                        }
                    }
		}
	      levelData[dit].mult(coef,0,0,1);
              if (has_src)
                {
                  (*basalThicknessSrc[lev])[dit].mult(coef,0,0,1);
                }
		
	    }
	    
	}
    }
}




void computeDischarge(Vector<LevelData<FArrayBox>* >& topography,
		      Vector<LevelData<FArrayBox>* >& thickness, 
		      Vector<Real>& dx, Vector<int>& ratio, 
		      Vector<std::string>& name, 
		      Vector<LevelData<FArrayBox>* >& data, 
		      Real a_iceDensity, Real a_waterDensity, Real a_gravity,
                      Real a_seaLevel)
{

  
  int numLevels = data.size();
  Vector<LevelData<FArrayBox>* > ccVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcThck(numLevels, NULL);
  Vector<LevelData<BaseFab<int> >* > ccMask(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDHCons(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccSMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccBMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      int comp = 0;
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccVel[lev] = new LevelData<FArrayBox>(grids,SpaceDim,IntVect::Unit);
      fcVel[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      fcThck[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      ccDH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccMask[lev] = new LevelData<BaseFab<int> >(grids,1,IntVect::Unit);
      ccDHCons[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccSMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccBMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  (*ccVel[lev])[dit].setVal(0.0);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      (*fcVel[lev])[dit][dir].setVal(0.0);
	    }
	}

      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
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
	  else if (name[j] == "dThickness/dt")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccDH[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "surfaceThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccSMB[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "basalThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccBMB[lev],Interval(0,0));
	      comp++;
	    }
	}
      CH_assert(comp == 5);

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

	  
	  const FArrayBox& thck = (*thickness[lev])[dit];
	  const FArrayBox& topg = (*topography[lev])[dit];
	  
	  FArrayBox usrf(topg.box(),1);

	  FORT_SURFACEHEIGHT(CHF_FRA1(usrf,0),
			     CHF_CONST_FRA1(thck,0),
			     CHF_CONST_FRA1(topg,0),
			     CHF_CONST_REAL(a_iceDensity),
			     CHF_CONST_REAL(a_waterDensity),
			     CHF_CONST_REAL(a_seaLevel),
			     CHF_BOX(topg.box()));

	  BaseFab<int>& mask = (*ccMask[lev])[dit];
	  int any = 0;

	  FORT_SETFLOATINGMASK(CHF_FIA1(mask,0),
			       CHF_CONST_FRA1(usrf,0),
			       CHF_CONST_FRA1(topg,0),
			       CHF_CONST_FRA1(thck,0),
			       CHF_INT(any),
			       CHF_REAL(a_iceDensity),
			       CHF_REAL(a_waterDensity),
			       CHF_REAL(a_seaLevel),
			       CHF_BOX(mask.box()));
    

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
		BaseFab<int> mask(grids[dit],1) ;

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
	  src.copy( (*ccSMB[lev])[dit]) ;
	  src.plus( (*ccBMB[lev])[dit]) ;
	  
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

      //work out (dh/dt)_c = smb - bmb + div(uh) and discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  FArrayBox& smb = (*ccSMB[lev])[dit];
    	  FArrayBox& bmb = (*ccBMB[lev])[dit];
    	  FArrayBox& dh = (*ccDH[lev])[dit];
	  FArrayBox& dhc = (*ccDHCons[lev])[dit];
    	  BaseFab<int>& mask = (*ccMask[lev])[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      dhc(iv) = smb(iv) + bmb(iv);
	    }

	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& vel = (*fcVel[lev])[dit][dir];
   
	      CH_assert(vel.norm(0) < 1.0e+12);

	      FArrayBox flux(vel.box(),1);
	      flux.copy((*fcVel[lev])[dit][dir]);
	      flux.mult((*fcThck[lev])[dit][dir]);
	     

	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  dhc(iv) += (flux(iv) - flux(iv + BASISV(dir)))/dx[lev]; // flux divergence
		  
		  Real epsThck = 10.0;// TODO fix magic number
		  if ((thck(iv) < epsThck) || (mask(iv) != GROUNDEDMASKVAL))
		    {
		      dh(iv) = 0.0;
		      dhc(iv) = 0.0;
		      smb(iv) = 0.0;
		      bmb(iv) = 0.0;
		      if (thck(iv + BASISV(dir)) > epsThck && (mask(iv + BASISV(dir)) == GROUNDEDMASKVAL) ) 
			{
			  discharge(iv) += -flux(iv + BASISV(dir)) / dx[lev];
			}
		      if (thck(iv - BASISV(dir)) > epsThck && (mask(iv - BASISV(dir)) == GROUNDEDMASKVAL) )
			{
			  discharge(iv) += flux(iv) / dx[lev];
			}
    		    } //end if (thck(iv) < epsThck) 
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

  Real sumDHC = computeSum(ccDHCons, ratio, dx[0], Interval(0,0), 0);
  pout() << " dhdtc = " << sumDH << " ";

  Real err = sumSMB + sumBMB - sumDHC -   sumDischarge;
  pout() << " (smb+bmb-dhtc-discharge) = " << err << " ";


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
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDH[lev] != NULL)
	{
	  delete ccDH[lev];ccDH[lev]= NULL;
	}
      if (ccMask[lev] != NULL)
	{
	  delete ccMask[lev];ccMask[lev]= NULL;
	}
      if (ccDHCons[lev] != NULL)
	{
	  delete ccDHCons[lev];ccDHCons[lev]= NULL;
	}
       if (ccSMB[lev] != NULL)
	{
	  delete ccSMB[lev];ccSMB[lev]= NULL;
	}
      if (ccBMB[lev] != NULL)
	{
	  delete ccBMB[lev];ccBMB[lev]= NULL;
	}
      
    }


}

void computeCalving(Vector<LevelData<FArrayBox>* >& topography,
                    Vector<LevelData<FArrayBox>* >& thickness, 
                    Vector<Real>& dx, Vector<int>& ratio, 
                    Vector<std::string>& name, 
                    Vector<LevelData<FArrayBox>* >& data, 
                    Real a_iceDensity, Real a_waterDensity,
                    Real a_gravity, Real a_seaLevel)
{

  
  int numLevels = data.size();
  Vector<LevelData<FArrayBox>* > ccVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcThck(numLevels, NULL);
  Vector<LevelData<BaseFab<int> >* > ccMask(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccCumCalvedIce(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccSMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccBMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccOMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccCMB(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > ccDischarge(numLevels, NULL);

  for (int lev = 0; lev < numLevels; lev++)
    {
      int comp = 0;
      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      ccVel[lev] = new LevelData<FArrayBox>(grids,SpaceDim,IntVect::Unit);
      fcVel[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      fcThck[lev] = new LevelData<FluxBox>(grids,1,IntVect::Unit);
      ccDH[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccMask[lev] = new LevelData<BaseFab<int> >(grids,1,IntVect::Unit);
      ccCumCalvedIce[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccSMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccBMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccOMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      ccCMB[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  (*ccVel[lev])[dit].setVal(0.0);
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      (*fcVel[lev])[dit][dir].setVal(0.0);
	    }
	}

      ccDischarge[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
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
	  else if (name[j] == "dThickness/dt")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccDH[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "surfaceThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccSMB[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "basalThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccBMB[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "calvedThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccCMB[lev],Interval(0,0));
	      comp++;
	    }
	  else if (name[j] == "calvedIceThickness")
	    {
	      data[lev]->copyTo(Interval(j,j),*ccCumCalvedIce[lev],Interval(0,0));
	      comp++;
	    }
	}
      CH_assert(comp == 7);

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

	  
	  const FArrayBox& thck = (*thickness[lev])[dit];
	  const FArrayBox& topg = (*topography[lev])[dit];
	  
	  FArrayBox usrf(topg.box(),1);

	  FORT_SURFACEHEIGHT(CHF_FRA1(usrf,0),
			     CHF_CONST_FRA1(thck,0),
			     CHF_CONST_FRA1(topg,0),
			     CHF_CONST_REAL(a_iceDensity),
			     CHF_CONST_REAL(a_waterDensity),
			     CHF_CONST_REAL(a_seaLevel),
			     CHF_BOX(topg.box()));

	  BaseFab<int>& mask = (*ccMask[lev])[dit];
	  int any = 0;

	  FORT_SETFLOATINGMASK(CHF_FIA1(mask,0),
			       CHF_CONST_FRA1(usrf,0),
			       CHF_CONST_FRA1(topg,0),
			       CHF_CONST_FRA1(thck,0),
			       CHF_INT(any),
			       CHF_REAL(a_iceDensity),
			       CHF_REAL(a_waterDensity),
			       CHF_REAL(a_seaLevel),
			       CHF_BOX(mask.box()));
    

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
		BaseFab<int> mask(grids[dit],1) ;

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
	  src.copy( (*ccSMB[lev])[dit]) ;
	  src.plus( (*ccBMB[lev])[dit]) ;
	  
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

      //work out (dh/dt)_c = smb - bmb + div(uh) and discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  FArrayBox& smb = (*ccSMB[lev])[dit];
    	  FArrayBox& bmb = (*ccBMB[lev])[dit];
    	  FArrayBox& cmb = (*ccCMB[lev])[dit];
    	  FArrayBox& omb = (*ccOMB[lev])[dit];
    	  FArrayBox& dh = (*ccDH[lev])[dit];
    	  BaseFab<int>& mask = (*ccMask[lev])[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      omb(iv) = bmb(iv) + cmb(iv);
	    }

	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& vel = (*fcVel[lev])[dit][dir];
   
	      CH_assert(vel.norm(0) < 1.0e+12);

	      FArrayBox flux(vel.box(),1);
	      flux.copy((*fcVel[lev])[dit][dir]);
	      flux.mult((*fcThck[lev])[dit][dir]);
	     

	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  
		  Real epsThck = 1.0e-10;// TODO fix magic number
		  if ((thck(iv) < epsThck) || (mask(iv) != GROUNDEDMASKVAL && mask(iv) != FLOATINGMASKVAL))
		    {
		      dh(iv) = 0.0;
		      smb(iv) = 0.0;
		      bmb(iv) = 0.0;
		      cmb(iv) = 0.0;
		      if (thck(iv + BASISV(dir)) > epsThck && (mask(iv + BASISV(dir)) == GROUNDEDMASKVAL || mask(iv + BASISV(dir)) == FLOATINGMASKVAL) ) 
			{
			  discharge(iv) += -flux(iv + BASISV(dir)) / dx[lev];
			}
		      if (thck(iv - BASISV(dir)) > epsThck && (mask(iv - BASISV(dir)) == GROUNDEDMASKVAL || mask(iv - BASISV(dir)) == FLOATINGMASKVAL) )
			{
			  discharge(iv) += flux(iv) / dx[lev];
			}
    		    } //end if (thck(iv) < epsThck) 
    		} //end loop over cells
    	    } // end loop over direction
    	} // end loop over grids
    } //end loop over levels

  pout() << " Calving:";
  Real sumDischarge = computeSum(ccDischarge, ratio, dx[0], Interval(0,0), 0);
  pout() << "dischargeFromEdge = " << sumDischarge << " ";
  
  Real sumDH = computeSum(ccDH, ratio, dx[0], Interval(0,0), 0);
  pout() << " dhdt = " << sumDH << " ";

  Real sumSMB = computeSum(ccSMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " smb = " << sumSMB << " ";

  Real sumBMB = computeSum(ccBMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " bmb = " << sumBMB << " ";

  Real sumOMB = computeSum(ccOMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " oceanMelt = " << sumOMB << " ";

  Real sumCMB = computeSum(ccCMB, ratio, dx[0], Interval(0,0), 0);
  pout() << " calvingMelt = " << sumCMB << " ";

  Real sumCumCalvedIce = computeSum(ccCumCalvedIce, ratio, dx[0], Interval(0,0), 0);
  pout() << " cumulatedCalvedIce = " << sumCumCalvedIce << " " << endl;


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
      if (ccDischarge[lev] != NULL)
	{
	  delete ccDischarge[lev];ccDischarge[lev]= NULL;
	}
      if (ccDH[lev] != NULL)
	{
	  delete ccDH[lev];ccDH[lev]= NULL;
	}
      if (ccMask[lev] != NULL)
	{
	  delete ccMask[lev];ccMask[lev]= NULL;
	}
      if (ccCumCalvedIce[lev] != NULL)
	{
	  delete ccCumCalvedIce[lev];ccCumCalvedIce[lev]= NULL;
	}
       if (ccSMB[lev] != NULL)
	{
	  delete ccSMB[lev];ccSMB[lev]= NULL;
	}
      if (ccBMB[lev] != NULL)
	{
	  delete ccBMB[lev];ccBMB[lev]= NULL;
	}
      if (ccOMB[lev] != NULL)
	{
	  delete ccOMB[lev];ccOMB[lev]= NULL;
	}
      if (ccCMB[lev] != NULL)
	{
	  delete ccCMB[lev];ccCMB[lev]= NULL;
	}
      
    }


}
		      
void computeStats(Vector<LevelData<FArrayBox>* >& topography, 
		  Vector<LevelData<FArrayBox>* >& thickness, 
		  Vector<LevelData<FArrayBox>* >& basalThicknessSrc, 
                  Vector<Real>& dx, Vector<int>& ratio,
		  Real iceDensity, Real waterDensity, Real gravity,
                  Real a_seaLevel)
{ 

  CH_TIMERS("computeStats");
  CH_TIMER("createSigmaCS",t1);
  CH_TIMER("integrateH",t2);
  CH_TIMER("integrateHab",t3);
  CH_TIMER("integrateGA",t4);
  CH_TIMER("integrateTA",t5);
  
  int numLevels = topography.size();

  bool has_src = (basalThicknessSrc[0] != NULL);  

  CH_START(t2);
   
  //Compute the total thickness
  // for (int lev=0; lev< numLevels; lev++)
  //  {
  //    coords[lev]->getH().copyTo(Interval(0,0),*tmp[lev],Interval(0,0));
  //  }
  Real iceVolumeAll = computeSum(thickness, ratio, dx[0], Interval(0,0), 0);
  Real totalMelt = 0.0;
  if (has_src) 
    {
      totalMelt = computeSum(basalThicknessSrc, ratio, dx[0], 
                             Interval(0,0), 0);
    }

  CH_STOP(t2);
  Real iceVolumeAbove = 0.0;
  Real groundedArea = 0.0;
  Real groundedPlusOpenLandArea = 0.0;
  Real floatingArea = 0.0;
  Real totalArea = 0.0;

  Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
  
  //Creating a LevelSigmaCS is expensive, so only do it if there is some ice
  if (iceVolumeAll > 1.0e-10)
    {
      CH_START(t1);
      
      IntVect sigmaCSGhost(2*IntVect::Unit);
       
      for (int lev = 0; lev < numLevels; lev++)
	{
	  const DisjointBoxLayout& levelGrids = topography[lev]->disjointBoxLayout();
	   
	  coords[lev] = RefCountedPtr<LevelSigmaCS> 
	    (new LevelSigmaCS(levelGrids, RealVect::Unit*dx[lev], sigmaCSGhost));
	  coords[lev]->setIceDensity(iceDensity);
	  coords[lev]->setWaterDensity(waterDensity);
	  coords[lev]->setGravity(gravity);
	  coords[lev]->setSeaLevel(a_seaLevel);          
	  topography[lev]->copyTo(Interval(0,0),coords[lev]->getTopography(),Interval(0,0));   
	  thickness[lev]->copyTo(Interval(0,0),coords[lev]->getH(),Interval(0,0));
	  coords[lev]->recomputeGeometry(NULL,0);
	   
	   
	}
       
      Vector<LevelData<FArrayBox>* > tmp(numLevels, NULL);
      Vector<LevelData<FArrayBox>* > tmp2(numLevels, NULL);
      Vector<LevelData<FArrayBox>* > tmp3(numLevels, NULL);
      for (int lev=0; lev< numLevels; lev++)
	{
	  const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
	  tmp[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);

	}
      CH_STOP(t1);
       
      CH_START(t3);
       
      {
	//Compute the total thickness above flotation;
	for (int lev=0; lev< numLevels; lev++)
	  {
	    coords[lev]->getThicknessOverFlotation().copyTo(Interval(0,0),*tmp[lev],Interval(0,0));
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
	    //tmp[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
	    tmp2[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
	    tmp3[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
	    for (DataIterator dit(grids);dit.ok();++dit)
	      {
		const BaseFab<int>& mask =  coords[lev]->getFloatingMask()[dit];
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
	    //tmp[lev] = new LevelData<FArrayBox>(grids,1,IntVect::Zero);
	    for (DataIterator dit(grids);dit.ok();++dit)
	      {
		const BaseFab<int>& mask =  coords[lev]->getFloatingMask()[dit];
		const Box& b = grids[dit];
		FArrayBox& a = (*tmp[lev])[dit];
		a.setVal(0.0);
		for (BoxIterator bit(b);bit.ok();++bit)
		  {
		    const IntVect& iv = bit();
		    if ((mask(iv) == GROUNDEDMASKVAL) || (mask(iv) == FLOATINGMASKVAL))
		      {
			a(iv) = 1.0;
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

  // calculate total amount of water in the domain (2 ways)
  Real bedrockBelowSeaLevel = 0;
  Real totalSeaWater = 0;
  Real totalBedrockBelowOcean = 0;
  {
    Vector<LevelData<FArrayBox>* > bedrockBelow(topography.size(), NULL);
    Vector<LevelData<FArrayBox>* > seaWater(topography.size(), NULL);
    Vector<LevelData<FArrayBox>* > bedrockBelowOcean(topography.size(), NULL);
    for (int lev=0; lev<topography.size(); lev++)
      {
        LevelData<FArrayBox>& levelTopo = *(topography[lev]);
        const LevelData<FArrayBox>& levelSurface = coords[lev]->getSurfaceHeight();
        const LevelData<BaseFab<int> >& mask = coords[lev]->getFloatingMask();
        bedrockBelow[lev] = new LevelData<FArrayBox>(levelTopo.getBoxes(),1,
                                                     IntVect::Zero);

        seaWater[lev] = new LevelData<FArrayBox>(levelTopo.getBoxes(),1,
                                                     IntVect::Zero);

        bedrockBelowOcean[lev] = new LevelData<FArrayBox>(levelTopo.getBoxes(),1,
                                                          IntVect::Zero);        
        
        
        DataIterator dit = levelTopo.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            const FArrayBox& thisTopo = levelTopo[dit];
            FArrayBox& thisBedrockBelow = (*bedrockBelow[lev])[dit];
            FArrayBox& thisBedrockBelowOcean = (*bedrockBelowOcean[lev])[dit];
            FArrayBox& thisSeawater = (*seaWater[lev])[dit];
            const FArrayBox& thisThickness = (*thickness[lev])[dit];
            const FArrayBox& thisSurface = levelSurface[dit];
            const BaseFab<int>& thisMask = mask[dit];
            thisBedrockBelow.setVal(0.0);
            thisBedrockBelowOcean.setVal(0.0);
            thisSeawater.setVal(0.0);
            BoxIterator bit = thisBedrockBelow.box();
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (thisTopo(iv,0) < a_seaLevel)
                  {
                    thisBedrockBelow(iv,0) = a_seaLevel - thisTopo(iv,0);
                  }
                if (thisThickness(iv,0) > 0)
                  {
                    // seawater is in cavity below lower surface and bedrock
                    thisSeawater(iv,0) = thisSurface(iv,0) - thisThickness(iv,0) - thisTopo(iv,0);
                  }
                else
                  {
                    // seawater is distance between topography and sealevel
                    thisSeawater(iv,0) = a_seaLevel - thisTopo(iv,0);
                  }
                // look only under ice shelves or open sea
                if ((thisMask(iv,0) == FLOATINGMASKVAL) || (thisMask(iv,0) == OPENSEAMASKVAL))
                  {
                    thisBedrockBelowOcean(iv,0) = a_seaLevel - thisTopo(iv,0);
                  }
              }
          }        
      }
    bedrockBelowSeaLevel = computeSum(bedrockBelow, ratio, dx[0], Interval(0,0), 0);
    totalSeaWater = computeSum(seaWater, ratio, dx[0], Interval(0,0), 0);

    totalBedrockBelowOcean = computeSum(bedrockBelowOcean, ratio, dx[0],
                                        Interval(0,0), 0);
  }
  
  // Approach 1: VaF*density ratio + integral(
  Real totalWaterVolume1 = bedrockBelowSeaLevel + iceDensity/waterDensity*iceVolumeAbove;
  // Approach 2: (total ice volume)*densityRatio + totalSeaWater
  Real totalWaterVolume2 = iceVolumeAll*iceDensity/waterDensity + totalSeaWater;
  
  
  
      
  pout() << " iceVolumeAll = " << iceVolumeAll << " ";
  pout() << " iceVolumeAbove = " << iceVolumeAbove << " ";
  pout() << " groundedArea = " << groundedArea << " ";
  pout() << " floatingArea = " << floatingArea << " ";
  pout() << " totalArea = " << totalArea << " ";
  pout() << " groundedPlusOpenLandArea = " << groundedPlusOpenLandArea << " ";
  pout() << " iceMassAll = " << iceVolumeAll*iceDensity << " ";
  pout() << " iceMassAbove = " << iceVolumeAbove*iceDensity << " ";
  pout() << " bedrockBelowSeaLevel = " << bedrockBelowSeaLevel << " ";
  pout() << " total seawater volume = " << totalSeaWater << " ";
  pout() << " totalWaterVolume = " << totalWaterVolume1 << " ";
  pout() << " totalWaterVolume2 = " << totalWaterVolume2 << " ";
  pout() << " bedrockBelowOcean = " << totalBedrockBelowOcean << "  ";
  
  if (has_src)
    {
      pout() << " Total Melt = " << totalMelt << "." << endl;
    }

  pout() << endl;
  
  
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
    if(argc < 6) 
      { 
	std::cerr << " usage: " << argv[0] << " <plot file> <ice_density> <water_density> <gravity> <sea level> [mask_file] [mask_no_start = 0] [mask_no_end = mask_no_start] " << std::endl; 
	exit(0); 
      }
    char* plotFile = argv[1];
    Real iceDensity = atof(argv[2]);
    Real waterDensity = atof(argv[3]);
    Real gravity = atof(argv[4]);
    Real seaLevel = atof(argv[5]);
    char* maskFile = (argc > 6)?argv[6]:NULL;
    int maskNoStart = 0;
    if (maskFile && argc > 7)
      {
    	maskNoStart = atoi(argv[7]);
      }
    int maskNoEnd = maskNoStart;
    if (maskFile && argc > 8)
      {
    	maskNoEnd = atoi(argv[8]);
      }

    Box domainBox;
    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    
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
	ReadAMRHierarchyHDF5(std::string(maskFile), mgrids, mdata, mname , 
			     mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);
      }
    
    CH_STOP(tp);

    Vector<LevelData<FArrayBox>* > thickness(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > topography(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > basalThicknessSrc(numLevels,NULL);
    int srcComp = -1;
    bool calving = false;
    for (int i=0; i< name.size(); i++)
      {
        if (name[i] == "basalThicknessSource") srcComp=i;
        if (name[i] == "calvedIceThickness") calving=true;
      }

    for (int lev=0;lev<numLevels;++lev)
      {
	thickness[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
	topography[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
        if (srcComp >=0)
          {
            basalThicknessSrc[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
          }
      }

    
    pout().setf(ios_base::scientific,ios_base::floatfield); 
    pout().precision(12);

   
    for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
      {
	createMaskedDEM(topography, thickness, basalThicknessSrc, mdata, 
                        name, data, ratio, dx, mcrseDx, maskNo, maskNo);
	pout() << " time = " << time  ;
	computeStats(topography, thickness, basalThicknessSrc, 
                     dx, ratio, iceDensity, waterDensity,gravity, seaLevel);
	if (maskFile)
	  {
	    pout() << " sector = " << maskNo;
	  }
#ifdef STATS_COMPUTE_DISCHARGE
	computeDischarge(topography, thickness, dx, ratio, name, data, iceDensity, waterDensity,gravity, seaLevel);
	if (calving)
	  {
	    pout() << endl;
	    computeCalving(topography, thickness, dx, ratio, name, data, iceDensity, waterDensity,gravity, seaLevel);
	  }
#endif
	pout() << endl;
      }

    // now compute stats for all selected sectors combined
    if (maskNoStart != maskNoEnd)
      {
	createMaskedDEM(topography, thickness, basalThicknessSrc, mdata, name, data, ratio, dx, mcrseDx, maskNoStart, maskNoEnd);
	pout() << " time = " << time  ;
	computeStats(topography, thickness, basalThicknessSrc, dx, ratio, iceDensity, waterDensity,gravity, seaLevel);
	if (maskFile)
	  {
	    pout() << " all sectors ";
	  }
#ifdef STATS_COMPUTE_DISCHARGE
	computeDischarge(topography, thickness, dx, ratio, name, data, iceDensity, waterDensity,gravity, seaLevel);
	if (calving)
	  {
	    pout() << endl;
	    computeCalving(topography, thickness, dx, ratio, name, data, iceDensity, waterDensity,gravity, seaLevel);
	  }
#endif
	pout() << endl;
      }
    
    for (int lev=0;lev<numLevels;++lev)
      {
	if (thickness[lev] != NULL) delete thickness[lev];
	if (topography[lev] != NULL) delete topography[lev];
	if (basalThicknessSrc[lev] != NULL) delete basalThicknessSrc[lev];
      }

		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
