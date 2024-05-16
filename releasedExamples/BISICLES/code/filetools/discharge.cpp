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
// discharge.cpp
// read in a bisicles plotfile (which must include thickness, Z_base, xVel, 
// yVel, dThickness/dt, surfaceThicknessSource and basalThicknessSource) 
// and a mask, and write out discharge from the ice sheet. 
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
		      Vector<std::string>& name, 
		      Vector<LevelData<FArrayBox>* >& data,
		      Vector<int>& ratio,
		      Vector<Real>& dx,
		      Real mcrseDx)
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
}




void computeDischarge(Vector<LevelData<FArrayBox>* >& topography,
		      Vector<LevelData<FArrayBox>* >& thickness, 
		      Vector<Real>& dx, Vector<int>& ratio, 
		      Vector<std::string>& name, 
		      Vector<LevelData<FArrayBox>* >& data, 
		      Vector<LevelData<FArrayBox>* >& sectorMask, 
		      Real a_iceDensity, Real a_waterDensity, Real a_gravity,
		      Real mcrseDx, 
		      int maskNo)
{

  
  int numLevels = data.size();
  Vector<LevelData<FArrayBox>* > ccVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcVel(numLevels, NULL);
  Vector<LevelData<FluxBox>* > fcThck(numLevels, NULL);
  Vector<LevelData<BaseFab<int> >* > ccMask(numLevels, NULL);
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
      ccMask[lev] = new LevelData<BaseFab<int> >(grids,1,IntVect::Unit);
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
      CH_assert(comp == 4);

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
	  Real seaLevel = 0.0;

	  FORT_SURFACEHEIGHT(CHF_FRA1(usrf,0),
			     CHF_CONST_FRA1(thck,0),
			     CHF_CONST_FRA1(topg,0),
			     CHF_CONST_REAL(a_iceDensity),
			     CHF_CONST_REAL(a_waterDensity),
			     CHF_CONST_REAL(seaLevel),
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
			       CHF_REAL(seaLevel),
			       CHF_BOX(mask.box()));
    

	    for (int dir = 0; dir < SpaceDim; ++dir)
	      {
		Box faceBox = grids[dit];
		faceBox.surroundingNodes(dir);
		FArrayBox& faceVel = (*fcVel[lev])[dit][dir];
	  Box grownFaceBox = faceBox;
	  CH_assert(faceVel.box().contains(grownFaceBox));
	  FArrayBox vface(faceBox,1);

		const FArrayBox& cellVel = (*ccVel[lev])[dit];
		BaseFab<int> mask(grids[dit],1) ;

		FORT_EXTRAPTOMARGIN(CHF_FRA1(faceVel,0),
                                    CHF_FRA1(vface,0),
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

      //work out discharge across boundary. 
      for (DataIterator dit(grids);dit.ok();++dit)
    	{
    	  FArrayBox& discharge = (*ccDischarge[lev])[dit];
    	  discharge.setVal(0.0);
    	  BaseFab<int>& mask = (*ccMask[lev])[dit];
    	  const FArrayBox& thck = (*thickness[lev])[dit];
    	  const Box& b = grids[dit];
	  
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
		  if ( std::abs ( (*sectorMask[lev])[dit](iv) - maskNo) < 1.0e-6 || maskNo == -1)
		    {
		  
		      Real epsThck = 10.0;// TODO fix magic number
		      if ((thck(iv) < epsThck) || (mask(iv) != GROUNDEDMASKVAL))
			{
			  if (thck(iv + BASISV(dir)) > epsThck && (mask(iv + BASISV(dir)) == GROUNDEDMASKVAL) ) 
			    {
			      discharge(iv) += -flux(iv + BASISV(dir)) / dx[lev];
			    }
			  if (thck(iv - BASISV(dir)) > epsThck && (mask(iv - BASISV(dir)) == GROUNDEDMASKVAL) )
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
      if (ccMask[lev] != NULL)
	{
	  delete ccMask[lev];ccMask[lev]= NULL;
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
	std::cerr << " usage: " << argv[0] << " <plot file> <ice_density> <water_density> <gravity> <mask_file> [mask_no_start = 0] [mask_no_end = mask_no_start] " << std::endl; 
	exit(0); 
      }
    char* plotFile = argv[1];
    Real iceDensity = atof(argv[2]);
    Real waterDensity = atof(argv[3]);
    Real gravity = atof(argv[4]);
    char* maskFile = argv[5];
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

    //load the sector mask
    Box mdomainBox;
    Vector<std::string> mname;
    Vector<LevelData<FArrayBox>* > mdata;
    Vector<DisjointBoxLayout > mgrids;
    Vector<int > mratio;
    int mnumLevels;
    Real mdt ,mcrseDx, mtime;

    ReadAMRHierarchyHDF5(std::string(maskFile), mgrids, mdata, mname , 
			     mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);
    
    CH_STOP(tp);

    Vector<LevelData<FArrayBox>* > sectorMask(numLevels,NULL);
    for (int lev=0;lev<numLevels;++lev)
      {
	sectorMask[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
	FillFromReference(*sectorMask[lev], *mdata[0], RealVect::Unit*dx[lev], RealVect::Unit*mcrseDx,true);
      }

    Vector<LevelData<FArrayBox>* > thickness(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > topography(numLevels,NULL);
    Vector<LevelData<FArrayBox>* > basalThicknessSrc(numLevels,NULL);
    int srcComp = -1;
    for (int i=0; i< name.size(); i++)
      {
        if (name[i] == "basalThicknessSource") srcComp=i;
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

    createMaskedDEM(topography, thickness, basalThicknessSrc, name, data, ratio, dx, mcrseDx);
   
    for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
      {
	pout() << " time = " << time  ;
	if (maskFile)
	  {
	    pout() << " sector = " << maskNo ;
	  }
	computeDischarge(topography, thickness, dx, ratio, name, data, sectorMask, iceDensity, waterDensity,gravity,mcrseDx,maskNo);
	pout() << endl;
      }

    // now compute stats for all selected sectors combined
    if (maskNoStart != maskNoEnd)
      {
	pout() << " time = " << time  ;
	if (maskFile)
	  {
	    pout() << " all sectors " << endl;
	  }
	computeDischarge(topography, thickness, dx, ratio, name, data, sectorMask, iceDensity, waterDensity,gravity,mcrseDx,-1);
	pout() << endl;
      }
    
    for (int lev=0;lev<numLevels;++lev)
      {
	if (thickness[lev] != NULL) delete thickness[lev];
	if (topography[lev] != NULL) delete topography[lev];
	if (basalThicknessSrc[lev] != NULL) delete basalThicknessSrc[lev];
	if (sectorMask[lev] != NULL) delete sectorMask[lev];
      }

		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
