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
// stats.cpp
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
		      Vector<LevelData<FArrayBox>* >& src, 
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

  for (int lev = 0; lev < numLevels; lev++)
    {

      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  (*topography[lev])[dit].setVal(0.0);
	  (*thickness[lev])[dit].setVal(0.0);
          (*src[lev])[dit].setVal(0.0);
	}

      bool activeSMB(false), activeBMB(false); // the 'active' SMB and BMB fields are more useful if available
      
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
	  else if (name[j] == "activeBasalThicknessSource")
	    {
	      activeBMB = true;
	      data[lev]->copyTo(Interval(j,j),*src[lev],Interval(0,0));
	    }
	  else if ((!activeBMB) && (name[j] == "basalThicknessSource"))
	    {
	      data[lev]->copyTo(Interval(j,j),*src[lev],Interval(0,0));
	    }
	  else if (name[j] == "activeSurfaceThicknessSource")
	    {
	      activeSMB = true;
	      data[lev]->copyTo(Interval(j,j),*src[lev],Interval(1,1));
	    }
	   else if ((!activeSMB) && (name[j] == "surfaceThicknessSource"))
	     {
	      data[lev]->copyTo(Interval(j,j),*src[lev],Interval(1,1));
	    }
	  else if (name[j] == "calvingFlux")
	    {
	      data[lev]->copyTo(Interval(j,j),*src[lev],Interval(2,2));
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
	      for (int i = 0; i < src[lev]->nComp(); i++)
		{
		  (*src[lev])[dit].mult(coef,0,i,1);
                }
		
	    }
	    
	}
    }
  
}


 
		      
void computeStats(Vector<LevelData<FArrayBox>* >& topography, 
		  Vector<LevelData<FArrayBox>* >& thickness, 
		  Vector<LevelData<FArrayBox>* >& src, 
                  Vector<Real>& dx, Vector<int>& ratio,
		  Real iceDensity, Real waterDensity, Real gravity)
{ 

  CH_TIMERS("computeStats");
  CH_TIMER("createSigmaCS",t1);
  CH_TIMER("integrateH",t2);
  CH_TIMER("integrateHab",t3);
  CH_TIMER("integrateGA",t4);
  CH_TIMER("integrateTA",t5);
  
  int numLevels = topography.size();

  // whole domain integrals
  CH_START(t2);
  Real iceVolumeAll = computeSum(thickness, ratio, dx[0], Interval(0,0), 0);
  Real totalBMB = computeSum(src, ratio, dx[0], Interval(0,0), 0);
  Real totalSMB = computeSum(src, ratio, dx[0], Interval(1,1), 0);
  Real totalCF = computeSum(src, ratio, dx[0], Interval(2,2), 0);
  CH_STOP(t2);

  // sub-domain integrals
  Real iceVolumeAbove(0.0), groundedArea(0.0), groundedSMB(0.0), groundedBMB(0.0);
  Real floatingArea(0.0), groundedPlusOpenLandArea(0.0);

  // need to create a LevelSigmaCS and compute mask
  CH_START(t1);
  Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
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
  //  CH_STOP(t1);
  
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
  
  //grounded and floating area, grounded SMB and BMB.
  // overwrites src (since totals have already been computed)
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
	    } // end for (BoxIterator bit(b);bit.ok();++bit)
	  (*src[lev])[dit].mult(a, 0, 0, 1); // mask BMB
	  (*src[lev])[dit].mult(a, 0, 1, 1); // mask SMB
	} // end for (DataIterator dit(grids);dit.ok();++dit)
    } // end   for (int lev=0; lev< numLevels; lev++)

  groundedArea = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
  floatingArea = computeSum(tmp2, ratio, dx[0], Interval(0,0), 0);
  groundedPlusOpenLandArea = computeSum(tmp3, ratio, dx[0], Interval(0,0), 0);
  groundedBMB = computeSum(src, ratio, dx[0], Interval(0,0), 0);
  groundedSMB = computeSum(src, ratio, dx[0], Interval(1,1), 0);
  CH_STOP(t4);

  // total area; overwrites tmp
  CH_START(t5);
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
  Real totalArea = computeSum(tmp, ratio, dx[0], Interval(0,0), 0);
  CH_STOP(t5);

  // clean up temporaries
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
  
    
  pout() << " iceVolumeAll = " << iceVolumeAll << " ";
  pout() << " iceVolumeAbove = " << iceVolumeAbove << " ";
  pout() << " groundedArea = " << groundedArea << " ";
  pout() << " floatingArea = " << floatingArea << " ";
  pout() << " totalArea = " << totalArea << " ";
  pout() << " groundedPlusOpenLandArea = " << groundedPlusOpenLandArea << " ";
  pout() << " iceMassAll = " << iceVolumeAll*iceDensity << " ";
  pout() << " iceMassAbove = " << iceVolumeAbove*iceDensity << " ";
  pout() << " totalBMBIn = " << totalBMB << " ";
  pout() << " totalSMBIn = " << totalSMB << " ";
  pout() << " totalCalvingFluxOut = " << totalCF << " ";
  pout() << " groundedBMBIn = " << groundedBMB << " " ;
  pout() << " groundedSMBIn = " << groundedSMB << " " ;
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
    Vector<LevelData<FArrayBox>* > src(numLevels,NULL);
  
    for (int lev=0;lev<numLevels;++lev)
      {
	thickness[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
	topography[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
	src[lev] = new LevelData<FArrayBox>(grids[lev],3,4*IntVect::Unit); //BMB, SMB, Calving Flux
      }

    
    pout().setf(ios_base::scientific,ios_base::floatfield); 
    pout().precision(12);

   
    for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
      {
	createMaskedDEM(topography, thickness, src, mdata, 
                        name, data, ratio, dx, mcrseDx, maskNo, maskNo);
	pout() << " time = " << time  ;
	computeStats(topography, thickness, src, 
                     dx, ratio, iceDensity, waterDensity,gravity);
	if (maskFile)
	  {
	    pout() << " sector = " << maskNo;
	  }
	pout() << endl;
      }

    // now compute stats for all selected sectors combined
    if (maskNoStart != maskNoEnd)
      {
	createMaskedDEM(topography, thickness, src, mdata, name, data, ratio, dx, mcrseDx, maskNoStart, maskNoEnd);
	pout() << " time = " << time  ;
	computeStats(topography, thickness, src, dx, ratio, iceDensity, waterDensity,gravity);
	if (maskFile)
	  {
	    pout() << " all sectors ";
	  }
      }
    
    for (int lev=0;lev<numLevels;++lev)
      {
	if (thickness[lev] != NULL) delete thickness[lev];
	if (topography[lev] != NULL) delete topography[lev];
	if (src[lev] != NULL) delete src[lev];
      }

		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
