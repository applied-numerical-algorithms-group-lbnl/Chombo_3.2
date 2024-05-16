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
// faces.cpp
// read in a bisicles plotfile and write face centered data.
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "JFNKSolver.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFrictionRelation.H"
#include "IceThicknessIBC.H"
#include "BasicThicknessIBC.H"
#include "VieliPayneIBC.H"
#include "MarineIBC.H"
#include "HumpIBC.H"
#include "LevelDataIBC.H"
#include "LevelSigmaCS.H"
#include "IceConstants.H"

inline std::string masktostr(int a_mask)
{
  std::string str;

  if (a_mask == GROUNDEDMASKVAL)
    {
      str="G";
    }
  else if (a_mask == FLOATINGMASKVAL)
    {
      str="F";
    }
  else if (a_mask == OPENSEAMASKVAL)
    {
      str="O";
    }
  else if (a_mask == OPENLANDMASKVAL)
    {
      str="L";
    }
  else
    {
      str="U";
    }
  
  return str;
}

void writeHeader()
{
  CH_assert(SpaceDim == 2);
  pout() << "face level facedir iminus jminus iplus jplus maskminus maskplus x y thk topg ub vb Tnx Tny" << std::endl;
}


//write SpaceDim lines of data centerd at a_iv - dir/2, 0 <= dir < SpaceDim
void writeFaces(const IntVect& a_iv, int a_lev, Real a_dx, 
		int a_faceDir,
		const BaseFab<int>& a_mask,
		const FArrayBox& a_vel, 
		const FArrayBox& a_topg,
		const FArrayBox& a_faceThck,
		const FArrayBox& a_faceVT)
{
  
  const IntVect& ivp = a_iv;
  //for (int fdir = 0; fdir < SpaceDim; ++fdir)
  int fdir = a_faceDir;
    {
      IntVect ivm = ivp-BASISV(fdir);
      pout() << "face " << a_lev << " " << fdir << " ";
      for (int dir = 0; dir < SpaceDim; dir++)
	  pout() << ivm[dir] << " ";
      for (int dir = 0; dir < SpaceDim; dir++)
	  pout() << ivp[dir] << " ";
      
      //MASKA
      //pout() << masktostr(a_mask(ivm)) << " " << masktostr(a_mask(ivp)) << " ";
      pout() << a_mask(ivm) << " " << a_mask(ivp) << " ";
      //(x,y) coords at the face
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  if (dir == fdir)
	    pout() << Real(ivp[dir])*a_dx << " ";
	  else
	    pout() << (Real(ivp[dir]) + 0.5)*a_dx << " ";
	}

      //thickness at face
      pout() << a_faceThck(ivp) << " ";

      //topography at face
      pout() << 0.5*(a_topg(ivp)+a_topg(ivm)) << " ";

      //velocities at face
      for (int dir = 0; dir < SpaceDim; dir++)
	pout() << 0.5*(a_vel(ivp,dir)+a_vel(ivm,dir)) << " ";

      //viscous tensor components at face;
      for (int dir = 0; dir < SpaceDim; dir++)
	pout() << a_faceVT(ivp,dir) << " ";

       pout() << std::endl;
    }

 

}
	       


int main(int argc, char* argv[]) {

#ifdef CH_MPI
#warning "Parallel faces might well not work";
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
   
    if(argc != 3) 
      { 
	std::cerr << " usage: " << argv[0] << " <plot_file> <config_file> " << std::endl; 
	exit(0); 
      }
    char* plot_file = argv[1];
    char* in_file = argv[2];

    ParmParse pp(0,0,NULL,in_file);
    ParmParse pp2("main");

    Box domainBox;
    

    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    Real dt ,crseDx, time;
    ReadAMRHierarchyHDF5(std::string(plot_file), grids, data, name , 
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

    //build valid region masks
    Vector<LevelData<BaseFab<int> >* > mask(numLevels,NULL);
    for (int lev = numLevels - 1; lev >= 0; lev--)
      {
	mask[lev] = new LevelData<BaseFab<int> >(grids[lev],1,IntVect::Unit);
	for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
	  {
	    (*mask[lev])[dit].setVal(1);
	    
	    if (lev < numLevels - 1)
	      {
		for (DataIterator fit(grids[lev+1]); fit.ok(); ++fit)
		  {
		    Box covered = grids[lev+1][fit];
		    covered.coarsen(ratio[lev]);
		    covered &= grids[lev][dit];
		    if (!covered.isEmpty())
		      {
			(*mask[lev])[dit].setVal(0,covered,0,1);
		      }
		  }
	      }
	  }	
      }

    //constitutive relation and rate factor
    ConstitutiveRelation* constRelPtr = ConstitutiveRelation::parse("main");
    if (constRelPtr == NULL)
      {
	MayDay::Error("undefined constitutiveRelation in inputs");
      }

    Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
    {
      ParmParse ppc("constants");
      ppc.query("seconds_per_unit_time",seconds_per_unit_time);
    }
    
    std::string rateFactorType = "constRate";
    pp2.query("rateFactor", rateFactorType);
    RateFactor* rateFactorPtr; 
    if (rateFactorType == "constRate")
      {
	ParmParse crPP("constRate");
	Real A = 9.2e-18  * seconds_per_unit_time/SECONDS_PER_TROPICAL_YEAR;
	crPP.query("A", A);
	ConstantRateFactor* crf = new ConstantRateFactor(A);
	rateFactorPtr = static_cast<RateFactor*>(crf);
	
      }
    else if (rateFactorType == "arrheniusRate")
      {
	ArrheniusRateFactor* arf = new ArrheniusRateFactor(seconds_per_unit_time);
	ParmParse arPP("ArrheniusRate");
	rateFactorPtr = static_cast<RateFactor*>(arf);
      }
    else if (rateFactorType == "patersonRate")
      {
	PatersonRateFactor* prf =  new PatersonRateFactor(seconds_per_unit_time);
	ParmParse arPP("PatersonRate");
	rateFactorPtr = static_cast<RateFactor*>(prf);
      }
    else if (rateFactorType == "zwingerRate")
      {
	ZwingerRateFactor* zrf = new ZwingerRateFactor(seconds_per_unit_time);
	ParmParse arPP("ZwingerRate");
	rateFactorPtr = static_cast<RateFactor*>(zrf);
      }


    BasalFrictionRelation* basalFrictionRelationPtr = BasalFrictionRelation::parse("main",0);
    if (!  basalFrictionRelationPtr )
      {
	MayDay::Error("undefined basalFrictionRelation in inputs");
      }
    
    //we only need this for the velocity boundary condition
    IceThicknessIBC* thicknessIBC;
    std::string problem_type;
    ParmParse geomPP("geometry");
    geomPP.get("problem_type", problem_type);
    if (problem_type == "basic")
      {
        thicknessIBC = new BasicThicknessIBC;
      }
    else if (problem_type == "VieliPayne")
      {
        VieliPayneIBC* ibcPtr = new VieliPayneIBC;     
        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
      }
     else if (problem_type == "marineIceSheet")
       {
        MarineIBC* ibcPtr = new MarineIBC;
	thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
       }
     else if (problem_type == "hump")
       {
         HumpIBC* ibcPtr = new HumpIBC;
         ParmParse humpPP("hump");
         thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
       }
     else if (problem_type == "LevelData")
       {
	 LevelDataIBC* ptr = new LevelDataIBC();
	 thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
       }
     else 
       {
         MayDay::Error("bad problem type");
       }

    ParmParse ppCon("constants");
    Real iceDensity;
    ppCon.query("ice_density",iceDensity);
    Real seaWaterDensity;
    ppCon.query("sea_water_density",seaWaterDensity);
    Real gravity;
    ppCon.query("gravity",gravity);

    ParmParse ppAmr("amr");
    Vector<int> ancells(3);
    ppAmr.getarr("num_cells", ancells, 0, ancells.size());
    
    //periodicity information is not stored in plot files
    bool is_periodic[SpaceDim];
    for (int dir=0; dir<SpaceDim; dir++)
      is_periodic[dir] = false;
    Vector<int> is_periodic_int(SpaceDim, 0);
    ppAmr.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
    for (int dir=0; dir<SpaceDim; dir++) 
    {
      is_periodic[dir] = (is_periodic_int[dir] == 1);
    }
   

    int nLayers = ancells[2];
    Vector<Real> faceSigma(nLayers+1);
    Real dsigma = 1.0 / Real(nLayers);
    for (unsigned int l = 0; l < faceSigma.size(); ++l)
      faceSigma[l] = dsigma * (Real(l));
    ppAmr.queryarr("sigma",faceSigma,0,faceSigma.size());

    //extract the topography, basal friction, thickness and velocity fields
    Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
    Vector<LevelData<FArrayBox>* > C(numLevels,NULL); // basal friction coefficient
    Vector<LevelData<FArrayBox>* > C0(numLevels,NULL); // basal friction coefficien
    Vector<LevelData<FArrayBox>* > vel(numLevels,NULL); // velocity field
    Vector<LevelData<FArrayBox>* > cellA(numLevels,NULL); // cell-centered rate factor
    Vector<LevelData<FluxBox>* > faceA(numLevels,NULL); // face-centered rate factor
    Vector<LevelData<FluxBox>* > faceVT(numLevels,NULL); // face-centered viscous tensor
 
    IntVect sigmaCSGhost(2*IntVect::Unit);

    for (int lev = 0; lev < numLevels; lev++)
      {
	coords[lev] = RefCountedPtr<LevelSigmaCS> 
	  (new LevelSigmaCS(grids[lev], vdx[lev], sigmaCSGhost));
	coords[lev]->setIceDensity(iceDensity);
	coords[lev]->setWaterDensity(seaWaterDensity);
	coords[lev]->setGravity(gravity);
	coords[lev]->setFaceSigma(faceSigma);

	C[lev] = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
	C0[lev] = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
	vel[lev] = new LevelData<FArrayBox>(grids[lev], SpaceDim, IntVect::Unit);
	cellA[lev] = new LevelData<FArrayBox>(grids[lev], faceSigma.size()-1, IntVect::Unit);
	faceA[lev] = new LevelData<FluxBox>(grids[lev], faceSigma.size()-1, IntVect::Zero);
	faceVT[lev] = new LevelData<FluxBox>(grids[lev], SpaceDim, IntVect::Zero);
	if (lev > 0)
	  {
	    coords[lev]->interpFromCoarse(*coords[lev-1],ratio[lev-1]);
	  }
	data[lev]->exchange();
	for (int j = 0; j < name.size(); j++)
	  {
	    if (name[j] == "Z_base")
	      {
		data[lev]->copyTo(Interval(j,j),coords[lev]->getTopography(),Interval(0,0));
	      }
	    else if (name[j] == "thickness")
	      {
		data[lev]->copyTo(Interval(j,j),coords[lev]->getH(),Interval(0,0));
	      }
	    else if (name[j] == "basal_friction")
	      {
		data[lev]->copyTo(Interval(j,j),*C[lev],Interval(0,0));
	      }
	    else if (name[j] == "xVel")
	      {
		data[lev]->copyTo(Interval(j,j),*vel[lev],Interval(0,0));
	      }
	    else if (name[j] == "yVel")
	      {
		data[lev]->copyTo(Interval(j,j),*vel[lev],Interval(1,1));
	      }

	    for (DataIterator dit(grids[lev]);dit.ok();++dit)
	      {
		(*cellA[lev])[dit].setVal(3.1536e-18);
		(*C0[lev])[dit].setVal(0.0);
	      }
	  }

	thicknessIBC->setGeometryBCs(*coords[lev],domain[lev],vdx[lev],time,dt);

	if (lev > 0)
	  {
	    coords[lev]->recomputeGeometry(coords[lev-1],ratio[lev-1]);
	  }
	else
	  {
	    coords[lev]->recomputeGeometry(NULL,0);
	  }
	CellToEdge(*cellA[lev],*faceA[lev]);
      }

    //these parameters don't matter because we don't solve anything. 
    Real vtopSafety = 1.0;
    int vtopRelaxMinIter = 4;
    Real vtopRelaxTol = 1.0;
    Real muMin = 0.0;
    Real muMax = 1.23456789e+300;
    IceNonlinearViscousTensor state(grids, ratio, domain, vdx, coords, vel, C, C0,
		       numLevels-1, *constRelPtr, *basalFrictionRelationPtr, 
		       *thicknessIBC, cellA, faceA, time, 
		       vtopSafety, vtopRelaxMinIter, vtopRelaxTol, 
		       muMin, muMax);
    state.setState(vel);
    state.computeViscousTensorFace(faceVT);



    writeHeader();
    for (int lev = 0; lev < numLevels; lev++)
      {
	
	for (DataIterator dit(grids[lev]);dit.ok();++dit)
	  {
	    //const FArrayBox& thisC = (*C[lev])[dit];

	    const BaseFab<int>& validMask = (*mask[lev])[dit];
	    const BaseFab<int>& floatingMask = coords[lev]->getFloatingMask()[dit];
	    //const FArrayBox& thk = coords[lev]->getH()[dit];
	    const FArrayBox& topg = coords[lev]->getTopography()[dit];
	    const FArrayBox& u = (*vel[lev])[dit];
	    //const Real& rhoi = coords[lev]->iceDensity();
	    //const Real& rhow = coords[lev]->waterDensity();
	    //const Real& gravity = coords[lev]->gravity();
	    //const Real flotationThicknessFactor = -rhow/rhoi;

	    for (int fdir = 0; fdir < SpaceDim; ++fdir)
	      {
		const FArrayBox& fvt =  (*faceVT[lev])[dit][fdir];
		
		const FArrayBox& fthk =  coords[lev]->getFaceH()[dit][fdir];

		Box b = grids[lev][dit];
		//b.grow(fdir,-1);
		b.growHi(fdir,1);
		for (BoxIterator bit(b);bit.ok();++bit)
		  {
		    const IntVect& iv = bit();
		    IntVect ivm = iv - BASISV(fdir);
		    if (validMask(iv) == 1 && validMask(ivm) ==1 )
		      {
			writeFaces(iv, lev , dx[lev], fdir, floatingMask, u, topg, fthk, fvt);
		      }
		  }
	      }
	  }
      }
  
    for (int lev = 0; lev < numLevels; lev++)
      {
	if (C[lev] != NULL)
	  {
	    delete C[lev];
	    C[lev] = NULL;
	  }

	if (C0[lev] != NULL)
	  {
	    delete C0[lev];
	    C0[lev] = NULL;
	  }
	
	if (vel[lev] != NULL)
	  {
	    delete C0[lev];
	    vel[lev] = NULL;
	  }

	if (cellA[lev] != NULL)
	  {
	    delete cellA[lev];
	    cellA[lev] = NULL;
	  }
	
	if (faceA[lev] != NULL)
	  {
	    delete faceA[lev];
	    faceA[lev] = NULL;
	  }
	
	if (faceVT[lev] != NULL)
	  {
	    delete faceVT[lev];
	    faceVT[lev] = NULL;
	  }

	if (mask[lev] != NULL)
	  {
	    delete mask[lev];
	    mask[lev] = NULL;
	  }
      }
		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
