//===========================================================================
// driver.cpp
//
//===========================================================================
#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"

#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "IceConstants.H"
#include "BasalFriction.H"
#include "BasalFrictionRelation.H"
#include "twistyStreamFriction.H"
#include "CopyBasalFriction.H"
#include "IceThicknessIBC.H"
#include "BasicThicknessIBC.H"
#include "VieliPayneIBC.H"
#include "MarineIBC.H"
#include "FortranInterfaceIBC.H"
#include "LevelDataIBC.H"
#include "PiecewiseLinearFillPatch.H"
#include "FineInterp.H"
#include "IceVelocitySolver.H"
//#include "PicardSolver.H"
#include "JFNKSolver.H"
#include "CoarseAverageFace.H"
#include "IceUtility.H"
#include "LevelSigmaCS.H"
#ifdef CH_USE_FAS
#include "FASIceSolver.H"
#endif

//#include "PetscSolver.H"
#include "PetscIceSolver.H"

#include "AMRIO.H"

/// tolerance for Real-Real comparison
Real realEps = 1.0e-10;

/// types of basal friction (beta) distributions
/** SinusoidalBeta is the one for exp C in Pattyn et al (2008)
 */
enum basalFrictionTypes {constantBeta = 0,
                         sinusoidalBeta,
                         sinusoidalBetay,
                         twistyStreamx,
                         NUM_BETA_TYPES};


enum velSolverTypes { Picard = 0,
                      JFNK = 1,
                      KnownVelocity = 2,
                      PetscNLSolver = 3,
                      FASMGAMR = 4,
                      NUM_SOLVER_TYPES};



//===========================================================================
// standalone ice-sheet velocity solver (for benchmarking)
//
//===========================================================================
int main(int argc, char* argv[]) {

  int ierr = 0;
  

#ifdef CH_USE_PETSC
 ierr = PetscInitialize(&argc, &argv,"./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif // end petsc conditional

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
    

    if(argc < 2) 
      { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);

    // this seems a bit silly, but the goal is to preserve 
    // the structure of the inputs files used in main code
    ParmParse pp2("main");
    ParmParse ppAmr("amr");
    ParmParse geomPP("geometry");
    


    int verbosity = 1;
    ppAmr.query("verbosity", verbosity);
    
    Real   velocity_solver_tolerance = 1e-10;
    ppAmr.query("velocity_solver_tolerance", velocity_solver_tolerance);

    // default is -1, which means use the solver's own defaults
    int maxSolverIterations = -1;
    ppAmr.query("max_solver_iterations",maxSolverIterations);

    string plotName;
    pp2.get("filename", plotName);

    // the combination of setVelToZero = false and interpFinestLevel = true
    // is to use the plotfile velocity field as the initial guess, but to 
    // replace the finest-level velocity field with one interpolated from 
    // the next coarser level. This is designed to (a) give the solver 
    // something to do, while still having a reasonable initial guess, and 
    // (b) mimic what happens after a regridding operation.

    bool setVelToZero = false;
    pp2.query("setVelToZero", setVelToZero);

    bool interpFinestLevel = true;
    pp2.query("interpFinestLevel", interpFinestLevel);

    bool doPlots = true;
    pp2.query("writePlotfile", doPlots);

    // do things this way so that by default we use the plotfile 
    // prefix used in the original run, but we can over-ride it with a 
    // benchmark-specific one.
    string plotRoot = "solverData";
    ppAmr.query("plot_prefix", plotRoot);
    pp2.query("plot_prefix", plotRoot);

    // get problem data from plotfile
    
    Vector<DisjointBoxLayout> vectGrids;
    Vector<LevelData<FArrayBox>* >  vectData;
    Vector<string> vectNames;
    Box domainBox;
    Real dx;
    Real dt;
    Real time;
    Vector<int> refRatio;
    int numLevels;
    
    int status;
    status = ReadAMRHierarchyHDF5(plotName,
                                  vectGrids,
                                  vectData,
                                  vectNames,
                                  domainBox,
                                  dx,
                                  dt,
                                  time,
                                  refRatio,
                                  numLevels);
    
    int finest_level = numLevels -1;    

    // write out number of cells per level in 
    if (verbosity > 2)
      {
        Box levelDomain(domainBox);
        
        pout() << "Number of cells per level:" << endl;
        for (int lev=0; lev < numLevels; lev++)
          {
            int numCells = 0;
            int numBoxes = 0;
            DataIterator levelDit = vectGrids[lev].dataIterator();
            for (levelDit.begin(); levelDit.ok(); ++levelDit)
              {
                numCells += vectGrids[lev][levelDit].numPts();
                ++numBoxes;
              }
            int domainCells = levelDomain.numPts();
            Real fillRatio = 100.0*numCells/domainCells;
            pout() << "     level " << lev << ": " 
                   << numBoxes << " boxes, " 
                   << numCells << " cells, "
                   << fillRatio << ", percent of domain" << endl;
            levelDomain.refine(refRatio[lev]);
          }
      }

    int max_level = finest_level;
    ppAmr.query("maxLevel", max_level);

    if (max_level > finest_level) max_level = finest_level;
    if (numLevels > (max_level +1)) numLevels = max_level +1;

    RealVect domainSize;
    Vector<Real> domSize(SpaceDim);
    pp2.getarr("domain_size", domSize, 0, SpaceDim);
    domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));

    // periodicity info required in inputs file
    bool is_periodic[SpaceDim];
    {
      Vector<int> is_periodic_int(SpaceDim, 0);
      ppAmr.getarr("is_periodic", is_periodic_int, 0, SpaceDim);

      for (int dir=0; dir<SpaceDim; dir++) 
        {
          is_periodic[dir] = (is_periodic_int[dir] == 1);
        }
    }

    // set up problem definition info from plotfile    
    Vector<ProblemDomain> amrDomains(numLevels);
    amrDomains[0] = ProblemDomain(domainBox, is_periodic);
    
    Vector<Real> amrDx(numLevels);
    RealVect computedDx;
    for (int dir=0; dir<SpaceDim; dir++)
      {
        computedDx[dir] = domSize[dir]/domainBox.size(dir);
      }

    // reality check
    if (abs(computedDx[0] - dx)/dx > realEps)
      {
        MayDay::Warning("dx in plotfile incompatible with computed dx");
        pout() << "computed dx = " << computedDx[0] 
               << ", dx in plotfile = " << dx << endl;        
      }

    amrDx[0] = computedDx[0];

    for (int lev=1; lev<numLevels; lev++)
      {
        amrDomains[lev] = amrDomains[lev-1];
        amrDomains[lev].refine(refRatio[lev-1]);
        
        amrDx[lev] = amrDx[lev-1]/refRatio[lev-1];
      }


    // need to incorporate periodicity info into DisjointBoxLayouts.
    // this gets a bit odd, but unfortunately periodicity info isn't stored
    // in plotfiles.
    // while we're here, extract what we need from plotData and copy it 
    Vector<DisjointBoxLayout> amrGrids(numLevels);
    Vector<RefCountedPtr<LevelSigmaCS > > vectCoordSys(numLevels);
    Vector<LevelData<FArrayBox>* > velocity(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > vectRhs(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > vectBeta(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > vectC0(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > vectA(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > vectSurface(numLevels, NULL); // L1L2 will need grad(S)
    Vector<LevelData<FArrayBox>* > calvedIceThickness(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > removedIceThickness(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > addedIceThickness(numLevels, NULL);
    Vector<LevelData<FArrayBox>* > vectMuCoef(max_level + 1,NULL);

    Real iceDensity = 910.0;
    Real seaWaterDensity = 1028.0;
    Real gravity = 9.81;

    // need to set things like ice density and gravity
    ParmParse ppCon("constants");
    ppCon.query("ice_density",iceDensity);
    ppCon.query("sea_water_density",seaWaterDensity);
    ppCon.query("gravity",gravity);        

    for (int lev=0; lev<numLevels; lev++)
      {
        const Vector<Box>& levelBoxes = vectGrids[lev].boxArray();
        const Vector<int>& levelProcs = vectGrids[lev].procIDs();
        
        amrGrids[lev].define(levelBoxes, levelProcs, amrDomains[lev]);

        // now define coordinate system objects
        IntVect sigmaCSGhost = 4*IntVect::Unit;
        RealVect levelDx = amrDx[lev]*RealVect::Unit;
        
        vectCoordSys[lev] = RefCountedPtr<LevelSigmaCS >
	  (new LevelSigmaCS(amrGrids[lev], levelDx, sigmaCSGhost));
        
        LevelSigmaCS& levelCS = *vectCoordSys[lev];

        // set physical constants
        levelCS.setIceDensity(iceDensity);
        levelCS.setWaterDensity(seaWaterDensity);
        levelCS.setGravity(gravity);

        velocity[lev] = new LevelData<FArrayBox>(amrGrids[lev],SpaceDim,
                                                 IntVect::Unit);
        vectRhs[lev] = new LevelData<FArrayBox>(amrGrids[lev],SpaceDim,
                                                 IntVect::Zero);
        vectBeta[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
                                                 IntVect::Unit);


        vectC0[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
                                                 IntVect::Unit);

        // default value for C0 is, fittingly enough, 0
        {
          DataIterator dit = amrGrids[lev].dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              (*vectC0[lev])[dit].setVal(0.0);
            }
        }
        calvedIceThickness[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
                                                 IntVect::Unit);
        removedIceThickness[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
                                                 IntVect::Unit);
        addedIceThickness[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
							  IntVect::Unit);
        {
          DataIterator dit = amrGrids[lev].dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              (*calvedIceThickness[lev])[dit].setVal(0.0);
              (*removedIceThickness[lev])[dit].setVal(0.0);
              (*addedIceThickness[lev])[dit].setVal(0.0);
            }
        }

	// default muCoef is 1.0
	vectMuCoef[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
						   IntVect::Unit);

	{
          DataIterator dit = amrGrids[lev].dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              (*vectMuCoef[lev])[dit].setVal(1.0);
	    }
	}

	vectSurface[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
						    2*IntVect::Unit);

        // temporaries for extracting data from plotData
        LevelData<FArrayBox>& levelThickness = levelCS.getH();
        LevelData<FArrayBox> levelZBase(amrGrids[lev], 1, 
                                        levelCS.getTopography().ghostVect());
        LevelData<FArrayBox> levelGrad(amrGrids[lev],SpaceDim,IntVect::Unit);
        LevelData<FluxBox> levelFaceGrad(amrGrids[lev],SpaceDim,IntVect::Zero);
	

        // copy data from plotData to new storage:
        LevelData<FArrayBox>& levelPlotData = *vectData[lev];

        // would have liked to use copyTo, but have to do this a 
        // bit differently to catch ghost cells
        DataIterator ditSrc = levelPlotData.dataIterator();
        DataIterator ditDest = amrGrids[lev].dataIterator();
        for (ditSrc.begin(), ditDest.begin(); ditSrc.ok(); ++ditDest, ++ditSrc)
          {
            const FArrayBox& plotDataFab = levelPlotData[ditSrc];
            // sanity check
            CH_assert(vectGrids[lev][ditSrc] == amrGrids[lev][ditDest]);
            
            // now loop over components and copy
            for (int comp=0; comp<levelPlotData.nComp(); comp++)
              {
                if (vectNames[comp] == "thickness")
                  {
                    levelThickness[ditDest].copy(plotDataFab, 
                                                 comp, 0, 1);
                  }
		if (vectNames[comp] == "Z_surface")
                  {
                    (*vectSurface[lev])[ditDest].copy(plotDataFab, 
                                                 comp, 0, 1);
                  }	
                else if (vectNames[comp] == "xVel")
                  {
                    (*velocity[lev])[ditDest].copy(plotDataFab,
                                                   comp, 0, 
                                                   SpaceDim);
                  }
                else if (vectNames[comp] == "Z_base")
                  {
                    levelZBase[ditDest].copy(plotDataFab,
                                             comp, 0, 1);
                  }
                else if (vectNames[comp] == "basal_friction")
                  {
                    (*vectBeta[lev])[ditDest].copy(plotDataFab,
                                                   comp, 0, 1);
                  }
                else if (vectNames[comp] == "C0")
                  {
                    (*vectC0[lev])[ditDest].copy(plotDataFab,
                                                 comp, 0, 1);
                  }                    
                else if (vectNames[comp] == "xRhs")
                  {
                    (*vectRhs[lev])[ditDest].copy(plotDataFab,
                                                  comp, 0, SpaceDim);
                  }
                else if (vectNames[comp] == "muCoef")
                  {
                    (*vectMuCoef[lev])[ditDest].copy(plotDataFab,
						     comp, 0, SpaceDim);
                  }		
              } // end loop over plot components
          } // end loop over boxes

        // handle background slope correctly -- set in LevelSigmaCS
        // and remove from topography 
        if (geomPP.contains("basalSlope") )
            {
              RealVect basalSlope;
              Vector<Real> basalSlopePP(SpaceDim, 0.0);
              geomPP.getarr("basalSlope", basalSlopePP, 0, SpaceDim);
              D_TERM(
                     basalSlope[0] = basalSlopePP[0];,
                     basalSlope[1] = basalSlopePP[1];,
                     basalSlope[2] = basalSlopePP[2];);

              levelCS.setBackgroundSlope(basalSlope);

              DataIterator dit = levelZBase.dataIterator();
              for (dit.begin(); dit.ok(); ++dit)
                {
                  FArrayBox& thisTopography = levelZBase[dit];
                  Real backgroundBase;
                  BoxIterator bit(thisTopography.box());
                  for (bit.begin(); bit.ok(); ++bit)
                    {
                      IntVect iv = bit();
                      RealVect loc(iv);
                      loc += 0.5*RealVect::Unit;
                      loc *= levelDx;
                      
                      backgroundBase = D_TERM(loc[0]*basalSlope[0],
                                              +loc[1]*basalSlope[1],
                                              +loc[2]*basalSlope[2]);
                      
                      thisTopography(iv,0) -= backgroundBase;
                    }
                      
                } // end loop over boxes in topography
            } // end if there is a background slope

        levelZBase.exchange();

        levelCS.setTopography(levelZBase);

        LevelSigmaCS* crseCSPtr = NULL;
        int nRefCrse = -1;
        if (lev > 0) 
          {
            crseCSPtr = &(*vectCoordSys[lev-1]);
            nRefCrse = refRatio[lev-1];

            // first, fill in thickness ghost cells
            PiecewiseLinearFillPatch  thicknessFiller
              (amrGrids[lev],amrGrids[lev-1],1, amrDomains[lev-1],
               nRefCrse, levelThickness.ghostVect()[0] );
            const LevelData<FArrayBox>& crseH = crseCSPtr->getH();
            thicknessFiller.fillInterp(levelThickness,crseH,crseH,
                                       0.0, 0, 0, 1);


            levelCS.interpFromCoarse(*crseCSPtr, nRefCrse,
                                     false, false, false);          
          }
        
        levelCS.recomputeGeometry(crseCSPtr, nRefCrse);
            

	//need to fill two ghost cells of vectSurface[lev];
        LevelData<FArrayBox>& levelZs = *vectSurface[lev];
	if (lev > 0)
          {
            int nGhost = levelZs.ghostVect()[0];
            PiecewiseLinearFillPatch surfaceFiller(amrGrids[lev],
                                                   amrGrids[lev-1],
                                                   1, 
                                                   amrDomains[lev-1],
                                                   refRatio[lev-1],
                                                   nGhost);
            
            // since we're not subcycling, don't need to interpolate in time
            Real time_interp_coeff = 0.0;
            surfaceFiller.fillInterp(levelZs, 
                                     *vectSurface[lev-1],
                                     *vectSurface[lev-1],
                                     time_interp_coeff,
                                     0, 0, 1);
            
            
          }
      
	levelZs.exchange();

	
        DataIterator dit = amrGrids[lev].dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            Box grownBox = amrGrids[lev][dit];
            grownBox.grow(sigmaCSGhost);
           
	    FArrayBox& thisZsurf = (*vectSurface[lev])[dit];

	    FArrayBox& grad = levelGrad[dit];
	    for (BoxIterator bit(grad.box()); bit.ok(); ++bit)
	      {
		IntVect iv = bit();
		for (int dir=0; dir<SpaceDim; dir++)
		  {
		    Real oneOnTwoDx = 0.5 / levelDx[dir];
		    
		    IntVect offsetVect = BASISV(dir);
		    IntVect hi = iv + offsetVect;
		    IntVect lo = iv - offsetVect;
		    
		    grad(iv,dir) = (thisZsurf(hi,0) - thisZsurf(lo,0)) * oneOnTwoDx;	   
		    
		  } // end loop over directions
	      } // end loop over cells
	    
	    
	    FluxBox& faceGrad = levelFaceGrad[dit];
            for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
	      {
		FArrayBox& thisFaceGrad = faceGrad[faceDir];
		const Box& faceBox = thisFaceGrad.box();
		BoxIterator fbit(faceBox);
		
		Real oneOnDx = 1.0 / levelDx[faceDir] ;
		
		for (fbit.begin(); fbit.ok(); ++fbit)
		  {
		    IntVect iv = fbit(); 
		    thisFaceGrad(iv, faceDir) =  
		      oneOnDx * (thisZsurf(iv,0) 
				 - thisZsurf(iv - BASISV(faceDir),0));
		  }

		if (SpaceDim == 2)
		  {
		    int transDir = (faceDir+1)%SpaceDim;
		    Real oneOnFourDx = 0.25 / levelDx[faceDir];
		    for (fbit.begin(); fbit.ok(); ++fbit)
		      {
			IntVect iv = fbit(); 
			thisFaceGrad(iv, transDir) = 
			  oneOnFourDx
			  * (thisZsurf(iv + BASISV(transDir),0) 
			     - thisZsurf(iv - BASISV(transDir),0)
			     + thisZsurf(iv + BASISV(transDir) 
					 - BASISV(faceDir),0) 
			     - thisZsurf(iv - BASISV(transDir) 
					 - BASISV(faceDir),0));
		      }
		  } 
	      }
	    
 
	    
        
          } // end loop over grids
        //levelCS.setGradSurface(levelGrad); 
        //levelCS.setGradSurfaceFace(levelFaceGrad);
        
      } // end loop over levels
    
  
 
    // ---------------------------------------------
    // set constitutive relation & rate factor 
    // (this info isn't in the plotfile)
    // ---------------------------------------------
    ConstitutiveRelation* constRelPtr = ConstitutiveRelation::parse("main");

    if (constRelPtr == NULL)
      {
	MayDay::Error("undefined constitutiveRelation in inputs");
      }

   
     
#if BISICLES_Z == BISICLES_LAYERED
    // num layers is the third component of num_cells in inputs file
    Vector<int> ancells(3);
    ppAmr.getarr("num_cells", ancells, 0, ancells.size());
    int nLayers = ancells[2];
    
    Vector<Real> faceSigma(nLayers+1);
    Real dsigma = 1.0 / Real(nLayers);
    for (unsigned int l = 0; l < faceSigma.size(); ++l)
      faceSigma[l] = dsigma * (Real(l));
    ppAmr.queryarr("sigma",faceSigma,0,faceSigma.size());

    for (int lev=0; lev<numLevels; lev++)
      {
	(*vectCoordSys[lev]).setFaceSigma(faceSigma);
	vectA[lev] = new LevelData<FArrayBox>(amrGrids[lev],nLayers,IntVect::Unit);
      }
#else
    for (int lev=0; lev<numLevels; lev++)
      {
	vectA[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,IntVect::Unit);
      }
#endif
    
    std::string rateFactorType = "constRate";
    pp2.query("rateFactor", rateFactorType);
    if (rateFactorType == "constRate")
      {
	Real constA = 9.2e-18;
	ParmParse crPP("constRate");
	crPP.query("A", constA);
	ConstantRateFactor rateFactor(constA);
	Real epsSqr0 = 1.0e-9;
	crPP.query("epsSqr0", epsSqr0);
 
	for (int lev=0; lev<numLevels; lev++)
	  {
	    for (DataIterator dit(amrGrids[lev]); dit.ok(); ++dit)
	      {
		(*vectA[lev])[dit].setVal(constA);	
	      }
	  }
      }
    else if (rateFactorType == "arrheniusRate")
      {
	ArrheniusRateFactor rateFactor(SECONDS_PER_TROPICAL_YEAR);
	ParmParse arPP("ArrheniusRate");
	Real epsSqr0 = 1.0e-9;
	arPP.query("epsSqr0", epsSqr0);

	//if we have an arrenhius rate, then we need to read temperature from
	//plot file and compute A.
	for (int lev=0; lev<numLevels; lev++)
	  {
	    const LevelSigmaCS& levelCS = *vectCoordSys[lev];
	    const LevelData<FArrayBox>& levelPlotData = *vectData[lev];
	    LevelData<FArrayBox> levelTemp(amrGrids[lev],levelCS.getSigma().size(),IntVect::Unit);
	    DataIterator ditSrc = levelPlotData.dataIterator();
	    DataIterator ditDest = amrGrids[lev].dataIterator();
	    for (ditSrc.begin(), ditDest.begin(); ditSrc.ok(); ++ditDest, ++ditSrc)
	      {
		const FArrayBox& plotDataFab = levelPlotData[ditSrc];
		// sanity check
		CH_assert(vectGrids[lev][ditSrc] == amrGrids[lev][ditDest]);
            
		// now loop over components and copy
		for (int comp=0; comp<levelPlotData.nComp(); comp++)
		  {
		    if (vectNames[comp] == "temperature0000")
		      {
			levelTemp[ditDest].copy(plotDataFab, comp, 0, levelCS.getSigma().size());
		      }
		  }
		
		IceUtility::computeA(*vectA[lev],levelCS.getSigma(),levelCS,&rateFactor,levelTemp);
	      }
	  }
      }
    // ---------------------------------------------
    // set basal friction coefficient and relation
    // ---------------------------------------------

    BasalFriction* basalFrictionPtr = NULL;

    std::string beta_type;
    geomPP.get("beta_type", beta_type);
    // read in type of beta^2 distribution
    
    if (beta_type == "constantBeta")
      {
        Real betaVal;
        geomPP.get("betaValue", betaVal);
        basalFrictionPtr = static_cast<BasalFriction*>(new constantFriction(betaVal));
      }
    else if (beta_type == "sinusoidalBeta")
      {
        Real betaVal, eps;
        RealVect omega(RealVect::Unit);
        Vector<Real> omegaVect(SpaceDim);
        geomPP.get("betaValue", betaVal);
        if (geomPP.contains("omega"))
          {
            geomPP.getarr("omega", omegaVect, 0, SpaceDim);
            omega = RealVect(D_DECL(omegaVect[0], omegaVect[1], omegaVect[2]));
          }
        geomPP.get("betaEps", eps);
        basalFrictionPtr = static_cast<BasalFriction*>(new sinusoidalFriction(betaVal, 
                                                                              omega, 
                                                                              eps,
                                                                              domainSize));
      }
    // keep this one around for backward compatibility, even if it
    // is a special case of sinusoidalBeta
    else if (beta_type == "sinusoidalBetay")
      {
        Real betaVal, eps, omegaVal;
        RealVect omega(RealVect::Zero);
        omega[1] = 1;
        
        geomPP.get("betaValue", betaVal);
        if (geomPP.contains("omega"))
          {
            geomPP.get("omega", omegaVal);
            omega[1] = omegaVal;
          }
        geomPP.get("betaEps", eps);
        basalFrictionPtr = static_cast<BasalFriction*>(new sinusoidalFriction(betaVal, 
                                                                              omega, 
                                                                              eps,
                                                                              domainSize));

        }
    else if (beta_type == "twistyStreamx")
      {
        Real betaVal, eps, magOffset;
        magOffset = 0.25;
        RealVect omega(RealVect::Unit);
        Vector<Real> omegaVect(SpaceDim);
        geomPP.get("betaValue", betaVal);
        if (geomPP.contains("omega"))
          {
            geomPP.getarr("omega", omegaVect, 0, SpaceDim);
            omega = RealVect(D_DECL(omegaVect[0], omegaVect[1], omegaVect[2]));
          }
        geomPP.query("magOffset", magOffset);
        geomPP.get("betaEps", eps);
        basalFrictionPtr = static_cast<BasalFriction*>(new twistyStreamFriction(betaVal, 
                                                                                omega, 
                                                                                magOffset, 
                                                                                eps,
                                                                                domainSize));          
      }
    else if (beta_type == "fortran" || beta_type == "LevelData")    
      {
        // this works a bit differently than the
        // fortran interface / LevelData basal friction in the main 
        // code. We will just read this in from the plotfile
        
        Vector<RefCountedPtr<LevelData<FArrayBox> > > vectBeta(numLevels);
        
        int betaComp = -1;
        for (int n = 0; n<vectNames.size(); n++)
          {
            if (vectNames[n] == "basal_friction") 
              {
                betaComp = n;
              }
          }
        
        if (betaComp < 0)
          {
            MayDay::Error("basal_friction not found in hdf5 file");
          }
        
        for (int lev=0; lev<numLevels; lev++)
          {
            const DisjointBoxLayout& levelGrids = amrGrids[lev];
            RefCountedPtr<LevelData<FArrayBox> > levelPtr(new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit) );
            vectBeta[lev] = levelPtr;            

            LevelData<FArrayBox>& levelBeta = *levelPtr;
            LevelData<FArrayBox>& levelPlotData = *vectData[lev];

            if (levelBeta.getBoxes() == levelPlotData.getBoxes())
              {
                // now copy from vectData -> vectBeta.
                // do a fab-by-fab copy rather than a copyTo
                DataIterator dit = levelGrids.dataIterator();
                for (dit.begin(); dit.ok(); ++dit)
                  {
                    FArrayBox& thisBeta = levelBeta[dit];
                    FArrayBox& thisPlotData = levelPlotData[dit];
                    thisBeta.copy(thisPlotData,betaComp, 0, 1);
                    //levelBeta[dit].copy(levelPlotData[dit],betaComp, 0, 1);
                  }
              } 
            else
              {
                // have to use copyTo here
                Interval srcComps(betaComp, betaComp);
                Interval destComps(0,0);
                levelPlotData.copyTo(srcComps, levelBeta, destComps);
              }
          }
               
        basalFrictionPtr = static_cast<BasalFriction*>(new CopyBasalFriction(vectBeta, amrDomains) );
        
      }
    else 
      {
        MayDay::Error("undefined beta_type in inputs");
      }

    BasalFrictionRelation* basalFrictionRelationPtr = BasalFrictionRelation::parse("main",0);
    if (!  basalFrictionRelationPtr )
      {
	MayDay::Error("undefined basalFrictionRelation in inputs");
      }

    // ---------------------------------------------
    // set IBC -- this includes initial ice thickness, 
    // and basal geometry
    // ---------------------------------------------

    
    IceThicknessIBC* thicknessIBCPtr;

    // this will contain the boxes in the index space of the 
    // original data in which the thickness will be cleared.
    Vector<Box> clearBoxes;
    int clearBoxesLevel = -1;



    std::string problem_type;
    geomPP.get("problem_type", problem_type);
    if (problem_type == "basic")
      {
        thicknessIBCPtr = new BasicThicknessIBC;
      }
    else if (problem_type == "VieliPayne")
      {
        VieliPayneIBC* ibcPtr = new VieliPayneIBC;
        ParmParse pvPP("vieliPayne");

        Real thickness, seaLevel, originElevation;
        RealVect basalSlope;
        pvPP.get("thickness", thickness);
        seaLevel = 0.0;
        pvPP.query("seaLevel", seaLevel);
        
        Vector<Real> vect(SpaceDim);
        pvPP.getarr("basal_slope", vect, 0, SpaceDim);
        basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));

        pvPP.get("originElevation", originElevation);

        ibcPtr->setParameters(thickness, basalSlope, 
                              originElevation, seaLevel);       

        thicknessIBCPtr = static_cast<IceThicknessIBC*>(ibcPtr);
      }
     else if (problem_type == "marineIceSheet")
      {
        MarineIBC* ibcPtr = new MarineIBC;
        ParmParse mPP("marineIceSheet");
	
	Real  seaLevel;
	seaLevel = 0.0;
	mPP.query("seaLevel", seaLevel);

	std::string thicknessType = "constant";
	mPP.query("thickness_type",thicknessType);
	RefCountedPtr<RealFunction<RealVect> > thicknessFunction;
	Vector<RefCountedPtr<RealFunction<RealVect> > > bedrockFunction(1);
	if (thicknessType == "constant")
	  {
	    Real thickness;
	    mPP.get("thickness", thickness);
	    RefCountedPtr<RealFunction<RealVect> > ptr(new ConstantRealFunction<RealVect>(thickness));
	    thicknessFunction =ptr;
	  }
        else if (thicknessType == "compactSupportConstant")
          {
	    Real thickness;
            Vector<Real> tmpIntVect(SpaceDim,0); 
	    mPP.get("thickness", thickness);
            mPP.getarr("loBound", tmpIntVect, 0, SpaceDim);
            RealVect loBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));
            mPP.getarr("hiBound", tmpIntVect, 0, SpaceDim);
            RealVect hiBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));

	    RefCountedPtr<RealFunction<RealVect> > ptr(new CompactSupportConstantRealFunction(thickness,loBound, hiBound));
	    thicknessFunction =ptr;

          }
        else if (thicknessType == "compactSupportInclinedPlane")
          {
	    Real originThickness;
            Vector<Real> tmpIntVect(SpaceDim,0); 
	    mPP.get("origin_thickness", originThickness);
            mPP.getarr("loBound", tmpIntVect, 0, SpaceDim);
            RealVect loBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));
            mPP.getarr("hiBound", tmpIntVect, 0, SpaceDim);
            RealVect hiBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));
            RealVect thicknessSlope;
            Vector<Real> vect(SpaceDim);
            mPP.getarr("thickness_slope", vect, 0, SpaceDim);
            thicknessSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
            RefCountedPtr<RealFunction<RealVect> > ptr(new CompactSupportInclinedPlaneFunction(originThickness, thicknessSlope,loBound, hiBound));
	    thicknessFunction =ptr;
            
          }
	else if (thicknessType == "step")
	  {
            int dir=0;
	    Real leftThickness, rightThickness, cutoff;
	    mPP.get("left_thickness", leftThickness);
	    mPP.get("right_thickness", rightThickness);
	    mPP.get("x_cutoff", cutoff);
            mPP.query("dir", dir);

	    RefCountedPtr<RealFunction<RealVect> > ptr(new StepRealFunction(leftThickness, rightThickness, cutoff, dir));
	    thicknessFunction =ptr;
	  }
	else if (thicknessType == "flowline")
	  {
	    Real dx; 
	    std::string file, set;
	    mPP.get("thickness_flowline_dx", dx);
	    mPP.get("thickness_flowline_file", file);
	    mPP.get("thickness_flowline_set", set);
	    RefCountedPtr<RealFunction<RealVect> > ptr(new ExtrudedPieceWiseLinearFlowline(file,set,dx));
	    thicknessFunction = ptr;
	  }
	else 
	  {
	    MayDay::Error("bad marineIceSheet.thicknessType");
	  }

	
	std::string geometry = "plane";
	mPP.query("geometry",geometry);

	if (geometry == "plane")
	  {
	    //inclined plane geometry
	    RealVect basalSlope;
	    Vector<Real> vect(SpaceDim);
	    mPP.getarr("basal_slope", vect, 0, SpaceDim);
	    basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
	    Real originElevation;
	    mPP.get("originElevation", originElevation);

	   
	    RefCountedPtr<RealFunction<RealVect> > ptr(new InclinedPlaneFunction(originElevation, basalSlope));
	    bedrockFunction[0] =  ptr;
	  }
	else if (geometry == "symmetricPlane")
	  {
	    //inclined plane geometry, symmetric about origin
	    RealVect basalSlope;
	    Vector<Real> vect(SpaceDim);
	    mPP.getarr("basal_slope", vect, 0, SpaceDim);
	    basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));

	    Real originElevation;
	    mPP.get("originElevation", originElevation);

            RealVect symmetryPoint(RealVect::Zero);
            mPP.getarr("symmetryPoint", vect, 0, SpaceDim);
            symmetryPoint = RealVect(D_DECL(vect[0], vect[1], vect[2]));
	   
	    RefCountedPtr<RealFunction<RealVect> > ptr(new SymmetricInclinedPlaneFunction(originElevation, basalSlope, symmetryPoint));
	    bedrockFunction[0] =  ptr;
	  }

        else if (geometry == "regroundingTest")
          {
	    //inclined plane geometry with a Gaussian bump
            bedrockFunction.resize(2);

	    RealVect basalSlope;
	    Vector<Real> vect(SpaceDim);
	    mPP.getarr("basal_slope", vect, 0, SpaceDim);
	    basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
	    Real originElevation;
	    mPP.get("originElevation", originElevation);

            // compose flat plane with Gaussian hump
	    RefCountedPtr<RealFunction<RealVect> > ptr1(new InclinedPlaneFunction(originElevation, basalSlope));
	    bedrockFunction[0] =  ptr1;
            
            Real bumpCenter;
            Real bumpRad;
            Real bumpMag;
            
            mPP.get("bumpCenter", bumpCenter);
            mPP.get("bumpRad", bumpRad);
            mPP.get("bumpMag",bumpMag);
            RefCountedPtr<RealFunction<RealVect> > ptr2(new GaussianFunctionX(bumpCenter, bumpRad, bumpMag));

	    bedrockFunction[1] =  ptr2;
          }
	else if (geometry == "Schoof")
	  {
	    //geometry of Schoof, 2007
	    Real originElevation;
	    mPP.get("originElevation", originElevation);
	    
	    Real lengthScaleFactor = 1.0;
	    mPP.query("schoofLengthScaleFactor",lengthScaleFactor);

	   
	    Real schoofCoeff2, schoofCoeff4, schoofCoeff6;
	    mPP.get("schoofCoeff2", schoofCoeff2);
	    mPP.get("schoofCoeff4", schoofCoeff4);
	    mPP.get("schoofCoeff6", schoofCoeff6);
	    
	    //RefCountedPtr<RealFunction<RealVect> > schoofBedrock
	    RefCountedPtr<RealFunction<RealVect> > ptr(new SchoofBedrockElevation(domainSize[SpaceDim-2] * lengthScaleFactor,
										  originElevation,
										  schoofCoeff2, schoofCoeff4, 
										  schoofCoeff6));

	    bedrockFunction[0] = ptr;
	   

	  }

	else if (geometry == "Katz")
	  {
	    //geometry of Katz and Worster, 2010
	    Real originElevation;
	    mPP.get("originElevation", originElevation);
	    
	    Real lengthScaleFactor = 1.0;
	    mPP.query("schoofLengthScaleFactor",lengthScaleFactor);

	    Real katzAlpha, katzSigma;
	    mPP.get("katzAlpha", katzAlpha);
	    mPP.get("katzSigma", katzSigma);

	    Real schoofCoeff2, schoofCoeff4, schoofCoeff6;
	    mPP.get("schoofCoeff2", schoofCoeff2);
	    mPP.get("schoofCoeff4", schoofCoeff4);
	    mPP.get("schoofCoeff6", schoofCoeff6);
	    
	    //RefCountedPtr<RealFunction<RealVect> > katzBedrock
	    RefCountedPtr<RealFunction<RealVect> > ptr(new KatzBedrockElevation(domainSize[SpaceDim-2],
										domainSize[SpaceDim-1],
										originElevation,
										katzAlpha, katzSigma,
										lengthScaleFactor,
										schoofCoeff2, schoofCoeff4, 
										schoofCoeff6));

	    bedrockFunction[0] = ptr;
	  

	  }
	else
	  {
	    MayDay::Error("bad marineIceSheet.geometry");
	  }
	
	ibcPtr->setParameters(thicknessFunction, bedrockFunction ,  seaLevel);
        thicknessIBCPtr = static_cast<IceThicknessIBC*>(ibcPtr);
      }

     else if (problem_type =="fortran")
       {
	 FortranInterfaceIBC* ibcPtr = new FortranInterfaceIBC;
	 // need to set thickness and topography

	 ParmParse interfacePP("glimmerInterface");
	 Vector<int> nGhost(SpaceDim, 0);
	 IntVect ghostVect;
	 interfacePP.queryarr("numGhost", nGhost, 0, SpaceDim);
	 {
	   ghostVect = IntVect(D_DECL(nGhost[0], nGhost[1], nGhost[2]));
	 }

        // this is about removing ice from regions which
        // don't affect the dynamics of the region, but which 
        // can cause our solvers problems. Siple Island comes to mind here
        // for now, store these regions in a separate file. There's probably 
        // a better way to do this.
        
        bool clearThicknessRegions = false;
        if (interfacePP.contains("clearThicknessRegionsFile"))
          {
            clearThicknessRegions = true;
            std::string clearFile;
            interfacePP.get("clearThicknessRegionsFile", clearFile);
            pp2.query("zeroH_level", clearBoxesLevel);

            if (procID() == uniqueProc(SerialTask::compute))
              {
                ifstream is(clearFile.c_str(), ios::in);
                if (is.fail())
                  {
                    MayDay::Error("Cannot open file with regions for thickness clearing");
                  }
                // format of file: number of boxes, then list of boxes.
                int numRegions;
                is >> numRegions;

                // advance pointer in file
                while (is.get() != '\n');

                clearBoxes.resize(numRegions);

                for (int i=0; i<numRegions; i++)
                  {
                    Box bx;
                    is >> bx;
                    while (is.get() != '\n');

                    clearBoxes[i] = bx;
                  }
                    
              } // end if serial proc
            // broadcast results
            broadcast(clearBoxes, uniqueProc(SerialTask::compute));
            
            ibcPtr->setThicknessClearRegions(clearBoxes);
          }
        

	 
	 thicknessIBCPtr = static_cast<IceThicknessIBC*>(ibcPtr);
	 
       }
     else if (problem_type =="LevelData")
       {
	 LevelDataIBC* ibcPtr = new LevelDataIBC;	 
	 thicknessIBCPtr = static_cast<IceThicknessIBC*>(ibcPtr);
       }
     else 
      {
        MayDay::Error("bad problem type");
      }

    IceVelocitySolver* solverPtr = NULL;
    
    // default is picard solver
    int solverType = JFNK;
    ppAmr.query("velocity_solver_type", solverType);

    if (solverType == Picard)
      {
	MayDay::Error("Picard solver was deprecated: use JFNK");
        //solverPtr = new PicardSolver;
      }
    else if (solverType == JFNK)
      {
        solverPtr = new JFNKSolver;
      }
#ifdef CH_USE_PETSC
  else if (solverType == PetscNLSolver)
    {
      solverPtr = new PetscIceSolver;
    }
#endif
#ifdef CH_USE_FAS
  else if (solverType == FASMGAMR)
    {
      FASIceSolver *solver = new FASIceSolver;
      solverPtr = solver;

      solver->setParameters("FASSolver");

    }
#endif
  else
    {
      MayDay::Error("Unknown Solver Type" );
    }

    RealVect dxCrse = amrDx[0]*RealVect::Unit;
    
    thicknessIBCPtr->setGridHierarchy(vectCoordSys, amrDomains);
    
    solverPtr->define(amrDomains[0],
                      constRelPtr,
                      basalFrictionRelationPtr,
                      amrGrids,
                      refRatio,
                      dxCrse,
                      thicknessIBCPtr,
                      numLevels);
    
    solverPtr->setVerbosity(verbosity);
    
    solverPtr->setTolerance(velocity_solver_tolerance);
    
    if (maxSolverIterations > 0)
      {
        solverPtr->setMaxIterations(maxSolverIterations);
	
      }
    
    
    if (setVelToZero)
      {
        for (int lev=0; lev<velocity.size(); lev++)
          {
            LevelData<FArrayBox>& levelVel = *velocity[lev];
            DataIterator dit = levelVel.dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                levelVel[dit].setVal(0.0);
              }
          }
      }
    else if (interpFinestLevel)
      {
        int finestLevel = numLevels -1;
        if (finestLevel > 0)
          {
            FineInterp interpolator(amrGrids[finestLevel], 
                                    SpaceDim, refRatio[finestLevel-1], 
                                    amrDomains[finestLevel]);
            
            interpolator.interpToFine(*velocity[finestLevel],
                                      *velocity[finestLevel-1]);
          }
      }
    // now do solve...
    // set convergence metric to zero here 
    Real convergenceMetric = 0;
    Real initialResidNorm, finalResidNorm;
    
    solverPtr->solve(velocity,
		     calvedIceThickness, addedIceThickness, removedIceThickness,
                     initialResidNorm, finalResidNorm,
                     convergenceMetric,
                     vectRhs, 
                     vectBeta, vectC0,
                     vectA,
		     vectMuCoef,
                     vectCoordSys,
                     time,
                     0, max_level);
    


    if (doPlots)
      {
        // simplest option is to take existing plot data and 
        // overwrite the velocity field into it
        for (int lev=0; lev<numLevels; lev++)
          {
            LevelData<FArrayBox>& levelPlotData = *vectData[lev];
            LevelData<FArrayBox>& levelVel = *velocity[lev];

            for (int comp=0; comp<levelPlotData.nComp(); comp++)
              {
                if (vectNames[comp] == "xVel")
                  {
                    Interval destComps(comp, comp+SpaceDim-1);
                    levelVel.copyTo(levelVel.interval(),
                                    levelPlotData, destComps);
                  }
              } // end loop over data components            

#if 0
            // do fab-by-fab copy in order to get ghost cells
            DataIterator dit = levelVel.dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                FArrayBox& plotDataFab = levelPlotData[dit];
                FArrayBox& velFab = levelVel[dit];

                for (int comp=0; comp<levelPlotData.nComp(); comp++)
                  {
                    if (vectNames[comp] == "xVel")
                      {
                        plotDataFab.copy(velFab, 0, comp, SpaceDim); 
                      }
                  } // end loop over data components
              } // end loop over grids
#endif
          } // end loop over levels

        // create filename
        string filename = "benchmark.";
        filename.append(plotName);
        if (verbosity > 4) 
          {
            pout() << "writing plotfile: " << filename << endl;
          }
        
        WriteAMRHierarchyHDF5(filename, amrGrids, vectData, vectNames, 
                              amrDomains[0].domainBox(), amrDx[0], dt, time,
                              refRatio, numLevels);
                              
        
      } // end if we're writing plots


    // clean up

    if (solverPtr != NULL)
      {
        delete solverPtr;
        solverPtr = NULL;
      }

    if (constRelPtr != NULL)
      {
        delete constRelPtr;
        constRelPtr = NULL;
      }

    if (thicknessIBCPtr != NULL)
      {
        delete thicknessIBCPtr;
        thicknessIBCPtr=NULL;
      }

    for (int lev=0; lev<numLevels; lev++)
      {
        if (velocity[lev] != NULL)
          {
            delete velocity[lev];
            velocity[lev] = NULL;
          }

        if (vectRhs[lev] != NULL)
          {
            delete vectRhs[lev];
            vectRhs[lev] = NULL;
          }

        if (vectBeta[lev] != NULL)
          {
            delete vectBeta[lev];
            vectBeta[lev] = NULL;
          }

        if (vectC0[lev] != NULL)
          {
            delete vectC0[lev];
            vectC0[lev] = NULL;
          }

        if (vectA[lev] != NULL)
          {
            delete vectA[lev];
            vectA[lev] = NULL;
          }	

        if (vectSurface[lev] != NULL)
          {
            delete vectSurface[lev];
            vectSurface[lev] = NULL;
          }

        if (vectData[lev] != NULL)
          {
            delete vectData[lev];
            vectData[lev] = NULL;
          }

        if (calvedIceThickness[lev] != NULL)
          {
            delete calvedIceThickness[lev];
            calvedIceThickness[lev] = NULL;
          }
        if (removedIceThickness[lev] != NULL)
          {
            delete removedIceThickness[lev];
            removedIceThickness[lev] = NULL;
          }
        if (addedIceThickness[lev] != NULL)
          {
            delete addedIceThickness[lev];
            addedIceThickness[lev] = NULL;
          }
      }

  }  // end nested scope
  

  CH_TIMER_REPORT();

#ifdef CH_MPI
#ifdef CH_USE_PETSC
  ierr = PetscFinalize(); CHKERRQ(ierr);
#else
  MPI_Finalize();
#endif // petsc conditional

#endif // mpi conditional
  
  return ierr;
}
