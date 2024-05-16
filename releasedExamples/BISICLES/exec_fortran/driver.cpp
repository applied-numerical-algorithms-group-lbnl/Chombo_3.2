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
// driver.cpp
//
//===========================================================================
#include <iostream>
#include <fstream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "AmrIce.H"
#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFriction.H"
#include "BasalFrictionRelation.H"
#include "MuCoefficient.H"
#include "twistyStreamFriction.H"
#include "singularStreamFriction.H"
#include "GaussianBumpFriction.H"
#include "IceThicknessIBC.H"
#include "BasicThicknessIBC.H"
#include "VieliPayneIBC.H"
#include "MarineIBC.H"
#include "HumpIBC.H"
#include "LevelDataIBC.H"
#include "MultiLevelDataIBC.H"
#include "IceInternalEnergyIBC.H"
#include "LevelDataTemperatureIBC.H"
#include "VerticalConductionInternalEnergyIBC.H"
#include "LevelDataBasalFriction.H"
#include "PiecewiseLinearFlux.H"
#include "SurfaceFlux.H"
#include "IceConstants.H"
#include "AMRDamage.H"
#include "AMRMelange.H"
#include "DamageConstitutiveRelation.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
//#include "LevelDataSurfaceFlux.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "FineInterp.H"
#include "ReadLevelData.H"
#include "memusage.H"
 
#ifdef CH_USE_PETSC
#include "petsc.h"
#endif 
#include "Regression.H"

/// types of basal friction (beta) distributions
/** SinusoidalBeta is the one for exp C in Pattyn et al (2008)
    guassianBump is used for the MISMIP3D perturbations tests.
 */
enum basalFrictionTypes {constantBeta = 0,
                         sinusoidalBeta,
                         sinusoidalBetay,
                         twistyStreamx,
			 gaussianBump,
			 singularStream,
                         NUM_BETA_TYPES};

/// main program for 2D ice sheet models
/**
   \callgraph
 */
int main(int argc, char* argv[]) {
  
  int ierr = 0;

#ifdef CH_USE_PETSC
  ierr = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 
#endif // end petsc conditional

  FineInterp::s_default_boundary_limit_type = 0;

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
    //std::string in_file_mangled(in_file);
    //in_file_mangled += "_mangled";
    //FileMangler(in_file,in_file_mangled.c_str());
    //ParmParse pp(argc-2,argv+2,NULL,in_file_mangled.c_str());

    

    
    ParmParse pp(argc-2,argv+2,NULL,in_file);
    ParmParse pp2("main");

    std::string poutBaseName = "pout";
    pp2.query("poutBaseName",poutBaseName);
    setPoutBaseName(poutBaseName);
    //make use of number_procs and rank for at least something.
    pout() << "number_procs = " << number_procs << ", rank = " << rank << std::endl;

    
    RealVect domainSize;
    Vector<Real> domSize(SpaceDim);
    pp2.getarr("domain_size", domSize, 0, SpaceDim);
    domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));


    AmrIce amrObject;
    // ---------------------------------------------
    // set constitutive relation & rate factor
    // ---------------------------------------------

    Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
    {
      ParmParse ppc("constants");
      ppc.query("seconds_per_unit_time",seconds_per_unit_time);
    }
    
    std::string rateFactorType = "constRate";
    pp2.query("rateFactor", rateFactorType);
    if (rateFactorType == "constRate")
      {
	ParmParse crPP("constRate");
	Real A = 9.2e-18 * seconds_per_unit_time/SECONDS_PER_TROPICAL_YEAR;
	crPP.query("A", A);
	ConstantRateFactor rateFactor(A);
	amrObject.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "arrheniusRate")
      {
	ArrheniusRateFactor rateFactor(seconds_per_unit_time);
	ParmParse arPP("ArrheniusRate");
	amrObject.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "patersonRate")
      {
	PatersonRateFactor rateFactor(seconds_per_unit_time);
	ParmParse arPP("PatersonRate");
	amrObject.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "zwingerRate")
      {
	ZwingerRateFactor rateFactor(seconds_per_unit_time);
	ParmParse arPP("ZwingerRate");
	amrObject.setRateFactor(&rateFactor);
      }

    ConstitutiveRelation* constRelPtr = ConstitutiveRelation::parse("main");

    if (constRelPtr == NULL)
      {
	MayDay::Error("undefined constitutiveRelation in inputs");
      }

    amrObject.setConstitutiveRelation(constRelPtr);
 
    std::string basalRateFactorType = "";
    pp2.query("basalRateFactor", basalRateFactorType);
    
    if (basalRateFactorType == "patersonRate")
      {
	PatersonRateFactor rateFactor(seconds_per_unit_time);
	rateFactor.setA0(1.0);
	amrObject.setBasalRateFactor(&rateFactor);
      }

    // ---------------------------------------------
    // set surface flux. 
    // ---------------------------------------------

    SurfaceFlux* surf_flux_ptr = SurfaceFlux::parse("surfaceFlux");
    if (surf_flux_ptr == NULL)
      {
	const std::string err("failed to parse surfaceFlux (maybe you have the old style surface_flux_type?");
	pout() << err << endl;
	MayDay::Error(err.c_str());
      }

    amrObject.setSurfaceFlux(surf_flux_ptr);
    //delete surf_flux_ptr;
    // ---------------------------------------------
    // set basal (lower surface) flux. 
    // ---------------------------------------------
    
    SurfaceFlux* basal_flux_ptr = SurfaceFlux::parse("basalFlux");
    if (basal_flux_ptr == NULL)
      {
	const std::string err("failed to parse basalFlux (maybe you have the old style basal_flux_type?");
	pout() << err << endl;
	MayDay::Error(err.c_str());
      }

    amrObject.setBasalFlux(basal_flux_ptr); 
    //delete basal_flux_ptr;
     // ---------------------------------------------
    // set topography (bedrock) flux. 
    // ---------------------------------------------
    
    SurfaceFlux* topg_flux_ptr = SurfaceFlux::parse("topographyFlux");
    if (topg_flux_ptr == NULL)
      {
	topg_flux_ptr = new zeroFlux();
      }
    amrObject.setTopographyFlux(topg_flux_ptr); 

    // ---------------------------------------------
    // set mu coefficient
    // ---------------------------------------------
    {
      MuCoefficient* muCoefPtr =  MuCoefficient::parseMuCoefficient("muCoefficient");

      if (muCoefPtr == NULL)
	{
	  const std::string err("failed to parse muCoefficient");
	  pout() << err << endl;
	  MayDay::Error(err.c_str());
	}
     
      amrObject.setMuCoefficient(muCoefPtr);
      delete muCoefPtr;
    }

    // ---------------------------------------------
    // set basal friction coefficient and relation
    // ---------------------------------------------

    ParmParse geomPP("geometry");
    
    BasalFriction* basalFrictionPtr 
      = BasalFriction::parse("geometry", domainSize);
    
    if (basalFrictionPtr == NULL)
      {
	MayDay::Error("undefined  geometry.beta_type in inputs");
      }
    
    amrObject.setBasalFriction(basalFrictionPtr);
    
    BasalFrictionRelation* basalFrictionRelationPtr 
      = BasalFrictionRelation::parse("main",0);

    amrObject.setBasalFrictionRelation(basalFrictionRelationPtr);

    // ---------------------------------------------
    // set IBC -- this includes initial ice thickness, 
    // and basal geometry
    // ---------------------------------------------

    
    IceThicknessIBC* thicknessIBC = NULL;

    std::string problem_type;
    geomPP.get("problem_type", problem_type);
    if (problem_type == "basic")
      {
        thicknessIBC = new BasicThicknessIBC;
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

        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
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
        else if (thicknessType == "circularSupportConstant")
          {
	    Real thickness, supportRadius;
            Vector<Real> tmpRealVect(SpaceDim,0); 
	    mPP.get("thickness", thickness);
            mPP.get("supportRadius", supportRadius);
            mPP.getarr("center", tmpRealVect, 0, SpaceDim);
            RealVect center(D_DECL(tmpRealVect[0],tmpRealVect[1],tmpRealVect[2]));

	    RefCountedPtr<RealFunction<RealVect> > ptr(new CircularSupportConstantRealFunction(thickness,center, supportRadius));
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
	else if (geometry == "gaussianHump")
	  {
	    //gaussian bedrock geometry, symmetric about origin
	    RealVect center;
	    Vector<Real> vect(SpaceDim);
	    mPP.getarr("center", vect, 0, SpaceDim);
	    center = RealVect(D_DECL(vect[0], vect[1], vect[2]));

	    RealVect radius;
	    mPP.getarr("radius", vect, 0, SpaceDim);
	    radius = RealVect(D_DECL(vect[0], vect[1], vect[2]));            

	    Real magnitude;
	    mPP.get("magnitude", magnitude);
            
            Real offset;
	    mPP.get("offset", offset);

	    RefCountedPtr<RealFunction<RealVect> > ptr(new GaussianFunction(center,
                                                                            radius,
                                                                            magnitude,
                                                                            offset));
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
        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
      }
     else if (problem_type == "hump")
       {
         HumpIBC* ibcPtr = new HumpIBC;
         ParmParse humpPP("hump");

         Real maxThickness, radSqr, baseElevation, minThickness, seaLevel;
         RealVect center, widthScale;
         
         // default values to be equivalent to hump in Glimmer-CISM
         radSqr = 0.125*domainSize[0]*domainSize[1];
         maxThickness = 2000.0*pow(radSqr,0.5);
         baseElevation = 0.0;
         minThickness = 0.0;
         widthScale = RealVect::Unit;
         // this just lowers the sea level so that it's not relevant...
         seaLevel = -10.0;
         center = 0.5*domainSize;

         humpPP.query("radSqr", radSqr);
         humpPP.query("maxThickness", maxThickness);
         humpPP.query("baseElevation", baseElevation);
         humpPP.query("minThickness", minThickness);
         if (humpPP.contains("center"))
           {
             Vector<Real> centerArr(SpaceDim);
             humpPP.getarr("center", centerArr, 0, SpaceDim);
             center = RealVect(D_DECL(centerArr[0], centerArr[1], 
                                      centerArr[2]));
           }

         if (humpPP.contains("widthScale"))
           {
             Vector<Real> factorArr(SpaceDim);
             humpPP.getarr("widthScale", factorArr, 0, SpaceDim);
             widthScale = RealVect(D_DECL(factorArr[0], factorArr[1], 
                                           factorArr[2]));
           }

         ibcPtr->setParameters(maxThickness,
                               radSqr,
                               baseElevation,
                               minThickness,
                               center,
                               seaLevel, 
                               widthScale);
         
         thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
       }
     else if (problem_type == "LevelData")
       {
	 //read geometry from an AMR Hierarchy, store in LevelDataIBC
	 ParmParse ildPP("inputLevelData");
	 std::string infile;
	 ildPP.get("geometryFile",infile);
	 std::string thicknessName = "thck";
	 ildPP.query("thicknessName",thicknessName);
	 std::string topographyName = "topg";
	 ildPP.query("topographyName",topographyName);

         // default values (only relevant when LevelData input doesn't
         // cover the domain
         Real defaultThickness = 0.0;
         ildPP.query("defaultThickness", defaultThickness);
         Real defaultTopography = -10000.0;
         ildPP.query("defaultTopography", defaultTopography);
         bool setDefaultValues = false;
         ildPP.query("setDefaultValues", setDefaultValues);
         
	 RefCountedPtr<LevelData<FArrayBox> > levelThck
	   (new LevelData<FArrayBox>());
	 RefCountedPtr<LevelData<FArrayBox> > levelTopg
	   (new LevelData<FArrayBox>());

	 Real dx;

	 Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
	 vectData.push_back(levelThck);
	 vectData.push_back(levelTopg);

	 Vector<std::string> names(2);
	 names[0] = thicknessName;
	 names[1] = topographyName;
	 readLevelData(vectData,dx,infile,names,1);

         // this is about removing ice from regions which
         // don't affect the dynamics of the region, but which 
         // can cause our solvers problems. 
         // for now, store these regions in a separate file. 
         // There's probably a better way to do this.
         
         // this will contain the boxes in the index space of the 
         // original LevelData in which the thickness will be cleared.
         
         if (ildPP.contains("clearThicknessRegionsFile"))
           {
             Vector<Box> clearBoxes;
             std::string clearFile;
             ildPP.get("clearThicknessRegionsFile", clearFile);
             
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
             
             // now loop over the thickness levelData and set intersections
             // with boxes to zero
             
             DataIterator dit = levelThck->dataIterator();
             for (dit.begin(); dit.ok(); ++dit)
               {
                 FArrayBox& thickFab = levelThck->operator[](dit);
                 const Box& fabBox = thickFab.box();
                 for (int boxno=0; boxno<clearBoxes.size(); boxno++)
                   {
                     Box intersectBox(fabBox);
                     intersectBox &= clearBoxes[boxno];
                     if (!intersectBox.isEmpty())
                       {
                         thickFab.setVal(0.0,intersectBox,0);
                       } // end if there's an intersection
                   } // end loop over clearboxes
               } // end loop over grids in thickness levelData

           } // end if we're setting thickness to zero
       
	 RealVect levelDx = RealVect::Unit * dx;
	 LevelDataIBC* ptr = new LevelDataIBC(levelThck,levelTopg,levelDx,
                                              defaultThickness,
                                              defaultTopography,
                                              setDefaultValues
                                              );
	 thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
       }
     else if (problem_type == "MultiLevelData")
       {
	 //read geometry from an AMR Hierarchy, store in MultiLevelDataIBC
	 ParmParse ildPP("inputLevelData");
	 std::string infile;
	 ildPP.get("geometryFile",infile);
	 std::string thicknessName = "thck";
	 ildPP.query("thicknessName",thicknessName);
	 std::string topographyName = "topg";
	 ildPP.query("topographyName",topographyName);
	
	
	 Real dx;
	 Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectData;
	

	 Vector<std::string> names(2);
	 names[0] = thicknessName;
	 names[1] = topographyName;
	 Vector<int> refRatio;
	 readMultiLevelData(vectData,dx,refRatio,infile,names,1);
	
	 RealVect crseDx = RealVect::Unit * dx;
	 MultiLevelDataIBC* ptr = new MultiLevelDataIBC
	   (vectData[0],vectData[1],crseDx,refRatio);
	 thicknessIBC = static_cast<IceThicknessIBC*>( ptr);

       }
#ifdef HAVE_PYTHON
     else if (problem_type == "Python")
       {
	 
	 ParmParse pyPP("PythonIBC");
	 std::string module;
	 pyPP.get("module",module);
	 std::string thckFuncName = "thickness";
	 pyPP.query("thicknessFunction",thckFuncName);
	 std::string topgFuncName = "topography";
	 pyPP.query("topographyFunction",topgFuncName);
	 std::string rhsFuncName = "";
	 pyPP.query("RHSFunction",rhsFuncName);
	 std::string faceVelFuncName = "";
	 pyPP.query("faceVelFunction",faceVelFuncName);
	 PythonInterface::PythonIBC* ptr = new PythonInterface::PythonIBC
	   (module, thckFuncName, topgFuncName, rhsFuncName,faceVelFuncName);
	 thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
       }
#endif
     else 
       {
         MayDay::Error("bad problem type");
       }

    amrObject.setThicknessBC(thicknessIBC);
    
    {
      // ---------------------------------------------
      // set surface heat boundary data 
      // ---------------------------------------------

      SurfaceFlux* surf_heat_boundary_data_ptr = SurfaceFlux::parse("surfaceHeatBoundaryData");
      ParmParse pps("surfaceHeatBoundaryData");
      bool diri = false; // flux boundary data by default
      pps.query("Dirichlett",diri);
      bool temp = true; //temperature boundary data by default
      pps.query("Temperature",temp);
      if (surf_heat_boundary_data_ptr == NULL)
	{
	  if (diri)
	    {
	      const std::string err("If surfaceHeatBoundaryData.Dirichlett = true, surfaceHeatBoundaryData.type must be set");
	      pout() << err << endl;
	      MayDay::Error(err.c_str());
	    }
	  else
	    {
	      const std::string warn("No surfaceHeatBoundaryData.type specified, so zero flux set. Only relevant for amr.isothermal = false");
	      pout() << warn << endl;
	      // only warn if we're on processor 0, otherwise we get 
	      // nproc copies of this warning to stderr
	      if (procID() == uniqueProc(SerialTask::compute))
		{
		  MayDay::Warning(warn.c_str());
		}
	      surf_heat_boundary_data_ptr = new zeroFlux();
	    }
	}

      amrObject.setSurfaceHeatBoundaryData(surf_heat_boundary_data_ptr, diri, temp);
      if (surf_heat_boundary_data_ptr != NULL)
      	{
      	  delete surf_heat_boundary_data_ptr;
      	  surf_heat_boundary_data_ptr=NULL;
      	}
    
      // ---------------------------------------------
      // set basal (lower surface) heat boundary data. 
      // ---------------------------------------------
      
      SurfaceFlux* basal_heat_boundary_data_ptr = SurfaceFlux::parse("basalHeatBoundaryData");
      if (basal_heat_boundary_data_ptr == NULL)
       	{
       	  basal_heat_boundary_data_ptr = new zeroFlux();
       	}
      
      amrObject.setBasalHeatBoundaryData(basal_heat_boundary_data_ptr);
      if (basal_heat_boundary_data_ptr != NULL)
	{
	  delete basal_heat_boundary_data_ptr;
	  basal_heat_boundary_data_ptr=NULL;
	}
      
    }
    
    {
      IceInternalEnergyIBC* internalEnergyIBC = NULL;
      ParmParse tempPP("temperature");
      std::string tempType("constant");
      tempPP.query("type",tempType);
      if (tempType == "constant")
	{
	  Real T = 258.0;
	  tempPP.query("value",T);
	  ConstantIceTemperatureIBC* ptr = new ConstantIceTemperatureIBC(T);
	  internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
	}
      else if (tempType == "LevelData")
	{
	  ParmParse ildPP("inputLevelData");
	  LevelDataTemperatureIBC* ptr = NULL;
	  ptr = LevelDataTemperatureIBC::parse(ildPP); CH_assert(ptr != NULL);
	  internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
	}
      else if (tempType == "VerticalConduction")
	{
	  VerticalConductionInternalEnergyIBC* ptr = NULL;
	  ptr = VerticalConductionInternalEnergyIBC::parse(tempPP); 
	  CH_assert(ptr != NULL);
	  internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
	}
#ifdef HAVE_PYTHON
      else if (tempType == "Python")
	{
	  ParmParse pyPP("PythonIceTemperatureIBC");
	  std::string module;
	  pyPP.get("module",module);
	  std::string funcName = "temperature";
	  pyPP.query("function",funcName);
	  internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>
	    (new PythonInterface::PythonIceTemperatureIBC(module, funcName));
	}
#endif
      else 
	{
	  MayDay::Error("bad temperature/internal energy type");
	}	
      amrObject.setInternalEnergyBC(internalEnergyIBC);
      if (internalEnergyIBC != NULL)
	{
	  delete internalEnergyIBC;
	}
    }

    amrObject.setDomainSize(domainSize);


    CalvingModel* calving_model_ptr = CalvingModel::parseCalvingModel("CalvingModel");
    if (calving_model_ptr == NULL)
      {
	calving_model_ptr = new NoCalvingModel;
      }
    amrObject.setCalvingModel(calving_model_ptr);
    
    { 
      //// TODO this is just a temporary means of initializing a 
      //// damage model observer etc. Once we have it working, it 
      //// probably needs to be bound up with MuCoefficient
      bool damage_model = false;
      pp2.query("damage_model",damage_model);
      if (damage_model)
	{
	  ///currently not deleting this, need to decide
	  ///how that will be done
	  DamageIceObserver* ptr = new DamageIceObserver();
	  amrObject.addObserver(ptr);

	  //whatever the constitutive relation was, wrap
	  //it up in a DamageConstitutiveRelation tied
	  //to the DamageIceObserver components
	  DamageConstitutiveRelation* dcrptr = 
	    new DamageConstitutiveRelation(constRelPtr, &ptr->damage());
	  amrObject.setConstitutiveRelation(dcrptr);

	  CalvingModel* d_calving_model_ptr = new DamageCalvingModel(calving_model_ptr, &ptr->damage());
	  amrObject.setCalvingModel(d_calving_model_ptr);
	  delete d_calving_model_ptr;
	  
	}
    }

    
    {
      /// initialize the melange model
      bool melange_model = false;
      pp2.query("melange_model",melange_model);
      if (melange_model)
	{
	  MelangeIceObserver* ptr = new MelangeIceObserver();
	  amrObject.addObserver(ptr);
	}
    }
    
    // set up initial grids, initialize data, etc.
    amrObject.initialize();
 

    int maxStep;
    Real maxTime;
    //Real startTime;
    pp2.get("maxTime", maxTime);
    pp2.get("maxStep", maxStep);
    
    amrObject.run(maxTime, maxStep);
    
    // clean up
    if (constRelPtr != NULL)
      {
        delete constRelPtr;
        constRelPtr = NULL;
      }

    if (surf_flux_ptr != NULL)
      {
        delete surf_flux_ptr;
        surf_flux_ptr = NULL;
      }
    
    if (basal_flux_ptr != NULL)
      {
        delete basal_flux_ptr;
        basal_flux_ptr = NULL;
      }

    if (topg_flux_ptr != NULL)
      {
        delete topg_flux_ptr;
        topg_flux_ptr = NULL;
      }

    if (calving_model_ptr != NULL)
      {
	delete calving_model_ptr;
	calving_model_ptr = NULL;
      }
    
    if (basalFrictionPtr != NULL)
      {
	delete basalFrictionPtr;
	basalFrictionPtr = NULL;
      }

    if (basalFrictionRelationPtr != NULL)
      {
	delete basalFrictionRelationPtr;
	basalFrictionRelationPtr = NULL;
      }

    if (thicknessIBC != NULL)
      {
        delete thicknessIBC;
        thicknessIBC=NULL;
      }

     
    
#ifdef CH_USE_HDF5
    {
      // finally, carry out an optional regression test
      
      ParmParse ppr("regression");
      
      std::string result_hdf5("");
      ppr.query("result_hdf5", result_hdf5);

      if (result_hdf5 != "")
	{
	  std::string reference_hdf5;
	  ppr.get("reference_hdf5", reference_hdf5);

	  Real tol = 1.0e-10;
	  ppr.query("tol",tol);
          Real norm = HDF5NormTest(result_hdf5, reference_hdf5);
      	  if (tol < norm)
	    {
	      ierr = 1;
              pout() << "FAILED HDF5NormTest!  norm = "
                     << norm 
                     << "   tolerance = " << tol 
                     << endl;
	    }
	}

    }
#endif

  }  // end nested scope
  

  // report memory usage and leaks
#ifdef REPORT_MEMORY
  dumpmemoryatexit();
#endif
  CH_TIMER_REPORT();

#ifdef HAVE_PYTHON
  Py_Finalize();
#endif

#ifdef CH_USE_PETSC
  ierr = PetscFinalize(); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Finalize();
#endif // mpi conditional
#endif // petsc conditional

  return ierr;
}


