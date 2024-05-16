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
// amrtocf.cpp
// Read data from hdf5 files containing Chombo block-structured AMR hierachies
// Write unstructurd data, plus the grid data needed to reconstruct
// the block structured data, to a CF compliant netcdf file
// Or the other way round
//===========================================================================

#include <iostream>
#include "AMRIO.H"
#include "UnstructuredIO.H"
#include "ParmParse.H"
#include "FieldNames.H"

void AMRtoCF(const std::string& ifile, const std::string& ofile,
	     const std::string& a_createdBy, 
	     const RealVect& a_origin,
	     const RealVect& a_loXYGrid,
	     const Real& a_timeOffset,
	     bool a_validOnly);

void CFtoAMR(const std::string& ifile, const std::string& ofile);




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
  
#ifdef CH_MPI
    std::string msg = argv[0] + std::string(" is serial only");
    MayDay::Error(msg.c_str());
#endif
    

    if(argc < 2) 
      { std::cerr << " usage: " << argv[0] << " <config_file> [additional key=value args]\n"; exit(0); }

    char* config_file = argv[1];

    ParmParse pp(argc-2,argv+2,NULL,config_file);

    pout() << "config_file read in " << config_file << std::endl;

    std::string ifile;
    pp.get("infile",ifile);

    std::string ofile;
    pp.get("outfile",ofile);
 
   
    
    if ( (ifile.size() >= 5) && (ifile.find(".hdf5") == ifile.size()-5))
      {
	bool validOnly = true;
	pp.query("valid_only",validOnly);

	Real timeOrigin = 0.0;
	pp.query("time_origin", timeOrigin);

	RealVect origin = RealVect::Zero;
	{
	  Vector<Real> t(SpaceDim,0.0);
	  pp.queryarr("origin",t,0,SpaceDim);
	  D_TERM(origin[0] = t[0];, origin[1] = t[1];, origin[2] = t[2];);
	}

	// Provide the x,y coordinates of the lower left hand corner of the grid to calculate the Polar Stereographic projection.
	// The grid in the plotfiles start from (0,0).
	// e.g. Origin for MCdataset-2015-03-03.nc is (-638000.0,-3349500.0, 0.0)
	RealVect loXYGrid = RealVect::Zero;
	{
	  Vector<Real> t(SpaceDim,0.0);
	  pp.queryarr("lo_XY_grid",t,0,SpaceDim);
	  D_TERM(loXYGrid[0] = t[0];, loXYGrid[1] = t[1];, loXYGrid[2] = t[2];);
	}

	std::string created;
	for (int i =0; i < argc; i++)
	  {
	    created += std::string(argv[i]) + " ";
	  }
	
	AMRtoCF(ifile,ofile,created,origin,loXYGrid,timeOrigin,validOnly);
      }
    else if ((ifile.size() >= 3) && (ifile.find(".nc") == ifile.size()-3))
      {
        CFtoAMR(ifile,ofile);
      }
    else
      {
	pout() << "unknown input file type [" <<  ifile << "]" << std::endl;
	MayDay::Error("unknown input file type");
      }


   

   
    
  }  // end nested scope
  CH_TIMER_REPORT();
  
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
void AMRtoCF(const std::string& ifile, const std::string& ofile,
	     const std::string& a_created, 
	     const RealVect& a_origin,
	     const RealVect& a_loXYGrid,
	     const Real& a_timeOrigin,
	     bool a_validOnly)
{
  

  //read the AMR data
  Vector<std::string> names;
  Vector<LevelData<FArrayBox>* > data;
  Vector<DisjointBoxLayout> grids;
  Vector<int> ratio;
  Real crseDx = 0.0, dt = 0.0, time = 0.0;
  Box domBox;
  int numLevels;
  int status = ReadAMRHierarchyHDF5
    (ifile,grids,data,names,domBox,crseDx,dt,time,
     ratio,numLevels);
  if (status != 0)
    {
      MayDay::Error("failed to read AMR hierarchy");
    }

 

  //extract valid data
  UnstructuredData usData(names.size(), RealVect::Unit*crseDx, domBox, ratio, time + a_timeOrigin,  a_origin);
  UnstructuredIO::BStoUS(usData, data, ratio, a_validOnly);
 
#ifdef HAVE_GDAL 

  int epsgval = 3413;  // ESPG number for Polar Stereographic North (70 degN, 45 degW) projection.

  Real x0 = a_loXYGrid[0];
  Real y0 = a_loXYGrid[1];

  gdalxytolatlon transformation
  (epsgval, x0, y0);

#else
  PolarStereographicCartesianToPolarTransformation transformation
    (0.08181922,  6.3781370e+6, 1.0 , 0.0, RealVect::Zero);
  //eccentricity, equatorial radius,
#endif 

  //write to CF file
  UnstructuredIO::writeCF ( ofile,  usData , names, a_created , transformation );

  for (int lev = 0; lev < numLevels; lev++)
    {
      if (data[lev] != NULL)
	{
	  delete data[lev];data[lev]=NULL;
	}
    }
}
void CFtoAMR(const std::string& ifile, const std::string& ofile)
{
  UnstructuredData usData;
  Vector<std::string> names;
  UnstructuredIO::readCF(usData, names, ifile);

  Vector<LevelData<FArrayBox> *> data;
  UnstructuredIO::USToBS(data, usData);

  Vector<DisjointBoxLayout> grids(data.size());
  for (int lev = 0; lev < grids.size(); lev++)
    {
      grids[lev] = data[lev]->disjointBoxLayout();
    }

  Real dt =0.0; Real time = 0.0;

  WriteAMRHierarchyHDF5
    (ofile,grids,data,names, usData.domain(0),usData.dx()[0][0],
     dt,time, usData.ratio() , usData.nLevel() );


  for (int lev = 0; lev < data.size(); lev++)
    {
      if (data[lev] != NULL)
	{
	  delete data[lev];data[lev]=NULL;
	}
    }
  
}


