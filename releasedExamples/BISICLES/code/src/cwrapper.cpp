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
// cwrapper.cpp
//===========================================================================
#include <iostream>
#include <fstream>
#include <cstring>
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
#include "LevelDataBasalFriction.H"
#include "PiecewiseLinearFlux.H"
#include "ComplexSurfaceFlux.H"
#include "IceConstants.H"
#include "AMRDamage.H"
#include "AMRMelange.H"
#include "DamageConstitutiveRelation.H"
#include "IceThermodynamics.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "FillFromReference.H"
#include "LevelDataSurfaceFlux.H"
#include "ReadLevelData.H"
#include "PetscSolver.H"
#include "cwrapper.H"
#include "cdriverconstants.h"

/// Associates an AmrIce with e.g LevelDataSurfaceFlux* objects intended for I/O with C/Fortran programs
/**
   The cdriver interface maintains a table of BISICLES instances in a map<int, BisiclesWrapper*>.
   That allows C/Fortran programs to refer to an in-memory BISICLES by an integer id.  BisiclesWrapper,
   therefore, just exists to support this.
 */
class BisiclesWrapper
{
public:
#ifdef CH_MPI
  MPI_Comm mpi_comm;
#endif
  AmrIce m_amrIce;
  AMRMelange* m_amrMelange;
  LevelDataSurfaceFlux* m_surface_flux;
  LevelDataSurfaceFlux* m_basal_flux;
  LevelDataSurfaceFlux* m_floating_ice_basal_flux;
  LevelDataSurfaceFlux* m_grounded_ice_basal_flux;
  LevelDataSurfaceFlux* m_surface_heat_boundary_data;
  bool m_surface_heat_boundary_dirichlett;
  LevelDataSurfaceFlux* m_basal_heat_boundary_data;

  //optional bedrock elevation change per unit time (GIA)
  LevelDataSurfaceFlux* m_topography_flux;

  //optional (initial) geometry data
  RefCountedPtr<LevelData<FArrayBox> >m_geometry_ice_thickness;
  RefCountedPtr<LevelData<FArrayBox> > m_geometry_bedrock_elevation;
  RealVect m_geometry_dx;
  ParmParse* m_pp;
  
  BisiclesWrapper();

  ~BisiclesWrapper();
};

BisiclesWrapper::BisiclesWrapper()
{
  
   m_surface_flux = NULL;
   m_basal_flux = NULL;
   m_floating_ice_basal_flux = NULL;
   m_grounded_ice_basal_flux = NULL;
   m_surface_heat_boundary_data = NULL;
   m_surface_heat_boundary_dirichlett = true;
   m_basal_heat_boundary_data = NULL;
   m_topography_flux = NULL;
   m_geometry_dx = -RealVect::Unit;
}

BisiclesWrapper::~BisiclesWrapper()
{
  if (m_surface_flux != NULL) 
    delete m_surface_flux;
  if (m_basal_flux != NULL) 
    delete m_basal_flux;
  if (m_floating_ice_basal_flux != NULL) 
    delete m_floating_ice_basal_flux;
  if (m_grounded_ice_basal_flux != NULL) 
    delete m_grounded_ice_basal_flux;
  if (m_surface_heat_boundary_data != NULL) 
    delete m_surface_heat_boundary_data;
  if (m_basal_heat_boundary_data != NULL) 
    delete m_basal_heat_boundary_data;
  if (m_topography_flux != NULL)
    delete m_topography_flux;
  if (m_pp != NULL)
    delete m_pp;
}


namespace bisicles_c_wrapper
{
  std::map<int, BisiclesWrapper*> instances;
}

//fortran wrappers
void f_bisicles_new_instance_(int *instance_key,  char *input_fname, const int *len_fname, const int *f_mpi_comm)
  {
    input_fname[*len_fname - 1] = 0; // null terminate the string
#ifdef CH_MPI
    MPI_Comm mpi_comm =  MPI_Comm_f2c(*f_mpi_comm);
    bisicles_new_instance(instance_key, input_fname, mpi_comm);
#else
    bisicles_new_instance(instance_key, input_fname, *f_mpi_comm);
#endif
   
  }

void f_bisicles_free_instance_(int *instance_key)
{
  bisicles_free_instance(instance_key);
}

void f_bisicles_write_checkpoint_(int *instance_key)
{
  bisicles_write_checkpoint(instance_key);
}

void f_bisicles_write_plot_(int *instance_key)
{
  bisicles_write_plot(instance_key);
}


void f_bisicles_read_checkpoint_(int *instance_key, char *checkpoint_fname, const int *len_fname)
{
  checkpoint_fname[*len_fname - 1] = 0; // null terminate the string
  bisicles_read_checkpoint(instance_key, checkpoint_fname);
}


void f_bisicles_set_2d_data_(int *instance_key,  double *data_ptr, const int *field, 
			   const double *dx, const int *dims, 
			   const int *boxlo, const int *boxhi)
{
  bisicles_set_2d_data(instance_key,  data_ptr, field, 
		       dx, dims, boxlo, boxhi);
}

void f_bisicles_set_2d_geometry_(int *instance_key,  double *thck_data_ptr, double *topg_data_ptr, 
			       const double *dx, const int *dims, 
			       const int *boxlo, const int *boxhi)
{
  bisicles_set_2d_geometry(instance_key,  thck_data_ptr, topg_data_ptr,dx, dims, boxlo, boxhi);
}

void f_bisicles_get_2d_data_(int *intance_id, double *data_ptr, const int *field,
			   const double *dx, const int *dims, 
			   const int *boxlo, const int *boxhi)
{
  bisicles_get_2d_data(intance_id, data_ptr, field,
		       dx, dims, boxlo, boxhi);
  
}

void f_bisicles_init_instance_(int *instance_key)
  {
    bisicles_init_instance(instance_key);
  }

void f_bisicles_advance_(int *instance_key, double *start_time, double *max_time, int *max_step)
  {
    std::cout << "f_bisicles_advance_" << *max_time << " " << *max_step << std::endl;
    bisicles_advance(instance_key, start_time, max_time, max_step);
  }


void f_bisicles_new_instance(int *instance_key,  char *input_fname, const int *len_fname, const int *f_mpi_comm)
  {
    input_fname[*len_fname - 1] = 0; // null terminate the string
#ifdef CH_MPI
    MPI_Comm mpi_comm =  MPI_Comm_f2c(*f_mpi_comm);
    bisicles_new_instance(instance_key, input_fname, mpi_comm);
#else
    bisicles_new_instance(instance_key, input_fname, *f_mpi_comm);
#endif
  }
void f_bisicles_free_instance(int *instance_key)
{
  bisicles_free_instance(instance_key);
}
void f_bisicles_set_2d_data(int *instance_key,  double *data_ptr, const int *field, 
			   const double *dx, const int *dims, 
			   const int *boxlo, const int *boxhi)
{
  bisicles_set_2d_data(instance_key,  data_ptr, field, 
		       dx, dims, boxlo, boxhi);
}

void f_bisicles_set_2d_geometry(int *instance_key,  double *thck_data_ptr, double *topg_data_ptr, 
			       const double *dx, const int *dims, 
			       const int *boxlo, const int *boxhi)
{
  bisicles_set_2d_geometry(instance_key,  thck_data_ptr, topg_data_ptr,dx, dims, boxlo, boxhi);
}

void f_bisicles_get_2d_data(int *intance_id, double *data_ptr, const int *field,
			   const double *dx, const int *dims, 
			   const int *boxlo, const int *boxhi)
{
  bisicles_get_2d_data(intance_id, data_ptr, field,
		       dx, dims, boxlo, boxhi);
  
}

void f_bisicles_init_instance(int *instance_key)
  {
    bisicles_init_instance(instance_key);
  }

void f_bisicles_advance(int *instance_key, double *start_time, double *max_time, int *max_step)
  {
    bisicles_advance(instance_key, start_time, max_time, max_step);
  }


void f_bisicles_write_checkpoint(int *instance_key)
{
  bisicles_write_checkpoint(instance_key);
}


void f_bisicles_write_plot(int *instance_key)
{
  bisicles_write_plot(instance_key);
}

void f_bisicles_read_checkpoint(int *instance_key, char *checkpoint_fname, const int *len_fname)
{
  checkpoint_fname[*len_fname - 1] = 0; // null terminate the string
  bisicles_read_checkpoint(instance_key, checkpoint_fname);
}

void bisicles_set_header_int(int *instance_key, const char* attr_key,  const int *val)
{
  bisicles_set_header(instance_key, attr_key, val);
}




void f_bisicles_set_header_int(int *instance_key, char* attr_key, const int *attr_key_len, const int *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_set_header_int(instance_key,attr_key,val);
}
void f_bisicles_set_header_dble(int *instance_key, char* attr_key, const int *attr_key_len, const double *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_set_header_dble(instance_key,attr_key,val);
}

void f_bisicles_set_header_char(int *instance_key, char* attr_key, const int *attr_key_len,  char *val, const int *len_val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  val[*len_val - 1] = 0;
  bisicles_set_header_char(instance_key,attr_key,val);
}


void f_bisicles_get_header_char(int *instance_key, char* attr_key, const int *attr_key_len,  int *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_get_header_int(instance_key,attr_key,val);
}
void f_bisicles_get_header_dble(int *instance_key, char* attr_key, const int *attr_key_len,  double *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_get_header_dble(instance_key,attr_key,val);
}
void f_bisicles_get_header_char(int *instance_key, char* attr_key, const int *attr_key_len,  char *val, const int *len_val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  val[*len_val - 1] = 0;
  bisicles_get_header_char(instance_key,attr_key,val);
}


void f_bisicles_set_header_int_(int *instance_key, char* attr_key, const int *attr_key_len, const int *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_set_header_int(instance_key,attr_key,val);
}

void f_bisicles_set_header_dble_(int *instance_key, char* attr_key, const int *attr_key_len, const double *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_set_header_dble(instance_key,attr_key,val);
}
void f_bisicles_set_header_char_(int *instance_key, char* attr_key, const int *attr_key_len,  char *val, const int *len_val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  val[*len_val - 1] = 0;
  bisicles_set_header_char(instance_key,attr_key,val);
}
void f_bisicles_get_header_char_(int *instance_key, char* attr_key, const int *attr_key_len,  int *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_get_header_int(instance_key,attr_key,val);
}

void f_bisicles_get_header_int_(int *instance_key, char* attr_key, const int *attr_key_len,  int *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_get_header_int(instance_key,attr_key,val);
}

void f_bisicles_get_header_dble_(int *instance_key, char* attr_key, const int *attr_key_len,  double *val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  bisicles_get_header_dble(instance_key,attr_key,val);
}
void f_bisicles_get_header_char_(int *instance_key, char* attr_key, const int *attr_key_len,  char *val, const int *len_val)
{
  attr_key[*attr_key_len - 1] = 0; // null terminate the string
  val[*len_val - 1] = 0;
  bisicles_get_header_char(instance_key,attr_key,val);
}



// initialize the AmrIce object in a wrapper. Any surface
// fluxes specified in the wrapper will be added to whatever
// is specified in the input file (a_innputfile) 
// computes the initial solve / load from checkpoint etc
// such that the object is ready to be used in timestepping
void init_bisicles_instance(BisiclesWrapper& a_wrapper)
{
#ifdef CH_MPI
  //Chombo_MPI::comm = a_wrapper.mpi_comm;
  MPI_Barrier(Chombo_MPI::comm);
#endif
//   int rank, number_procs;
// #ifdef CH_MPI
//   MPI_Comm_rank(Chombo_MPI::comm, &rank);
//   MPI_Comm_size(Chombo_MPI::comm, &number_procs);
// #else
//   rank=0;
//   number_procs=1;
// #endif
  
  
  AmrIce& amrObject = a_wrapper.m_amrIce;
  
  ParmParse pp2("main");

  ///\todo Check how pout will work with multiple instances
  std::string poutBaseName = "pout";
  pp2.query("poutBaseName",poutBaseName);
  setPoutBaseName(poutBaseName);

  RealVect domainSize;
  Vector<Real> domSize(SpaceDim);
  pp2.getarr("domain_size", domSize, 0, SpaceDim);
  domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));

 
  Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
  {
    ParmParse ppc("constants");
    ppc.query("seconds_per_unit_time",seconds_per_unit_time);
  }
   
  

  // ---------------------------------------------
  // set constitutive relation & rate factor
  // ---------------------------------------------
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
 
  // -------------------------------------------------
  // set surface flux.
  // --------------------------------------------------

  SurfaceFlux* surf_flux_ptr = SurfaceFlux::parse("surfaceFlux");
  if ( surf_flux_ptr == NULL )
    {
      const std::string err("failed to parse surfaceFlux (maybe you have the old style surface_flux_type?");
      pout() << err << endl;
      MayDay::Error(err.c_str());
    }

  if (a_wrapper.m_surface_flux)
    {
      AxbyFlux* ptr = new AxbyFlux(1.0, a_wrapper.m_surface_flux, 1.0, surf_flux_ptr);	
      amrObject.setSurfaceFlux(ptr);
      delete(ptr);
    }
  else
    {
      amrObject.setSurfaceFlux(surf_flux_ptr);
    }
  
 

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

  if ( a_wrapper.m_floating_ice_basal_flux || a_wrapper.m_grounded_ice_basal_flux)
    {
      //need to specify floating and grounded ice fluxes etc
      /// \todo fix this mess...
      SurfaceFlux*   f; 
      if (a_wrapper.m_floating_ice_basal_flux == NULL)
	{
	  f = new zeroFlux();
	}
      else
	{
	   f = a_wrapper.m_floating_ice_basal_flux;
	}
      /// \todo fix this mess...
      SurfaceFlux*   g;
      if (a_wrapper.m_grounded_ice_basal_flux == NULL)
	{
	  g = new zeroFlux();
	}
      else
	{
	  g = a_wrapper.m_grounded_ice_basal_flux;
	}
      MaskedFlux* m = new MaskedFlux(g,f,f,g);
      AxbyFlux* ptr = new AxbyFlux(1.0, m , 1.0, basal_flux_ptr);
      amrObject.setBasalFlux(ptr); 
      delete(ptr);
    }
  else if (a_wrapper.m_basal_flux)
    {
      AxbyFlux* ptr = new AxbyFlux(1.0, a_wrapper.m_basal_flux, 1.0, basal_flux_ptr);
      amrObject.setBasalFlux(ptr);
      delete(ptr);
    }
  else
    {
      amrObject.setBasalFlux(basal_flux_ptr);
    }


  // ---------------------------------------------
  // set  topography (bedrock) flux
  // ---------------------------------------------
  
  SurfaceFlux* topg_flux_ptr = SurfaceFlux::parse("topographyFlux");
  if (topg_flux_ptr == NULL)
    {
      topg_flux_ptr = new zeroFlux();
    }

   if (a_wrapper.m_topography_flux)
     {
       AxbyFlux* ptr = new AxbyFlux(1.0, a_wrapper.m_topography_flux, 1.0, topg_flux_ptr);	
       amrObject.setTopographyFlux(ptr);
       delete(ptr);
     }
   else
     {
      amrObject.setTopographyFlux(topg_flux_ptr);
}

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
  else if (problem_type == "MemoryLevelData")
    {
      //read geometry from somewhere in memory and store in a LevelDataIBC
      CH_assert(!a_wrapper.m_geometry_ice_thickness.isNull());
      CH_assert(!a_wrapper.m_geometry_bedrock_elevation.isNull());
      LevelDataIBC* ptr = new LevelDataIBC
	(a_wrapper.m_geometry_ice_thickness, a_wrapper.m_geometry_bedrock_elevation, a_wrapper.m_geometry_dx);
      thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
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
      LevelDataIBC* ptr = new LevelDataIBC(levelThck,levelTopg,levelDx);
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

  IceInternalEnergyIBC* temperatureIBC = NULL;
  ParmParse tempPP("temperature");
  std::string tempType("constant");
  tempPP.query("type",tempType);
  if (tempType == "constant")
    {
      Real T = 258.0;
      tempPP.query("value",T);
      ConstantIceTemperatureIBC* ptr = new ConstantIceTemperatureIBC(T);
      temperatureIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
    }
  else if (tempType == "LevelData")
    {
      ParmParse ildPP("inputLevelData");
      LevelDataTemperatureIBC* ptr = NULL;
      ptr = LevelDataTemperatureIBC::parse(ildPP); CH_assert(ptr != NULL);
      temperatureIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
    }
#ifdef HAVE_PYTHON
  else if (tempType == "Python")
    {
      ParmParse pyPP("PythonIceTemperatureIBC");
      std::string module;
      pyPP.get("module",module);
      std::string funcName = "temperature";
      pyPP.query("function",funcName);
      temperatureIBC  = static_cast<IceInternalEnergyIBC*>
	(new PythonInterface::PythonIceTemperatureIBC(module, funcName));
    }
#endif
  else 
    {
      MayDay::Error("bad temperature type");
    }	
	
  amrObject.setInternalEnergyBC(temperatureIBC);
    
  {
    // ---------------------------------------------
    // set surface heat boundary data 
    // ---------------------------------------------
        
    ParmParse pps("surfaceHeatBoundaryData");
    std::string surface_heat_type = "";
    pps.query("type",surface_heat_type);
    
    if (surface_heat_type == "MemoryLevelData")
      {
	//first try the in-memory interface data
	CH_assert(a_wrapper.m_surface_heat_boundary_data != NULL);
	amrObject.setSurfaceHeatBoundaryData(a_wrapper.m_surface_heat_boundary_data, 
					     a_wrapper.m_surface_heat_boundary_dirichlett, true);
      }
    else
      {
	//try the standalone BISICLES boundary data options
	bool diri = true; //Dirichlett boundary data by default
	pps.query("Dirichlett",diri);
	bool temp = true; //Dirichlett boundary data is temperature by default
	pps.query("Dirichlett",temp);
	SurfaceFlux* surf_heat_boundary_data_ptr =
	  SurfaceFlux::parse("surfaceHeatBoundaryData");
	if ( surf_heat_boundary_data_ptr == NULL)
	  {
	    if (!diri)
	      {
		surf_heat_boundary_data_ptr = new zeroFlux();
	      }
	  }
	amrObject.setSurfaceHeatBoundaryData(surf_heat_boundary_data_ptr, diri, temp);
	if (surf_heat_boundary_data_ptr != NULL)
	  {
	    delete surf_heat_boundary_data_ptr;
	    surf_heat_boundary_data_ptr=NULL;
	  }
      }
  
    
      // ---------------------------------------------
      // set basal (lower surface) heat boundary data. 
      // ---------------------------------------------
    ParmParse ppb("basalHeatBoundaryData");
    std::string basal_heat_type = "";
    ppb.query("type",basal_heat_type);
    if (basal_heat_type == "MemoryLevelData")
      {
	CH_assert(a_wrapper.m_basal_heat_boundary_data != NULL);
	amrObject.setBasalHeatBoundaryData(a_wrapper.m_basal_heat_boundary_data);//first try the in-memory interface data
      }
    else
      {
	//try the standalone BISICLES boundary data options
	SurfaceFlux* basal_heat_boundary_data_ptr = 
	  SurfaceFlux::parse("basalHeatBoundaryData");
	if ( basal_heat_boundary_data_ptr == NULL)
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
  }

     
  amrObject.setDomainSize(domainSize);

  CalvingModel* calving_model_ptr = CalvingModel::parseCalvingModel("CalvingModel");
  if (calving_model_ptr == NULL)
    {
      calving_model_ptr = new NoCalvingModel;
    }
  amrObject.setCalvingModel(calving_model_ptr);

  
  {
    /// initialize the damge model
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
	a_wrapper.m_amrMelange = & ptr->melange();
	amrObject.addObserver(ptr);
      }
  }


  // set up initial grids, initialize data, etc. 
  amrObject.initialize();

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

  if (temperatureIBC != NULL)
    {  
      delete temperatureIBC;
      temperatureIBC=NULL;
    }

   if (calving_model_ptr != NULL)
      {
	delete calving_model_ptr;
	calving_model_ptr = NULL;
      }

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif


}  
void advance_bisicles_instance(BisiclesWrapper* wrapper_ptr,  double start_time, double max_time, int max_step  )
{
  if (wrapper_ptr != NULL)
    {
      wrapper_ptr->m_amrIce.setTime(start_time);
      wrapper_ptr->m_amrIce.run(max_time,max_step);
    }
}

namespace bisicles_c_wrapper
{

  /// On-the-fly definition of a DisjointBoxLayout which assumes level data arranged as one  rectangular box per processor
  /** assume that each process contributes one rectangular block of data and do it will call this function once per data field.
      If a processor has no data, then specify a box outside of the problem domain 
      we will need to think again if we want anything flash
  */
  void defineDBL(DisjointBoxLayout& a_dbl,  const int *dims, const int *boxlo, const int *boxhi)
  {
    IntVect lo; D_TERM(lo[0] = boxlo[0];, lo[1] = boxlo[1];, lo[2] = boxlo[2]);
    IntVect hi; D_TERM(hi[0] = boxhi[0];, hi[1] = boxhi[1];, hi[2] = boxhi[2]);

    Box domBox;
    if (procID() == uniqueProc(SerialTask::compute))
      {
	IntVect dlo; D_TERM(dlo[0] = 0;, dlo[1] = 0;, dlo[2] = 0);
	IntVect dhi; D_TERM(dhi[0] = dims[0] - 1;, dhi[1] = dims[1] - 1;, dhi[2] = dims[2]-1);
	domBox.define(Box(dlo,dhi));
      }
    broadcast(domBox,uniqueProc(SerialTask::compute));
   
    //first off, accumulate the boxes from all the processors
    Vector<Box> tboxes;
    gather(tboxes, Box(lo,hi),  uniqueProc(SerialTask::compute));
    Vector<int> tprocIDs;
    gather(tprocIDs, procID(), uniqueProc(SerialTask::compute));
    

    //now discard boxes outside the domain 
    Vector<Box> boxes;
    Vector<int> procIDs;
    if (procID() == uniqueProc(SerialTask::compute))
      {
	for (int i = 0; i < tboxes.size(); i++)
	  {
	    if (tboxes[i].intersects(domBox))
	      {
		boxes.push_back(tboxes[i]);
		procIDs.push_back(tprocIDs[i]);
	      }
	  }
      }

   
      
    broadcast(boxes,uniqueProc(SerialTask::compute));
    broadcast(procIDs,uniqueProc(SerialTask::compute));
    a_dbl.define(boxes, procIDs);

    
    pout() << "proc: " << procID() << " " << domBox << " ";
    for (int i = 0; i < procIDs.size(); i++)
      {
	pout() <<  procIDs[i] << " " << boxes[i] << endl;
      }

  }
}


void bisicles_set_2d_geometry
(BisiclesWrapper* wrapper_ptr,   double *thck_data_ptr, double *topg_data_ptr, 
 const double *dx, const int *dims, const int *boxlo, const int *boxhi)
{
  if (wrapper_ptr != NULL)
    {
      DisjointBoxLayout dbl;
      bisicles_c_wrapper::defineDBL(dbl,dims,boxlo,boxhi);
      RefCountedPtr<LevelData<FArrayBox> > thck_ptr(new LevelData<FArrayBox> (dbl, 1, IntVect::Zero));
      RefCountedPtr<LevelData<FArrayBox> > topg_ptr(new LevelData<FArrayBox> (dbl, 1, IntVect::Zero));
      DataIterator dit(dbl);
      dit.reset();

      if (dit.ok())
	{
	 
	  (*thck_ptr)[dit].define(dbl[dit], 1, thck_data_ptr) ;
	 
	  (*topg_ptr)[dit].define(dbl[dit], 1, topg_data_ptr) ;
	}
      RealVect dxv; D_TERM(dxv[0] = dx[0];, dxv[1] = dx[1];, dxv[2] = dx[2]);

      wrapper_ptr->m_geometry_ice_thickness = thck_ptr;
      wrapper_ptr->m_geometry_bedrock_elevation = topg_ptr;
      wrapper_ptr->m_geometry_dx = dxv;
    }
}

void bisicles_set_2d_data
(BisiclesWrapper* wrapper_ptr,   double *data_ptr, const int *field, 
 const double *dx, const int *dims, const int *boxlo, const int *boxhi)
{
  if (wrapper_ptr != NULL)
    {
      DisjointBoxLayout dbl;
      bisicles_c_wrapper::defineDBL(dbl,dims,boxlo,boxhi);
     
      RefCountedPtr<LevelData<FArrayBox> > ptr(new LevelData<FArrayBox> (dbl, 1, IntVect::Zero));
      DataIterator dit(dbl);
      dit.reset();

      if (dit.ok())
	{
	  
	  (*ptr)[dit].define(dbl[dit], 1, data_ptr) ;
	}

      RealVect dxv; D_TERM(dxv[0] = dx[0];, dxv[1] = dx[1];, dxv[2] = dx[2]);

      switch (*field)
	{
	case BISICLES_FIELD_SURFACE_FLUX:
	  wrapper_ptr->m_surface_flux = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  break;
	case BISICLES_FIELD_BASAL_FLUX:
	  wrapper_ptr->m_basal_flux = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  break;
	case BISICLES_FIELD_FLOATING_ICE_BASAL_FLUX:
	  wrapper_ptr->m_floating_ice_basal_flux = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  break;
	case BISICLES_FIELD_GROUNDED_ICE_BASAL_FLUX:
	  wrapper_ptr->m_grounded_ice_basal_flux = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  break;
	case BISICLES_FIELD_SURFACE_TEMPERATURE:
	  wrapper_ptr->m_surface_heat_boundary_data = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  wrapper_ptr->m_surface_heat_boundary_dirichlett = true;
	  break;
	case BISICLES_FIELD_SURFACE_HEAT_FLUX:
	  wrapper_ptr->m_surface_heat_boundary_data = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  wrapper_ptr->m_surface_heat_boundary_dirichlett = false;
	  break;
	case BISICLES_FIELD_BASAL_HEAT_FLUX:
	  wrapper_ptr->m_basal_heat_boundary_data = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  break;
	case BISICLES_FIELD_TOPOGRAPHY_FLUX:
	  wrapper_ptr->m_topography_flux = new LevelDataSurfaceFlux(ptr, dxv, 0.0);
	  break;
	

	default: 
	  MayDay::Error("bisicles_set_2d_data: unknown (or unimplemented) field");
	}
      

    }
}




void bisicles_get_2d_data
(BisiclesWrapper* wrapper_ptr,   double *data_ptr, const int *field, 
 const double *dx, const int *dims, const int *boxlo, const int *boxhi)
{
  if (wrapper_ptr != NULL)
    {
      DisjointBoxLayout dbl;
      bisicles_c_wrapper::defineDBL(dbl,dims,boxlo,boxhi);
      RefCountedPtr<LevelData<FArrayBox> > ptr(new LevelData<FArrayBox> (dbl, 1, IntVect::Zero));

      DataIterator dit(dbl);
      dit.reset();
      if (dit.ok())
	{
	  
	  (*ptr)[dit].define(dbl[dit], 1, data_ptr) ;
	}

      RealVect dxv; D_TERM(dxv[0] = dx[0];, dxv[1] = dx[1];, dxv[2] = dx[2]);
      AmrIce& amrIce = wrapper_ptr->m_amrIce;
      int n = amrIce.finestLevel() + 1;
      Vector<LevelData<FArrayBox>* > data(n);
      Vector<RealVect> amrDx(n);

      switch (*field)
	{
	case BISICLES_FIELD_SURFACE_ELEVATION:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(&(amrIce.geometry(lev)->getSurfaceHeight()));
	      amrDx[lev] = amrIce.dx(lev);
	    }

	  flattenCellData(*ptr,dxv,data,amrDx,true);
	  break;
	  
	case BISICLES_FIELD_BEDROCK_ELEVATION:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(&(amrIce.geometry(lev)->getTopography()));
	      amrDx[lev] = amrIce.dx(lev);
	    }
	  
	  flattenCellData(*ptr,dxv,data,amrDx,true);	
	  
	  break;

	case BISICLES_FIELD_ICE_THICKNESS:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(&(amrIce.geometry(lev)->getH()));
	      amrDx[lev] = amrIce.dx(lev);
	    }
	  
	  flattenCellData(*ptr,dxv,data,amrDx,true);	
	  
	  break;

	case BISICLES_FIELD_ICE_THICKNESS_HARMONIC:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(&(amrIce.geometry(lev)->getH()));
	      amrDx[lev] = amrIce.dx(lev);
	    }
	  
	  flattenCellData(*ptr,dxv,data,amrDx,true, CoarseAverage::harmonic);	
	  
	  break;

	 case BISICLES_FIELD_MELANGE_THICKNESS:

	   if (wrapper_ptr->m_amrMelange)
	     {
	       for (int lev = 0; lev < n ; lev++)
		 {
		   data[lev] = const_cast<LevelData<FArrayBox>* >
		     (  wrapper_ptr->m_amrMelange->melangeThickness(lev) );
		   
	       amrDx[lev] = amrIce.dx(lev);
		 }
	  
	       flattenCellData(*ptr,dxv,data,amrDx,true);	
	  
	       
	     }
	   else
	     {
	       MayDay::Warning("bisicles_get_2d_data: no melange model");
	     }

	   break;
	   
	case BISICLES_FIELD_SURFACE_TEMPERATURE:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(amrIce.surfaceInternalEnergy()[lev]);
	      amrDx[lev] = amrIce.dx(lev);
	    }
	  
	  flattenCellData(*ptr,dxv,data,amrDx,true);	
	  //convert internal energy to temperature
	  for (DataIterator dit( (*ptr).disjointBoxLayout()); dit.ok(); ++dit)
	    {
	      FArrayBox& T = (*ptr)[dit];
	      FArrayBox E(T.box(),1);
	      FArrayBox w(T.box(),1);
	      FArrayBox p(T.box(),1);
	      p.setVal(0.0);
	      E.copy(T);
	      IceThermodynamics::decomposeInternalEnergy(T,w,E,p,T.box());
	    }

	  break;
	  
	case BISICLES_FIELD_SURFACE_HEAT_FLUX:
	  
	  for (int lev = 0; lev < n ; lev++)
	    {
	      data[lev] = const_cast<LevelData<FArrayBox>* >(amrIce.surfaceHeatFlux()[lev]);
	      amrDx[lev] = amrIce.dx(lev);
	    }
	  
	  flattenCellData(*ptr,dxv,data,amrDx,true);	
	  
	  break;


	default: 
	  MayDay::Error("bisicles_get_2d_data: unknown (or unimplemented) field");
	}
     

    }
}




void bisicles_new_instance(int *instance_key, const char *input_fname, MPI_Comm mpi_comm)
{
  //fake argc
#define NARGS 2

  BisiclesWrapper* ptr = new BisiclesWrapper;
#ifdef CH_MPI
  CH_assert(mpi_comm);
  ptr->mpi_comm = mpi_comm;
  Chombo_MPI::comm = mpi_comm;
  
  int rank, number_procs;
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);

#ifdef CH_USE_PETSC
  PETSC_COMM_WORLD = Chombo_MPI::comm;
  PetscInitialize(PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
#endif


  MPI_Barrier(Chombo_MPI::comm);
#endif
  CH_assert(ptr != NULL);

  #define NARGS 2
  int argc = NARGS;
  char *argv[NARGS];
  char argv0[] = "cwrapper";
  argv[0] = argv0;
  char argv1[] = "drivel";
  argv[1] = argv1;
  ptr->m_pp = new ParmParse(argc-2,argv+2,NULL,input_fname);

  if  (bisicles_c_wrapper::instances.size() == 0)
    {
      *instance_key = 0;
    }
  else
    {
      *instance_key  = (--bisicles_c_wrapper::instances.end())->first + 1;
    }
  bisicles_c_wrapper::instances[*instance_key] = ptr;

  CH_assert(bisicles_c_wrapper::instances.size() == 1); //just one instance for now. need at least to fix petsc init otherwise
  
}


void bisicles_init_instance(int *instance_key)
{
  if (instance_key)
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL)
	    {
	      init_bisicles_instance( *(i->second) );
	    }
	}
    } 
}

void bisicles_free_instance(int *instance_key)
{
  
  if (instance_key)
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL)
	    {
	      delete i->second;
	      i->second = NULL;
	      bisicles_c_wrapper::instances.erase(i->first);
	    }
	}
    }
#ifdef CH_USE_PETSC
  int ierr = PetscFinalize(); 
#endif

}


void bisicles_advance(int *instance_key, double *start_time,  double *max_time, int *max_step)
{
  if (instance_key && max_time && max_step && start_time)
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  advance_bisicles_instance(i->second,  *start_time, *max_time, *max_step);
	}
    }
}

void bisicles_set_2d_geometry(int *instance_key,  double *thck_data_ptr, double *topg_data_ptr, 
			      const double *dx, const int *dims, 
			      const int *boxlo, const int *boxhi)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  bisicles_set_2d_geometry(i->second, thck_data_ptr, topg_data_ptr, dx, dims, boxlo, boxhi);
	}
    }
}


void bisicles_set_2d_data(int *instance_key,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  bisicles_set_2d_data(i->second, data_ptr, field, dx, dims, boxlo, boxhi);
	}
    }
}

void bisicles_get_2d_data(int *instance_key,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  bisicles_get_2d_data(i->second, data_ptr, field, dx, dims, boxlo, boxhi);
	}
    }
}

///write a checkpoint file 
void bisicles_write_checkpoint(int *instance_key)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL)
	    {
	      AmrIce& amrIce = i->second->m_amrIce;
	      amrIce.writeCheckpointFile();
	    }
	}
    }
}

///write a plot file 
void bisicles_write_plot(int *instance_key)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL)
	    {
	      AmrIce& amrIce = i->second->m_amrIce;
	      amrIce.writePlotFile();
	    }
	}
    }
}



///read a checkpoint file 
void bisicles_read_checkpoint(int *instance_key, const char *checkpoint_fname)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL && checkpoint_fname != NULL)
	    {
	      AmrIce& amrIce = i->second->m_amrIce;
	      std::string s = checkpoint_fname;
	      amrIce.restart(std::string(s));
	    }
	}
    }
}

///set header data
template<typename T>
void bisicles_set_header(int *instance_key, const char *attr_key, const T *val)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL && attr_key != NULL && val != NULL)
	    {
	      AmrIce& amrIce = i->second->m_amrIce;
	      std::string s(attr_key);
	      amrIce.setHeader(attr_key, *val);
	    }
	}
    }
}

///get header data
template<typename T>
void bisicles_get_header(int *instance_key, const char *attr_key, T *val)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL && attr_key != NULL && val != NULL)
	    {
	      AmrIce& amrIce = i->second->m_amrIce;
	      amrIce.getHeader(attr_key, *val);
	    }
	}
    }
}
 
void bisicles_get_header(int *instance_key, const char *attr_key, char *val)
{
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  if (i->second != NULL && attr_key != NULL && val != NULL)
	    {
	      AmrIce& amrIce = i->second->m_amrIce;
	      std::string v(val);
	      int l = v.size();
	      amrIce.getHeader(attr_key, v);
	      strncpy(val,v.c_str(),l); //should pad val with 0 if l > v.size()
	    }
	}
    }
}


 
void bisicles_set_header_dble(int *instance_key, const char* attr_key, const double *val)
{
  bisicles_set_header(instance_key, attr_key, val);
}
void bisicles_set_header_char(int *instance_key, const char* attr_key, const char *val)
{
  bisicles_set_header(instance_key, attr_key, val);
}
  
void bisicles_get_header_int(int *instance_key, const char* attr_key, int *val)
{
  bisicles_get_header(instance_key, attr_key, val);
}  

void bisicles_get_header_dble(int *instance_key, const char* attr_key, double *val)
{
  bisicles_get_header(instance_key, attr_key, val);
}

void bisicles_get_header_char(int *instance_key, const char* attr_key, char *val)
{
  bisicles_get_header(instance_key, attr_key, val);
}

BisiclesWrapper* bisicles_instance(int *instance_key)
{
  BisiclesWrapper* rc = NULL;
  if (instance_key) ///\todo : check all pointers
    {
      std::map<int, BisiclesWrapper*>::iterator i 
	= bisicles_c_wrapper::instances.find(*instance_key) ;
      if (i != bisicles_c_wrapper::instances.end())
	{
	  rc = i->second;
	}
    }
  return rc;
}

void bisicles_push_pop_thin_ice(BisiclesWrapper* wrapper_ptr, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				 const int *boxlo, const int *boxhi, bool pop)

{
  if (wrapper_ptr)
    {
      pout() << "bisicles_push_pop_thin_ice!" << endl;

      
      AmrIce& amrIce = wrapper_ptr->m_amrIce;
      
      int finest_level = amrIce.finestLevel();
      
      Vector<LevelData<FArrayBox>* > delta_h(finest_level + 1);
      Vector<RealVect> amrDx(finest_level + 1);
      // the uniform mesh data storted at data_ptr provides an upper bound on the thickness to be added
      DisjointBoxLayout dbl;
      bisicles_c_wrapper::defineDBL(dbl,dims,boxlo,boxhi);
      LevelData<FArrayBox> max_extra_thk(dbl, 1, IntVect::Zero);
      DataIterator dit(dbl);
      dit.reset();
      CH_assert(dit.ok());
      pout() <<  *data_ptr << endl;
      max_extra_thk[dit].define(dbl[dit], 1, data_ptr);
      pout() << "bisicles_push_pop_thin_ice  max_extra_thk[dit] : " << max_extra_thk[dit].max() << "," << max_extra_thk[dit].max() << endl; 
	

      //need a vector version of dx
      RealVect dxv; D_TERM(dxv[0] = dx[0];, dxv[1] = dx[1];, dxv[2] = dx[2]);

      // compute the thickness increment, store in new_thck
      for (int lev = 0; lev <= finest_level; lev++)
	{
	  const DisjointBoxLayout& grids = amrIce.grids(lev);
	  LevelData<FArrayBox> level_max_extra_thk(grids,1,IntVect::Zero);
	  const LevelData<FArrayBox>& thk = amrIce.geometry(lev)->getH();
	  const LevelData<BaseFab<int> >& mask = amrIce.geometry(lev)->getFloatingMask();
	  delta_h[lev] = new LevelData<FArrayBox>(grids,thk.nComp(),thk.ghostVect());
	  amrDx[lev] = amrIce.dx(lev);
	  FillFromReference(level_max_extra_thk, max_extra_thk, amrDx[lev] , dxv, true);

	  for (DataIterator dit(grids); dit.ok(); ++dit)
	    {
	      FArrayBox& dh = (*delta_h[lev])[dit];
	      dh.setVal(0.0);
	      const FArrayBox& h =  (thk)[dit];
	      FArrayBox& dh_max =  level_max_extra_thk[dit];
	      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		  bool ground = (mask[dit](iv) == FLOATINGMASKVAL) ||  (mask[dit](iv) == GROUNDEDMASKVAL);
		  if ( h(iv) < *thin_ice_limit && ground)
		    {
		      if (pop)
			{
			  dh(iv) = -h(iv);
			}
		      else
			{
			  dh(iv) = dh_max(iv);
			}
		    }
		}
	    }
	}
      
      wrapper_ptr->m_amrIce.incrementIceThickness(delta_h);

      //finally, flatten dh into max_extra_thk, and, hence, the input data_ptr
      flattenCellData(max_extra_thk,dxv,delta_h,amrDx,true);
      
      //clean up
      for (int lev = 0; lev <= finest_level; lev++)
	{
	  if (delta_h[lev])
	    {
	      delete delta_h[lev]; delta_h[lev] = NULL;
	    }
	}
      
    }

  pout() << "!bisicles_push_pop_thin_ice" << endl;
}


void bisicles_push_thin_ice(BisiclesWrapper* wrapper_ptr, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)
{
  bisicles_push_pop_thin_ice(wrapper_ptr, data_ptr, thin_ice_limit, dx, dims, boxlo, boxhi, false);
}


void bisicles_push_thin_ice(int *instance_key, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)

{

  BisiclesWrapper* instance = bisicles_instance(instance_key);
  if (instance) 
    {
      bisicles_push_thin_ice(instance, data_ptr, thin_ice_limit, dx, dims, boxlo, boxhi);
    }
}

void bisicles_pop_thin_ice(BisiclesWrapper* wrapper_ptr, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)

{
  bisicles_push_pop_thin_ice(wrapper_ptr, data_ptr, thin_ice_limit, dx, dims, boxlo, boxhi, true);
}

void bisicles_pop_thin_ice(int *instance_key, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)
{
  
  BisiclesWrapper* instance = bisicles_instance(instance_key);
  if (instance) 
    {
      bisicles_pop_thin_ice(instance, data_ptr, thin_ice_limit, dx, dims, boxlo, boxhi);
    }
  
}

void f_bisicles_push_thin_ice(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)
{
  bisicles_push_thin_ice(intance_id, data_ptr, thin_ice_limit,dx,dims, boxlo, boxhi);
}
  
void f_bisicles_pop_thin_ice(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)
{
  bisicles_pop_thin_ice(intance_id, data_ptr, thin_ice_limit,dx,dims, boxlo, boxhi);
}

void f_bisicles_push_thin_ice_(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)
{
  bisicles_push_thin_ice(intance_id, data_ptr, thin_ice_limit,dx,dims, boxlo, boxhi);
}
  
void f_bisicles_pop_thin_ice_(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi)
{
  bisicles_pop_thin_ice(intance_id, data_ptr, thin_ice_limit,dx,dims, boxlo, boxhi);
}
