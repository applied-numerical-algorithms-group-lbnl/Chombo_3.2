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
// ncio.cpp
// read/write FABs from netcdf files
//===========================================================================

#include "fabncio.H"
#include "netcdf.h"

#ifdef HAVE_GDAL
#include "ogr_spatialref.h"
#endif

#include "NamespaceHeader.H"


#ifdef HAVE_NETCDF
/// utility functions needed by writeFAB/readFAB
namespace NCIO
{
  /// write CF compliant spatial ref (ie projection) attributes
  void writeCFSpatialRef(int ncID, int a_epsg);

  /// define netcdf dimensions and variables needed to store a_dd
  void defineCF(int ncID, const DomainDiagnosticData& a_dd);

  /// write a_dd to netcdf
  void writeCF(int ncID, const DomainDiagnosticData& a_dd);
}

/// attempt to write CF projection data to a netcdf group, built from an epsg code
void NCIO::writeCFSpatialRef(int ncID, int a_epsg)
{

  int rc, crsID;
  rc = nc_def_var(ncID, "crs", NC_INT, 0, 0, &crsID);
  
  
  if (a_epsg <= 0)
    {
      rc = nc_put_att_text(ncID, crsID, "crs_wkt", 7, "unknown");
    }
  else
    {
     
      // write epsg, in case we don't have gdal
      rc = nc_put_att_int(ncID, crsID, "EPSG", NC_INT, 1, &a_epsg);
      
#ifdef HAVE_GDAL
      //projection
      
      OGRSpatialReference crs;
      crs.importFromEPSG( a_epsg );
     
      // Allowed to include WKT in CF, not sure how useful that is
      char* wkt;
      crs.exportToWkt( &wkt );
      rc = nc_put_att_text(ncID, crsID, "crs_wkt", strlen(wkt), wkt);

      // gdal's netcdf driver apparently can work with wkt (as well as CF)
      rc = nc_put_att_text(ncID, NC_GLOBAL, "spatial_ref", strlen(wkt), wkt);

      // attempt to construct CF grid_mapping
      const char *projection = crs.GetAttrValue("PROJECTION");
      if (projection)
	{
	  /// we might want some other projections at some point
	  if (EQUAL(projection, SRS_PT_POLAR_STEREOGRAPHIC))
	    {
	      rc = nc_put_att_text(ncID, crsID, "grid_mapping_name", 19, "polar_stereographic");

	      // cf seems to want latitude_of_projection_origin
	      //but it isn't really part of the projection definition
	      double lat_0 = crs.GetProjParm(SRS_PP_LATITUDE_OF_ORIGIN);
	      lat_0 = (lat_0 > 0.0)?90.0:-90.0;
	      int rc = nc_put_att_double(ncID, crsID, "latitude_of_projection_origin", NC_DOUBLE, 1, &lat_0);

	      auto fw = [ncID,crsID,&crs](const char* a_gdal_name,
					  const char* a_cf_name ){
		double d = crs.GetProjParm(a_gdal_name,0.0);
		int rc = nc_put_att_double(ncID, crsID, a_cf_name, NC_DOUBLE, 1, &d);
	      };
	      fw(SRS_PP_CENTRAL_MERIDIAN,SRS_PP_CENTRAL_MERIDIAN);
	      fw(SRS_PP_SCALE_FACTOR,SRS_PP_SCALE_FACTOR);
	      fw(SRS_PP_LATITUDE_OF_ORIGIN,"standard_parallel");
	      fw(SRS_PP_FALSE_EASTING,SRS_PP_FALSE_EASTING);
	      fw(SRS_PP_FALSE_NORTHING,SRS_PP_FALSE_NORTHING);	       
	     }
	}
#endif
      
    } // end if a_epsg > 0


}

void NCIO::defineCF(int ncID, const DomainDiagnosticData& a_dd)
{

  int rc;

  //(timestep) dimension
  const Vector<Real>* time = a_dd.m_cf_stuff[0].data;
  if (time && time->size() > 0)
    {
  
      int tdimID[1];
      int one = 1;
      rc = nc_def_dim(ncID, "ts", time->size(),  &tdimID[0]);

      for (int i = 0; i < a_dd.m_cf_stuff.size(); i++)
	{
	  pout() << "Diagnostic info: " << a_dd.m_cf_stuff[i].cf_name 
		 << ", " << a_dd.m_cf_stuff[i].units 
		 << ", " << a_dd.m_cf_stuff[i].long_name 
		 << endl;
	}

      for (int i = 0; i < a_dd.m_cf_stuff.size(); ++i)
	{
	  int varID;
	  rc = nc_def_var(ncID, a_dd.m_cf_stuff[i].cf_name.c_str(), NC_DOUBLE, one ,
       			  &tdimID[0], &varID);
	  const char* c = a_dd.m_cf_stuff[i].cf_name.c_str();
	  rc = nc_put_att_text(ncID, varID, "standard_name", strlen(c), c); 
	  c = a_dd.m_cf_stuff[i].units.c_str();
	  rc = nc_put_att_text(ncID, varID, "units", strlen(c), c); 
	  c = a_dd.m_cf_stuff[i].long_name.c_str();
	  rc = nc_put_att_text(ncID, varID, "long_name", strlen(c), c); 

	}
      
    }
}

void NCIO::writeCF(int ncID, const DomainDiagnosticData& a_dd)
{


  int rc;
  
  for (int i = 0; i < a_dd.m_cf_stuff.size(); ++i)
    {
      int varID;
      rc = nc_inq_varid(ncID, a_dd.m_cf_stuff[i].cf_name.c_str(), &varID);
	
      const Vector<Real>* v = a_dd.m_cf_stuff[i].data;
      if (v && v->size() > 0)
	{
	  size_t start = 0;
	  size_t count = v->size();
	  rc = nc_put_vara_double(ncID, varID, &start, &count, &(*v)[0]);
	}
    }
 
}

void NCIO::writeFAB(const std::string& a_file,
		    const Vector<std::string>& a_names,
		    const Vector<std::string>& a_cf_standard_names,
		    const Vector<std::string>& a_cf_units,
		    const Vector<std::string>& a_cf_long_names,
		    const FArrayBox& a_fab, 
		    const Real& a_dx,
		    const RealVect& a_x0,
		    int a_epsg,
		    const DomainDiagnosticData& a_dd,
		    const std::string& a_flattenInfo,
		    const HDF5HeaderData& a_file_header)
{
  int rc; int ncID;
  if ( (rc = nc_create(a_file.c_str(), NC_CLOBBER, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }


  //global attributes
  rc = nc_put_att_text(ncID, NC_GLOBAL, "Conventions", 6, "CF-1.7");
  for (auto i = a_file_header.m_string.begin(); i != a_file_header.m_string.end(); ++i)
    {
      string sname = i->first;
      if (sname.find("component_00") == string::npos)
	{
	  rc = nc_put_att_text(ncID, NC_GLOBAL, i->first.c_str(), i->second.size(), i->second.c_str());
	}
    }
  for (auto i = a_file_header.m_int.begin(); i != a_file_header.m_int.end(); ++i)
    {
      int iatt = i->second;
      string iname = i->first;
      if (iname.find("num_comps") == string::npos)
	{
	  rc = nc_put_att_int(ncID, NC_GLOBAL, i->first.c_str(), NC_INT, 1, &iatt);
	}
    }
  for (auto i = a_file_header.m_real.begin(); i != a_file_header.m_real.end(); ++i)
    {
      Real ratt = i->second;
      rc = nc_put_att_double(ncID, NC_GLOBAL, i->first.c_str(), NC_DOUBLE, 1, &ratt);
    }
  rc = nc_put_att_double(ncID, NC_GLOBAL, "dx", NC_DOUBLE, 1, &a_dx);

  rc = nc_put_att_text(ncID, NC_GLOBAL, "Conversion_history", a_flattenInfo.size(), a_flattenInfo.c_str());

  //spatial reference
  writeCFSpatialRef(ncID, a_epsg);

  //define diagnostic data
  defineCF(ncID, a_dd);
  
  // fab data
  // dimensions
  std::string xname[SpaceDim] = {D_DECL("x","y","z")};
  int dimID[SpaceDim];

  int bufSize = 1;
  size_t start[SpaceDim] = {D_DECL(0,0,0)};
  size_t count[SpaceDim];
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      int n = 1 + a_fab.box().bigEnd()[dir] - a_fab.box().smallEnd()[dir];
      count[dir] = n;
      bufSize*=n;
      if ( (rc = nc_def_dim(ncID, xname[dir].c_str(), n, &dimID[dir])) != NC_NOERR)
	{
	  MayDay::Error("failed to define dimension");
	}
    }

  // x,y,z data
  for (int dir = 0; dir < SpaceDim; ++dir)
  {
    
    int varID;
    int one = 1;
    if ( (rc = nc_def_var(ncID, xname[dir].c_str(), NC_DOUBLE,
			  one, &dimID[dir], &varID)) != NC_NOERR)
      {
	MayDay::Error("failed to define x");
      }

    string s = "projection_"; s.append(xname[dir]); s.append("_coordinate");
    const char* c = s.c_str();
    rc = nc_put_att_text(ncID, varID, "standard_name", strlen(c), c);
    rc = nc_put_att_text(ncID, varID, "units", 1, "m"); // BISICLES x unit is m.
    
    
  }
  int FORTRAN_dimID[SpaceDim];
  size_t FORTRAN_start[SpaceDim];
  size_t FORTRAN_count[SpaceDim];
  for (int i = 0 ; i < SpaceDim ; i++)
    {
      FORTRAN_dimID[i] = dimID[SpaceDim-i-1];
      FORTRAN_count[i] = count[SpaceDim-i-1];
      FORTRAN_start[i] = start[SpaceDim-i-1];
    }
  for (int i =0; i < a_names.size(); i++)
    {
      int varID;
      if ( (rc = nc_def_var(ncID, a_names[i].c_str(), NC_DOUBLE,
			    SpaceDim, FORTRAN_dimID, &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to define variable");
	}
       rc = nc_put_att_text(ncID, varID, "grid_mapping", 3, "crs" );

       if (a_cf_standard_names[i] != "")
	 {
	   const char* c = a_cf_standard_names[i].c_str();
	   rc = nc_put_att_text(ncID, varID, "standard_name", strlen(c), c); 
	   c = a_cf_units[i].c_str();
	   rc = nc_put_att_text(ncID, varID, "units", strlen(c), c); 
	   c = a_cf_long_names[i].c_str();
	   rc = nc_put_att_text(ncID, varID, "long_name", strlen(c), c); 
	 }
       else
	 {
	   const char* c = a_names[i].c_str();
	   rc = nc_put_att_text(ncID, varID, "standard_name", strlen(c), c);
	 }
       
    }

  if ( (rc = nc_enddef(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to define netcdf file");
    }
  
  //write diagnostic data
  writeCF(ncID, a_dd);

  // write fab data
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      int varID;
      if ( (rc = nc_inq_varid(ncID, xname[dir].c_str(), &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable id");
	}
      double *xptr = new double[count[dir]];

      xptr[0] =  a_x0[dir] + (Real(a_fab.box().smallEnd()[dir])+0.5)*a_dx;
      for (int i = 1; i < count[dir]; ++i)
	{
	  xptr[i] = xptr[i-1] + a_dx; 
	}
      if ( (rc = nc_put_vara_double(ncID, varID, &start[dir], &count[dir], xptr)) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}
      
      delete xptr;
    }

  double* dptr = new double[bufSize];
  for (int i =0; i < a_names.size(); i++)
    {
      int varID;
      if ( (rc = nc_inq_varid(ncID, a_names[i].c_str(), &varID)) != NC_NOERR)
	{
	  MayDay::Error("failed to find variable id");
	}
       
      Interval ivl(i,i);
      a_fab.linearOut((void*)dptr,a_fab.box(),ivl);
      if ( (rc = nc_put_vara_double(ncID, varID, FORTRAN_start, FORTRAN_count, dptr)) != NC_NOERR)
	{
	  MayDay::Error("failed to write data");
	}

    }

  if ( (rc = nc_close(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to close netcdf file");
    }
  if (dptr != NULL)
    delete[] dptr;
 
}

//construct a fab with n=a_var.size() components,
//filled with data loaded from a netcdf file
//with a one-dimensional variable x and 
//SpaceDim-dimensional variables var[0]-var[n-1]
void NCIO::readFAB(const std::string& a_file,
		const Vector<std::string>& a_var,
		FArrayBox& a_fab,
		Real& a_dx)
{

  int rc; int ncID;
  if ( (rc = nc_open(a_file.c_str(), NC_NOWRITE, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  {
    //determine a_dx
    int xID;
    if ( (rc = nc_inq_varid(ncID, "x", &xID) ) != NC_NOERR)
      {
	MayDay::Warning("nc_inq_varid failed to find id for x : setting dx = 5.0km");
	a_dx = 5.0e+3;
      }
    else
      {
	int nxdims;
	if ( (rc = nc_inq_varndims(ncID, xID, &nxdims) )  != NC_NOERR)
	  {
	    MayDay::Error("nc_inq_varndims failed to find ndims for x");
	  }
	if (nxdims != 1)
	  {
	    MayDay::Error("wrong dimensions in x");
	  }
	double x[2];
	size_t start[2] = {0,0};
	size_t count[2] = {2,2};
	if ( ( rc =  nc_get_vara_double(ncID, xID, start, count, x) ) != NC_NOERR)
	  {
	    MayDay::Error("nc_get_vara_double failed to read first two values of x");
	  }
	a_dx = x[1] - x[0];
      }
    if (a_dx <= 0.0)
      MayDay::Error("a_dx <= 0.0");
  }

  double* dptr = NULL;
  Vector<int> varID(a_var.size());
  size_t dimLength[SpaceDim];
  for (int i =0; i < a_var.size(); ++i)
    {
      if ( (rc = nc_inq_varid(ncID, a_var[i].c_str(), &varID[i]) ) != NC_NOERR)
	{
	  MayDay::Error("nc_inq_varid failed");
	}

      int ndims;
      if ( (rc = nc_inq_varndims(ncID, varID[i], &ndims ) ) != NC_NOERR)
	{
	  MayDay::Error("nc_inq_varndims failed");
	}
      if (ndims != SpaceDim)
	{
	  MayDay::Error("wrong dimensions in variable");
	}

      if (i == 0)
	{
	  int dimID[SpaceDim];

	  if ( (rc = nc_inq_vardimid(ncID, varID[i], dimID ) ) != NC_NOERR)
	    {
	      MayDay::Error("nc_inq_vardimid  failed");
	    }
	  
	  
	  
	  IntVect hi; int bufSize = 1;
	  for (int dir = 0; dir < SpaceDim; ++dir)
	    {
	      
	      if ( (rc =  nc_inq_dimlen  (ncID, dimID[dir], &dimLength[dir]) ) != NC_NOERR)
		{
		  MayDay::Error("nc_inq_dimlen failed");
		}
	     
	      bufSize *= dimLength[dir];
	    }
	  if (SpaceDim == 2)
	    {
	      hi[0] = dimLength[1]-1;
	      hi[1] = dimLength[0]-1;
	    }
	  else
	    {
	      MayDay::Error("2d only for now");
	    }
	  a_fab.define(Box(IntVect::Zero,hi), a_var.size());
	  
	  dptr = new double[bufSize];

	 

	} //end if (i == 0)
      size_t start[SpaceDim] = {D_DECL(0,0,0)};
      if ( ( rc =  nc_get_vara_double(ncID, varID[i], start, dimLength, dptr) ) != NC_NOERR)
	{
	  MayDay::Error("nc_get_vara_double failed");
	}
      a_fab.linearIn(dptr,a_fab.box(),Interval(i,i));

    }

  if (dptr != NULL)
    delete[] dptr;

 if ( (rc = nc_close(ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to close netcdf file");
    }
}

#endif
 

#include "NamespaceFooter.H"
