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
// UnstructuredIO.cpp
// Operations related to reading and writing data from/to valid regions
// of AMR Hierarchies
//===========================================================================
#include "UnstructuredIO.H"
#include "FieldNames.H"
#include "IntVectSet.H"
#include "CoarseAverage.H"
#include "NamespaceHeader.H"


UnstructuredData::UnstructuredData
(int a_nComp,  const RealVect& a_crseDx, const Box& a_crseDomain,
 const Vector<int>& a_ratio, const Real& a_time, const RealVect& a_x0)
{
  define(a_nComp, a_crseDx, a_crseDomain, a_ratio, a_time, a_x0);
}

void UnstructuredData::define
(int a_nComp,  const RealVect& a_crseDx, const Box& a_crseDomain,
 const Vector<int>& a_ratio, const Real& a_time, const RealVect& a_x0)
{
  m_nComp = a_nComp;
  m_crseDomain = a_crseDomain;
  m_ratio = a_ratio;
  m_time = a_time;
  m_x0 = a_x0;
  m_nLevel = m_ratio.size();
  m_dx.resize(m_nLevel);
  m_dx[0] = a_crseDx;
  for (int lev = 1; lev < m_nLevel; lev++)
    m_dx[lev] = m_dx[lev-1] / Real(m_ratio[lev-1]);
    
  m_field.resize(m_nComp);
  m_x.resize(SpaceDim);
  m_iv.resize(SpaceDim);
  m_isDefined = true;
}

void UnstructuredData::resize(int a_size)
{
  CH_TIME("UnstructuredData::resize");
  CH_assert(m_isDefined);

  m_level.resize(a_size);
 

  for (int dir = 0; dir < SpaceDim; dir++)
    {
      m_iv[dir].resize(a_size);
      m_x[dir].resize(a_size);
    }
  
  for (int ic = 0; ic < nComp() ; ic++)
    {
      m_field[ic].resize(a_size);
    }

}


void UnstructuredData::append
(int a_lev, const Box& a_box, const IntVect& a_iv,  const Vector<Real>& a_data)
{
  CH_TIME("UnstructuredData::append");

  CH_assert(m_isDefined);
  CH_assert(a_lev < m_nLevel);
  CH_assert(a_data.size() == m_field.size());

    m_level.push_back(a_lev);
    m_levelBoxSet.insert(std::make_pair(a_lev,a_box)); 

    for (int dir = 0; dir < SpaceDim; dir++)
      {
	m_iv[dir].push_back(a_iv[dir]);  
      }
   
    for (int ic = 0; ic < m_nComp; ic++)
      {
	m_field[ic].push_back(a_data[ic]);
      }

    for (int dir = 0; dir < SpaceDim; dir++)
      {
	m_x[dir].push_back(m_x0[dir] + m_dx[a_lev][dir]*(Real(a_iv[dir]) + 0.5));  
      }
  }

///fill  a_nodeCoord with the a_dir co-ordinate of node a_node. 
/**
   \param Vector<Real>& a_nodeCoord will be resized to this.size(). 
   \param a_node[dir] = 0 implies lo side,   1 implies high side
                        
*/
void UnstructuredData::computeNodeCoord(Vector<Real>& a_nodeCoord, int a_dir, const IntVect& a_node) const
{
  
  CH_TIME("UnstructuredData::computeNodeCoord");
  CH_assert(m_isDefined);
  if (a_nodeCoord.size() != nCell())
    {
      a_nodeCoord.resize(nCell());
    }

  CH_assert( a_node[a_dir] == 0 || a_node[a_dir]  == 1);

  Real f =  (a_node[a_dir] == 1)?0.5:-0.5;
  for (int i =0; i <  a_nodeCoord.size(); i++)
    {
      a_nodeCoord[i] = m_x[a_dir][i] + f * m_dx[m_level[i]][a_dir];
    }

}

//given unstructured data  produce block structured data
void UnstructuredIO::USToBS 
( Vector<LevelData<FArrayBox>*>& a_bsData, 
  const UnstructuredData& a_usData)
{
  CH_TIME("UnstructuredIO::UnstructuredtoBS");

  //need to start from scratch
  CH_assert(a_bsData.size() == 0);
  a_bsData.resize(a_usData.nLevel());

 
  for (int lev = 0; lev < a_usData.nLevel(); lev++)
    {
      //need to extract a disjoint box layout  
      Vector<Box> boxes;
      const UnstructuredData::LevelBoxSet& lbs = a_usData.levelBoxSet();
      for (UnstructuredData::LevelBoxSet::const_iterator it = lbs.begin(); it != lbs.end(); it++)
	{
	  if (it->first == lev)
	    {
	      boxes.push_back(it->second);
	    }
	}
     
      Vector<int> procID(0);
      DisjointBoxLayout grids(boxes, procID, ProblemDomain(a_usData.domain(lev)));
      
      a_bsData[lev] = new LevelData<FArrayBox>(grids, a_usData.nComp(), IntVect::Unit); }

  for (int lev =  a_usData.nLevel() - 1 ; lev >=0; lev--)
    {
      //start from the finest level, so that if there *is* data stored for the 
      //invalid regions, it overwrites the coarse average we carry out
      const DisjointBoxLayout& grids = a_bsData[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
      	{
	   //\todo : optimization. at the momemt, we are going through 
	  // all the data for every FAB.
	  FArrayBox& data = (*a_bsData[lev])[dit];
      	  for (int i =0; i <  a_usData.nCell(); i++)
      	    {
	      if (a_usData.level()[i] == lev)
		{
		  IntVect iv;
		  for (int dir =0; dir < SpaceDim; dir++)
		    {
		      iv[dir] = a_usData.iv(dir)[i];
		    }
		  if (data.box().contains(iv))
		    {
		      for (int comp = 0; comp < data.nComp(); comp++)
			{
			  data(iv,comp) = a_usData.field(comp)[i];
			}
		    }
		}
      	    }
      	}

      if (lev > 0)
	{
	  CoarseAverage av(grids , a_bsData[lev]->nComp(), a_usData.ratio()[lev-1]);
	  av.averageToCoarse(*a_bsData[lev-1],   *a_bsData[lev] );
	}

    } // end loop over levels
}


//given block structured data  produce unstructured data 
void UnstructuredIO::BStoUS
( UnstructuredData& a_usData , 
  const Vector<LevelData<FArrayBox>*>& a_bsData, 
  const Vector<int>& a_ratio, 
  bool a_validonly)
{
  CH_TIME("UnstructuredIO::BStoUnstructured");

  int numLevels = a_bsData.size();
  Vector<LevelData<BaseFab<int> >* > mask(numLevels,NULL);
  if (a_validonly)
    {
      //build valid data mask
      for (int lev = numLevels - 1; lev >= 0; lev--)
	{
	  const DisjointBoxLayout& grids = a_bsData[lev]->disjointBoxLayout();
	  mask[lev] = new LevelData<BaseFab<int> >(grids,1,IntVect::Zero);
	  for (DataIterator dit(grids); dit.ok(); ++dit)
	    {
	      (*mask[lev])[dit].setVal(1);
	      
	      if (lev < numLevels - 1)
		{
		  const DisjointBoxLayout& fgrids = a_bsData[lev+1]->disjointBoxLayout();
		  for (DataIterator fit(fgrids) ; fit.ok(); ++fit)
		    {
		      Box covered = fgrids[fit];
		      covered.coarsen(a_ratio[lev]);
		      covered &= grids[dit];
		      if (!covered.isEmpty())
			{
			  (*mask[lev])[dit].setVal(0,covered,0,1);
			}
		    }
		}
	    }
	  
	}
    }

  for (int lev = 0; lev < a_bsData.size(); lev++)
    {
      LevelData<FArrayBox>& levelData = *a_bsData[lev];
      for (DataIterator dit = levelData.dataIterator(); dit.ok(); ++dit)
	{
	  const Box& noGhostBox = levelData.disjointBoxLayout()[dit];
	  const Box& box = (a_validonly)?(noGhostBox):(levelData[dit].box());
	   
	  for (BoxIterator bit(box);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	     
	      if ( (!a_validonly) || ( (*mask[lev])[dit](iv) == 1))
		{
		  Vector<Real> v(levelData.nComp());
		  for (int ic = 0; ic < v.size(); ic++)
		    v[ic] = levelData[dit](iv,ic);
		  a_usData.append(lev, noGhostBox , iv,v);
		}
	    }
	}
    }
  
  //clean up mask
  for (int lev = 0; lev < mask.size(); lev++)
    {
      if (mask[lev] != NULL)
	{
	  delete mask[lev];mask[lev]=NULL;
	}
    }


}

/// read unstructureddata from a NetCDF-CF compliant file
/// along with the mesh data needed to reconstruct 
/// a Chombo AMR hierarchy  from it

void UnstructuredIO::readCF ( UnstructuredData& a_usData, 
		       Vector<std::string>& a_names, 
		       const std::string& a_file)
{
  
  CH_TIME("UnstructuredIO::readCF");
  int rc, ncID; 
  std::string ivname[SpaceDim] = {D_DECL("i","j","k")};
 

  //open file
  if ( (rc = nc_open(a_file.c_str(), NC_NOWRITE, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  // number of amr levels
  int levelDimID;
  if ( (rc = nc_inq_dimid(ncID, "level", &levelDimID) ) != NC_NOERR)
    {
      MayDay::Error("failed to determine level dimension id ");
    }
    
  size_t nLevel;
  if ( (rc = nc_inq_dimlen(ncID, levelDimID, &nLevel) ) != NC_NOERR)
    {
      MayDay::Error("failed to determine number of levels ");
    }

  // id of cell dimension (cell-centered data has a cell dimension only)
  int cellDimID;
  if ( (rc = nc_inq_dimid(ncID, "cell", &cellDimID) ) != NC_NOERR)
    {
      MayDay::Error("failed to determine cell dimension id ");
    }

  //read cell count;
  size_t nCell;
  if ( (rc = nc_inq_dimlen(ncID, cellDimID, &nCell) ) != NC_NOERR) 
    {
      MayDay::Error("failed to determine  number of cells");
    }

  // id of box dimension
  int boxDimID;
  if ( (rc = nc_inq_dimid(ncID, "box", &boxDimID) ) != NC_NOERR)
    {
      MayDay::Error("failed to determine box dimension id ");
    }

  //read box count;
  size_t nBox;
  if ( (rc = nc_inq_dimlen(ncID, boxDimID, &nBox) ) != NC_NOERR) 
    {
      MayDay::Error("failed to determine  number of boxes");
    }

  

  //refinement ratio
  Vector<int> ratio(nLevel); 
  readCFVar(ncID, "amr_ratio", ratio);

  //coarse level mesh spacing
  Vector<Real> dx(SpaceDim);
  readCFVar(ncID, "amr_crse_dx", dx);
  RealVect crseDx(dx);
   
  //coarse domain
  Box crseDomain;
  {
    Vector<int> corner(SpaceDim*SpaceDim);
    IntVect lo;
    IntVect hi;
    readCFVar(ncID, "amr_crse_domain", corner);
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	lo[dir] = corner[dir];
	hi[dir] = corner[dir+SpaceDim];
      }
    crseDomain.define(lo,hi);
  }
  

 
  // total number of variables in the netcdf file
  int nVar; 
  if ( ( rc = nc_inq_nvars    (ncID, &nVar) ) != NC_NOERR)
    {
      MayDay::Error("nc_inq_nvars failed");
    }
  
  //number and IDs of all netcdf variables
  Vector<int> ncVarID(nVar);

  if ( ( rc = nc_inq_varids    (ncID, &nVar, &ncVarID[0]) ) != NC_NOERR)
    {
      MayDay::Error("nc_inq_varids failed");
    }
  
  //number, name and IDs of cell-centered component fields
  Vector<int> ccVarID;
  
  for (int i =0; i < ncVarID.size(); i++)
    {
      int nDims;
      if ( ( rc = nc_inq_varndims  (ncID, ncVarID[i], &nDims ) ) != NC_NOERR)
	{
	  MayDay::Error("nc_inq_varndims failed");
	}
      //variables with one dim are candidates
      if (nDims == 1)
	{
	  int dimID;
	  if ( ( rc = nc_inq_vardimid  (ncID, ncVarID[i], &dimID ) ) != NC_NOERR)
	    {
	      MayDay::Error("nc_inq_vardimid failed");
	    }
	  
	  if (dimID == cellDimID)
	    {
	      //check that the amr_name attribute exists 
	      size_t l;
	      nc_type t;
	      if ( ( rc = nc_inq_att  (ncID, ncVarID[i], "amr_name", &t , &l  ) ) == NC_NOERR)
		{
		  if (l > 0 && t == NC_CHAR)
		    {
		      
		      char* c = new char[l+1];
		      if ( ( rc = nc_get_att_text  (ncID, ncVarID[i], "amr_name", c  ) ) != NC_NOERR)
			{
			  MayDay::Error("nc_get_att_text failed");
			}
		      c[l] = 0;
		      a_names.push_back(c);
		      delete c;
		      ccVarID.push_back(ncVarID[i]);
		    }
		}
	      }
	}
    }
  
  RealVect origin = RealVect::Zero;// we don't care about the space origin
  Real time = 0.0; //\todo : compute the time from the CF file 
  a_usData.define(ccVarID.size(), crseDx ,crseDomain, ratio, time, origin); 
  a_names.resize(ccVarID.size());
  a_usData.resize(nCell);
  
  //read in level and grid vector data
  {
    readCFVar(ncID, "amr_level", a_usData.level());
    
    std::string ivname[SpaceDim] = {D_DECL("i","j","k")};
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	std::string s = "amr_levelgrid_" + ivname[dir] + "_index"; 
	readCFVar(ncID,s,a_usData.iv(dir));	
      }
  }
  
  //field data
  //need to look up the scale factors
  Vector<FieldNames::CFRecord> cfRecord;
  FieldNames::CFLookup(cfRecord,  a_names);
  for (int i =0; i < ccVarID.size(); i++)
    {
      Vector<Real>& field = a_usData.field(i);
      readCFVar(ncID, ccVarID[i], field);
      if (cfRecord[i].scale() != 1.0)
	{
	  Real oneOnScale = 1.0 / cfRecord[i].scale();
	  for (int j = 0; j < field.size(); j++) 
	    field[j] *= oneOnScale;
	}
    }

  //box data
  {
    Vector<int> boxLevel(nBox);
    Vector<Vector<int> > lo(SpaceDim);
    Vector<Vector<int> > hi(SpaceDim);
    readCFVar(ncID, "amr_box_level", boxLevel);
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	std::string s = "amr_box_level_lo_" + ivname[dir] + "_index";
	lo[dir].resize(nBox);
	readCFVar(ncID, s, lo[dir]);
	s = "amr_box_level_hi_" + ivname[dir] + "_index"; 
	hi[dir].resize(nBox);
	readCFVar(ncID, s, hi[dir]);
      }

    UnstructuredData::LevelBoxSet& lbs = a_usData.levelBoxSet();
    for (int i = 0; i < boxLevel.size(); i++)
      {
	const int& lev = boxLevel[i];
	IntVect ivlo, ivhi;
	for (int dir = 0; dir < SpaceDim; dir++)
	  {
	    ivlo[dir] = lo[dir][i];
	    ivhi[dir] = hi[dir][i];  
	  }
	Box box(ivlo,ivhi);
	lbs.insert(std::make_pair( lev, box));
      }
  }

  


  if ( (rc = nc_close(ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to close netcdf file");
    }
}


/// write unstructured data to a NetCDF-CF compliant file, 
/// along with the mesh data needed to reconstruct 
/// a Chombo AMR hierarchy  from it
void UnstructuredIO::writeCF ( const std::string& a_file,  
			       const UnstructuredData& a_usData , 
			       const Vector<std::string>& a_names, 
			       const std::string& a_created,
			       const Transformation& a_latlonTransformation)
{

  CH_TIME("UnstructuredIO::writeCF");

  int rc; int ncID; int varID;

  std::string ivname[SpaceDim] = {D_DECL("i","j","k")};
  std::string xname[SpaceDim] = {D_DECL("x","y","z")};

  //create new file
  if ( (rc = nc_create(a_file.c_str(), NC_CLOBBER, &ncID) ) != NC_NOERR) 
    {
      MayDay::Error("failed to open netcdf file");
    }

  //define the netCDF dimensions 
  size_t nCell = a_usData.nCell();
  int cellDimID = defineCFDimension(ncID, nCell, "cell");

  size_t nVertex = 4; //all cells are quads
  int vertexDimID = defineCFDimension(ncID, nVertex, "vertex");

  size_t nLevel = a_usData.nLevel();
  int levelDimID =  defineCFDimension(ncID, nLevel, "level");

  int spaceDimID =  defineCFDimension(ncID, SpaceDim, "spacedim");

  size_t nBox = a_usData.levelBoxSet().size();
  int boxDimID =  defineCFDimension(ncID, nBox , "box");

  //set global attributes
  {
    std::string gs("unstructured");
    if ( (rc =  nc_put_att_text (ncID, NC_GLOBAL , "gridtype", 
				 gs.size(), gs.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add gridtype attribute");
      }

    std::string s("CF-1.6");
    if ( (rc =  nc_put_att_text (ncID, NC_GLOBAL , "Conventions", 
				 s.size(), s.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add Conventions attribute");
      }

    if ( (rc =  nc_put_att_text (ncID, NC_GLOBAL , "CreatedBy", 
				 a_created.size(), a_created.c_str())) != NC_NOERR)
      {
	MayDay::Error("failed to add CreatedBy attribute");
      }

  }

  


  //definition of data related to block structured mesh : 
  //(ratio, dx, coarse domain, level, box_level, box_lo_i, etc , IntVect components)
  {
    std::string s = "amr_ratio" ;
    defineCFVar(ncID, 1, &levelDimID, NC_INT,  s, "1","",s);
  }
  
  {
    std::string s = "amr_crse_dx" ;
    defineCFVar(ncID, 1, &spaceDimID, NC_DOUBLE, s, "1","",s); 
  }
  
  {
    std::string s = "amr_crse_domain" ;
    defineCFVar(ncID, 1, &vertexDimID, NC_INT, s, "1","",s); 
  }

  {
    std::string s = "amr_level" ;
    defineCFVar(ncID, 1, &cellDimID, NC_INT, s, "1","",s); 
  }
  
  {
    std::string s = "amr_box_level" ;
    defineCFVar(ncID, 1, &boxDimID, NC_INT, s, "1","",s); 
  }
  
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      std::string s = "amr_box_level_lo_" + ivname[dir] + "_index"; 
      defineCFVar(ncID, 1, &boxDimID, NC_INT, s, "1","",s);
      s = "amr_box_level_hi_" + ivname[dir] + "_index"; 
      defineCFVar(ncID, 1, &boxDimID, NC_INT, s, "1","",s);
    }
  
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      std::string s = "amr_levelgrid_" + ivname[dir] + "_index"; 
      defineCFVar(ncID, 1, &cellDimID, NC_INT, s, "1","",s);
    }
  
  //definition of cell center x,y,[z]
  
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      const std::string& name =  xname[dir];
      std::string stdname = "projection_"  + xname[dir] + "_coordinate";
      int varID = defineCFVar(ncID, 1, &cellDimID, NC_DOUBLE, name , "m", stdname ,"");
      
      std::string s = xname[dir] + "_vertices";
      if ( (rc =  nc_put_att_text (ncID, varID, "bounds", 
				   s.size(), s.c_str())) != NC_NOERR)
	{
	  MayDay::Error("failed to add bounds attribute to x");
	}
    }

  //definition of cell-centered lat and lon variables
  {
    int varID = defineCFVar(ncID, 1, &cellDimID, NC_DOUBLE, "lat" ,
			    "degrees_north", "" ,"latitude");

    std::string s = "lat_vertices";

    if ( (rc =  nc_put_att_text (ncID, varID, "bounds", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add bounds attribute to lat");
      }
    
  }
  
  {
   
    int varID = defineCFVar(ncID, 1, &cellDimID, NC_DOUBLE, "lon" ,
			    "degrees_east", "" ,"longitude");
    std::string s = "lon_vertices";

    if ( (rc =  nc_put_att_text (ncID, varID, "bounds", 
    				 s.size(), s.c_str())) != NC_NOERR)
      {
    	MayDay::Error("failed to add bounds attribute to lon");
      }
  }
  
  //definition of node-centered x/y values
  {
    int dimID[2] = {cellDimID, vertexDimID};
    std::string s("x_vertices");
    defineCFVar(ncID, 2, dimID, NC_DOUBLE, s , "m", "" ,s);
    s = "y_vertices";
    defineCFVar(ncID, 2, dimID, NC_DOUBLE, s , "m", "" ,s);
  }

  //definition of node-centered lat/lon values
  {
    int dimID[2] = {cellDimID, vertexDimID};
    std::string s("lat_vertices");
    defineCFVar(ncID, 2, dimID, NC_DOUBLE, s , "degrees_north", "" ,s);
    s = "lon_vertices";
    defineCFVar(ncID, 2, dimID, NC_DOUBLE, s , "degrees_east", "" ,s);
  }

  //convert the input names to CF names
  Vector<FieldNames::CFRecord> cfRecord;
  FieldNames::CFLookup(cfRecord,  a_names);
  

  //definition of field variables
  for (int ic = 0; ic < a_usData.nComp(); ic++)
    {
      int varID = defineCFVar(ncID, 1, &cellDimID, NC_DOUBLE, 
			      cfRecord[ic].name() , 
			      cfRecord[ic].units(), 
			      cfRecord[ic].standardName(), 
			      cfRecord[ic].longName());
      std::string s = "lat lon x y";
      
      if ( (rc =  nc_put_att_text (ncID, varID, "coordinates", 
				   s.size(), s.c_str())) != NC_NOERR)
	{
	  MayDay::Error("failed to add coordinates attribute");
	}
      
      s = a_names[ic];
      if (s.size() > 0)
	{
	  if ( (rc =  nc_put_att_text (ncID, varID, "amr_name", 
				       s.size(), s.c_str())) != NC_NOERR)
	    {
	      MayDay::Error("failed to add amr_name attribute");
	    }
	}
    }
  

  if ( (rc = nc_enddef(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to define netcdf file");
    }

 
  //write data

  //amr data / ratio
  writeCFVar(ncID, "amr_ratio", a_usData.ratio());

  //amr data / coarse dx
  Vector<Real> crseDx(SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      crseDx[dir] = a_usData.dx()[0][dir];
    }
  writeCFVar(ncID, "amr_crse_dx", crseDx);


  //coarse domain
  Vector<int> corner(SpaceDim*SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      corner[dir] = a_usData.domain(0).smallEnd()[dir];
      corner[dir+SpaceDim] = a_usData.domain(0).bigEnd()[dir];
    }
  writeCFVar(ncID, "amr_crse_domain", corner);

  // amr data / cell level
  writeCFVar(ncID, "amr_level", a_usData.level());
    
  // amr data / box level
  {
    Vector<int> boxLevel;
    Vector<Vector<int> > lo(SpaceDim);
    Vector<Vector<int> > hi(SpaceDim);
    const UnstructuredData::LevelBoxSet& lbs = a_usData.levelBoxSet();
    for (UnstructuredData::LevelBoxSet::const_iterator it 
	   = lbs.begin(); it != lbs.end(); it++)
      {
	boxLevel.push_back(it->first);
	for (int dir = 0; dir < SpaceDim; dir++)
	  {
	    lo[dir].push_back(it->second.smallEnd()[dir]);
	    hi[dir].push_back(it->second.bigEnd()[dir]);
	  }
      }
    writeCFVar(ncID, "amr_box_level", boxLevel);
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	std::string s = "amr_box_level_lo_" + ivname[dir] + "_index";
	writeCFVar(ncID, s, lo[dir]);
	s = "amr_box_level_hi_" + ivname[dir] + "_index"; 
	writeCFVar(ncID, s, hi[dir]);
      }
  }

  // amr data / grid vectors 
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      std::string s = "amr_levelgrid_" + ivname[dir] + "_index"; 
      writeCFVar(ncID, s , a_usData.iv(dir));
    }


  //x,y
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      writeCFVar(ncID, xname[dir] , a_usData.x(dir));
    }

  {
    // cell-center latitude and longitude data
    D_TERM(const Vector<Real>& x = a_usData.x(0);,
	   const Vector<Real>& y = a_usData.x(1);,
	   const Vector<Real>& z = a_usData.x(2););
    
    Vector<Real> lat(nCell);
    Vector<Real> lon(nCell);
    
    for (int i = 0; i < x.size(); i++)
      {
	RealVect X(D_DECL(x[i],y[i],z[i]));
	RealVect L = a_latlonTransformation.transform(X);
	lat[i] = L[0] / M_PI * 180.0;
	lon[i] = L[1] / M_PI * 180.0;
      }
    
    writeCFVar(ncID, "lat" , lat);
    writeCFVar(ncID, "lon" , lon);

  }

  {
    // cell node x,y, latitude and longitude data
    Vector< Vector<Real> > nodeX(SpaceDim) ;
    Vector<Real> nodeLat(nCell);
    Vector<Real> nodeLon(nCell);

    Box unit(IntVect::Zero,IntVect::Unit);
    //for rectangular cells, iterating across every node
    //is a bit wasteful as e.g x(0,1) = x(0,0), but maybe
    //one day we might have mapped co-oords. Presumably the
    //I/O will dominate anyway
    size_t start[2] = {0,0};
    size_t count[2] = {nCell,1};

    Vector<IntVect> ivs;
    ivs.push_back(IntVect(0,0));
    ivs.push_back(IntVect(1,0));
    ivs.push_back(IntVect(1,1));
    ivs.push_back(IntVect(0,1));

    for (int j=0; j < ivs.size(); j++)
      {
	const IntVect& iv = ivs[j];
	
	for (int dir = 0; dir < SpaceDim; dir++)
	  {
	    a_usData.computeNodeCoord( nodeX[dir], dir, iv);
	  }

	for (int i = 0; i < nCell; i++)
	  {
	    RealVect X(D_DECL(nodeX[0][i],nodeX[1][i],nodeX[2][i]));
	    RealVect L = a_latlonTransformation.transform(X);
	    nodeLat[i] = L[0] / M_PI * 180.0;
	    nodeLon[i] = L[1] / M_PI * 180.0;
	    


	  }

	if ( (rc = nc_inq_varid(ncID,"x_vertices", &varID)) != NC_NOERR)
	  {
	    MayDay::Error("failed to find variable id");
	  }
	
	if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeX[0][0]))) != NC_NOERR)
	  {
	    MayDay::Error("failed to write x_vertices data");
	  }
	
	if ( (rc = nc_inq_varid(ncID,"y_vertices", &varID)) != NC_NOERR)
	  {
	    MayDay::Error("failed to find variable id");
	  }
	
	if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeX[1][0]))) != NC_NOERR)
	  {
	    MayDay::Error("failed to write y_vertices data");
	  }
	


	if ( (rc = nc_inq_varid(ncID,"lat_vertices", &varID)) != NC_NOERR)
	  {
	    MayDay::Error("failed to find variable id");
	  }
	
	if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeLat[0]))) != NC_NOERR)
	  {
	    MayDay::Error("failed to write lat_vertices data");
	  }
	
	
	if ( (rc = nc_inq_varid(ncID,"lon_vertices", &varID)) != NC_NOERR)
	      {
		MayDay::Error("failed to find variable id");
	      }
	    
	    if ( (rc = nc_put_vara_double(ncID, varID, start, count, &(nodeLon[0]))) != NC_NOERR)
	      {
		MayDay::Error("failed to write lon_vertices data");
	      }
	  

	start[1] += 1;
      }
  
  }

  //scalar fields
  for (int ic = 0; ic < a_usData.nComp(); ic++)
    {
      Vector<Real> scaledField = a_usData.field(ic);
      for (int i =0; i < scaledField.size(); i++)
	{
	  scaledField[i] *= cfRecord[ic].scale();
	}

      writeCFVar(ncID,cfRecord[ic].name(),scaledField );
    }
      
  //close file
  if ( (rc = nc_close(ncID) ) != NC_NOERR)
    {
      MayDay::Error("failed to close netcdf file");
    }

}

int UnstructuredIO::defineCFDimension
(int a_ncID, 
 size_t a_len, 
 const std::string& a_name)
{

  int rc,dimID;
  if ( (rc = nc_def_dim(a_ncID, a_name.c_str() , a_len, &dimID)) != NC_NOERR)
    {
      std::string msg = "failed to define " + a_name + " dimension";
      MayDay::Error(msg.c_str());
    }
  return dimID;
}


int UnstructuredIO::defineCFVar
(int a_ncID,  int a_nDim, int* a_dimID, nc_type a_type,
 const std::string& a_name,
 const std::string& a_unit,
 const std::string& a_stdName, 
 const std::string& a_longName)
{

  int rc, varID;

  if ( (rc = nc_def_var(a_ncID, a_name.c_str(), a_type,
			a_nDim, a_dimID, &varID)) != NC_NOERR)
    {
      std::string msg = "failed to define " + a_name;
      MayDay::Error(msg.c_str());
    }
  
  if (a_stdName != "")
    {
      if ( (rc =  nc_put_att_text (a_ncID, varID, "standard_name", 
			       a_stdName.size(), a_stdName.c_str())) != NC_NOERR)
	{
	  std::string msg = "failed to add " + a_name + "standard_name attribute" ;
	  MayDay::Error(msg.c_str());
	}
    }
    
  if (a_longName != "")
    {
      if ( (rc =  nc_put_att_text (a_ncID, varID, "long_name", 
				   a_longName.size(), a_longName.c_str())) != NC_NOERR)
	{
	  std::string msg = "failed to add " + a_name + "long_name attribute" ;
	  MayDay::Error(msg.c_str());
	}
    }
      
 if ( (rc =  nc_put_att_text (a_ncID, varID, "units", 
			      a_unit.size(), a_unit.c_str())) != NC_NOERR)
   {
     std::string msg = "failed to add " + a_name + "unitsattribute" ;
     MayDay::Error(msg.c_str());
   }  


 return varID;

}

void UnstructuredIO::readCFVar
(int a_ncID, const std::string& a_name, Vector<int>& a_data)
{
  int varID, rc;
  if ( (rc = nc_inq_varid(a_ncID, a_name.c_str(), &varID)) != NC_NOERR)
    {
      std::string msg = "failed to find variable " + a_name;
      MayDay::Error(msg.c_str());
    }
  readCFVar(a_ncID,varID,a_data);
}

void UnstructuredIO::readCFVar
(int a_ncID, int  a_varID , Vector<int>& a_data)
{
  int rc;
  size_t start = 0;
  size_t count = a_data.size();
  if ( (rc = nc_get_vara_int(a_ncID, a_varID, &start, &count, &a_data[0])) != NC_NOERR)
    {
      std::string msg = "failed to read from var";
      MayDay::Error(msg.c_str() );
    }   
}


void UnstructuredIO::readCFVar
(int a_ncID, const std::string& a_name, Vector<Real>& a_data)
{

  int varID, rc;
  if ( (rc = nc_inq_varid(a_ncID, a_name.c_str(), &varID)) != NC_NOERR)
    {
      std::string msg = "failed to find variable " + a_name;
      MayDay::Error(msg.c_str());
    }
  readCFVar(a_ncID,varID,a_data);
  
}

void UnstructuredIO::readCFVar
(int a_ncID, int  a_varID , Vector<Real>& a_data)
{
  int rc;
  size_t start = 0;
  size_t count = a_data.size();
  if ( (rc = nc_get_vara_double(a_ncID, a_varID, &start, &count, &a_data[0])) != NC_NOERR)
    {
      std::string msg = "failed to read from var";
      MayDay::Error(msg.c_str() );
    }   
}

void UnstructuredIO::writeCFVar
(int a_ncID, const std::string& a_name, const Vector<int>& a_data)
{

  int varID, rc;
  if ( (rc = nc_inq_varid(a_ncID, a_name.c_str(), &varID)) != NC_NOERR)
    {
      std::string msg = "failed to find variable " + a_name;
      MayDay::Error(msg.c_str());
    }
  
  size_t start = 0;
  size_t count = a_data.size();
  if ( (rc = nc_put_vara_int(a_ncID, varID, &start, &count, &a_data[0])) != NC_NOERR)
    {
      std::string msg = "failed to write to" +  a_name;
      MayDay::Error(msg.c_str() );
    }   
}

void UnstructuredIO::writeCFVar
(int a_ncID, const std::string& a_name, const Vector<Real>& a_data)
{

  int varID, rc;
  if ( (rc = nc_inq_varid(a_ncID, a_name.c_str(), &varID)) != NC_NOERR)
    {
      std::string msg = "failed to find variable " + a_name;
      MayDay::Error(msg.c_str());
    }
  
  size_t start = 0;
  size_t count = a_data.size();
  if ( (rc = nc_put_vara_double(a_ncID, varID, &start, &count, &a_data[0])) != NC_NOERR)
    {
      std::string msg = "failed to write to" +  a_name;
      MayDay::Error(msg.c_str() );
    }   
}




#include "NamespaceFooter.H"
