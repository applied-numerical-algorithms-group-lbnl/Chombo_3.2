 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

/*===========================================================================
 libamrfile.cpp
 C-interface functions to read Chombo's amr files, read and modify their data, 
 and write them 
==========================================================================*/
#include "libamrfile.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include <fstream>
#include "AMRIO.H"


class AMRHierarchy
{
  Vector<LevelData<FArrayBox>* > m_data;
  Vector<DisjointBoxLayout> m_grids;
  Vector<std::string> m_names;
  Vector<int> m_ratio;
  Vector<Real> m_dx;
  Real m_dt;
  Real m_time;
  int m_nLevel;
  bool m_ok;
  Box m_crseBox;
  std::map<std::string,int> m_nameCompMap;
  Vector<Real> m_sigma;
  
  AMRHierarchy();

  void updateNameCompMap()
  {
    m_nameCompMap.clear();
    for (int i = 0; i < m_names.size(); i++)
      {
	m_nameCompMap[m_names[i]] = i;
      }
  }

public:

  AMRHierarchy(const std::string& file)
  {
  
    Real crseDx;

    //check the file access first, because ReadAMRHierarchyHDF5 will 
    //throw an exception otherwise
    int status = LIBAMRFILE_ERR_READ_FAILED ;
    std::ifstream ifs(file.c_str());
    bool access = ifs.good();
    ifs.close();
   
    if (access)
      {
	 status = ReadAMRHierarchyHDF5(file, m_grids, m_data, m_names, 
					  m_crseBox, crseDx, m_dt,
					  m_time, m_ratio, m_nLevel);
      
	 m_dx.resize(m_nLevel,crseDx);
	 for (int lev = 1; lev < m_nLevel; lev++)
	   {
	     m_dx[lev] = m_dx[lev-1]/Real(m_ratio[lev-1]);
	     
	   } 
	 for (int lev = 0; lev < m_nLevel; lev++)
	   {
	     m_data[lev]->exchange();
	   }

	 for (int lev = m_nLevel - 1; lev > 0; lev--)
	   {
	     //coarse average to lower levels
	     CoarseAverage av(m_grids[lev],m_names.size(), m_ratio[lev-1]);
	     av.averageToCoarse(*m_data[lev-1],*m_data[lev]);
	   }


	 updateNameCompMap();
      }
    m_ok = status == 0;
  }

  //construct single level amr hierarchy
  AMRHierarchy(int nx, int ny, double dx, int n_comp, int n_ghost)
  {
    //one box for now
    IntVect lim; lim[0] = nx - 1 ; lim[1] = ny -1 ;
    
    m_crseBox = Box(IntVect::Zero,lim);
    DisjointBoxLayout dbl(Vector<Box>(1,m_crseBox),Vector<int>(1,0));
    m_grids.push_back(dbl);
    m_dx.push_back(dx);
    m_data.push_back(new LevelData<FArrayBox>(m_grids[0], n_comp, IntVect::Unit * n_ghost));
    m_time = 0.0;
    m_dt = 1.0;
    m_nLevel = 1;
    m_ratio.push_back(2);
    m_names.resize(n_comp);
    for (int i =0; i < n_comp; i++)
      {
	char name[16];
	sprintf(name,"component%i",i);
	m_names[i] = name; 
      }
    updateNameCompMap();
    m_ok = true;
  }
  


  bool ok() const {return m_ok;}
  const Vector<LevelData<FArrayBox>* >& data() const {return m_data;}
  const Vector<DisjointBoxLayout>& grids() const {return m_grids;}
  const Vector<Real>& dx() const  {return m_dx;}
  const Vector<std::string>& names() const  {return m_names;}
  const Vector<int>& ratio() const {return m_ratio;}
  const Vector<Real>& sigma() const {return m_sigma;}
  const std::map<std::string,int>& nameCompMap() const {return m_nameCompMap;}
  int nLevel() const  {return m_nLevel;}
  Real time() const {return m_time;}

  void setName(int comp, const std::string& name)
  {
    if (comp >= 0 && m_names.size() > comp )
      {
	m_names[comp] = name;
	updateNameCompMap();
      }
  }

  int write(const std::string& file)
  {
    int status = 0;
    if (m_ok)
      {
	WriteAMRHierarchyHDF5(file, m_grids, m_data, 
			      m_names, m_crseBox, m_dx[0], m_dt,
			      m_time, m_ratio, m_nLevel);
      }
    return status;
  }

  ~AMRHierarchy()
  {
    for (int lev = 0; lev < m_data.size(); lev++)
      {
	if (m_data[lev] != NULL)
	  {
	    delete m_data[lev];
	    m_data[lev] = NULL;
	  }
      }
  }

};

namespace libamrfile
{
  std::map<int, AMRHierarchy*> g_store; 
}

void amr_read_file(int *status, int *amr_id, const char *file)
{

  if (!status)
    return;

  if (libamrfile::g_store.size() > LIBAMRFILE_MAX_AMR_HIERARCHIES )
    {
      *status = LIBAMRFILE_ERR_TOO_MANY_AMR_HIERARCHIES ;
      return;
    }
	    
  AMRHierarchy* h = new AMRHierarchy(file);
 
  if (h->ok())
    {
      if (libamrfile::g_store.size() == 0)
	*amr_id = 0;
      else
	*amr_id  = (--libamrfile::g_store.end())->first + 1;
      libamrfile::g_store[*amr_id] = h;
      *status = 0;
    }
  else
    {
      *status = LIBAMRFILE_ERR_READ_FAILED;
    }
}
void amr_read_file_R(int *status, int *amr_id, char **file)
{
  amr_read_file(status,amr_id,*file);
}

void amr_create_coarse_2d(int *status, int *amr_id, 
		       const int *nx, const int *ny,  
		       const double* dx, const int* n_comp, const int* n_ghost)
{

  if (!status)
    return;

  if (!(amr_id && D_TERM(nx, && ny, && nz) && n_comp && n_ghost))
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }

  if (libamrfile::g_store.size() > LIBAMRFILE_MAX_AMR_HIERARCHIES )
    {
      *status = LIBAMRFILE_ERR_TOO_MANY_AMR_HIERARCHIES ;
      return;
    }

  AMRHierarchy* h = new AMRHierarchy(*nx, *ny, *dx, *n_comp, *n_ghost);
 
  if (h->ok())
    {
      if (libamrfile::g_store.size() == 0)
	*amr_id = 0;
      else
	*amr_id  = (--libamrfile::g_store.end())->first;
      libamrfile::g_store[*amr_id] = h;
      *status = 0;
    }
  else
    {
      *status = LIBAMRFILE_ERR_CREATE_FAILED;
    }
  

}


void amr_write_file(int *status, int *amr_id, const char *file)
{

  if (!status)
    return;

  if (!amr_id)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }

  std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
  if (i == libamrfile::g_store.end())
    {
      *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
      return;
    }

  AMRHierarchy* h = i->second;
  if (!h || !(h->ok()))
    {
      *status = LIBAMRFILE_ERR_BAD_AMR_HIERARCHY;
      return;
    }

  *status = h->write(file);

}


void amr_write_file_R(int *status, int *amr_id, char **file)
{
  amr_write_file(status,amr_id,*file);
}


void amr_free_all()
{
  for (std::map<int, AMRHierarchy*>::iterator i = libamrfile::g_store.begin();
       i != libamrfile::g_store.end(); ++i)
    {
      if (i->second != NULL)
	{
	  delete i->second;
	  i->second = NULL;
	  libamrfile::g_store.erase(i->first);
	}
    }
  
}

void amr_free(int *status, int *amr_id)
{

  if (!status)
    return;

  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  if (i->second != NULL)
	    {
	      delete i->second;
	      i->second = NULL;
	      libamrfile::g_store.erase(i->first);
	    }
	  *status = 0;
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
    }

}

void amr_query_comp_name(int *status, char *name, const int* amr_id, const int* comp, const int* namelen)
{

  if (!status)
    return;

  if (amr_id && comp && namelen && name)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  if (i->second)
	    {
	      const AMRHierarchy& h = *i->second;
	      if (*comp < h.names().size())
		{
		  int c = *comp;
		  int l = *namelen;
		  strncpy( name, h.names()[c].c_str(), l);
		  *status = 0;
		}
	      else
		{
		  *status = LIBAMRFILE_ERR_NO_SUCH_COMP;
		}
	    }
	  else
	    {
	      *status = LIBAMRFILE_ERR_NULL_POINTER;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }

}

void amr_query_comp_id(int *status, int *comp, 
		       const int* amr_id, const char* name , const int* namelen)
{

  if (!status)
    return;

  if (amr_id && comp && namelen && name && namelen)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  if (i->second)
	    {
	      const AMRHierarchy& h = *i->second;
	      std::map<std::string,int>::const_iterator j = h.nameCompMap().find(std::string(name, *namelen));
	      if (j == h.nameCompMap().end())
		{
		  *comp = -1;
		  *status = LIBAMRFILE_ERR_NO_SUCH_COMP;
		}
	      else
		{
		  *comp = j->second;
		  *status = 0;
		}
	      
	    }
	  else
	    {
	      *status = LIBAMRFILE_ERR_NULL_POINTER;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }

}


void amr_query_comp_name_R(int *status, char **name, const int* amr_id, const int* comp, const int* namelen)
{
  amr_query_comp_name(status, *name, amr_id, comp, namelen);
}

void amr_query_comp_id_R(int *status, int *comp, const int* amr_id, const char** name, const int* namelen)
{
  amr_query_comp_id(status, comp, amr_id, *name, namelen);
}

void amr_set_comp_name(int *status, const char *name, const int* amr_id, const int* comp)
{
  
  if (!status)
    return;

  if (!(amr_id && comp ))
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }


  std::map<int, AMRHierarchy*>::iterator i = libamrfile::g_store.find(*amr_id);
  if (i != libamrfile::g_store.end())
    {
      if (i->second)
	{
	  AMRHierarchy& h = *i->second;
	  if (*comp < h.names().size())
	    {
	      std::string s = name;
	      h.setName(*comp,s);
	    }
	    else
	      {
		*status = LIBAMRFILE_ERR_NO_SUCH_COMP;
	      }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NULL_POINTER;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
    }
}



void amr_set_comp_name_R(int *status, const char **name, const int* amr_id, const int* comp)
{
  if (!name)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }
  
  amr_set_comp_name(status, *name, amr_id, comp);
}


void amr_query_n_level(int *status, int *n_level, const int *amr_id)
{
  
  if (!status)
    return;

  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  *n_level = i->second->nLevel();
	  *status = 0;
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}
 

void amr_query_time(int *status, double *time, const int *amr_id)
{
  
  if (!status)
    return;

  if (amr_id && time)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  *time = i->second->time();
	  *status = 0;
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}


void amr_query_domain_corners(int *status, int *lo, int* hi, const int *amr_id, const int *level)
{

  if (!status)
    return;

  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  if (lo && hi && level)
	    {
	      if ( (*level < 0) || (*level >= i->second->nLevel()))
		{
		  *status = LIBAMRFILE_ERR_NO_SUCH_LEVEL ;
		}
	      else
		{
		  const Box& dombox = i->second->grids()[*level].physDomain().domainBox();
		  for (int dir = 0; dir < SpaceDim; dir++)
		    {
		      lo[dir] = dombox.smallEnd()[dir];
		      hi[dir] = dombox.bigEnd()[dir];
		    }
		  *status = 0;
		}
	    }
	  else
	    {
	      *status = LIBAMRFILE_ERR_NULL_POINTER ;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}


//
void bisicles_read_sigma(int *status, double *sigma, const int* n_sigma, const int* amr_id)
{
  if (!status)
    return;

  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  AMRHierarchy* h = i->second;
	  
	    
	}
    }
}


void amr_query_n_comp(int *status, int* n_comp, const int* amr_id)
{

  if (!status)
    return;

  if (amr_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  AMRHierarchy* h = i->second;
	  if (h)
	    {
	      *n_comp = h->data()[0]->nComp();
	      *status = 0;
	    }
	  else
	    {
	      *status = LIBAMRFILE_ERR_NULL_POINTER;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}


void amr_query_n_fab(int *status, int *n_fab, const int *amr_id, const int *level_id)
{

  if (amr_id && level_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  int nLevel =  i->second->nLevel();
	  if (*level_id < nLevel)
	    {
	      *n_fab = i->second->grids()[*level_id].size();
	      *status = 0;
	    }
	  else
	    {
	      *status =  LIBAMRFILE_ERR_NO_SUCH_LEVEL;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }

}

void amr_query_fab_dimensions_2d(int *status, 
			      int *nx, int *ny,
			      int *ncomp, 
			      const int *amr_id, 
			      const int *level_id, 
			      const int* fab_id)
{
  
  if (!status)
    return;

  if (!nx)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
      return;
    }

  if (!ny)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
      return;
    }


  if (ncomp && amr_id && level_id && fab_id)
    {
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  int nLevel =  i->second->nLevel();

	  if (*level_id < nLevel)
	    {
	      const LevelData<FArrayBox>& ldf = *i->second->data()[*level_id];
	      DataIterator dit(i->second->grids()[*level_id]);
	      for (int j=0; j < *fab_id; j++)
		{
		  ++dit;
		}
		if (dit.ok())
		  {
		    const Box& b = i->second->grids()[*level_id][dit];
		    *nx = b.bigEnd()[0] - b.smallEnd()[0] + 1;
		    *ny = b.bigEnd()[1] - b.smallEnd()[1] + 1;
		    *ncomp = ldf.nComp();
		    *status = 0;
		  }
		else
		  {
		    *status = LIBAMRFILE_ERR_NO_SUCH_FAB;
		  }
	    }
	  else
	    {
	      *status =  LIBAMRFILE_ERR_NO_SUCH_LEVEL;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
    }
}

void amr_read_fab_data_2d(int *status, 
		       double *fab_data, 
		       double *x_data, 
		       double *y_data,
		       const int *amr_id, 
		       const int *level_id, 
		       const int* fab_id,
		       const int* comp_id,
		       const int* nghost)
{

  
  if (!status)
    return;

  if (!x_data)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER  ;
      return;
    }

  if (!y_data)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }

  
  if (fab_data  && amr_id && level_id && fab_id && comp_id && nghost)
    {

     

      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  int nLevel =  i->second->nLevel();

	  if (*level_id < nLevel)
	    {
	      const LevelData<FArrayBox>& ldf = *i->second->data()[*level_id];
	      int maxnghost=std::min(ldf.ghostVect()[0],ldf.ghostVect()[1]);

		if ( *nghost < 0 || ( *nghost >  maxnghost ))
		{
		  *status = LIBAMRFILE_ERR_BAD_NGHOST;
		  return;
		}
		  
	      if (*comp_id < 0 || *comp_id >= ldf.nComp())
		{
		  *status = LIBAMRFILE_ERR_NO_SUCH_COMP;
		  return;
		}


	      Real dx = i->second->dx()[*level_id];
	      DataIterator dit(i->second->grids()[*level_id]);
	      for (int j=0; j < *fab_id; j++)
		{
		  ++dit;
		}
		if (dit.ok())
		  {
		    Box b = i->second->grids()[*level_id][dit];
		    b.grow(*nghost);


		    const FArrayBox& fab = ldf[dit];
		    fab.linearOut((void*)fab_data,b,Interval(*comp_id,*comp_id));

		    double *xptr = x_data;
		    for (int ix = b.smallEnd()[0] ; ix <= b.bigEnd()[0]; ix++)
		      {
			*xptr++ = (double(ix)+0.5)*dx;
		      }

		    double *yptr = y_data;
		    for (int ix = b.smallEnd()[1] ; ix <= b.bigEnd()[1]; ix++)
		      {
			*yptr++ = (double(ix)+0.5)*dx;
		      }

		    *status = 0;
		  }
		else
		  {
		    *status = LIBAMRFILE_ERR_NO_SUCH_FAB;
		  }
	    }
	  else
	    {
	      *status =   LIBAMRFILE_ERR_NO_SUCH_LEVEL;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}

/**fill a_destData with piecewise constant data from a_srcData, 
   when a_srcData is lower resolution, and the refinement ratio is a_nRef. 
 */
int upFill(LevelData<FArrayBox>& a_destData, const LevelData<FArrayBox>&  a_srcData, int a_nRef, int a_order)
{
  int err = 0;
  const DisjointBoxLayout& srcGrids = a_srcData.disjointBoxLayout();
  const DisjointBoxLayout& destGrids = a_destData.disjointBoxLayout();
  bool coarsenable = destGrids.coarsenable(a_nRef);
  if ( coarsenable) 
    {
 
      const ProblemDomain& fineDomain = destGrids.physDomain();
      FineInterp interpolator(destGrids, a_destData.nComp(), a_nRef, fineDomain);
      if (a_order == 0)
	{
	  interpolator.pwcinterpToFine(a_destData, a_srcData, true);
	}
      else if (a_order == 1)
	{
	  interpolator.interpToFine(a_destData, a_srcData, true);
	}
    }
  else if (a_nRef%2 == 0)
    {
      //interpolate coarse data recursively, refining the coarser level by a factor of 2 on each recursion.
      DisjointBoxLayout stepGrids;
      refine(stepGrids,srcGrids,2);
      LevelData<FArrayBox> stepData(stepGrids,a_srcData.nComp(), a_srcData.ghostVect());
      int err = upFill(stepData,  a_srcData,  2, a_order );
      if (err == 0)
	err = upFill(a_destData, stepData, a_nRef/2, a_order);
    }
  else
    {
      err =  LIBAMRFILE_ERR_BAD_REFINEMENT_RATIO;
    }
  return err;
}


void amr_read_box_data_2d(int *status, 
			      double *comp_data, 
			      double *x_data, 
			      double *y_data,
			      const int *amr_id, 
			      const int *level_id, 
			      const int *lo,
			      const int *hi,
			      const int* comp_id,
			      const int* interp_order)
{

  if (!status)
    return;

  if (!x_data)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER  ;
      return;
    }

  if (!y_data)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER;
      return;
    }

  if (!interp_order)
    {
      *status =  LIBAMRFILE_ERR_BAD_INTERPOLATION_ORDER;
      return;
    }
  //interpolation order must be 0 or 1
  if (*interp_order < 0)
    {
      *status =  LIBAMRFILE_ERR_BAD_INTERPOLATION_ORDER;
      return;
    }
  if (*interp_order > 1)
    {
      *status =  LIBAMRFILE_ERR_BAD_INTERPOLATION_ORDER;
      return;
    }
  if (!comp_id)
    {
      *status =  LIBAMRFILE_ERR_NULL_POINTER  ;
      return;
    }


  if (comp_data && x_data && y_data && amr_id && level_id && lo && hi && comp_id )
    {

      
      std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
      if (i != libamrfile::g_store.end())
	{
	  int nLevel =  i->second->nLevel();
	  
	  if (*level_id < nLevel)
	    {
	      const LevelData<FArrayBox>& ldf = *i->second->data()[*level_id];
	      
	      if (*comp_id < 0 || *comp_id >= ldf.nComp())
		{
		  *status = LIBAMRFILE_ERR_NO_SUCH_COMP;
		  return;
		}

	      Real dx = i->second->dx()[*level_id];
	      
	      for (int dir = 0; dir < SpaceDim; dir++)
		{
		  if (hi[dir] < lo[dir])
		    {
		      *status = LIBAMRFILE_ERR_BAD_BOX;
		      return;
		    }
		}
	      Box box(IntVect(lo[0],lo[1]),IntVect(hi[0],hi[1]));
	      if (!box.intersects(ldf.disjointBoxLayout().physDomain().domainBox()))
		{
		  *status = LIBAMRFILE_ERR_BAD_BOX;
		  return;
		}

	      double *xptr = x_data;
	      for (int ix = box.smallEnd()[0] ; ix <= box.bigEnd()[0]; ix++)
		{
		  *xptr++ = (double(ix)+0.5)*dx;
		}

	      double *yptr = y_data;
	      for (int ix = box.smallEnd()[1] ; ix <= box.bigEnd()[1]; ix++)
		{
		  *yptr++ = (double(ix)+0.5)*dx;
		}
	      
	      LevelData<FArrayBox> destData(DisjointBoxLayout(Vector<Box>(1,box),Vector<int>(1,0),
							      ldf.disjointBoxLayout().physDomain())
					    ,1,IntVect::Zero);
	      *status = 0;
	      
	      for (int lev = 0; lev < *level_id; lev++)
		{
		  //coarse-to fine interpolation 
		  

		  //work out the refinement ration from level lev to level *level_id
		  int nRef = 1;
		  for (int l = lev; l < *level_id; l++)
		    {
		      nRef *= i->second->ratio()[l];
		    }

		  LevelData<FArrayBox> alias; 
		  aliasLevelData(alias, i->second->data()[lev] ,Interval(*comp_id,*comp_id));
		  *status = upFill(destData, alias  , nRef , *interp_order);
		  //const DisjointBoxLayout& fineGrids = destData.getBoxes();
		  //FineInterp fi(fineGrids,1,nRef,fineGrids.physDomain());
		  //fi.interpToFine(destData,alias, true);
		}
	      if (*status == 0)
		{
		  //same level copy
		  ldf.copyTo(Interval(*comp_id,*comp_id), destData, Interval(0,0));
		  
		  //fine-to-coarse averaging not needed? Assumes it has been done to ldf already
		  
		  //copy to argument buffer (hopefully allocated correctly)
		  DataIterator dit(destData.disjointBoxLayout());
		  dit.reset();
		  destData[dit].linearOut((void*)comp_data,box,Interval(0,0));
		}	 
	    }
	  else
	    {
	      *status =   LIBAMRFILE_ERR_NO_SUCH_LEVEL;
	    }
	}
      else
	{
	  *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
	}
    }
  else
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER;
    }
}



void amr_write_fab_data_2d(int *status, 
			   double *fab_data, 
			   int *nx, int *ny,
			   const int *amr_id, 
			   const int *level_id, 
			   const int* fab_id,
			   const int* comp_id,
			   const int* nghost)
{
  if (!status)
    return;
  
   if (!nx)
    {
      *status = LIBAMRFILE_ERR_NULL_POINTER ;
      return;
    }

   if (!ny)
     {
       *status = LIBAMRFILE_ERR_NULL_POINTER ;
       return;
     }



   if ( !(fab_data && amr_id && level_id && fab_id && comp_id && nghost))
     {
       *status = LIBAMRFILE_ERR_NULL_POINTER ;
       return;
     }
     
   std::map<int, AMRHierarchy*>::const_iterator i = libamrfile::g_store.find(*amr_id);
   if (i == libamrfile::g_store.end())
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_AMR_ID;
       return;
     }

   AMRHierarchy* h = i->second;
   if (!h)
     {
       *status = LIBAMRFILE_ERR_NULL_POINTER;
       return;
     }

   int nLevel =  i->second->nLevel();
   if (*level_id >= nLevel)
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_LEVEL;
       return;
     }
   
   LevelData<FArrayBox>& ldf = *i->second->data()[*level_id];
   int maxnghost=std::min(ldf.ghostVect()[0],ldf.ghostVect()[1]);
   
   if (*nghost < 0 || ( *nghost >  maxnghost ))
     {
       *status = LIBAMRFILE_ERR_BAD_NGHOST;
       return;
     }
   
   if (*comp_id < 0 || *comp_id >= ldf.nComp())
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_COMP;
       return;
     }

   DataIterator dit(i->second->grids()[*level_id]);

   for (int j=0; j < *fab_id; j++)
     {
       ++dit;
     }
   if (! dit.ok())
     {
       *status = LIBAMRFILE_ERR_NO_SUCH_FAB;
       return;
     }

   Box b = i->second->grids()[*level_id][dit];
   b.grow(*nghost);
   FArrayBox& fab = ldf[dit];
   fab.linearIn((void*)fab_data,b,Interval(*comp_id,*comp_id));

   status = 0;
		    
}


