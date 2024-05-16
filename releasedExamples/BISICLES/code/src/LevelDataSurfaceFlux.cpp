#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// concrete class encapsulating surface fluxes determined  
// by copying, coarsening, or interpolating a LevelData<FArrayBox>
// covering an entire domain. Objects can be defined either
// by specifying a RefCountedPtr<LevelData<FArrayBox> > , or by specifying
// a std::map<Real,string> mapping time ti to an hdf5 file f. In the
// latter case, the flux at time t is found by linear interploation
// between time max(ti <= t) and min(ti > t) 

// \todo replace the std::map<Real,string> mechanism with
// a suitable abstraction

#include "LevelDataSurfaceFlux.H"
#include "ReadLevelData.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "AmrIceBase.H"
#include "NamespaceHeader.H"

void LevelDataSurfaceFlux::surfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
    
  const Real& time = a_amrIce.time();
  const RealVect& levelDx = a_amrIce.dx(a_level);

  if (m_timeFileMap != NULL)
    {
      std::map<Real,std::string>::const_iterator start,end;
      if (m_timeFileMap->size() > 1)
	{
	  end = m_timeFileMap->begin();
	  while (end != m_timeFileMap->end() && end->first <= time )
	    {
	      ++end;
	    }
	  start = end;
	  if (start != m_timeFileMap->begin())
	    --start;
	  if (end == m_timeFileMap->end())
	    end = start;
	}
      else
	{
	  start = end = m_timeFileMap->begin();
	}

      Vector<std::string> name(1,m_name);
     
      if (start->first != m_startTime)
	{
	  //load m_startFlux
	  pout() << " LevelDataSurfaceFlux::surfaceThicknessFlux loading start time data " << start->second << std::endl;
	  Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_startFlux);
	  Real dx;
	  readLevelData(data,dx,start->second,name,1);
	  for (int dir=0;dir<SpaceDim;dir++)
	    m_dx[dir] = dx;
	  m_startTime = start->first;
	}
      if (end->first != m_endTime)
	{
	  if ((start != end) & (m_linearInterp))
	    {
	      //load m_endFlux
	      pout() << " LevelDataSurfaceFlux::surfaceThicknessFlux loading end time data for linear interpolation " << end->second << std::endl;
	      Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_endFlux);
	      Real dx;
	      readLevelData(data, dx , end->second,name,1);
	      for (int dir=0;dir<SpaceDim;dir++)
		CH_assert(m_dx[dir] = dx);
	      m_endTime = end->first;
	    }
	  else
	    {
	      pout() << " LevelDataSurfaceFlux::surfaceThicknessFlux piecewise const " << endl;
	    }
	  
	  m_endTime = end->first;
	}
    }
  
  for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(m_defaultValue);
    }

  FillFromReference(a_flux, *m_startFlux, levelDx ,m_dx,m_verbose);

  if (m_linearInterp)
    {
      // Linear interpolation of surface fluxes in time
      if (time > m_startTime && m_startTime < m_endTime)
	{
     	  
	  Real w = std::min(1.0 , (time - m_startTime) / (m_endTime - m_startTime)); 
	  LevelData<FArrayBox> tmp; tmp.define(a_flux);
	  for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
	    {
	      tmp[dit].setVal(m_defaultValue);
	    }

          
	  FillFromReference(tmp, *m_endFlux, levelDx ,m_dx,m_verbose);
	  for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
	    {
	      tmp[dit] *=w;
	      a_flux[dit] *= (1.0-w);
	      a_flux[dit] += tmp[dit];
	    }
	}
    }
  
  const ProblemDomain& domain = a_flux.disjointBoxLayout().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_flux, domain, dir, Side::Lo);
	ReflectGhostCells(a_flux, domain, dir, Side::Hi);
      }
    }
  a_flux.exchange();
  
  

}

void MultiLevelDataSurfaceFlux::surfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
    
  const Real& time = a_amrIce.time();
  const RealVect& levelDx = a_amrIce.dx(a_level);

  if (m_timeFileMap != NULL)
    {
      std::map<Real,std::string>::const_iterator start,end;
      if (m_timeFileMap->size() > 1)
	{
	  end = m_timeFileMap->begin();
	  while (end != m_timeFileMap->end() && end->first <= time )
	    {
	      ++end;
	    }
	  start = end;
	  if (start != m_timeFileMap->begin())
	    --start;
	  if (end == m_timeFileMap->end())
	    end = start;
	}
      else
	{
	  start = end = m_timeFileMap->begin();
	}

      Vector<std::string> name(1,m_name);
      
      if (start->first != m_startTime)
	{
	  //load m_startFlux
	  pout() << " LevelDataSurfaceFlux::surfaceThicknessFlux loading start time data " << start->second << std::endl;
	  if (m_linearInterp)
	    {
	      pout() << " LevelDataSurfaceFlux::surfaceThicknessFlux linearly interpolating between time data " << endl;
	    }
	  else
	    {
	      pout() << " LevelDataSurfaceFlux::surfaceThicknessFlux piecewise const " << endl;
	    }
	  Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > data;
	  Real dxCrse;
	  readMultiLevelData(data,dxCrse,m_ratio,start->second,name,1);
	  for (int dir=0;dir<SpaceDim;dir++)
	    m_dxCrse[dir] = dxCrse;
	  m_startFlux.resize(data[0].size());
	  for (int lev = 0; lev < m_startFlux.size(); lev++)
	    {
	      m_startFlux[lev] = data[0][lev];
	    }


	  m_startTime = start->first;
	}
      if (end->first != m_endTime)
	{
	  if (start == end)
	    {
	      m_endFlux = m_startFlux;
	    }
	  else
	    {
	      //load m_endFlux
	      pout() << " LevelDataSurfaceFlux::surfaceThicknessFlux loading end time data " << end->second << std::endl;
	      Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > data;
	      Real dxCrse;
	      Vector<int> ratio;
	      readMultiLevelData(data,dxCrse,ratio,start->second,name,1);
	      for (int dir=0;dir<SpaceDim;dir++)
		CH_assert(m_dxCrse[dir] = dxCrse);
              m_endFlux.resize(data[0].size());
	      for (int lev = 0; lev < m_endFlux.size(); lev++)
		{
		  CH_assert(ratio[lev] = m_ratio[lev]);
		  m_endFlux[lev] = data[0][lev];
		}

	      m_endTime = end->first;
	    }
	  m_endTime = end->first;
	}
    }
  
  for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(m_defaultValue);
    }
  
  RealVect dx(m_dxCrse);
  for (int refDataLev = 0; refDataLev < m_startFlux.size(); refDataLev++) 
    {
      FillFromReference(a_flux, *m_startFlux[refDataLev], levelDx ,dx,m_verbose);
      dx /= Real(m_ratio[refDataLev]);
    }

  if (time > m_startTime && m_startTime < m_endTime)
    {
      
      Real w;
      if (m_linearInterp)
	{
	  // Linear interpolation between surface fluxes
	  w = std::min(1.0 , (time - m_startTime) / (m_endTime - m_startTime)); 
	}
      else
	{
	  // Piecewise constant
	  w=0.0;
	} 
     

      LevelData<FArrayBox> tmp; tmp.define(a_flux);
      for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
	{
	  tmp[dit].setVal(m_defaultValue);
	}

        
      dx = m_dxCrse;
      for (int refDataLev = 0; refDataLev < m_endFlux.size(); refDataLev++) 
	{
	  FillFromReference(tmp, *m_endFlux[refDataLev], levelDx ,dx,m_verbose);
	  dx /= Real(m_ratio[refDataLev]);
	}
      for (DataIterator dit= a_flux.dataIterator(); dit.ok(); ++dit)
  	{
  	  tmp[dit] *=w;
  	  a_flux[dit] *= (1.0-w);
  	  a_flux[dit] += tmp[dit];
  	}
    }
  
  const ProblemDomain& domain = a_flux.disjointBoxLayout().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_flux, domain, dir, Side::Lo);
	ReflectGhostCells(a_flux, domain, dir, Side::Hi);
      }
    }
  a_flux.exchange();
  
  

}



#include "NamespaceFooter.H"
