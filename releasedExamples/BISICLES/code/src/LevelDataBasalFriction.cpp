#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LevelDataBasalFriction.H"
#include "FillFromReference.H"
#include "ReadLevelData.H"
#include "NamespaceHeader.H"

// define basal friction coefficient C and place in a_betaSqr
/* Either copy data from a single stored level, or load files as needed
    depending on a_time
*/
void LevelDataBasalFriction::setBasalFriction(LevelData<FArrayBox>& a_C,
					      LevelSigmaCS& a_coordSys,
					      Real a_time,
					      Real a_dt)
{
  if (m_timeFileMap == NULL)
    {
      // no time-file map, just a stored LevelData<FArrayBox> 
      FillFromReference(a_C, *m_C, a_coordSys.dx(),m_dx,m_verbose);
    }
  else
    {
      //work out which of the times that map to a file  bracket the current time
      std::map<Real,std::string>::const_iterator start,end;
      if (m_timeFileMap->size() > 1)
	{
	  end = m_timeFileMap->begin();
	  while (end != m_timeFileMap->end() && end->first <= a_time )
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
	  //load m_C
	  pout() << " LevelDataBasalFriction::setBasalFriction loading start time data " << start->second << std::endl;
	  Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_C);
	  Real dx;
	  readLevelData(data,dx,start->second,name,1);
	  for (int dir=0;dir<SpaceDim;dir++)
	    m_dx[dir] = dx;
	  m_startTime = start->first;
	}
      if (end->first != m_endTime)
	{
	  m_endTime = end->first;
	  if (start == end)
	    {
	      m_endC = m_C;
	      }
	  else
	    {
	      pout() << " LevelDataBasalFriction::setBasalFriction loading end time data " << end->second << std::endl;
	      Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_endC);
	      Real dx;
	      readLevelData(data, dx , end->second,name,1);
	      for (int dir=0;dir<SpaceDim;dir++)
		CH_assert(m_dx[dir] = dx);
	    }
	} //end bracket calculation
      FillFromReference(a_C, *m_C, a_coordSys.dx(),m_dx,m_verbose);
      if (a_time > m_startTime && m_startTime < m_endTime)
	{
	  Real w = std::min(1.0 , (a_time - m_startTime) / (m_endTime - m_startTime));
	  LevelData<FArrayBox> tmp; tmp.define(a_C);
	  FillFromReference(tmp, *m_endC, a_coordSys.dx() ,m_dx,m_verbose);
	  for (DataIterator dit=a_C.dataIterator(); dit.ok(); ++dit)
	    {
	      tmp[dit] *=w;
	      a_C[dit] *= (1.0-w);
	      a_C[dit] += tmp[dit];
	    }
	}
    } // end if (m_timeFileMap != NULL)
}

// define basal friction coefficient C and place in a_betaSqr
/* Copy data from a single stored AMR HIerarchy
 */
void MultiLevelDataBasalFriction::setBasalFriction(LevelData<FArrayBox>& a_C, 
						   RealVect a_dx,
						   Real a_time,
						   Real a_dt)
{
  // rather than doing single-level-at-a-time approach, switch to calling 
  // flattenCellData function, which fills each level from the AMR hierarchy 
  // in a way that places no restricions on whether the fine grids are 
  // coarseable down to all of the the coarser-level grids (which is a 
  // side-effect requirement of the original approach)
#if 0
  RealVect dx(m_dxCrse);
  for (int refDataLev = 0; refDataLev < m_C.size(); refDataLev++)
    {
      
      if (a_dx[0] < dx[0] && refDataLev > 0)
	{
	  //We need a LevelData<FArrayBox> whose DisjointBoxLayout covers
	  //a_C's, but m_C[refDataLev] doesn't necessarily provide one. Since we secretly
	  //want to be lisp programmers, we turn to recursion to build one,
	  //throwing performance worries to the wind.
	  //\todo consider modifying the BasalFriction interface to avoid recompuing levels 1 to n-1
	  const int& nRef = int(dx[0]/a_dx[0]);
	  DisjointBoxLayout crseDBL;
	  coarsen(crseDBL, a_C.disjointBoxLayout(), nRef);
	  LevelData<FArrayBox> crseC(crseDBL, 1 , IntVect::Zero);
	  this->setBasalFriction(crseC , dx, a_time, a_dt);
	  FillFromReference(a_C, crseC , a_dx, dx ,m_verbose);
	}
      else
	{
	  FillFromReference(a_C, *m_C[refDataLev], a_dx, dx ,m_verbose);
	}
      dx /= Real(m_ratio[refDataLev]);
    }
# endif

  Vector<RealVect> refDx (m_C.size(), -1.0*RealVect::Unit);
  refDx[0] = m_dxCrse;
  for (int lev=1; lev<refDx.size(); lev++)
    {
      refDx[lev] = refDx[lev-1]/m_ratio[lev-1];
    }
  
  flattenCellData(a_C,
                  a_dx,
                  m_C,
                  refDx,
                  m_verbose);

}
#include "NamespaceFooter.H"
