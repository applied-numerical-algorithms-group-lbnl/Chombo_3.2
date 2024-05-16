#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "MuCoefficient.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "ParmParse.H"
#include "ReadLevelData.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "NamespaceHeader.H"

///assemble a MuCoefficient  object from ParmParse input, return pointer
MuCoefficient* MuCoefficient::parseMuCoefficient(const char* a_prefix)
{
  MuCoefficient* ptr = NULL;
  ParmParse muPP(a_prefix);
  std::string muCoefType = "unit";
  muPP.query("type",muCoefType );
  if (muCoefType == "unit")
    {
      ptr = static_cast<MuCoefficient*>(new UnitMuCoefficient());
    }
  else if (muCoefType == "LevelData")
    {
      //read a one level muCoef from an AMR Hierarchy, and  store it in a LevelDataMuCoeffcient
      ParmParse ildPP("inputLevelData");
      std::string muCoefName = "muCoef";
      ildPP.query("muCoefName",muCoefName);
      std::string infileFormat = "";
      ildPP.query("muCoefFileFormat",infileFormat);

      if (infileFormat == "")
	//basic case : one file is specified. this might be unnecessary since the multiple file case
	  //includes the one file case, but for now we are holding onto the older code
	{
	  std::string infile;
	  ildPP.get("muCoefFile",infile);
	  
	  RefCountedPtr<LevelData<FArrayBox> > levelMuCoef (new LevelData<FArrayBox>());
	  Vector<RefCountedPtr<LevelData<FArrayBox> > > vectMuCoef;
	  vectMuCoef.push_back(levelMuCoef);
	  Vector<std::string> names(1);
	  names[0] = muCoefName;
	  Real dx;
	  readLevelData(vectMuCoef,dx,infile,names,1);
	  RealVect levelDx = RealVect::Unit * dx;
	  ptr = static_cast<MuCoefficient*>
	    (new LevelDataMuCoefficient(levelMuCoef,levelDx));
	}
      else
	{
	  int n = 1;
	  ildPP.query("muCoefFileSteps",n);
	  
	  int startTime = 0, timeStep = 1;
	  ildPP.query("muCoefFileStartTime", startTime);
	  ildPP.query("muCoefFileTimeStep", timeStep);
	  RefCountedPtr<std::map<Real,std::string> > tf
	    (new std::map<Real,std::string>);
	  for (int i =0; i < n; i++)
	    {
	      char* file = new char[infileFormat.length()+32];
	      int j = startTime + i * timeStep;
	      sprintf(file, infileFormat.c_str(),j);
	      tf->insert(make_pair(Real(j), file));
	      delete file;
	    }
	  
	  ptr = static_cast<MuCoefficient*>
	    (new LevelDataMuCoefficient(tf,muCoefName));
	}
	  
      }
    else if (muCoefType == "MultiLevelData")
      {
	//read a multi level muCoef from an AMR Hierarchy, and  store it in a MultiLevelDataMuCoeffcient
	 ParmParse ildPP("inputLevelData");
	 std::string infile;
	 ildPP.get("muCoefFile",infile);
	 std::string muCoefName = "muCoef";
	 ildPP.query("muCoefName",muCoefName);
	 
	 Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectMuCoef;
	 Vector<std::string> names(1);
	 names[0] = muCoefName;
	 Real dx;
	 Vector<int> ratio;
	 readMultiLevelData(vectMuCoef,dx,ratio,infile,names,1);
	 RealVect dxCrse = RealVect::Unit * dx;
	 ptr = static_cast<MuCoefficient*>
	   (new MultiLevelDataMuCoefficient(vectMuCoef[0],dxCrse,ratio));
	 
      }
#ifdef HAVE_PYTHON
    else if (muCoefType == "Python")
      {
	ParmParse pyPP("PythonMuCoefficient");
	std::string module;
	pyPP.get("module",module);
	std::string funcName = "muCoefficient";
	pyPP.query("function",funcName);
	ptr =  static_cast<MuCoefficient*>
	  (new PythonInterface::PythonMuCoefficient(module, funcName));
      }
#endif
    else
      {
	MayDay::Error("undefined MuCoefficient in inputs");
      }
  
  return ptr;
  
}


MuCoefficient* UnitMuCoefficient::new_muCoefficient() const
{
  UnitMuCoefficient* ptr = new UnitMuCoefficient();
  return static_cast<MuCoefficient*>(ptr);

}



void UnitMuCoefficient::setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
					 LevelSigmaCS& a_coordSys,
					 Real a_time,
					 Real a_dt)
{
  for (DataIterator dit = a_cellMuCoef.dataIterator(); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].setVal(1.0);
    }
}


MuCoefficient* AxbyMuCoefficient::new_muCoefficient() const
{
  return static_cast<MuCoefficient*>
    (new AxbyMuCoefficient(m_a, m_x,m_b, m_y) );
}

AxbyMuCoefficient::AxbyMuCoefficient
(const Real& a_a, MuCoefficient* a_x, 
 const Real& a_b, MuCoefficient* a_y)
{

  m_a = a_a;
  m_b = a_b;
  
  CH_assert(a_x != NULL);
  CH_assert(a_y != NULL);
  m_x = a_x->new_muCoefficient();
  m_y = a_y->new_muCoefficient();
  CH_assert(m_x != NULL);
  CH_assert(m_y != NULL);

}

AxbyMuCoefficient::~AxbyMuCoefficient()
{
  if (m_x != NULL)
    {
      delete m_x; m_x = NULL;
    }
  if (m_y != NULL)
    {
      delete m_y; m_y = NULL;
    }
}

void AxbyMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_cellMuCoef,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{

  LevelData<FArrayBox> y_mu(a_cellMuCoef.disjointBoxLayout(),1,a_cellMuCoef.ghostVect());
  m_x->setMuCoefficient(a_cellMuCoef, a_coordSys, a_time, a_dt);
  m_y->setMuCoefficient(y_mu, a_coordSys, a_time, a_dt);
  for (DataIterator dit(a_cellMuCoef.disjointBoxLayout()); dit.ok(); ++dit)
    {
      a_cellMuCoef[dit].axby(a_cellMuCoef[dit],y_mu[dit],m_a,m_b);
    }

}



MuCoefficient* LevelDataMuCoefficient::new_muCoefficient() const
{
  LevelDataMuCoefficient* ptr;
  if (m_timeFileMap == NULL)
    {
      ptr = new LevelDataMuCoefficient(m_muCoef,m_dx);
    }
  else
    {
      RefCountedPtr<std::map<double, std::basic_string<char> > > r = m_timeFileMap;
      ptr = new LevelDataMuCoefficient(r,m_name);
    }
  
  return static_cast<MuCoefficient*>(ptr);

}
void LevelDataMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_muCoef,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{
  if (m_timeFileMap == NULL)
    {
      // no time-file map, just a stored LevelData<FArrayBox> 
      FillFromReference(a_muCoef, *m_muCoef, a_coordSys.dx(),m_dx,m_verbose);
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
	  //load m_muCoef
	  pout() << " LevelDataMuCoefficient::setMuCoefficient loading start time data " << start->second << std::endl;
	  Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_muCoef);
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
	      m_endMuCoef = m_muCoef;
	      }
	  else
	    {
	      pout() << " LevelDataMuCoefficient::setMuCoefficient loading end time data " << end->second << std::endl;
	      Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,m_endMuCoef);
	      Real dx;
	      readLevelData(data, dx , end->second,name,1);
	      for (int dir=0;dir<SpaceDim;dir++)
		CH_assert(m_dx[dir] = dx);
	    }
	} //end bracket calculation
      FillFromReference(a_muCoef, *m_muCoef, a_coordSys.dx(),m_dx,m_verbose);
      if (a_time > m_startTime && m_startTime < m_endTime)
	{
	  Real w = std::min(1.0 , (a_time - m_startTime) / (m_endTime - m_startTime));
	  LevelData<FArrayBox> tmp; tmp.define(a_muCoef);
	  FillFromReference(tmp, *m_endMuCoef, a_coordSys.dx() ,m_dx,m_verbose);
	  for (DataIterator dit=a_muCoef.dataIterator(); dit.ok(); ++dit)
	    {
	      tmp[dit] *=w;
	      a_muCoef[dit] *= (1.0-w);
	      a_muCoef[dit] += tmp[dit];
	    }
	}
    } // end if (m_timeFileMap != NULL)

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      const ProblemDomain domain = a_muCoef.disjointBoxLayout().physDomain();
      if (!(domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_muCoef, domain, dir, Side::Lo);
	  ReflectGhostCells(a_muCoef, domain, dir, Side::Hi);
	}
    }
   
}

MuCoefficient* MultiLevelDataMuCoefficient::new_muCoefficient() const
{
  MultiLevelDataMuCoefficient* ptr = new MultiLevelDataMuCoefficient(m_muCoef, m_dxCrse, m_ratio);
  return static_cast<MuCoefficient*>(ptr);

}


void MultiLevelDataMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_cellMuCoef,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt) 
{
    
  
    setMuCoefficient(a_cellMuCoef,a_coordSys.dx(),a_time,a_dt); 
  
}

void MultiLevelDataMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_cellMuCoef,
 RealVect a_dx,
 Real a_time,
 Real a_dt)
{

  
  Vector<RealVect> refDx (m_muCoef.size(), -1.0*RealVect::Unit);
  refDx[0] = m_dxCrse;
  for (int lev=1; lev<refDx.size(); lev++)
    {
      refDx[lev] = refDx[lev-1]/m_ratio[lev-1];
    }

  flattenCellData(a_cellMuCoef,
                  a_dx,
                  m_muCoef,
                  refDx,
                  m_verbose);
  

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      const ProblemDomain domain = a_cellMuCoef.disjointBoxLayout().physDomain();
      if (!(domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_cellMuCoef, domain, dir, Side::Lo);
	  ReflectGhostCells(a_cellMuCoef, domain, dir, Side::Hi);
	}
    }

  
}


#include "NamespaceFooter.H"
