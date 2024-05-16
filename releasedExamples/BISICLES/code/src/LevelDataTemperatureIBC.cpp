#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  LevelDataTemperatureIBC.cpp
// ============
//
// PhysIBC-derived class which stores initial temperature  data
// and imposes either periodic or reflection boundary conditions

#include "LevelDataTemperatureIBC.H"
#include "IceConstants.H"
#include "IceThermodynamics.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "ReadLevelData.H"
#include "AmrIceBase.H"
#include "NamespaceHeader.H"


LevelDataTemperatureIBC* 
LevelDataTemperatureIBC::parse(ParmParse& a_pp)
{
  bool dataIsInternalEnergy(false);
  a_pp.query("readInternalEnergy", dataIsInternalEnergy);

  std::string infile;
  std::string layer0Name(""); 
  std::string surfaceName("");
  
  if (dataIsInternalEnergy)
    {
      a_pp.get("internalEnergyFile",infile);
      layer0Name = "internalEnergy000000";
      a_pp.query("internalEnergyName",layer0Name);
      a_pp.query("surfaceInternalEnergyName",surfaceName);
    }
  else
    {
      a_pp.get("temperatureFile",infile);
      layer0Name = "temp000000";
      a_pp.query("temperatureName",layer0Name);
      a_pp.query("surfaceTemperatureName",surfaceName);
    }
  	
  Real defaultTemperature = 258.0; // 
  a_pp.query("defaultTemperature", defaultTemperature);
	
  RefCountedPtr<LevelData<FArrayBox> > levelBulkData
      (new LevelData<FArrayBox>());
  Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
  vectData.push_back(levelBulkData);
  Vector<std::string> names(1);
  names[0] = layer0Name;
  Real dx;
  ParmParse ppAmr ("amr");
  Vector<int> ancells(3,0); 
  ppAmr.queryarr("num_cells", ancells, 0, ancells.size());
  if (ancells[0] == 0)
    {
      ParmParse ppGeo ("geometry");
      ppGeo.getarr("num_cells", ancells, 0, ancells.size());
    }
  readLevelData(vectData,dx,infile,names,ancells[2]);
  RealVect levelDx = RealVect::Unit * dx;
  
  RefCountedPtr<LevelData<FArrayBox> > levelSurfaceData(new LevelData<FArrayBox>());
    
  if (surfaceName == "")
    {
      //not ideal, but in this case copy the top layer data to the surface.
      levelSurfaceData->define( levelBulkData->disjointBoxLayout(), 1, levelBulkData->ghostVect());
      levelBulkData->copyTo(Interval(0,0),*levelSurfaceData,Interval(0,0));
    }
  else
    {
      names.resize(1);
      names[0] = surfaceName;
      vectData[0] = levelSurfaceData;
      readLevelData(vectData,dx,infile,names,1);
      ; 
      if (dx != levelDx[0])
	{
	  pout() << "surface temperature dx = " << dx << " but bulk temperature mesh dx = " << levelDx[0] << endl;
	  CH_assert(dx == levelDx[0]);
	  MayDay::Error("dx != levelDx[0]");
	}
    }
  
  RefCountedPtr<LevelData<FArrayBox> > levelBasalHeatFlux(new LevelData<FArrayBox>());
  std::string basalHeatFluxName = "";
  a_pp.query("basalHeatFluxName",basalHeatFluxName);
  if (basalHeatFluxName == "")
    {
      //if no basal heat flux is given, assume zero flux
      levelBasalHeatFlux->define( levelBulkData->disjointBoxLayout(), 1, levelBulkData->ghostVect());
      for (DataIterator dit = levelBasalHeatFlux->dataIterator(); dit.ok(); ++dit)
	{
	  (*levelBasalHeatFlux)[dit].setVal(0.0);
	}
    }
  else
    {
      names.resize(1);
      names[0] = basalHeatFluxName;
      vectData[0] = levelBasalHeatFlux;
      readLevelData(vectData,dx,infile,names,1);
      if (dx != levelDx[0])
	{
	  pout() << "basal temperature dx = " << dx << " but bulk temperature mesh dx = " << levelDx[0] << endl;
	  CH_assert(dx == levelDx[0]);
	  MayDay::Error("dx != levelDx[0]");
	}
    }
  return new LevelDataTemperatureIBC
    (levelBulkData,levelSurfaceData,levelBasalHeatFlux,levelDx, defaultTemperature, dataIsInternalEnergy);

}

LevelDataTemperatureIBC::LevelDataTemperatureIBC
(RefCountedPtr<LevelData<FArrayBox> > a_bulkData, 
 RefCountedPtr<LevelData<FArrayBox> > a_surfaceData,
 RefCountedPtr<LevelData<FArrayBox> > a_basalHeatFlux,
 const RealVect& a_dx, const Real& a_defaultTemperature,
 const bool& a_dataIsInternalEnergy)
{
  m_bulkData = a_bulkData;
  m_surfaceData = a_surfaceData;
  m_basalHeatFlux = a_basalHeatFlux;
  m_defaultTemperature =  a_defaultTemperature;
  m_defaultE = IceThermodynamics::m_ice_heat_capacity * m_defaultTemperature ;
  m_dx = a_dx;
  m_dataIsInternalEnergy = a_dataIsInternalEnergy;
}

LevelDataTemperatureIBC::~LevelDataTemperatureIBC()
{
  
}

void LevelDataTemperatureIBC::define(const ProblemDomain& a_domain,
				     const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

void LevelDataTemperatureIBC::basalHeatFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  FillFromReference(a_flux,*m_basalHeatFlux,a_amrIce.dx(a_level),m_dx,true);
  const ProblemDomain& domain = a_amrIce.grids(a_level).physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_flux, domain, dir, Side::Lo);
	ReflectGhostCells(a_flux, domain, dir, Side::Hi);
      }
    }

}

void LevelDataTemperatureIBC::initializeIceInternalEnergy(LevelData<FArrayBox>& a_E,
							  LevelData<FArrayBox>& a_tillWaterDepth,
							  LevelData<FArrayBox>& a_surfaceE, 
							  LevelData<FArrayBox>& a_basalE,
							  const AmrIceBase& a_amrIce, 
							  int a_level, Real a_dt)
{

  if (true)
    {
      pout() << " LevelDataIBC::initializeIceInternalEnergy" << endl;
    }

  
  const LevelSigmaCS& coordSys = *a_amrIce.geometry(a_level); 
  const DisjointBoxLayout dbl = coordSys.grids();

  if (m_dataIsInternalEnergy)
    {
      //simple case
      for (DataIterator dit(dbl);dit.ok();++dit)
	{
	  a_E[dit].setVal(m_defaultE);
	  a_surfaceE[dit].setVal(m_defaultE);
	}
      FillFromReference(a_E,*m_bulkData,coordSys.dx(),m_dx,true);
      FillFromReference(a_surfaceE,*m_surfaceData,coordSys.dx(),m_dx,true);
    }
  else
    {
      // need to convert to internal energy
      LevelData<FArrayBox> T(dbl, a_E.nComp(), a_E.ghostVect());
      LevelData<FArrayBox> sT(dbl, 1, a_E.ghostVect());
      
      // set default tempeature - applies to regions of the domain not
      // covered by the data
      for (DataIterator dit(dbl);dit.ok();++dit)
	{
	  T[dit].setVal(m_defaultTemperature);
	  sT[dit].setVal(m_defaultTemperature);
	}
      
      FillFromReference(T,*m_bulkData,coordSys.dx(),m_dx,true);
      FillFromReference(sT,*m_surfaceData,coordSys.dx(),m_dx,true);
      {
	const ProblemDomain& domain = coordSys.grids().physDomain();
	for (int dir = 0; dir < SpaceDim; ++dir)
	  {
	    if (!(domain.isPeriodic(dir))){
	      ReflectGhostCells(T, domain, dir, Side::Lo);
	      ReflectGhostCells(T, domain, dir, Side::Hi);
	      ReflectGhostCells(sT, domain, dir, Side::Lo);
	      ReflectGhostCells(sT, domain, dir, Side::Hi);
	    }
	  }
      }
      
      for (DataIterator dit(dbl);dit.ok();++dit)
	{
	  FArrayBox w( a_E[dit].box(), a_E[dit].nComp()); //water fraction, set to zero for now
	  w.setVal(0.0);
	  IceThermodynamics::composeInternalEnergy(a_E[dit],T[dit],w, a_E[dit].box() );
	}
      
      
      for (DataIterator dit(dbl);dit.ok();++dit)
	{
	  FArrayBox w( a_surfaceE[dit].box(), 1); //water fraction, set to zero for now 
	  w.setVal(0.0);
	  IceThermodynamics::composeInternalEnergy(a_surfaceE[dit],sT[dit],w,a_surfaceE[dit].box() );
	}
    } // end conversion
  
  for (DataIterator dit(dbl);dit.ok();++dit)
    {
      //todo : read this also
      a_tillWaterDepth[dit].setVal(0.0);
    }
   
  const ProblemDomain& domain = coordSys.grids().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_E, domain, dir, Side::Lo);
	ReflectGhostCells(a_E, domain, dir, Side::Hi);
	ReflectGhostCells(a_surfaceE, domain, dir, Side::Lo);
	ReflectGhostCells(a_surfaceE, domain, dir, Side::Hi);
      }
    }

  a_E.exchange();
  a_surfaceE.exchange();

}


LevelDataTemperatureIBC* 
LevelDataTemperatureIBC::new_internalEnergyIBC()
{
  return new LevelDataTemperatureIBC(m_bulkData,m_surfaceData,m_basalHeatFlux,
				     m_dx,m_defaultTemperature,m_dataIsInternalEnergy);
}




#include "NamespaceFooter.H"
