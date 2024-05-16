 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "IceInternalEnergyIBC.H"
#include "IceThermodynamics.H"
#include "IceConstants.H"
#include "BoxIterator.H"
#include "ExtrapGhostCells.H"
#include "ReflectGhostCells.H"
#include "NamespaceHeader.H"

ConstantIceTemperatureIBC* 
ConstantIceTemperatureIBC::new_internalEnergyIBC()
{
  return new ConstantIceTemperatureIBC(m_T);
}

/// set a basal heat flux to zero. units are Joules / Year
void ConstantIceTemperatureIBC::basalHeatFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  for (DataIterator dit = a_flux.dataIterator();dit.ok();++dit)
    {
      a_flux[dit].setVal(0.0);
    }
}

#if BISICLES_Z == BISICLES_LAYERED
void 
ConstantIceTemperatureIBC::initializeIceInternalEnergy
(LevelData<FArrayBox>& a_E,
 LevelData<FArrayBox>& a_tillWaterDepth,
 LevelData<FArrayBox>& a_surfaceE, 
 LevelData<FArrayBox>& a_basalE,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  
  for (DataIterator dit(a_E.disjointBoxLayout()); dit.ok() ; ++dit)
    {
      FArrayBox T(a_E[dit].box(),1);
      T.setVal(m_T);
      FArrayBox w(a_E[dit].box(),1);
      w.setVal(0.0);
      FArrayBox E(a_E[dit].box(),1);
      IceThermodynamics::composeInternalEnergy(E,T,w,a_E[dit].box());
      for (int l = 0; l < a_E[dit].nComp(); ++l)
	{
	  a_E[dit].copy(E,0,l,1);
	}	    
      a_surfaceE[dit].copy(E);
      a_basalE[dit].copy(E);
      a_tillWaterDepth[dit].setVal(0.0);
    }
}

void 
ConstantIceTemperatureIBC::setIceInternalEnergyBC
(LevelData<FArrayBox>& a_E,
 LevelData<FArrayBox>& a_tillWaterDepth,
 LevelData<FArrayBox>& a_surfaceE, 
 LevelData<FArrayBox>& a_basalE,
 const LevelSigmaCS& a_coordSys)
{
  const ProblemDomain& domain = a_coordSys.grids().physDomain();
  ExtrapGhostCells( a_E, domain);
  ExtrapGhostCells( a_basalE, domain);
  ExtrapGhostCells( a_surfaceE, domain);
  ExtrapGhostCells( a_tillWaterDepth, domain);
}


#elif BISICLES_Z == BISICLES_FULLZ
#error BISICLES_FULLZ not implemented
#endif

void 
IceInternalEnergyIBC::initialize(LevelData<FArrayBox>& a_U)
{
  //we shouldn't get here
  MayDay::Error("IceInternalEnergyIBC::initialize not implemented");
}

void  IceInternalEnergyIBC::primBC
(FArrayBox&            a_WGdnv,
 const FArrayBox&      a_Wextrap,
 const FArrayBox&      a_W,
 const int&            a_dir,
 const Side::LoHiSide& a_side,
 const Real&           a_time)
{

  // do nothing in periodic case
  if (!m_domain.isPeriodic(a_dir))
    {
      
#if CH_SPACEDIM == 2
      // 2D case : We are at either an ice divide (where u = 0)
      // or an outflow (where extrapolation is fine)
      // \todo : support an inflow with a known internalEnergy?
   int lohisign;
   
   //find the strip of cells just inside the domain. DFM notes that 
   //a_WGdnv might have more than one layer of ghost cells, so this
   //is slightly more complicated than just testing tmp is at the domain
   //edge
   Box tmp = a_WGdnv.box();
   lohisign = sign(a_side);
   tmp.shiftHalf(a_dir,lohisign);
   const Box& dbox = m_domain.domainBox();
   Box ghostBox = (a_side == Side::Lo)?adjCellLo(dbox,a_dir,1):adjCellHi(dbox,a_dir,1);
   ghostBox &= tmp;
 
   if (!ghostBox.isEmpty() && !m_domain.contains(tmp))
     {
       tmp &= m_domain;
       Box boundaryBox = (a_side == Side::Lo)?bdryLo(tmp,a_dir):bdryHi(tmp,a_dir);
       BoxIterator bit(boundaryBox);
       for (bit.begin(); bit.ok(); ++bit){
	 const IntVect& i = bit();
	 a_WGdnv(i,0) = a_Wextrap(i,0);
	 
       }
     }

#elif CH_SPACEDIM == 3
  MayDay::Error("IceThicknessIBC::primBC DIM = 3 not yet impplemented");
#else
  MayDay::Error("IceThicknessIBC::primBC only supports DIM = 2 and DIM = 3");
#endif

    }
}

void ReflectionIceInternalEnergyIBC::setIceInternalEnergyBC
(LevelData<FArrayBox>& a_E, 
 LevelData<FArrayBox>& a_tillWaterDepth,
 LevelData<FArrayBox>& a_surfaceE, 
 LevelData<FArrayBox>& a_basalE,
 const LevelSigmaCS& a_coordSys)
{
  const ProblemDomain& domain = a_coordSys.grids().physDomain();
  
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_E, domain, dir, Side::Lo);
	  ReflectGhostCells(a_E, domain, dir, Side::Hi);
	  ReflectGhostCells(a_tillWaterDepth, domain, dir, Side::Lo);
	  ReflectGhostCells(a_tillWaterDepth, domain, dir, Side::Hi);
	  ReflectGhostCells(a_surfaceE, domain, dir, Side::Lo);
	  ReflectGhostCells(a_surfaceE, domain, dir, Side::Hi);
	  ReflectGhostCells(a_basalE, domain, dir, Side::Lo);
	  ReflectGhostCells(a_basalE, domain, dir, Side::Hi);
	}
    }
}

#include "NamespaceFooter.H"
