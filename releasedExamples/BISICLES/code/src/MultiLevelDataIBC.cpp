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
// MultiLevelDataIBC.cpp
// ============
//
// PhysIBC-derived class which stores initial topography and thickness data
// and imposes either periodic or reflection boundary conditions


#include "MultiLevelDataIBC.H"
#include "ReflectGhostCells.H"
#include "FillFromReference.H"
#include "IceConstants.H"
#include "NamespaceHeader.H"

MultiLevelDataIBC::MultiLevelDataIBC
(const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_thck,
 const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_topg, 
 const RealVect& a_dxCrse, const Vector<int> & a_refRatio)
  :m_isBCsetUp(false),m_verbose(true),m_thck(a_thck),
   m_topg(a_topg),m_dxCrse(a_dxCrse),m_refRatio(a_refRatio)
{
  m_nLevel = a_thck.size();
  CH_assert(a_topg.size() == m_nLevel);
  CH_assert(a_refRatio.size() >= m_nLevel);
}

MultiLevelDataIBC::~MultiLevelDataIBC()
{
  
}

void MultiLevelDataIBC::define(const ProblemDomain& a_domain,
			       const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}


void  reflectionBC_MLDIBC(FArrayBox& a_vel,
			  const Box& a_valid,
			  const ProblemDomain& a_domain,
			  Real a_dx,
			  bool a_homogeneous)
{
  
  const IntVect ghostVect = IntVect::Unit;
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_vel, a_domain,ghostVect, dir, Side::Lo);
	  ReflectGhostCells(a_vel, a_domain,ghostVect, dir, Side::Hi);
	  
	  Box ghostBoxLo = adjCellBox(a_valid, dir, Side::Lo, 1);
              
	  if(!a_domain.domainBox().contains(ghostBoxLo))
	    {
	      ghostBoxLo &= a_vel.box();
	      a_vel.mult(-1.0,ghostBoxLo,dir,1);
		  
	    }
	  Box ghostBoxHi = adjCellBox(a_valid, dir, Side::Hi, 1);
	  if(!a_domain.domainBox().contains(ghostBoxHi))
	    {
	      ghostBoxHi &= a_vel.box();
	      a_vel.mult(-1.0,ghostBoxHi,dir,1);
	    }

	}
    }
}


RefCountedPtr<CompGridVTOBC> 
MultiLevelDataIBC::velocitySolveBC()
{
  if (!m_isBCsetUp)
    {
      m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(reflectionBC_MLDIBC));
      m_isBCsetUp = true;
    }
  return m_velBCs;
}

void MultiLevelDataIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
					 const RealVect& a_dx,
					 const RealVect& a_domainSize,
					 const Real& a_time, 
					 const LevelSigmaCS* a_crseCoords,
					 const int a_refRatio)
{

  CH_TIME("MultiLevelDataIBC::initializeIceGeometry");
  if (m_verbose)
    {
      pout() << " MultiLevelDataIBC::initializeIceGeometry" << endl;
    }
  
  Real refRatio = m_dxCrse[0]/a_dx[0];
  //sanity check : a_dx must be an even integer multiple or divisor of m_dxCrse
  Real testDx = a_dx[0]*refRatio;
 
  if ( ((int(refRatio)%2 != 0)
	&& (int(refRatio != 1))
	&& (int(1.0/refRatio)%2 != 0))
       || (abs(testDx - m_dxCrse[0])/m_dxCrse[0] > TINY_NORM))
    {
      MayDay::Error("MultiLevelDataIBC::initializeIceGeometry incompatible a_dx and m_dxCrse");
    }

  //interpolate level from coarser data, or copy.
  int nRef = int(refRatio);
  RealVect dx = m_dxCrse;
  for (int lev = 0; lev < m_nLevel; lev++)
    {
      
      if (nRef > 1)
	{
	  //required level is finer than the data level, so interpolate
	  if ( (a_crseCoords != NULL) && a_refRatio <= nRef)
	    {
	      //we have a valid coarse level LevelSigmaCS, so we might as well interpolate from that
	      if (m_verbose)
		{
		  pout() << "...interpolating from coarse LevelSigmaCS rather than reference data" << endl;
		}
	      a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
	    }
	  else
	    {
	      const DisjointBoxLayout& crseGrids = m_thck[lev]->disjointBoxLayout();
	      LevelSigmaCS crseCoords(crseGrids, dx, a_coords.ghostVect());
	      for (DataIterator dit(crseGrids); dit.ok(); ++dit)
		{
		  crseCoords.getH()[dit].copy ( (*m_thck[lev])[dit]);
		  crseCoords.getTopography()[dit].copy ( (*m_topg[lev])[dit]);
		}
	      crseCoords.recomputeGeometry(NULL,0);
	      a_coords.interpFromCoarse(crseCoords, nRef);
	    }
	  
	}
     
      nRef /= m_refRatio[lev];
      dx /= Real(m_refRatio[lev]);
      
    }

  //coarse average level from finer data, or simple copies
  Real oneOnRef = refRatio;
  dx = m_dxCrse;
  for (int lev = 0; lev < m_nLevel; lev++)
    {
      nRef = int(1.0/oneOnRef + 1.0e-6);
      if (nRef >= 1)
	{
	  //required level is coarser than the current data level; coarsen from there
	  FillFromReference(a_coords.getH(),*m_thck[lev],a_dx,dx,m_verbose);
	  FillFromReference(a_coords.getTopography(),*m_topg[lev],a_dx,dx,m_verbose);
	  //\todo reground if neeeded
	  
	}
      oneOnRef /= Real(m_refRatio[lev]);
      dx /= Real(m_refRatio[lev]);
    } 

  Real dt = 0.0;
  setGeometryBCs(a_coords, a_coords.grids().physDomain(), a_dx, a_time, dt);

}

bool MultiLevelDataIBC::regridIceGeometry(LevelSigmaCS& a_coords,
				     const RealVect& a_dx,
				     const RealVect& a_domainSize,
				     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio)
{
  CH_TIME("MultiLevelDataIBC::regridIceGeometry");
  if (m_verbose)
    {
      pout() << " MultiLevelDataIBC::regridIceGeometry" << endl;
    }
  
  Real refRatio = m_dxCrse[0]/a_dx[0];
  //sanity check : a_dx must be an even integer multiple or divisor of m_dxCrse
  Real testDx = a_dx[0]*refRatio;
 
  if ( ((int(refRatio)%2 != 0)
	&& (int(refRatio != 1))
	&& (int(1.0/refRatio)%2 != 0))
       || (abs(testDx - m_dxCrse[0])/m_dxCrse[0] > TINY_NORM))
    {
      MayDay::Error("MultiLevelDataIBC::regridIceGeometry incompatible a_dx and m_dxCrse");
    }




  //interpolate level from coarser data, or copy.
  int nRef = int(refRatio);

  if ( (a_crseCoords != NULL) && a_refRatio <= nRef)
    {
      return false;
    }
  
  RealVect dx = m_dxCrse;
  for (int lev = 0; lev < m_nLevel; lev++)
    {
      if (nRef > 1)
	{
	 //required level is finer than the data level, so interpolate.
	  if ( (a_crseCoords != NULL) && a_refRatio <= nRef)
	    {
	      //we have a valid coarse level LevelSigmaCS, so we might as well interpolate from that
	      //thi function is deprecated.
	      MayDay::Error("MultiLevelDataIBC::regridIceGeometry interploation deprecated");
	      a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
	    }
	  else
	    {
	      FillFromReference(a_coords.getTopography(),*m_topg[lev],a_dx,dx,m_verbose);
	      if (a_crseCoords!= NULL)
		{
		  // fill ghost cells
		  a_coords.interpFromCoarse(*a_crseCoords, a_refRatio, 
				     false, false, false);
		}
	      
	    }
	  
	}
      nRef /= m_refRatio[lev];
      dx /= Real(m_refRatio[lev]);
      
    }

  //coarse average level from finer data
  Real oneOnRef = refRatio;
  dx = m_dxCrse;
  for (int lev = 0; lev < m_nLevel; lev++)
    {
      nRef = int(1.0/oneOnRef + 1.0e-6);
      if (nRef >= 1)
	{
	  //required level is coarser than the current data level, 
	  FillFromReference(a_coords.getTopography(),*m_topg[lev],a_dx,dx,m_verbose);
	  //\todo reground if neeeded
	  
	}
      oneOnRef /= Real(m_refRatio[lev]);
      dx /= Real(m_refRatio[lev]);
    } 

  return true;



}

IceThicknessIBC* MultiLevelDataIBC::new_thicknessIBC()
{
 
  MultiLevelDataIBC* retval = new MultiLevelDataIBC(m_thck,m_topg,m_dxCrse,m_refRatio) ;
  return static_cast<IceThicknessIBC*>(retval);
}


/// set non-periodic ghost cells for surface height z_s. 
void MultiLevelDataIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
				       LevelSigmaCS& a_coords,
				       const ProblemDomain& a_domain,
				       const RealVect& a_dx,
				       Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Hi);
	}
    }
}

/// set non-periodic ghost cells for thickness & topography
void MultiLevelDataIBC::setGeometryBCs(LevelSigmaCS& a_coords,
				  const ProblemDomain& a_domain,
				  const RealVect& a_dx,
				  Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir))){
	ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Lo);
	ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Hi);
	ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Lo);
	ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Hi);
      }
    }
}

void MultiLevelDataIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("MultiLevelDataIBC::initialize not implemented");
}

// set boundary fluxes
void MultiLevelDataIBC::primBC(FArrayBox&            a_WGdnv,
			   const FArrayBox&      a_Wextrap,
			   const FArrayBox&      a_W,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const Real&           a_time)
{
 if (!m_domain.isPeriodic(a_dir))
   {
     int lohisign;
     Box tmp = a_WGdnv.box();
     // Determine which side and thus shifting directions
     lohisign = sign(a_side);
     tmp.shiftHalf(a_dir,lohisign);
     // (DFM - 5/28/10) this little dance with the ghostBox is a bit 
     // of a kluge to handle the case where a_WGdnv has more than one layer   
     // of ghosting, in which case just testing agains tmp isn't 
     // sufficient to determine whether you're up against the domain edge
     Box ghostBox = (a_side == Side::Lo)?
       adjCellLo(m_domain.domainBox(),a_dir, 1):
       adjCellHi(m_domain.domainBox(),a_dir, 1);
     ghostBox &= tmp;
     // Is there a domain boundary next to this grid
     if (!ghostBox.isEmpty() && !m_domain.contains(tmp))
       {
	 tmp &= m_domain;
	 Box boundaryBox = (a_side == Side::Lo)?
	   bdryLo(tmp,a_dir):bdryHi(tmp,a_dir);
	 BoxIterator bit(boundaryBox);
	 for (bit.begin(); bit.ok(); ++bit){
	   const IntVect& i = bit();
	   a_WGdnv(i,0) = std::max(0.0,a_Wextrap(i,0));
	 }
       }
   }
}
#include "NamespaceFooter.H"
