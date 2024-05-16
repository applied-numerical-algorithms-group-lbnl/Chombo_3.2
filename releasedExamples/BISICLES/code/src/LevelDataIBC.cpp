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
//  LevelDataIBC.cpp
// ============
//
// PhysIBC-derived class which stores initial topography and thickness data
// and imposes either periodic or reflection boundary conditions

#include "LevelDataIBC.H"
#include "ReflectGhostCells.H"
#include "FillFromReference.H"
#include "CoarseAverage.H"
#include "IceConstants.H"
#include "NamespaceHeader.H"

LevelDataIBC::LevelDataIBC( RefCountedPtr<LevelData<FArrayBox> > a_thck, 
			    RefCountedPtr<LevelData<FArrayBox> > a_topg,
			    const RealVect& a_dx,
                            Real a_defaultThickness,
                            Real a_defaultTopography,
                            bool a_setDefaultValues
                            )
                            :m_thck(a_thck),m_topg(a_topg),m_dx(a_dx),m_default_thickness(a_defaultThickness), m_default_topography(a_defaultTopography), m_set_default_values(a_setDefaultValues), m_isBCsetUp(false),m_verbose(true)
{
  // Construction means nothing to me. It's a lot like Vienna.
  if (!m_thck.isNull())
    {
      for (DataIterator dit( m_thck->disjointBoxLayout());dit.ok();++dit)
        {
          //CH_assert( (*m_thck)[dit].min() >= 0.0);
          FArrayBox& thck = (*m_thck)[dit];
          for (BoxIterator bit(thck.box());bit.ok();++bit)
            {
              if (thck(bit()) < 0.0)
                {
                  thck(bit()) = 0.0;
                }
            }
          
        }
    }
}

LevelDataIBC::LevelDataIBC()
  :m_isBCsetUp(false),m_verbose(true)
{

}


LevelDataIBC::~LevelDataIBC()
{
  
}

void LevelDataIBC::define(const ProblemDomain& a_domain,
			  const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
  //m_dx = a_dx * RealVect::Unit;
}

void LevelDataIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("LevelDataIBC::initialize not implemented");
}

// set boundary fluxes
void LevelDataIBC::primBC(FArrayBox&            a_WGdnv,
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

void  LevelDataIBC::setBdrySlopes(FArrayBox&       a_dW,
				  const FArrayBox& a_W,
				  const int&       a_dir,
				  const Real&      a_time)
{
  // one-sided differences sounds fine with me, so do nothing...
}


void LevelDataIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here
  MayDay::Error("LevelDataIBC::artViscBC not implemented");
}

/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC> LevelDataIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  return m_velBCs; 
}

//set all components of u to zero in ghost regions
void iceDivideBC_LDBC(FArrayBox& a_state,
		      const Box& a_valid,
		      const ProblemDomain& a_domain,
		      Real a_dx,
		      bool a_homogeneous)
{
if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
		  ghostBoxLo &= a_state.box();
		  a_state.setVal(0.0,  ghostBoxLo, 0, a_state.nComp());
		  
                }
	      Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
		  ghostBoxHi &= a_state.box();
		  a_state.setVal(0.0,  ghostBoxHi, 0, a_state.nComp());
                }

            } // end if is not periodic in ith direction
        }
    }

}

void  iceDivideIIBC_LDBC(FArrayBox& a_vel,
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


void LevelDataIBC::setupBCs()
{
  

  m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(iceDivideIIBC_LDBC));
  //m_velBCs = iceDivideBC_LDBC;
  m_isBCsetUp = true;
}

/// set non-periodic ghost cells for surface height z_s. 
void LevelDataIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
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
void LevelDataIBC::setGeometryBCs(LevelSigmaCS& a_coords,
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

void LevelDataIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
					 const RealVect& a_dx,
					 const RealVect& a_domainSize,
					 const Real& a_time, 
					 const LevelSigmaCS* a_crseCoords,
					 const int a_refRatio)
{

  if (m_verbose)
    {
      pout() << " LevelDataIBC::initializeIceGeometry" << endl;
    }

  Real tolerance = 1.0e-6;
  Real refRatio = m_dx[0]/a_dx[0];
  //sanity check
  Real testDx = a_dx[0]*refRatio;
  if (Abs(testDx - m_dx[0])/m_dx[0] > TINY_NORM)
    {
      MayDay::Error("LevelDataIBC::initializeIceGeometry incompatible a_dx and m_dx");
    }

  // if desired, fill with default values first
  if (m_set_default_values)
    {
      LevelData<FArrayBox>& thickness = a_coords.getH();
      LevelData<FArrayBox>& topography = a_coords.getTopography();
      DataIterator dit = thickness.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          thickness[dit].setVal(m_default_thickness);
          topography[dit].setVal(m_default_topography);
        }
    }
  
  if (refRatio > 1.0 + tolerance)
    {
      //fill by interpolation
      if ( (a_crseCoords != NULL) &&
	   (a_refRatio * a_dx[0] <= m_dx[0] * (1.0 + tolerance)))
	{
	  //we have a valid coarse LevelSigmaCS, so we might as well interpolate from that
	  if (m_verbose)
	    {
	      pout() << " ...interpolating data from coarse LevelSigmaCS with refinement ratio = " 
		 << a_refRatio << endl;
	    }

	  a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
	}
      else
	{
	  int nRef = (int)(refRatio + tolerance);
	  if (m_verbose)
	    {
	      pout() << " ...interpolating data from coarse LevelData with refinement ratio = " 
		     << nRef << endl;
	    }
	  //maybe a bit inefficient, but at least we get the same interpolation scheme
	  const DisjointBoxLayout& crseGrids = m_thck->disjointBoxLayout();
	  LevelSigmaCS crseCoords(crseGrids,m_dx, a_coords.ghostVect());
	  for (DataIterator dit(crseGrids); dit.ok(); ++dit)
	    {
	      crseCoords.getH()[dit].copy ( (*m_thck)[dit]);
	      crseCoords.getTopography()[dit].copy ( (*m_topg)[dit]);
	    }
	  crseCoords.recomputeGeometry(NULL,0);
	  a_coords.interpFromCoarse(crseCoords, nRef);
	}
    }
  else 
    {
      FillFromReference(a_coords.getH(),*m_thck,a_dx,m_dx,m_verbose);
      FillFromReference(a_coords.getTopography(),*m_topg,a_dx,m_dx,m_verbose);
      
      if ( refRatio < 1.0 + tolerance)
	{ 
	  int nRef = (int)(1.0/refRatio + tolerance);
	  groundCoarse(a_coords, nRef);
	}
    }
  Real dt = 0.0;
  setGeometryBCs(a_coords, a_coords.grids().physDomain(), a_dx, a_time, dt);  
}

bool LevelDataIBC::regridIceGeometry(LevelSigmaCS& a_coords,
				     const RealVect& a_dx,
				     const RealVect& a_domainSize,
				     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio)
{

    
  Real tolerance = 1.0e-6;
  if (a_dx[0] + tolerance < m_dx[0])
    return false; // if the requested grid is finer than the stored DEM, the best approach is interpolation

  Real refRatio = m_dx[0]/a_dx[0];
  //sanity check
  Real testDx = a_dx[0]*refRatio;

  if (Abs(testDx - m_dx[0])/m_dx[0] > TINY_NORM)
    {
      MayDay::Error("LevelDataIBC::regridIceGeometry incompatibe a_dx and m_dx");
    }

  if (a_crseCoords != NULL &&
      a_refRatio * a_dx[0] <= m_dx[0] * (1.0 + tolerance))
     {
       // in this (common) case, interpolation from a_crseCoords is as good as it gets
       // But we should not get here now.
       MayDay::Error("LevelDataIBC::regridIceGeometry interpolation from a_crseCoords deprecated");
       if (m_verbose)
	 {
	   pout() << " ...interpolating data from coarse LevelSigmaCS with refinement ratio = " 
		  << a_refRatio << endl;
	 }
       a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
     }
   else
     {
       FillFromReference(a_coords.getTopography(),*m_topg,a_dx,m_dx,m_verbose);

       if (a_crseCoords!= NULL)
	 {
	   // fill ghost cells
	   a_coords.interpFromCoarse(*a_crseCoords, a_refRatio, 
				     false, false, false);
	 }

       if ( refRatio < 1.0 + tolerance)
	{
	  int nRef = (int)(1.0/refRatio + tolerance);
	  groundCoarse(a_coords, nRef);
	}

     }

  return true;
}


IceThicknessIBC* LevelDataIBC::new_thicknessIBC()
{
  LevelDataIBC* retval = new LevelDataIBC(m_thck,m_topg,m_dx,
                                          m_default_thickness,
                                          m_default_topography,
                                          m_set_default_values);
  return static_cast<IceThicknessIBC*>(retval);
}


//If there were any grounded cells on the data level but none on the
//computed level, raise the topography by just enough to ground the
//entire cell. 
void LevelDataIBC::groundCoarse(LevelSigmaCS& a_coords, int a_nRef)
{
  
  if (m_verbose)
    {
      pout() << " building data level mask "<< endl;
    }
  
  const DisjointBoxLayout& refGrids = m_thck->disjointBoxLayout();
  LevelData<FArrayBox> refMask(refGrids,1,IntVect::Zero);
  Real rhoi = a_coords.iceDensity();
  Real rhoo = a_coords.waterDensity();
  Real l0 = a_coords.seaLevel();
  Real r = rhoi/rhoo;
  Real f = (1.0-r);
  for (DataIterator dit(refGrids); dit.ok(); ++dit)
    {
      FArrayBox& mask = refMask[dit];
      const FArrayBox& topg = (*m_topg)[dit];
      const FArrayBox& thck = (*m_thck)[dit];
      for (BoxIterator bit(refGrids[dit]);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  mask(iv) = ( (thck(iv)+topg(iv)) >= f * thck(iv))?1.0:0.0; ;
	}
    }
  
  const DisjointBoxLayout& crseGrids = a_coords.grids();
  LevelData<FArrayBox> crseMask(crseGrids,1,IntVect::Zero);
  
  CoarseAverage avg(refGrids,crseGrids,1,a_nRef,IntVect::Zero);
  avg.averageToCoarse(crseMask, refMask);
  
  if (m_verbose)
    {
      pout() << " ...raising topography "<< endl;
    }
  
  Real tol = 0.24;
  for (DataIterator dit(crseGrids); dit.ok(); ++dit)
    {
      const FArrayBox& mask = crseMask[dit];
      FArrayBox& topg = a_coords.getTopography()[dit];
      const FArrayBox& thck = a_coords.getH()[dit];
      for (BoxIterator bit(crseGrids[dit]);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  if (mask(iv) > tol)
	    {
	      topg(iv) = std::max(topg(iv),tiny_thickness + l0 - r*thck(iv)) ;
	    }
	}
    }
}



#include "NamespaceFooter.H"


