#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "SurfaceFlux.H"
#include "ComplexSurfaceFlux.H"
#include "LevelDataSurfaceFlux.H"
#include "GroundingLineLocalizedFlux.H"
#include "HotspotFlux.H"
#include <map>
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "IceConstants.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "BisiclesF_F.H"
#include "ParmParse.H"
#include "AmrIceBase.H"
#include "FortranInterfaceIBC.H"
#include "FillFromReference.H"

#include "NamespaceHeader.H"


// constructor
ProductSurfaceFlux::ProductSurfaceFlux  (SurfaceFlux* a_flux1, SurfaceFlux* a_flux2)
{
  m_flux1 = a_flux1;
  m_flux2 = a_flux2;
}


/// destructor
ProductSurfaceFlux::~ProductSurfaceFlux()
{
  // I think we should be deleting m_flux1 and m_flux2 here
}

/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* 
ProductSurfaceFlux::new_surfaceFlux()
{
  SurfaceFlux* f1 = m_flux1->new_surfaceFlux();
  SurfaceFlux* f2 = m_flux2->new_surfaceFlux();
  return static_cast<SurfaceFlux*>(new ProductSurfaceFlux(f1,f2));
}

/// define source term for thickness evolution and place it in flux
/** dt is included in case one needs integrals or averages over a
    timestep. flux should be defined in meters/second in the current 
    implementation. 
*/
void 
ProductSurfaceFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					 const AmrIceBase& a_amrIce, 
					 int a_level, Real a_dt)
{
  LevelData<FArrayBox> f2(a_flux.getBoxes(), a_flux.nComp(), a_flux.ghostVect());
  // compute flux1, put in a_flux, compute flux2 in f2, then multiply
  m_flux1->surfaceThicknessFlux(a_flux, a_amrIce, a_level, a_dt);
  m_flux2->surfaceThicknessFlux(f2, a_amrIce, a_level, a_dt);
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_flux[dit].mult(f2[dit]);
    }
}




/// constructor
MaskedFlux::MaskedFlux(SurfaceFlux* a_groundedIceFlux, SurfaceFlux* a_floatingIceFlux,
		       SurfaceFlux* a_openSeaFlux, SurfaceFlux* a_openLandFlux)
  :m_groundedIceFlux(a_groundedIceFlux),m_floatingIceFlux(a_floatingIceFlux),
   m_openSeaFlux(a_openSeaFlux),m_openLandFlux(a_openLandFlux)
{
  CH_assert(a_groundedIceFlux);
  CH_assert(a_floatingIceFlux);
  CH_assert(a_openSeaFlux);
  CH_assert(a_openLandFlux);
}
/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* MaskedFlux::new_surfaceFlux()
{
  SurfaceFlux* f = m_floatingIceFlux->new_surfaceFlux();
  SurfaceFlux* g = m_groundedIceFlux->new_surfaceFlux();
  SurfaceFlux* s = m_openSeaFlux->new_surfaceFlux();
  SurfaceFlux* l = m_openLandFlux->new_surfaceFlux();
  return static_cast<SurfaceFlux*>(new MaskedFlux(g,f,s,l));
}

void MaskedFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				      const AmrIceBase& a_amrIce, 
				      int a_level, Real a_dt)
{

  //somewhat ineffcient, because we compute all fluxes everywhere. 
  //At some point, come back and only compute (say) grounded ice flux
  //in boxes where at least some of the ice is grounded.

  //first, grounded ice values
  m_groundedIceFlux->surfaceThicknessFlux(a_flux,a_amrIce,a_level,a_dt);

  //floating,open sea,open land ice values
  std::map<int,SurfaceFlux*> mask_flux;
  mask_flux[FLOATINGMASKVAL] = m_floatingIceFlux;
  mask_flux[OPENSEAMASKVAL] =  m_openSeaFlux ;
  mask_flux[OPENLANDMASKVAL] = m_openLandFlux;
  LevelData<FArrayBox> tmpFlux;
  tmpFlux.define(a_flux);
  for (std::map<int,SurfaceFlux*>::iterator i = mask_flux.begin(); i != mask_flux.end(); ++i)
    {
      i->second->surfaceThicknessFlux(tmpFlux, a_amrIce,a_level,a_dt);
      for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask =  a_amrIce.geometry(a_level)->getFloatingMask()[dit];
      
	  Box box = mask.box();
	  box &= a_flux[dit].box();

	  int m = i->first;
	  FORT_MASKEDREPLACE(CHF_FRA1(a_flux[dit],0),
			     CHF_CONST_FRA1(tmpFlux[dit],0),
			     CHF_CONST_FIA1(mask,0),
			     CHF_CONST_INT(m),
			     CHF_BOX(box));
	}
    }

}

SurfaceFlux* AxbyFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new AxbyFlux(m_a, m_x,m_b, m_y) );
}

AxbyFlux::AxbyFlux(const Real& a_a, SurfaceFlux* a_x, 
		   const Real& a_b, SurfaceFlux* a_y)
{

  m_a = a_a;
  m_b = a_b;
  
  CH_assert(a_x != NULL);
  CH_assert(a_y != NULL);
  m_x = a_x->new_surfaceFlux();
  m_y = a_y->new_surfaceFlux();
  CH_assert(m_x != NULL);
  CH_assert(m_y != NULL);

}

AxbyFlux::~AxbyFlux()
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

void AxbyFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  const AmrIceBase& a_amrIce, 
					  int a_level, Real a_dt)
{

  LevelData<FArrayBox> y_flux(a_flux.disjointBoxLayout(),1,a_flux.ghostVect());
  m_x->surfaceThicknessFlux(a_flux, a_amrIce, a_level,a_dt );
  m_y->surfaceThicknessFlux(y_flux, a_amrIce, a_level,a_dt );
  for (DataIterator dit(a_flux.disjointBoxLayout()); dit.ok(); ++dit)
    {
      a_flux[dit].axby(a_flux[dit],y_flux[dit],m_a,m_b);
    }
  

}


/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* CompositeFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new CompositeFlux(m_fluxes));
}

CompositeFlux::CompositeFlux(const Vector<SurfaceFlux*>& a_fluxes)
{
  m_fluxes.resize(a_fluxes.size());
  for (int i =0; i < a_fluxes.size(); i++)
    {
      CH_assert(a_fluxes[i] != NULL);
      m_fluxes[i] =  a_fluxes[i]->new_surfaceFlux();
    }
}

CompositeFlux::~CompositeFlux()
{
  for (int i =0; i < m_fluxes.size(); i++)
    {
      if (m_fluxes[i] != NULL)
	{
	  delete m_fluxes[i];
	  m_fluxes[i] = NULL;
	}
    }
}

void CompositeFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  const AmrIceBase& a_amrIce, 
					  int a_level, Real a_dt)
{
  m_fluxes[0]->surfaceThicknessFlux(a_flux, a_amrIce, a_level,a_dt );

  // this is hardly effcient... but it is convenient
  LevelData<FArrayBox> tmpFlux(a_flux.disjointBoxLayout(),1,a_flux.ghostVect());
  for (int i = 1; i <  m_fluxes.size(); i++)
    {
      m_fluxes[i]->surfaceThicknessFlux(tmpFlux, a_amrIce, a_level,a_dt );
      for (DataIterator dit(a_flux.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  a_flux[dit] += tmpFlux[dit];
	}
    }

}



SurfaceFlux* BoxBoundedFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>( new BoxBoundedFlux(m_lo, m_hi,m_startTime,m_endTime,m_fluxPtr));
}


void BoxBoundedFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  const AmrIceBase& a_amrIce, 
					  int a_level, Real a_dt)
{

  Real time = a_amrIce.time();
  

  for (DataIterator dit(a_flux.disjointBoxLayout()); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(0.0);
    }

  if (time >= m_startTime && time < m_endTime)
    {
      // this is hardly efficient... but it is convenient
      LevelData<FArrayBox> tmpFlux(a_flux.disjointBoxLayout(),1,a_flux.ghostVect());
      m_fluxPtr->surfaceThicknessFlux(tmpFlux, a_amrIce, a_level,a_dt);
      const RealVect& dx = a_amrIce.dx(a_level);
      
      IntVect ilo,ihi;
      for (int dir =0; dir < SpaceDim; dir++)
	{
	  ilo[dir] = int(m_lo[dir]/dx[dir] - 0.5);
	  ihi[dir] = int(m_hi[dir]/dx[dir] - 0.5);
	}
      
    
      for (DataIterator dit(a_flux.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  const Box& b = a_flux[dit].box();
	  if (b.intersects(Box(ilo,ihi)))
	    { 
	      Box sub(max(b.smallEnd(),ilo), min(b.bigEnd(),ihi));
	      a_flux[dit].plus(tmpFlux[dit],sub,0,0,1);
	    }
	}
    }

}

// --------------------------------------------------------------
// fortran interface surface flux
// --------------------------------------------------------------

/// class which takes an input fortran array 
/** averages or interpolates as necessary to fill the flux
 */

/// constructor
fortranInterfaceFlux::fortranInterfaceFlux()
  : m_isValSet(false)
{
}

/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* 
fortranInterfaceFlux::new_surfaceFlux()
{
  if (m_verbose)
    {
      pout() << "in fortranInterfaceFlux::new_surfaceFlux" << endl;
    }

  fortranInterfaceFlux* newPtr = new fortranInterfaceFlux;

  newPtr->m_grids = m_grids;
  newPtr->m_gridsSet = m_gridsSet;
    // keep these as aliases, if they're actually defined
  if (!m_inputFlux.box().isEmpty())
    {
      newPtr->m_inputFlux.define(m_inputFlux.box(), 
                                 m_inputFlux.nComp(),
                                 m_inputFlux.dataPtr());
    }

  if (!m_ccInputFlux.box().isEmpty())
    {
      newPtr->m_ccInputFlux.define(m_ccInputFlux.box(), 
                                   m_ccInputFlux.nComp(),
                                   m_ccInputFlux.dataPtr());
    }      
  
  newPtr->m_inputFluxLDF = m_inputFluxLDF;
  
  newPtr->m_fluxGhost = m_fluxGhost;
  newPtr->m_inputFluxDx = m_inputFluxDx;
  newPtr->m_grids = m_grids;
  newPtr->m_gridsSet = m_gridsSet;
  
  newPtr->m_verbose = m_verbose;

  newPtr->m_isValSet = m_isValSet;
  return static_cast<SurfaceFlux*>(newPtr);  
}

/// define source term for thickness evolution and place it in flux
/** dt is included in case one needs integrals or averages over a
    timestep. flux should be defined in meters/second in the current 
    implementation. 
*/
void 
fortranInterfaceFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					   const AmrIceBase& a_amrIce, 
					   int a_level, Real a_dt)
{
  CH_assert(m_isValSet);

  // this looks a lot like the code in FortranInterfaceIBC

  DisjointBoxLayout levelGrids = m_grids;
  RealVect dx = a_amrIce.dx(a_level);

  FillFromReference(a_flux,
                    *m_inputFluxLDF,
                    dx, m_inputFluxDx,
                    m_verbose);
#if 0
  // refinement ratio for flux
  Real refRatio = m_inputFluxDx[0]/dx[0];
 
  Real tolerance = 1.0e-6;

  if (refRatio > 1 + tolerance)
    {
      // importFlux coarser than what we want, have to interpolate
      //int nRef = (int)(refRatio + tolerance);

    }
  else if (refRatio < 1 + tolerance)
    {
      // importFlux finer than what we want, have to average
      //int nRef = (int)(1.0/refRatio + tolerance);
    }
  else
    {
      // same size, just copy  
      m_inputFluxLDF->copyTo(a_flux);
    }
  
#endif

}

/// set fortran array-valued surface flux
void
fortranInterfaceFlux::setFluxVal(Real* a_data_ptr,
                                 const int* a_dimInfo,
                                 const int* a_boxlo, const int* a_boxhi, 
                                 const Real* a_dew, const Real* a_dns,
                                 const IntVect& a_offset,
                                 const IntVect& a_nGhost,
                                 const ProblemDomain& a_domain,
                                 const bool a_nodal)

{

  m_fluxGhost = a_nGhost;
  m_nodalFlux = a_nodal;
  m_domain = a_domain;

  // dimInfo is (SPACEDIM, nz, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering
  //cout << "a_dimonfo" << a_dimInfo[0] << a_dimInfo[1] << endl;  

  if (m_verbose)
    {
      pout() << "In FortranInterfaceIBC::setFlux:" << endl;
      pout() << " -- entering setFAB..." << endl;
    }
  
  FortranInterfaceIBC::setFAB(a_data_ptr, a_dimInfo,a_boxlo, a_boxhi,
                              a_dew,a_dns,a_offset,a_nGhost,
                              m_inputFlux, m_ccInputFlux, a_nodal);

  if (m_verbose)
    {
      pout() << "... done" << endl;
    }

  // if we haven't already set the grids, do it now
  if (!gridsSet())
    {
      if (m_verbose) 
        {
          pout() << " -- entering setGrids" << endl;
        }
      Box gridBox(m_ccInputFlux.box());
      gridBox.grow(-a_nGhost);
      FortranInterfaceIBC::setGrids(m_grids, gridBox, m_domain, m_verbose);
      m_gridsSet = true;
      if (m_verbose)
        {
          pout() << " -- out of setGrids" << endl;
        }
    }
  


  m_inputFluxDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // now define LevelData and copy from FAB->LevelData 
  // (at some point will likely change this to be an  aliased 
  // constructor for the LevelData, but this should be fine for now....
  
  // if nodal, we'd like at least one ghost cell for the LDF
  // (since we'll eventually have to average back to nodes)
  IntVect LDFghost = m_fluxGhost;
  if (a_nodal && (LDFghost[0] == 0))
    {
      LDFghost += IntVect::Unit;
    }
      
  RefCountedPtr<LevelData<FArrayBox> > localLDFPtr(new LevelData<FArrayBox>(m_grids, 1, LDFghost) );

  m_inputFluxLDF = localLDFPtr;
  // fundamental assumption that there is no more than one box/ processor 
  // don't do anything if there is no data on this processor
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box copyBox = (*m_inputFluxLDF)[dit].box();
      copyBox &= m_inputFlux.box();
      (*m_inputFluxLDF)[dit].copy(m_inputFlux, copyBox);
      
    } // end DataIterator loop

  m_isValSet = true;
}

PiecewiseLinearFlux::PiecewiseLinearFlux(const Vector<Real>& a_abscissae, 
					 const Vector<Real>& a_ordinates, 
					 Real a_minWaterDepth)
  :m_abscissae(a_abscissae),m_ordinates(a_ordinates),
   m_minWaterDepth(a_minWaterDepth)
{
  CH_assert(m_abscissae.size() == m_ordinates.size());
}


SurfaceFlux* PiecewiseLinearFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new PiecewiseLinearFlux(m_abscissae,m_ordinates,m_minWaterDepth));
}

void PiecewiseLinearFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					       const AmrIceBase& a_amrIce, 
					       int a_level, Real a_dt)
{
  Vector<Real> dx(m_abscissae.size());
  Vector<Real> db(m_abscissae.size());

  const LevelData<FArrayBox>& levelH = a_amrIce.geometry(a_level)->getH();
  const LevelData<FArrayBox>& levelS = a_amrIce.geometry(a_level)->getSurfaceHeight();
  const LevelData<FArrayBox>& levelR = a_amrIce.geometry(a_level)->getTopography();
  for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit)
    {

      FORT_PWLFILL(CHF_FRA1(a_flux[dit],0),
		   CHF_CONST_FRA1(levelH[dit],0),
		   CHF_CONST_VR(m_abscissae),
		   CHF_CONST_VR(m_ordinates),
		   CHF_VR(dx),CHF_VR(db),
		   CHF_BOX(a_flux[dit].box()));
       
      if (m_minWaterDepth > 0.0)
	{
	  
	  FArrayBox D(a_flux[dit].box(),1);
	  FORT_WATERDEPTH(CHF_FRA1(D,0),
			  CHF_CONST_FRA1(levelH[dit],0),
			  CHF_CONST_FRA1(levelS[dit],0),
			  CHF_CONST_FRA1(levelR[dit],0),
			  CHF_BOX(a_flux[dit].box()));
  
	  
	  FORT_ZEROIFLESS(CHF_FRA1(a_flux[dit],0),
			  CHF_CONST_FRA1(D,0),
			  CHF_CONST_REAL(m_minWaterDepth),
			  CHF_BOX(a_flux[dit].box()));

	}

    }
}





/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* NormalizedFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new NormalizedFlux(m_direction, m_amplitude));
}

NormalizedFlux::NormalizedFlux(SurfaceFlux* a_direction, const Real& a_amplitude)
{
  m_direction = a_direction->new_surfaceFlux();
  m_amplitude = a_amplitude;
}

NormalizedFlux::~NormalizedFlux()
{
  if (m_direction != NULL)
    {
      delete m_direction; m_direction = NULL;
    }
}

void NormalizedFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					  const AmrIceBase& a_amrIce, 
					  int a_level, Real a_dt)
{
  // Need to compute the norm over *all* the levels, which will
  // mean that under typical circumstances the flux and its norm is computed n_level times.
  // alternatives would involve risky assumptions or interface redesign

  Vector<LevelData<FArrayBox>*>  flux;
  for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
    {
      flux.push_back(new LevelData<FArrayBox>(a_amrIce.grids(lev), a_flux.nComp(), a_flux.ghostVect()));
      m_direction->surfaceThicknessFlux(*flux[lev], a_amrIce, lev, a_dt );
    }
  Real norm = computeNorm(flux, a_amrIce.refRatios(), a_amrIce.dx(0)[0],  Interval(0,0), 1);
  for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
    {
      if (flux[lev] != NULL)
	{
	  delete flux[lev]; flux[lev] = NULL;
	}
    }			  
  m_direction->surfaceThicknessFlux(a_flux, a_amrIce, a_level, a_dt );
  Real factor = m_amplitude/norm;
  for (DataIterator dit (a_flux.disjointBoxLayout()); dit.ok(); ++dit)
    {
      a_flux[dit] *= factor;
    }
}


/// factory method
/** return a pointer to a new SurfaceFlux object
 */
SurfaceFlux* TargetThicknessFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new TargetThicknessFlux(m_target, m_timescale));
}

TargetThicknessFlux::TargetThicknessFlux(SurfaceFlux* a_target, const Real& a_timescale)
{
  m_target = a_target->new_surfaceFlux();
  m_timescale = a_timescale;
}

TargetThicknessFlux::~TargetThicknessFlux()
{
  if (m_target != NULL)
    {
      delete m_target; m_target = NULL;
    }
}

void TargetThicknessFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
					       const AmrIceBase& a_amrIce, 
					       int a_level, Real a_dt)
{
  
  m_target->surfaceThicknessFlux(a_flux, a_amrIce, a_level, a_dt );
  Real factor = std::max(1.0e-10,1.0/m_timescale);
  for (DataIterator dit (a_flux.disjointBoxLayout()); dit.ok(); ++dit)
    {
      a_flux[dit] -= a_amrIce.geometry(a_level)->getH()[dit];
      a_flux[dit] *= factor;
    }
}




SurfaceFlux*  FloatingDivUHLocalizedFlux::new_surfaceFlux()
{
  return static_cast<SurfaceFlux*>(new  FloatingDivUHLocalizedFlux(m_flux, m_mesh_spacing));
}

FloatingDivUHLocalizedFlux::FloatingDivUHLocalizedFlux(SurfaceFlux* a_flux, const Real& a_mesh_spacing)
{
  m_flux = a_flux->new_surfaceFlux();
  m_mesh_spacing = a_mesh_spacing;
}

FloatingDivUHLocalizedFlux::~FloatingDivUHLocalizedFlux()
{
  if (m_flux != NULL)
    {
      delete m_flux; m_flux = NULL;
    }
}

void FloatingDivUHLocalizedFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux_data,
					  const AmrIceBase& a_amrIce, 
					  int a_level, Real a_dt)
{

  // I can't see any reason to use this class when m_mesh_spacing < a_amrIce.dx(0)[0],
  // and doing so is ineffcient given the implementation, so for now, disallow
  CH_assert( m_mesh_spacing >=  a_amrIce.dx(0)[0]);


  //evaluate the input flux on level 0
  DisjointBoxLayout baseDBL =  a_amrIce.grids(0);
  Real baseDx = a_amrIce.dx(0)[0];
  LevelData<FArrayBox> baseFlux(baseDBL,a_flux_data.nComp() ,a_flux_data.ghostVect());
  m_flux->evaluate(baseFlux,a_amrIce, 0 , a_dt);

  // construct uniform mesh with spacing m_mesh_spacing
  // (first check for compatible grids)
  Real ratio = (m_mesh_spacing / baseDx);
  CH_assert( Abs(Real(int(ratio)) - ratio) < TINY_NORM);
  CH_assert( int(ratio)% 2 == 0);
  int nRef = int(ratio);
  DisjointBoxLayout crseDBL;
  coarsen(crseDBL,baseDBL,nRef);
  LevelData<FArrayBox> crseFlux(crseDBL,a_flux_data.nComp() ,a_flux_data.ghostVect());
  
  // coarsen input flux to crseDBL, ie. integrate over the coarse cells  
  CoarseAverage av(baseDBL, crseDBL, 1, nRef, baseFlux.ghostVect());
  av.averageToCoarse(crseFlux, baseFlux);

  // compute w and w div(uh) 
  Vector<LevelData<FArrayBox>* > w;
  Vector<LevelData<FArrayBox>* > divuh;
  Vector<LevelData<FArrayBox>* > wdivuh;
  Vector<RealVect> amrDx;
  for (int lev = 0; lev <=  a_amrIce.finestLevel(); lev++)
    {
      Real dxl = a_amrIce.dx(lev)[0];
      amrDx.push_back(a_amrIce.dx(lev));
      const DisjointBoxLayout& dbl =  a_amrIce.grids(lev);
      w.push_back(new LevelData<FArrayBox>(dbl, 1, IntVect::Unit));
      divuh.push_back(new LevelData<FArrayBox>(dbl, 1, IntVect::Unit));
      wdivuh.push_back(new LevelData<FArrayBox>(dbl, 1, IntVect::Unit));
      const LevelData<FArrayBox>& h = a_amrIce.geometry(lev)->getH();
      const LevelData<FArrayBox>& u = *a_amrIce.velocity(lev);
      const LevelData<BaseFab<int> >& mask = a_amrIce.geometry(lev)->getFloatingMask();
      for (DataIterator dit(dbl); dit.ok(); ++dit)
	{ 
	  for (BoxIterator bit(dbl[dit]); bit.ok(); ++bit)
	    {
	      // hoping to get away with this scheme for div(uh)...
	      const IntVect& iv = bit();
	      (*divuh[lev])[dit](iv) = 0.0;
	      for (int dir = 0; dir < SpaceDim; ++dir)
		{
		  (*divuh[lev])[dit](iv) += 1.0 / dxl *
		    (  h[dit](iv+BASISV(dir))* u[dit](iv+BASISV(dir)) -
		       h[dit](iv-BASISV(dir))* u[dit](iv-BASISV(dir)) );
		}
	      (*divuh[lev])[dit](iv) = std::max(-100.0,std::min(0.0,  (*divuh[lev])[dit](iv) ));
	      (*w[lev])[dit](iv) = (mask[dit](iv) == FLOATINGMASKVAL)?1.0:0.0;

	      
	      (*wdivuh[lev])[dit](iv) = (*w[lev])[dit](iv) * (*divuh[lev])[dit](iv);
	    }
	}
    }

  // coarsen w and w div(uv) ie. integrate over the coarse cells
  LevelData<FArrayBox> crseW(crseDBL,a_flux_data.nComp() ,a_flux_data.ghostVect());
  RealVect crseDx = RealVect::Unit *  m_mesh_spacing;
  flattenCellData( crseW, crseDx , w, amrDx , false);
  LevelData<FArrayBox> crseWdivuh(crseDBL,a_flux_data.nComp() ,a_flux_data.ghostVect());
  flattenCellData( crseWdivuh, crseDx , wdivuh, amrDx, false);
  
  // compute Q on the coarse grid
  LevelData<FArrayBox> crseQ(crseDBL,a_flux_data.nComp() ,a_flux_data.ghostVect());
  for (DataIterator dit(crseDBL); dit.ok(); ++dit)
    {
      for (BoxIterator bit(crseDBL[dit]); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (crseW[dit](iv) > 1.0e-10)
	    {
	      crseQ[dit](iv) = - (crseWdivuh[dit](iv) - crseFlux[dit](iv))/(crseW[dit](iv));
	    }
	  else
	    {
	      crseQ[dit](iv) = 0.0;
	    }
	}
    }
  // compute Q on the base grid

  // compute Q on the current level

  //FineInterp interpolator(dbl, 1, nRef, dbl.physDomain());
  LevelData<FArrayBox> Q(a_amrIce.grids(a_level), a_flux_data.nComp() ,a_flux_data.ghostVect());
  FillFromReference(Q,crseQ, amrDx[a_level], crseDx,false);
  
  // finally, f = w * (Q + div(uh))
  
  for (DataIterator dit(a_amrIce.grids(a_level)); dit.ok(); ++dit)
    {
      FArrayBox t(a_amrIce.grids(a_level)[dit],1);
      t.copy(Q[dit]);
      t += (*divuh[a_level])[dit];
      t *= (*w[a_level])[dit];
      a_flux_data[dit].setVal(0.0);
      a_flux_data[dit].copy(t);
    }
  
  //clean up.
  for (int lev = 0; lev <=  a_amrIce.finestLevel(); lev++)
    {
      delete w[lev]; w[lev] = NULL;
      delete wdivuh[lev]; wdivuh[lev] = NULL;
      delete divuh[lev]; wdivuh[lev] = NULL;
    }
}



#ifdef HAVE_PYTHON
#include "signal.h"
#endif
#include "NamespaceFooter.H"
