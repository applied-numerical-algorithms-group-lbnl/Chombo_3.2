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
// by ISMIP6 ocean forcing 

#include "ISMIP6OceanForcing.H"
#include "ReadLevelData.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "AmrIceBase.H"
#include "IceConstants.H"
#include "computeSum.H"
#include "NamespaceHeader.H"

SurfaceFlux* ISMIP6OceanForcing::new_surfaceFlux() 
{
  // relies on the default copy constructor
  ISMIP6OceanForcing* ptr = new ISMIP6OceanForcing(*this);
  return static_cast<SurfaceFlux*>(ptr);
}


ISMIP6OceanForcing::ISMIP6OceanForcing(ParmParse& a_pp)
{

  std::string file_format;
  a_pp.get("file_format",file_format);

  m_start_year = 1995;
  a_pp.query("start_year",m_start_year);
  m_uniform_source_year = m_start_year - 1;
  
  m_end_year = 2100;
  a_pp.query("end_year",m_end_year);
 
  m_name = "thermal_forcing_0000";
  a_pp.query("name", m_name);

  m_basin_file = "";
  a_pp.query("basin_file",m_basin_file);
  
  m_basin_var_name = "basinNumber";
  a_pp.query("basin_var_name",m_basin_var_name);

  m_n_basin = 16; // imbie2 basins
  a_pp.query("n_basin",m_n_basin);
  
  m_deltaT_file = "";
  a_pp.get("deltaT_file",m_deltaT_file);
  
  m_deltaT_var_name = "deltaT_basin";
  a_pp.query("deltaT_var_name",m_deltaT_var_name);

  m_anomaly = false;
  a_pp.query("anomaly",m_anomaly);

  m_n_layer = 30; // IMSIP6 default
  a_pp.query("n_layer",m_n_layer);
  
  m_dz = 60.0; // ISMIP6 default (metres)
  a_pp.query("dz",m_dz);

  m_basin_mean_min_thickness = 0.0;
  a_pp.query("basin_mean_min_thickness", m_basin_mean_min_thickness);

  
  m_local = true;
  a_pp.query("local", m_local);

  {
    // default value for factor in the formula
    // src = factor * |Tf + delta_T|(Tf + delta_T)
    //ISMIP6 local/non-local defaults for gamma
    Real gamma0 = (m_local)?11075.4506451341:14477.3367602277;
    a_pp.query("gamma0", gamma0);
    // \todo something sensible with constants
    Real rhoi = 918.0;
    Real rhoo = 1028.0;
    Real L = 3.34e+5;
    Real Cp = 3974.0;
    m_factor = - gamma0 * std::pow( (rhoo * Cp) / (rhoi * L), 2);
  }

  // allow user to set factor directly
  a_pp.query("factor", m_factor);
  
  /// populate the time -> file map
  for (int year = m_start_year; year <= m_end_year; year++)
    {
      char* file = new char[file_format.length()+32];
      sprintf(file, file_format.c_str(),year);
      m_year_file.insert(make_pair(year, file));
      delete[] file;
    }


  m_tf_is_source = false; // false: normal ISMIP 6 meaning of Tf
  a_pp.query("tf_is_source", m_tf_is_source);
  
  m_uniform_source = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
}

void ISMIP6OceanForcing::computeSource
(LevelData<FArrayBox>& a_source,
 LevelData<FArrayBox>& a_TFb,
 LevelData<FArrayBox>& a_TFb_basin_mean,
 LevelData<FArrayBox>& a_deltaT,
 Real a_factor)
{

  
  const DisjointBoxLayout& grids = a_source.disjointBoxLayout(); 			      
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  a_source[dit](iv) = a_factor
	    * (a_TFb[dit](iv) + a_deltaT[dit](iv))
	    * Abs(a_TFb_basin_mean[dit](iv) + a_deltaT[dit](iv));
	}
    }
}


void ISMIP6OceanForcing::computeTFb(LevelData<FArrayBox>& a_TFb,
				    const LevelData<FArrayBox>& a_TF,
				    const AmrIceBase& a_amrIce)
{

  CH_TIME("ISMIP6OceanForcing::updateThermalForcing");
  pout() << " ISMIP6OceanForcing::updateThermalForcing() " << std::endl;

  const DisjointBoxLayout& grids = a_TF.disjointBoxLayout();
  a_TFb.define(grids, 1, IntVect::Zero);
  CH_assert(a_TF.nComp() == m_n_layer);
  
  // which mesh level in AmrIce corresponds to TF/Tfb?
  int lev = 0;
  // for the ISMIP6 exercise, lev == 0 (or should be).
  // we will be lazy and only support that simple case,
  // if we want to be more general we need to coarsen/refine...
  Real tol = 1.0e-6;
  CH_assert(Abs(a_amrIce.dx(lev)[0] - m_dx[0]) < tol);
  
  // need s,h to compute the ice shelf draft zb
  // although TF and the geometry are on the same uniform mesh, they
  // are not on the same DisjointBoxLayout
  LevelData<FArrayBox> s(grids,1,IntVect::Zero);
  a_amrIce.geometry(lev)->getSurfaceHeight().copyTo(Interval(0,0), s, Interval(0,0));
  LevelData<FArrayBox> h(grids,1,IntVect::Zero);
  a_amrIce.geometry(lev)->getH().copyTo(Interval(0,0), h , Interval(0,0));


  Real z_min = m_dz * 0.5;
  Real z_max = m_dz * (0.5 + Real(m_n_layer));
    
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real zb = - s[dit](iv) + h[dit](iv); //z downward

	  int km = std::max(0,std::min(m_n_layer-1, int(zb/m_dz - 0.5))); // layer closer to surface
	  int kp = std::max(0,std::min(m_n_layer-1, int(zb/m_dz + 0.5))); // layer closer to bed

	  //if ((iv[0] == 185) && (iv[1] == 350))
	  //  {
	  //    /// near the GL of PIG.
	  //    pout() << " PIG " << iv << " " << km << " " << kp << " " << zb << std::endl;
	  //  }
	  
	  Real TFb = a_TF[dit](iv,0);
	  if (kp > 0){
	    // below midpoint of top layer
	    Real zp = Real(kp + 0.5) * m_dz;
	    if (zb < zp)
	      {
		//between midpoint of the top layer and bottom layer
		Real wm = (zp - zb)/(m_dz);
		TFb = wm * a_TF[dit](iv,km) + (1.0 - wm)*a_TF[dit](iv,kp);
	      }
	    else
	      {
		// below the midpoint of the bottom layer
		TFb = a_TF[dit](iv,kp);
	      }
	  }
	  a_TFb[dit](iv) = TFb;
	} // end loop over cells
    } // end loop over grids
  
}  

/// Read a field TF(x,y,z) and compute TFb(x,y) = TF(x,y,z = b)
/**
   readUniformSource is intended to compute source = - melt 
   on a uniform mesh that covers the entire domain. It should only 
   need to be called once per year (when new TF is needed). It 
   could be called more
   often if the ice shelf geometry changes rapidly. However, 
   because it computes source even in regions outside the ice
   shelf, it does not necessairly need to be called simply 
   because the GL or CF have moved
*/
void ISMIP6OceanForcing::readUniformSource
( RefCountedPtr<LevelData<FArrayBox> >& a_source,
  RealVect& a_dx,
  int a_year,
  const AmrIceBase& a_amrIce)
  
{

  // Read TFb(x,y,z)
  const std::string& file_name = m_year_file[a_year];
  pout() << "    reading "
	 << file_name << ":" << m_name << "..." << std::endl;
  Real dx;
  Vector<std::string> names(1);
  names[0] = m_name;
  RefCountedPtr<LevelData<FArrayBox> > TF(new LevelData<FArrayBox>);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData(1, TF);
  readLevelData(vectData,dx,file_name,names,m_n_layer);
  a_dx = RealVect::Unit * dx;

  /// work out TFb(x,y,z) = TF(x,y,z=b)
  RefCountedPtr<LevelData<FArrayBox> > TFb(new LevelData<FArrayBox>);
  computeTFb(*TFb, *TF, a_amrIce);

  // read deltaT if needed
  if ( m_deltaT.isNull() )
    {
      m_deltaT = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
      vectData[0] = m_deltaT;
      names[0] = m_deltaT_var_name;
      readLevelData(vectData,dx,m_deltaT_file,names,1);
    }
  // m_deltaT is not on the same DBL as TFb...
  LevelData<FArrayBox> deltaT(TF->disjointBoxLayout(), 1, IntVect::Zero);
  m_deltaT->copyTo(Interval(0,0), deltaT, Interval(0,0));

  /// local formula.
  RefCountedPtr<LevelData<FArrayBox> > TFb_basin_mean = TFb;
  if (!m_local)
    {
      //read basin_mask if needed
      if ( m_basin_mask.isNull())
	{
	  m_basin_mask = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
	  vectData[0] = m_basin_mask;
	  names[0] = m_basin_var_name;
	  readLevelData(vectData,dx,m_basin_file,names,m_n_basin); 
	}
      // m_basin_mask is not on the same DBL as TFb...
      LevelData<FArrayBox> basin_mask(TFb->disjointBoxLayout(), m_n_basin, IntVect::Zero);
      m_basin_mask->copyTo(Interval(0, m_n_basin-1), basin_mask, Interval(0,m_n_basin-1 ));
      TFb_basin_mean = RefCountedPtr<LevelData<FArrayBox> >
	(new LevelData<FArrayBox>(TFb->disjointBoxLayout(), 1, IntVect::Zero));
      computeBasinMeans(*TFb_basin_mean, *TFb, basin_mask, a_amrIce);
    }
  
  /// compute the source term
  a_source->define(TF->disjointBoxLayout(), 1, IntVect::Unit);
  if (m_tf_is_source) // option for alanna
    {
      TF->copyTo(Interval(0,0), *a_source, Interval(0,0));
    }
  else
    {
      computeSource(*a_source, *TFb, *TFb_basin_mean, deltaT,  m_factor);
    }
  
}

void ISMIP6OceanForcing::computeBasinMeans(LevelData<FArrayBox>&a_TFb_basin_mean,
					   LevelData<FArrayBox>&a_TFb,
					   LevelData<FArrayBox>&a_basin_mask,
					   const AmrIceBase& a_amrIce)
{
  CH_TIME("ISMIP6OceanForcing::computeBasinMeans");
  pout() << " ISMIP6OceanForcing::computeBasinMeans() " << std::endl;

  // which mesh level in AmrIce corresponds to TFb?
  int lev = 0;
  // for the ISMIP6 exercise, lev == 0 (or should be).
  // we will be lazy and only support that simple case,
  // if we want to be more general we need to coarsen/refine...
  Real tol = 1.0e-6;
  CH_assert(Abs(a_amrIce.dx(lev)[0] - m_dx[0] < tol));

  const DisjointBoxLayout& grids = a_TFb.disjointBoxLayout();
  //set a_TFb_basin_mean = 0 since we are going to accumulate upon it
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      a_TFb_basin_mean[dit].setVal(0.0);
    }
  
  
  // need s,h,b to compute an ice shelf mask
  // although TF and the geometry are on the same uniform mesh, they
  // are not on the same DisjointBoxLayout
  LevelData<FArrayBox> s(grids,1,IntVect::Zero);
  a_amrIce.geometry(lev)->getSurfaceHeight().copyTo(Interval(0,0), s, Interval(0,0));
  LevelData<FArrayBox> h(grids,1,IntVect::Zero);
  a_amrIce.geometry(lev)->getH().copyTo(Interval(0,0), h , Interval(0,0));
  LevelData<FArrayBox> b(grids,1,IntVect::Zero);
  a_amrIce.geometry(lev)->getTopography().copyTo(Interval(0,0), b , Interval(0,0));

  LevelData<FArrayBox> TFb_masked(grids,1,IntVect::Zero);
  LevelData<FArrayBox> basin_shelf_mask(grids,1,IntVect::Zero);
  for (int i_basin = 0; i_basin < a_basin_mask.nComp(); i_basin++)
    {
      
      //re-compute Tfb_masked, shelf_mask for basin i_basin
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  TFb_masked[dit].setVal(0.0);
	  basin_shelf_mask[dit].setVal(0.0);
	  const Box& bx = grids[dit];
	  for (BoxIterator bit(bx); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      if ( (h[dit](iv) > m_basin_mean_min_thickness)
		   && (s[dit](iv) > h[dit](iv) + b[dit](iv) ) )
		{
		  basin_shelf_mask[dit](iv) = a_basin_mask[dit](iv,i_basin);
		  TFb_masked[dit](iv) = a_TFb[dit](iv) * a_basin_mask[dit](iv,i_basin);
		}
	    }
	} // end loop over grids

      // integrals
      Real shelf_area = std::max(TINY_NORM,computeSum(basin_shelf_mask, NULL, 0, m_dx[0], Interval(0,0)));
      Real Tfb_mean = computeSum(TFb_masked, NULL, 0, m_dx[0], Interval(0,0)) / shelf_area;

      // accumulate result
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  FArrayBox t(grids[dit],1);
	  t.copy(basin_shelf_mask[dit]);
	  t *=  Tfb_mean;
	  a_TFb_basin_mean[dit] += t;
	}
    } // end loop over basins
  
}



void ISMIP6OceanForcing::updateUniformSource(Real a_time,
					     const AmrIceBase& a_amrIce)
{

  CH_TIME("ISMIP6OceanForcing::updateThermalForcing");
  pout() << " ISMIP6OceanForcing::updateThermalForcing() " << std::endl;
  m_uniform_source_year = int(a_time);

  readUniformSource(m_uniform_source, m_dx, m_uniform_source_year, a_amrIce); 

  if (m_anomaly)
    {
      // we do need to re-read the initial TF, becasue the
      // ice shelf has changed.
      RealVect dx;
      RefCountedPtr<LevelData<FArrayBox> > m_source_0 (new LevelData<FArrayBox>);
      readUniformSource(m_source_0, dx, m_start_year, a_amrIce); 
      CH_assert(dx == m_dx);

      // \todo I seem to have created different DBL... revisit this
      LevelData<FArrayBox> t(m_uniform_source->disjointBoxLayout(), 1, IntVect::Zero);
      m_source_0->copyTo(Interval(0,0), t, Interval(0,0));
      for (DataIterator dit(m_uniform_source->disjointBoxLayout()); dit.ok(); ++dit)
	{
	  (*m_uniform_source)[dit] -= t[dit];
	}
    }
  
}

void ISMIP6OceanForcing::surfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  CH_TIME("ISMIP6OceanForcing::surfaceThicknessFlux");
  pout() << "ISMIP6OceanForcing::surfaceThicknessFlux " << std::endl;
  
  const Real& time = a_amrIce.time();
  if (time < m_uniform_source_year)
    {
      MayDay::Error("Time out of range for ISMIP6 Ocean forcing");
    }
  //update TF once per year
  //pout() << " update? " << time << " " << m_uniform_source_year << " " << (time >= Real(m_uniform_source_year + 1)) << std::endl;
  if (time >= Real(m_uniform_source_year + 1))
    {
      updateUniformSource(time, a_amrIce);
    }
  
  const RealVect& levelDx = a_amrIce.dx(a_level);
  FillFromReference(a_flux, *m_uniform_source, levelDx ,m_dx, true);
  
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
