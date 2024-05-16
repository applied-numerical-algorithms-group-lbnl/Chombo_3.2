#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DomainDiagnosticData.H"
#include "AmrIce.H"
#include "IceConstants.H"
#include "computeSum.H"
#include "NamespaceHeader.H"

/// set up a struct of cf names to storage
void DomainDiagnosticData::setCFdata()
{

  cfDiagnostic cf_info;

  cf_info.short_name = CFIO_DIAGNOSTIC_TIME_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_TIME_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_TIME_NAME;
  cf_info.units = "yr";
  cf_info.data = &m_time;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_ICE_VOLUME_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_ICE_VOLUME_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_ICE_VOLUME_LONG_NAME;
  cf_info.units = "kg";
  cf_info.data = &m_ice_volume;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_ICE_VAF_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_ICE_VAF_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_ICE_VAF_LONG_NAME;
  cf_info.units = "kg";
  cf_info.data = &m_ice_vaf;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_LONG_NAME;
  cf_info.units = "m^2";
  cf_info.data = &m_ice_grounded_area;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_LONG_NAME;
  cf_info.units = "m^2";
  cf_info.data = &m_ice_floating_area;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_total_smb;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_total_bmb;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_floating_total_bmb;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_CALVING_FLUX_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_CALVING_FLUX_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_CALVING_FLUX_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_total_calving_flux;
  m_cf_stuff.push_back(cf_info);

  cf_info.short_name = CFIO_DIAGNOSTIC_ICE_FRONT_CALVING_AND_MELTING_SHORT_NAME;
  cf_info.cf_name = CFIO_DIAGNOSTIC_ICE_FRONT_CALVING_AND_MELTING_CF_NAME;
  cf_info.long_name = CFIO_DIAGNOSTIC_ICE_FRONT_CALVING_AND_MELTING_LONG_NAME;
  cf_info.units = "kg yr^-1";
  cf_info.data = &m_ice_total_calving_and_ice_front_melting_flux;
  m_cf_stuff.push_back(cf_info);

}

/// diagnostic functions -- integrates variables over domain
Real
DomainDiagnosticData::computeTotalIce
(const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const Vector<int>& a_refRatio, Real a_crseDx,
 int a_finestLevel) const
{
  Vector<LevelData<FArrayBox>* > thickness(a_finestLevel+1, NULL);
  for (int lev=0; lev<=a_finestLevel; lev++)
    {
      const LevelSigmaCS& levelCoords = *a_coordSys[lev];
      // need a const_cast to make things all line up right
      // (but still essentially const)
      thickness[lev] = const_cast<LevelData<FArrayBox>* >(&levelCoords.getH());
    }

  Interval thicknessInt(0,0);
  Real totalIce = computeSum(thickness, a_refRatio, a_crseDx, thicknessInt, 0);

  return totalIce;

}

Real
DomainDiagnosticData::computeVolumeAboveFlotation
  (const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   const Vector<int>& a_refRatio, Real a_crseDx,
   int a_finestLevel) const
{

  //Compute the total thickness above flotation
  Vector<LevelData<FArrayBox>* > thk(a_finestLevel+1, NULL);
  for (int lev=0; lev <= a_finestLevel ; lev++)
    {
      const LevelSigmaCS& levelCoords = *a_coordSys[lev];
      // need a const_cast to make things all line up right
      // (but still essentially const)
      thk[lev] = const_cast<LevelData<FArrayBox>*>(&levelCoords.getThicknessOverFlotation());
    }
  Real VAF = computeSum(thk, a_refRatio, a_crseDx, Interval(0,0), 0);
  return VAF;
}

Real 
DomainDiagnosticData::computeTotalGroundedIce
  (const Vector<DisjointBoxLayout>& a_grids,
   const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
   const Vector<int>& a_refRatio, Real a_crseDx,
   int a_finestLevel) const
{
  
  Real totalGroundedIce = 0;

  Vector<LevelData<FArrayBox>* > vectGroundedThickness(a_finestLevel+1, NULL);

  for (int lev=0; lev<=a_finestLevel; lev++)
    {
      const LevelData<FArrayBox>& levelThickness = a_coordSys[lev]->getH();
      // temporary with only ungrounded ice
      vectGroundedThickness[lev] = new LevelData<FArrayBox>(a_grids[lev],1,
							    IntVect::Zero);

      LevelData<FArrayBox>& levelGroundedThickness = *vectGroundedThickness[lev];
      // now copy thickness to       
      levelThickness.copyTo(levelGroundedThickness);

      const LevelData<BaseFab<int> >& levelMask = a_coordSys[lev]->getFloatingMask();
      // now loop through and set to zero where we don't have grounded ice.
      // do this the slow way, for now
      DataIterator dit=levelGroundedThickness.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  const BaseFab<int>& thisMask = levelMask[dit];
	  FArrayBox& thisThick = levelGroundedThickness[dit];
	  BoxIterator bit(thisThick.box());
	  for (bit.begin(); bit.ok(); ++bit)
	    {
	      IntVect iv = bit();
	      if (thisMask(iv,0) != GROUNDEDMASKVAL)
		{
		  thisThick(iv,0) = 0.0;
		}
	    }
	}
    

    }

  // now compute sum
  Interval thicknessInt(0,0);
  totalGroundedIce = computeSum(vectGroundedThickness, a_refRatio,
				a_crseDx, thicknessInt, 0);

  
  // clean up temp storage
  for (int lev=0; lev<vectGroundedThickness.size(); lev++)
    {
      if (vectGroundedThickness[lev] != NULL)
	{
	  delete vectGroundedThickness[lev];
	  vectGroundedThickness[lev] = NULL;
	}
    }

  return totalGroundedIce;

}

void 
DomainDiagnosticData::computeAreaFraction(LevelData<FArrayBox>& a_area,
					  const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
					  int a_maskVal,
					  int a_level,
					  int a_finestLevel) const
{
  CH_assert(a_level <= a_finestLevel)
    {
      const LevelData<BaseFab<int> >& levelMask = a_coordSys[a_level]->getFloatingMask();
      // set a_area = 1.0 where we have grounded ice
      DataIterator dit=a_area.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  const BaseFab<int>& thisMask = levelMask[dit];
	  FArrayBox& thisIce = a_area[dit];
	  thisIce.setVal(0.0);
	  BoxIterator bit(thisIce.box());
	  for (bit.begin(); bit.ok(); ++bit)
	    {
	      IntVect iv = bit();
	      if (thisMask(iv,0) == a_maskVal)
		{
		  thisIce(iv,0) = 1.0;
		}
	    }
	}
    }
}

Real 
DomainDiagnosticData::computeArea(int a_maskVal,
				  const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
				  const Vector<DisjointBoxLayout>& a_grids,
				  const Vector<int>& a_refRatio, Real a_crseDx,
				  int a_finestLevel) const
{
  Real sum_area = 0.0;

  //compute per-cell area fraction
  Vector<LevelData<FArrayBox>* > area_fraction(a_finestLevel+1, NULL);
  for (int lev=0; lev<=a_finestLevel; lev++)
    {
      area_fraction[lev] = new LevelData<FArrayBox>(a_grids[lev],1,IntVect::Zero);
      computeAreaFraction(*area_fraction[lev], a_coordSys, a_maskVal,lev, a_finestLevel );
    }

  // integrate over domain
  sum_area = computeSum(area_fraction, a_refRatio, a_crseDx, Interval(0,0), 0);

  // clean up temp storage
  for (int lev=0; lev<area_fraction.size(); lev++)
    {
      if (area_fraction[lev] != NULL)
	{
	  delete area_fraction[lev];
	  area_fraction[lev] = NULL;
	}
    }
  return sum_area;
}

Real 
DomainDiagnosticData::computeGroundedArea(const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
					  const Vector<DisjointBoxLayout>& a_grids,
					  const Vector<int>& a_refRatio, Real a_crseDx,
					  int a_finestLevel) const
{
  return computeArea(GROUNDEDMASKVAL, a_coordSys, a_grids, a_refRatio, a_crseDx, 
		     a_finestLevel);
}

Real 
DomainDiagnosticData::computeFloatingArea(const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
					  const Vector<DisjointBoxLayout>& a_grids,
					  const Vector<int>& a_refRatio, Real a_crseDx,
					  int a_finestLevel) const
{
  return computeArea(FLOATINGMASKVAL, a_coordSys, a_grids, a_refRatio, a_crseDx,
		     a_finestLevel);
}

Real 
DomainDiagnosticData::computeFluxOverIce(const Vector<LevelData<FArrayBox>* > a_flux,
					 const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
					 const Vector<DisjointBoxLayout>& a_grids,
					 const Vector<int>& a_refRatio, Real a_crseDx,
					 int a_finestLevel) const
{

  //compute sum of a flux component over ice
  //construct fluxOverIce
  Vector<LevelData<FArrayBox>* > fluxOverIce ( a_finestLevel+1, NULL);
  for (int lev = 0; lev <= a_finestLevel ; lev++)
    {
      fluxOverIce[lev] = new
	LevelData<FArrayBox>(a_grids[lev],1, IntVect::Zero);
      const LevelData<FArrayBox>& thk = a_coordSys[lev]->getH();
      //const LevelData<FArrayBox>* flux = a_flux[lev];
       
      for (DataIterator dit(a_grids[lev]); dit.ok(); ++dit)
	{
	  const Box& box =  a_grids[lev][dit];
	  const FArrayBox& source = (*a_flux[lev])[dit];
	  const FArrayBox& dit_thck = thk[dit];
	  FArrayBox& dit_fluxOverIce = (*fluxOverIce[lev])[dit];
	     
	  for (BoxIterator bit(box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      // set fluxOverIce to source if thck > 0
	      if (dit_thck(iv) < 1e-10)
		{
		  dit_fluxOverIce(iv) = 0.0;
		}
	      else
		{
		  dit_fluxOverIce(iv) = source(iv);
		}
	    }
	    
	}
    }
  // compute sum
  Real tot_per_year = computeSum(fluxOverIce, a_refRatio,a_crseDx,
				Interval(0,0), 0);

  //free storage
  for (int lev = 0; lev <= a_finestLevel ; lev++)
    {
      if (fluxOverIce[lev] != NULL)
	{
	  delete fluxOverIce[lev]; fluxOverIce[lev] = NULL;


	}
    }

  return tot_per_year;
}

Real 
DomainDiagnosticData::computeFluxOverMaskedIce(int a_maskVal,
					       const Vector<LevelData<FArrayBox>* > a_flux,
					       const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
					       const Vector<DisjointBoxLayout>& a_grids,
					       const Vector<int>& a_refRatio, Real a_crseDx,
					       int a_finestLevel) const
{
  Real tot_per_year = 0.0;

  //compute per-cell area fraction
  Vector<LevelData<FArrayBox>* > area_fraction(a_finestLevel+1, NULL);
  for (int lev=0; lev<=a_finestLevel; lev++)
    {
      area_fraction[lev] = new LevelData<FArrayBox>(a_grids[lev],1,IntVect::Zero);
      computeAreaFraction(*area_fraction[lev], a_coordSys, a_maskVal,lev, a_finestLevel );
      for (DataIterator dit(a_grids[lev]); dit.ok(); ++dit)
	{

	  (*area_fraction[lev])[dit].mult((*a_flux[lev])[dit]);
	}
    }
  

  // integrate over domain
  tot_per_year = computeSum(area_fraction, a_refRatio, a_crseDx, Interval(0,0), 0);

  // clean up temp storage
  for (int lev=0; lev<area_fraction.size(); lev++)
    {
      if (area_fraction[lev] != NULL)
	{
	  delete area_fraction[lev];
	  area_fraction[lev] = NULL;
	}
    }

  return tot_per_year;
}

Real 
DomainDiagnosticData::computeDeltaVolumeOverIce(const Vector<LevelData<FArrayBox>* >& a_old_thickness,
						const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
						const Vector<DisjointBoxLayout>& a_grids,
						const Vector<int>& a_refRatio, Real a_crseDx,
						int a_finestLevel) const
{

  //compute sum of a flux component over ice
  //construct deltaVolOverIce
  Vector<LevelData<FArrayBox>* > deltaVolOverIce ( a_finestLevel+1, NULL);
  for (int lev = 0; lev <= a_finestLevel ; lev++)
    {
      deltaVolOverIce[lev] = new
	LevelData<FArrayBox>(a_grids[lev],1, IntVect::Zero);
      const LevelData<FArrayBox>& thk = a_coordSys[lev]->getH();
      //const LevelData<FArrayBox>* flux = a_flux[lev];
       
      for (DataIterator dit(a_grids[lev]); dit.ok(); ++dit)
	{
	  const Box& box =  a_grids[lev][dit];
	  const FArrayBox& oldH = (*a_old_thickness[lev])[dit];
	  const FArrayBox& dit_thck = thk[dit];
	  FArrayBox& dit_deltaVolOverIce = (*deltaVolOverIce[lev])[dit];
	     
	  for (BoxIterator bit(box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      // set deltaVolOverIce to source if thck > 0
	      if (dit_thck(iv) < 1e-10)
		{
		  dit_deltaVolOverIce(iv) = 0.0;
		}
	      else
		{
		  dit_deltaVolOverIce(iv) = dit_thck(iv)-oldH(iv);
		}
	    }
	    
	}
    }
  // compute sum
  Real tot_per_year = computeSum(deltaVolOverIce, a_refRatio, a_crseDx,
				Interval(0,0), 0);

  //free storage
  for (int lev = 0; lev <= a_finestLevel ; lev++)
    {
      if (deltaVolOverIce[lev] != NULL)
	{
	  delete deltaVolOverIce[lev]; deltaVolOverIce[lev] = NULL;


	}
    }

  return tot_per_year;
}

Real 
DomainDiagnosticData::computeTotalFlux(const Vector<LevelData<FArrayBox>* > a_flux,
				       const Vector<int>& a_refRatio, 
				       Real a_crseDx) const
{

  //compute sum of a flux for whole domain  
  Real tot_per_year = computeSum(a_flux, a_refRatio, a_crseDx, Interval(0,0), 0);

  return tot_per_year;
}

void
DomainDiagnosticData::initDiagnostics
(AmrIce& a_amrIce,
 const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatio, Real a_crseDx,
 Real a_time, int a_finestLevel)
{

  m_report_grounded_ice = false;
  m_report_area = false;
  m_report_total_flux = false;
  m_report_calving = false;

  ParmParse pp("amr");

  pp.query("report_sum_grounded_ice", m_report_grounded_ice);
  pp.query("report_ice_area", m_report_area);
  pp.query("report_total_flux", m_report_total_flux);
  pp.query("report_calving", m_report_calving);

  m_initialSumIce = computeTotalIce(a_coordSys, a_refRatio, a_crseDx, a_finestLevel);
  m_lastSumIce = m_initialSumIce;

  m_initialSumGroundedIce = computeTotalGroundedIce
    (a_grids, a_coordSys, a_refRatio, a_crseDx, a_finestLevel);
  m_lastSumGroundedIce = m_initialSumGroundedIce;
  m_initialVolumeAboveFlotation = computeVolumeAboveFlotation
    (a_coordSys, a_refRatio, a_crseDx, a_finestLevel);
  m_lastVolumeAboveFlotation = m_initialVolumeAboveFlotation; 

  m_diagnostic_values[0]=a_time;
  m_diagnostic_values[1]=m_initialSumIce;
  m_diagnostic_values[2]=computeVolumeAboveFlotation
    (a_coordSys, a_refRatio, a_crseDx, a_finestLevel);
  m_diagnostic_values[3]=computeGroundedArea
    (a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
  m_diagnostic_values[4]=computeFloatingArea
    (a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
  m_diagnostic_values[5]=NAN;
  m_diagnostic_values[6]=NAN;
  m_diagnostic_values[7]=NAN;

  record(a_amrIce);

}

// DIAGNOSTICS
// Diagnostic routine -- compute discharge

void 
DomainDiagnosticData::computeDischarge
(const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys, 
 const Vector<LevelData<FluxBox>* >& a_vectFluxes,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<Real>& a_dx,
 const Vector<int>& a_refRatio,
 Real a_time, Real a_offsetTime,
 int a_cur_step, int a_finestLevel, int a_verbosity)
{

  Real sumDischarge = 0.0;
  Real sumGroundedDischarge = 0.0;
  Real sumDischargeToOcean = 0.0;

  Vector<LevelData<FArrayBox>* > vectDischarge ( a_finestLevel+1, NULL);
  Vector<LevelData<FArrayBox>* > vectGroundedDischarge ( a_finestLevel+1, NULL);
  Vector<LevelData<FArrayBox>* > vectDischargeToOcean ( a_finestLevel+1, NULL);

  for (int lev=0; lev<=a_finestLevel; lev++)
    {
      vectDischarge[lev] = new LevelData<FArrayBox>(a_grids[lev],1,
							    IntVect::Zero);
      LevelData<FArrayBox>& levelDischarge = *vectDischarge[lev];
      vectGroundedDischarge[lev] = new LevelData<FArrayBox>(a_grids[lev],1,
							    IntVect::Zero);
      LevelData<FArrayBox>& levelGroundedDischarge = *vectGroundedDischarge[lev];
      vectDischargeToOcean[lev] = new LevelData<FArrayBox>(a_grids[lev],1,
							    IntVect::Zero);
      LevelData<FArrayBox>& levelDischargeToOcean = *vectDischargeToOcean[lev];

      const LevelData<FArrayBox>& levelThickness =  a_coordSys[lev]->getH();
      const LevelData<BaseFab<int> >& levelMask = a_coordSys[lev]->getFloatingMask();

      DataIterator dit=levelDischarge.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  const FluxBox& vflux = (*a_vectFluxes[lev])[dit];
	  const BaseFab<int>& mask = levelMask[dit];
	  const FArrayBox& thk = levelThickness[dit];

	  FArrayBox& discharge = levelDischarge[dit];
	  FArrayBox& groundedDischarge = levelGroundedDischarge[dit];
	  FArrayBox& dischargeToOcean = levelDischargeToOcean[dit];
	  discharge.setVal(0.0);
	  groundedDischarge.setVal(0.0);
	  dischargeToOcean.setVal(0.0);

	  for (int dir=0; dir<SpaceDim; dir++)
	    {

	      const FArrayBox& flux = vflux[dir];
	      BoxIterator bit(discharge.box());
	      for (bit.begin(); bit.ok(); ++bit)
		{
		  IntVect iv = bit();
		  Real smallThk = 10.0;
		  if ((thk(iv) < smallThk) || (mask(iv) != GROUNDEDMASKVAL))
		    {
		      if (thk(iv + BASISV(dir)) > smallThk && (mask(iv + BASISV(dir)) == GROUNDEDMASKVAL) ) 
			{
			  groundedDischarge(iv) += -flux(iv + BASISV(dir)) / a_dx[lev];
			}
		      if (thk(iv - BASISV(dir)) > smallThk && (mask(iv - BASISV(dir)) == GROUNDEDMASKVAL) )
			{
			  groundedDischarge(iv) += flux(iv) / a_dx[lev];
			}

		    }		  
		  if (thk(iv) < tiny_thickness) 
		    {
		      if (thk(iv + BASISV(dir)) > tiny_thickness)
			{
			  discharge(iv) += -flux(iv + BASISV(dir)) / a_dx[lev];
			}
		      if (thk(iv - BASISV(dir)) > tiny_thickness)
			{
			  discharge(iv) += flux(iv) / a_dx[lev];
			}

		    }
		  if ((thk(iv) < tiny_thickness) && (mask(iv) == OPENSEAMASKVAL)) 
		    {
		      if (thk(iv + BASISV(dir)) > tiny_thickness)
			{
			  dischargeToOcean(iv) += -flux(iv + BASISV(dir)) / a_dx[lev];
			}
		      if (thk(iv - BASISV(dir)) > tiny_thickness)
			{
			  dischargeToOcean(iv) += flux(iv) / a_dx[lev];
			}

		    }

		}
	    } // end direction 
	}

    } // end loop over levels
  
  // now compute sum
    sumDischarge = computeSum(vectDischarge, a_refRatio,
  				a_dx[0], Interval(0,0), 0);
    sumGroundedDischarge = computeSum(vectGroundedDischarge, a_refRatio,
  				a_dx[0], Interval(0,0), 0);
    sumDischargeToOcean = computeSum(vectDischargeToOcean, a_refRatio,
  				a_dx[0], Interval(0,0), 0);

  if (a_verbosity > 0) 
    {
      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
	     << ": DischargeFromIceEdge = " << sumDischarge << " m3/y " << endl;

      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
	     << ": DischargeFromGroundedIce = " << sumGroundedDischarge << " m3/y " << endl;
      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
	     << ": DischargeToOcean = " << sumDischargeToOcean << " m3/y " << endl;


    }  

  // clean up temp storage
  for (int lev=0; lev<vectDischarge.size(); lev++)
    {
      if (vectDischarge[lev] != NULL)
	{
	  delete vectDischarge[lev];
	  vectDischarge[lev] = NULL;
	}
    }
  for (int lev=0; lev<vectGroundedDischarge.size(); lev++)
    {
      if (vectGroundedDischarge[lev] != NULL)
	{
	  delete vectGroundedDischarge[lev];
	  vectGroundedDischarge[lev] = NULL;
	}
    }
  for (int lev=0; lev<vectDischargeToOcean.size(); lev++)
    {
      if (vectDischargeToOcean[lev] != NULL)
	{
	  delete vectDischargeToOcean[lev];
	  vectDischargeToOcean[lev] = NULL;
	}
    }

}

void 
DomainDiagnosticData::endTimestepDiagnostics
(const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const Vector<LevelData<FArrayBox>* > a_old_thickness,
 const Vector<LevelData<FArrayBox>* > a_divThicknessFlux, 
 const Vector<LevelData<FArrayBox>* > a_basalThicknessSource,
 const Vector<LevelData<FArrayBox>* > a_surfaceThicknessSource,
 const Vector<LevelData<FArrayBox>* > a_volumeThicknessSource,
 const Vector<LevelData<FArrayBox>* > a_calvedIceThickness,
 const Vector<LevelData<FArrayBox>* > a_addedIceThickness,
 const Vector<LevelData<FArrayBox>* > a_removedIceThickness,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatio, Real a_crseDx,
 Real a_time, Real a_offsetTime, Real a_dt,
 int a_cur_step, int a_finestLevel, int a_verbosity)
{

      Real sumIce = computeTotalIce(a_coordSys, a_refRatio, a_crseDx, a_finestLevel);

      m_diagnostic_values[0] = a_time;
      m_diagnostic_values[1] = sumIce;

      Real VAF=0.0;
      Real groundedArea = 0.0, floatingArea = 0.0;
      Real sumBasalFlux = 0.0, sumBasalFluxOverFloatingIce = 0.0;
      Real sumSurfaceFlux = 0.0;
      Real sumCalvedIce = 0.0, cflux = 0.0;
      Real sumVolFlux = 0.0;

      VAF = computeVolumeAboveFlotation(a_coordSys, a_refRatio, a_crseDx, a_finestLevel);
      groundedArea = computeGroundedArea(a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
      floatingArea = computeFloatingArea(a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
      sumBasalFlux = computeTotalFlux(a_basalThicknessSource, a_refRatio, a_crseDx);
      sumBasalFluxOverFloatingIce = computeFluxOverMaskedIce(FLOATINGMASKVAL, 
	     a_basalThicknessSource, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
      sumSurfaceFlux = computeTotalFlux(a_surfaceThicknessSource, a_refRatio, a_crseDx);
      sumCalvedIce = computeSum(a_calvedIceThickness, a_refRatio,a_crseDx,
					Interval(0,0), 0);
      sumVolFlux = computeTotalFlux(a_volumeThicknessSource, a_refRatio, a_crseDx);

      if (a_dt > 0)
	{
	  cflux=sumCalvedIce/a_dt;
	}

      m_diagnostic_values[2] = VAF;
      m_diagnostic_values[3] = groundedArea;
      m_diagnostic_values[4] = floatingArea;
      m_diagnostic_values[5] = sumSurfaceFlux;
      m_diagnostic_values[6] = sumBasalFlux;	    
      m_diagnostic_values[7] = sumBasalFluxOverFloatingIce;	    
      m_diagnostic_values[8] = cflux;
      m_diagnostic_values[9] = cflux+sumVolFlux;

      if (a_verbosity > 0) 
	{
	  Real diffSum = sumIce - m_lastSumIce;
	  Real totalDiffSum = sumIce - m_initialSumIce;
  
	  Real sumGroundedIce = 0.0, diffSumGrounded = 0.0, totalDiffGrounded = 0.0;
	  Real diffVAF = 0.0, totalDiffVAF = 0.0;
	  Real sumBasalFluxOverIce = 0.0;
	  Real sumRemovedIce = 0.0, sumAddedIce = 0.0;
	  Real sumAccumCalvedIce = 0.0;
	  Real diffAccumCalvedIce;
	  Real totalLostIce = 0.0;
	  //Real totalLostOverIce = 0.0;
	  Real sumDeltaVolumeOverIce = 0.0;
	  Real sumSurfaceFluxOverIce = 0.0;
	  Real sumDivThckFluxOverIce = 0.0, sumDivThckFlux = 0.0;
	  Real sumVolFluxOverIce = 0.0;
	  //Real sumCalvedOverIce = 0.0;
	  Real sumRemovedOverIce = 0.0, sumAddedOverIce = 0.0;
	  if (m_report_grounded_ice)
	    {
	      sumGroundedIce = computeTotalGroundedIce(a_grids, a_coordSys, a_refRatio, a_crseDx, a_finestLevel);
	      diffSumGrounded = sumGroundedIce - m_lastSumGroundedIce;
	      totalDiffGrounded = sumGroundedIce - m_initialSumGroundedIce;      
	      m_lastSumGroundedIce = sumGroundedIce;
      
	      diffVAF = VAF -  m_lastVolumeAboveFlotation;
	      totalDiffVAF = VAF - m_initialVolumeAboveFlotation;

	      m_lastVolumeAboveFlotation = VAF;


	    }
	  
	  if (m_report_total_flux)

	    {
	      sumDeltaVolumeOverIce = computeDeltaVolumeOverIce
		(a_old_thickness, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
	      sumBasalFluxOverIce = computeFluxOverIce
		(a_basalThicknessSource, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
	      sumSurfaceFluxOverIce = computeFluxOverIce
		(a_surfaceThicknessSource, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
	      sumDivThckFluxOverIce = computeFluxOverIce
		(a_divThicknessFlux, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
	      sumDivThckFlux = computeTotalFlux(a_divThicknessFlux, a_refRatio, a_crseDx);

	      sumVolFluxOverIce = computeFluxOverIce
		(a_volumeThicknessSource, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);

	    }
	  
	  if (m_report_calving)

	    {
	      // sumAccumCalvedIce = computeSum(m_melangeThickness, a_refRatio,a_crseDx,
	      // 				Interval(0,0), 0);
	      diffAccumCalvedIce=sumAccumCalvedIce-m_lastSumCalvedIce;
	      sumRemovedIce = computeSum(a_removedIceThickness, a_refRatio,a_crseDx,
					 Interval(0,0), 0);
	      sumAddedIce = computeSum(a_addedIceThickness, a_refRatio,a_crseDx,
					 Interval(0,0), 0);
	      //sumCalvedOverIce = computeFluxOverIce(a_calvedIceThickness);
	      sumRemovedOverIce = computeFluxOverIce
		(a_removedIceThickness, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
	      sumAddedOverIce = computeFluxOverIce
		(a_addedIceThickness, a_coordSys, a_grids, a_refRatio, a_crseDx, a_finestLevel);
	      totalLostIce = sumCalvedIce+sumRemovedIce+sumAddedIce;
	      //totalLostOverIce = sumCalvedOverIce+sumRemovedOverIce+sumAddedOverIce;

	      m_lastSumCalvedIce = sumAccumCalvedIce;

	    }


	  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) " 
		 << ": sum(ice) = " << sumIce 
		 << " ( " << diffSum
		 << " " << totalDiffSum
		 << " )" << endl;
      
	  if (m_report_grounded_ice)
	    {
	      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
		     << ": sum(grounded ice) = " << sumGroundedIce 
		     << " ( " << diffSumGrounded
		     << " " << totalDiffGrounded
		     << " )" << endl;

	      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
		     << ": VolumeAboveFlotation = " << VAF
		     << " ( " << diffVAF
		     << " " << totalDiffVAF
		     << " )" << endl;
	    } 
	  if (m_report_area)
	    {
	      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
		     << ": GroundedArea = " << groundedArea << " m2 " << endl;

	      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
		     << ": FloatingArea = " << floatingArea << " m2 " << endl;

	    } 

	  if (m_report_total_flux)
	    {
	      if (a_dt > 0)
		{
		  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			 << ": BasalFlux = " << sumBasalFluxOverIce << " m3/yr " 
			 << " ( " << sumBasalFlux 
			 << "  " << sumBasalFlux-sumBasalFluxOverIce
			 << " )"
			 << endl;

		  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			 << ": SurfaceFlux = " << sumSurfaceFluxOverIce << " m3/yr  " 
			 << " ( " << sumSurfaceFlux 
			 << "  " << sumSurfaceFlux-sumSurfaceFluxOverIce 
			 << " )"
			 << endl;

		  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			 << ": DivergenceThicknessFlux = " << sumDivThckFluxOverIce << " m3/yr " 
			 << " ( " << sumDivThckFlux 
			 << "  " << sumDivThckFlux-sumDivThckFluxOverIce
			 << " )"
			 << endl;
		  if (abs(sumVolFlux) > 1.0e-6)
		    // Print the sums of volumeThickessnessSource if they are not zero. 
		    // \todo: If these qualities are not zero their contribution to the mass convervation equation, Domain error and Ice sheet error below, need be included.
		    {
		      pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			     << ": VolumeFlux = " << sumVolFluxOverIce << " m3/yr  " 
			     << " ( " << sumVolFlux 
			     << "  " << sumVolFlux-sumVolFluxOverIce 
			     << " )"
			     << endl;
		    }
		}
	    }



	  if (m_report_calving)
	    {
	      if (a_dt > 0)
		{
		  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			 << ": AccumCalvedIce = " << sumAccumCalvedIce << " m3 " 
			 << " ( " << diffAccumCalvedIce << "  " << diffAccumCalvedIce - totalLostIce << " ) " << endl;
		  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			 << ": CalvedIce = " << sumCalvedIce << " m3 " << " RemovedIce = " << sumRemovedIce << " m3 " << " AddedIce = " << sumAddedIce << " m3 Sum " << totalLostIce << " m3 " << endl;
		}
	    }

	  if (m_report_calving && m_report_total_flux)
	    {
	      if (a_dt > 0)
		{
		  Real adjflux=(sumRemovedIce+sumAddedIce)/a_dt;
		  Real calvingerr=sumSurfaceFlux+sumBasalFlux-(cflux+diffSum/a_dt+adjflux);
		  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			 << ": Domain error = " << calvingerr << " m3/yr"
			 << " ( dV/dt = " << diffSum/a_dt 
			 << " calving flux = " << cflux
			 << " SMB = " << sumSurfaceFlux
			 << " BMB = " << sumBasalFlux
			 << " adjustment flux to maintain front = " << adjflux
			 << " )"  << endl;
	      
		  adjflux=(sumRemovedOverIce+sumAddedOverIce)/a_dt;
		  Real err=sumSurfaceFluxOverIce+sumBasalFluxOverIce-
		    (sumDivThckFluxOverIce+sumDeltaVolumeOverIce/a_dt+adjflux);
		  pout() << "Step " << a_cur_step << ", time = " << a_time << " ( " << a_offsetTime << " ) "
			 << ": Ice sheet error = " << err << " m3/yr"
			 << " ( dV/dt = " << sumDeltaVolumeOverIce/a_dt 
			 << " flux = " << sumDivThckFluxOverIce
			 << " smb = " << sumSurfaceFluxOverIce
			 << " bmb = " << sumBasalFluxOverIce
			 << " adjustment flux to maintain front = " << adjflux
			 << " )" << endl;

		}
	    }
	}

      m_lastSumIce = sumIce;

}

void DomainDiagnosticData::record(AmrIce& a_amrIce)
{

  Real dt = 1.2345678e+300;
  size_t ithis = m_time.size();
  size_t iprev = ithis - 1;
  if (ithis > 0)
    {
      dt = a_amrIce.time() - m_time[iprev]; 
    }

  Real rhoi = a_amrIce.m_iceDensity;
  if (dt > 0)
    {

      CH_assert(a_amrIce.time() == m_diagnostic_values[0]);

      m_time.push_back( a_amrIce.time() );
      m_ice_volume.push_back(m_diagnostic_values[1]*rhoi);
      m_ice_vaf.push_back(m_diagnostic_values[2]*rhoi);
      m_ice_grounded_area.push_back(m_diagnostic_values[3]);
      m_ice_floating_area.push_back(m_diagnostic_values[4]);
      m_ice_total_smb.push_back(m_diagnostic_values[5]*rhoi);
      m_ice_total_bmb.push_back(m_diagnostic_values[6]*rhoi);
      m_ice_floating_total_bmb.push_back(m_diagnostic_values[7]*rhoi);
      m_ice_total_calving_flux.push_back(m_diagnostic_values[8]*rhoi);
      m_ice_total_calving_and_ice_front_melting_flux.push_back(m_diagnostic_values[9]*rhoi);

    }
}

inline 
void last_to_first_resize(Vector<Real>& a_v)
{
  a_v[0] = a_v[a_v.size()-1];
  a_v.resize(1);
}


void DomainDiagnosticData::reset()
{

 if (m_time.size() > 0)
   {
     /// copy the last items from the previous set of records
     last_to_first_resize(m_time);
     last_to_first_resize(m_ice_volume);
     last_to_first_resize(m_ice_vaf);
     last_to_first_resize(m_ice_grounded_area);
     last_to_first_resize(m_ice_floating_area);
     last_to_first_resize(m_ice_total_smb);
     last_to_first_resize(m_ice_total_bmb);
     last_to_first_resize(m_ice_floating_total_bmb);
     last_to_first_resize(m_ice_total_calving_flux);
     last_to_first_resize(m_ice_total_calving_and_ice_front_melting_flux);
   }
 else
   {
     m_time.resize(0);
     m_ice_volume.resize(0);
     m_ice_vaf.resize(0);
     m_ice_grounded_area.resize(0);
     m_ice_floating_area.resize(0);
     m_ice_total_smb.resize(0);
     m_ice_total_bmb.resize(0);
     m_ice_floating_total_bmb.resize(0);
     m_ice_total_calving_flux.resize(0);
     m_ice_total_calving_and_ice_front_melting_flux.resize(0);
   }
}

#ifdef CH_USE_HDF5

void readV(HDF5Handle& a_handle, Vector<Real>& a_data,
	   const std::string& a_name)
{

  H5E_auto_t efunc; void* edata;
  // need these to turn auto error messaging off then back 
  
#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  hid_t set = H5Dopen(a_handle.groupID(), a_name.c_str());
  H5Eset_auto(efunc, edata);
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  hid_t set = H5Dopen2(a_handle.groupID(), a_name.c_str(), H5P_DEFAULT);
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
#endif
  if (set >= 0)
    {
      hid_t space = H5Dget_space(set);
      if (space >= 0)
	{
	  hsize_t dim,maxdim;
	  H5Sget_simple_extent_dims(space, &dim, &maxdim);
	  a_data.resize(dim);
	  readDataset(set, space, &a_data[0], 0,  a_data.size());
	  H5Sclose(space);
	}
      H5Dclose(set);
    }
}

void readAttr(hid_t loc_id, const std::string& a_name, std::string& a_value)
{
  herr_t ret = 0;
  char messg[1024];

#ifdef H516
  hid_t attr   = H5Aopen_name(loc_id, a_name.c_str());
#else
  hid_t attr   = H5Aopen_by_name(loc_id,".", a_name.c_str(), H5P_DEFAULT, H5P_DEFAULT );
#endif
  hid_t atype  = H5Aget_type(attr);
  hid_t aclass = H5Tget_class(atype);
  char* buf = NULL;  size_t size = 0;

  size = H5Tget_size(atype);
  buf = new char[size+1];
  ret = H5Aread(attr, atype, buf);
  if (ret < 0) 
    {
      sprintf(messg,"DomainDiagnosticData:: Problem reading attribute %s",a_name.c_str());
      MayDay::Warning(messg);
    }

  buf[size] = 0; // for some reason HDF5 is not null terminating strings correctly
  a_value = std::string(buf);

  delete[] buf;
  H5Tclose(atype);
  H5Aclose(attr);

}

void readStruct(HDF5Handle& a_handle, cfDiagnostic& a_cf_info)
{

  H5E_auto_t efunc; void* edata;
  // need these to turn auto error messaging off then back 

  Vector<Real>& data = *a_cf_info.data;
  std:: string data_name = a_cf_info.short_name;
  
#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  hid_t set = H5Dopen(a_handle.groupID(), data_name.c_str());
  H5Eset_auto(efunc, edata);
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  hid_t set = H5Dopen2(a_handle.groupID(), data_name.c_str(), H5P_DEFAULT);
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
#endif
  if (set >= 0)
    {
      hid_t space = H5Dget_space(set);
      if (space >= 0)
	{
	  hsize_t dim,maxdim;
	  H5Sget_simple_extent_dims(space, &dim, &maxdim);
	  data.resize(dim);
	  readDataset(set, space, &data[0], 0,  data.size());
	  H5Sclose(space);
	}

      H5Dclose(set);
    }
}

void writeV(HDF5Handle& a_handle, Vector<Real>& a_data,
	    const std::string& a_name)
{
  hid_t set,space,type;
  createDataset(set, space, a_handle, a_name, &a_data[0], a_data.size());
  if (procID() == 0)
    {
      writeDataset(set, space, &a_data[0], 0,  a_data.size());
    }
  H5Sclose(space);
  H5Dclose(set);
}

void writeAttr(hid_t loc_id, const std::string& a_name, const std::string& a_value)
{
  H5E_auto_t efunc; void* edata;
#ifdef H516
  H5Eget_auto(&efunc, &edata);
#else
  H5Eget_auto2(H5E_DEFAULT, &efunc, &edata);
#endif
  herr_t ret = 0;
  char messg[1024];

  hid_t s_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(s_type, a_value.length());
  hid_t aid  = H5Screate(H5S_SCALAR);
#ifdef H516
  H5Eset_auto(NULL, NULL);
  hid_t attr = H5Acreate(loc_id, a_name.c_str(), s_type,
			 aid, H5P_DEFAULT);
#else
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  hid_t attr = H5Acreate2(loc_id, a_name.c_str(), s_type,
			  aid, H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0)
    {
      H5Adelete(loc_id, a_name.c_str());
#ifdef H516
      attr = H5Acreate(loc_id, a_name.c_str(), s_type,
		       aid, H5P_DEFAULT);
#else
      attr = H5Acreate2(loc_id, a_name.c_str(), s_type,
			aid, H5P_DEFAULT,H5P_DEFAULT);
#endif
      if (attr < 0)
	{
	  sprintf(messg,"DomainDiagnosticData:: Problem creating attribute %s",a_name.c_str());
	  MayDay::Warning(messg);
	}
    }
#ifdef H516
  H5Eset_auto(efunc, edata);
#else
  H5Eset_auto2(H5E_DEFAULT, efunc, edata);
#endif
  char* tmp = (char*)a_value.c_str();
  ret = H5Awrite(attr, s_type, tmp);
  if (ret < 0) 
    {
      sprintf(messg,"DomainDiagnosticData:: Problem writing attribute %s",a_name.c_str());
      MayDay::Warning(messg);
    }

  H5Sclose(aid);
  H5Aclose(attr);
  H5Tclose(s_type);

}

void writeStruct(HDF5Handle& a_handle, const cfDiagnostic& a_cf_info)
{
  hid_t set,space,type;

  Vector<Real>& data = *a_cf_info.data;
  std:: string data_name = a_cf_info.short_name;

  createDataset(set, space, a_handle, data_name, &data[0], data.size());
  if (procID() == 0)
    {
      writeDataset(set, space, &data[0], 0,  data.size());
    }

  writeAttr(set, "Short name", a_cf_info.short_name);
  writeAttr(set, "Long name", a_cf_info.long_name);
  writeAttr(set, "Units", a_cf_info.units);
  writeAttr(set, "Standard name", a_cf_info.cf_name);

  H5Sclose(space);
  H5Dclose(set);
}

void DomainDiagnosticData::write(HDF5Handle& a_handle)
{  
  if (a_handle.pushGroup(HDF5_SUBGROUP_NAME) == 0)
    {     
      for (int i = 0; i < m_cf_stuff.size(); ++i)
	{
	  cfDiagnostic cf_info = m_cf_stuff[i];
	  writeStruct(a_handle,cf_info);
	}
      a_handle.popGroup();
    }
}

void DomainDiagnosticData::read(HDF5Handle& a_handle)
{
  if (a_handle.pushGroup(HDF5_SUBGROUP_NAME) == 0)
    {
      for (int i = 0; i < m_cf_stuff.size(); ++i)
	{
	  cfDiagnostic cf_info = m_cf_stuff[i];
	  readStruct(a_handle,cf_info);
	}
      a_handle.popGroup();
    }
}

DomainDiagnosticData::DomainDiagnosticData()
{
  setCFdata();
}

DomainDiagnosticData::DomainDiagnosticData(const DomainDiagnosticData& a)
{
  setCFdata();
}

DomainDiagnosticData&  DomainDiagnosticData::operator=(const DomainDiagnosticData& a)
{
  setCFdata();
  return *this;
}



#endif
#include "NamespaceFooter.H"
