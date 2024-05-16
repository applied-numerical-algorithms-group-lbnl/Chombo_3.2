#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "CalvingModel.H"
#include "MaskedCalvingModel.H"
#include "CrevasseCalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

/// a default implementation
/**
   most models provide a criterion, rather than a rate.
 */
void
CalvingModel::getCalvingRate(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  DataIterator dit = a_calvingRate.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_calvingRate[dit].setVal(0.0);
    }
}

void
CalvingModel::getWaterDepth(LevelData<FArrayBox>& a_waterDepth, const AmrIce& a_amrIce,int a_level)
{
  DataIterator dit = a_waterDepth.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_waterDepth[dit].setVal(0.0);
    }
}

void
VariableRateCalvingModel::getCalvingRate(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  m_calvingRate->evaluate(a_calvingRate, a_amrIce, a_level, 0.0);
}


void 
DeglaciationCalvingModelA::applyCriterion
(LevelData<FArrayBox>& a_thickness, 
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce, 
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();

      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }
	      
	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}

void DomainEdgeCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const DisjointBoxLayout& grids = levelCoords.grids();
  const ProblemDomain domain = grids.physDomain();
  const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
  const IntVect ghost = a_thickness.ghostVect();
  //const LevelData<FArrayBox>& vt  = *a_amrIce.viscousTensor(a_level);
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      //const Box& gridBox = grids[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  if (!domain.isPeriodic(dir))
	    {

	      if (m_frontLo[dir] > 0)
		{
		  Box loBox = adjCellLo(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  loBox.grow(transverseVect);
		  loBox &= a_thickness[dit].box();
		  for (BoxIterator bit(loBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv + BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;
		      if (a_iceFrac[dit].box().contains(iv))
			a_iceFrac[dit](iv) = 0.0;
		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
					  a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

		    }
		}
	      
	      if (m_frontHi[dir] > 0)
		{
		  Box hiBox = adjCellHi(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  hiBox.grow(transverseVect);
		  hiBox &= a_thickness[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv - BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;
		      if (a_iceFrac[dit].box().contains(iv))
			a_iceFrac[dit](iv) = 0.0;
		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
				      a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = levelMask[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (m_preserveSea && mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (m_preserveLand && mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  thck(iv) = std::max(thck(iv),0.0);

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}

    } // end loop over boxes

}

void ProximityCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  Real time = a_amrIce.time();
  bool calvingActive = (time >= m_startTime && time < m_endTime);
  calvingActive = false;
  pout() << " time = " << time 
	 << " m_startTime = " <<  m_startTime
	 << " m_endTime = " <<  m_endTime
	 << "calvingActive = " << calvingActive
	 << std::endl;
  if (true || calvingActive)
    {
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& proximity = *a_amrIce.groundingLineProximity(a_level);
      const LevelData<FArrayBox>& velocity = *a_amrIce.velocity(a_level);
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& thck = a_thickness[dit];
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  const FArrayBox& prox = proximity[dit];
	  const FArrayBox& vel = velocity[dit];
	  Box b = thck.box();b &= prox.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real prevThck = thck(iv);
	      Real vmod = std::sqrt(vel(iv,0)*vel(iv,0) + vel(iv,1)*vel(iv,1));
	      if (prox(iv) < m_proximity && calvingActive && vmod > m_velocity)
		{
		  //thck(iv) *= 0.5; thck(iv) = max(thck(iv),10.0);
		  thck(iv) = 0.0;
		}
	      if (mask(iv) == OPENSEAMASKVAL)
		{
		   thck(iv) = 0.0;
		}
	      if (mask(iv) == FLOATINGMASKVAL)
		{
		  thck(iv) = max(thck(iv),1.0);
		}

	      // Record gain/loss of ice
	      if (calved.box().contains(iv))
		{
		  updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		}

	    }
	}
    }
}



CalvingModel* CalvingModel::parseCalvingModel(const char* a_prefix)
{

  CalvingModel* ptr = NULL;
  std::string type = "";
  ParmParse pp(a_prefix);
  pp.query("type",type);
  
  if (type == "NoCalvingModel")
    {
      ptr = new NoCalvingModel;
    }
  else if (type == "DomainEdgeCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new DomainEdgeCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "FixedFrontCalvingModel")
    {
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      ptr = new DeglaciationCalvingModelA
	(0.0,  1.0e+10, minThickness, -1.2345678e+300, 1.2345678e+300);
    }
  else if (type == "DeglaciationCalvingModelA")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.get("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelA
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "DeglaciationCalvingModelB")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelB
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "ProximityCalvingModel")
    {
      Real proximity = 0.0;
      pp.get("proximity", proximity );
      Real velocity = 0.0;
      pp.query("velocity", velocity );
      Real startTime = -1.2345678e+300;
      pp.get("startTime",  startTime);
      Real endTime = 1.2345678e+300;
      pp.get("endTime",  endTime);
      ptr = new ProximityCalvingModel(proximity,velocity, startTime, endTime);
    }
  else if (type == "FlotationCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new FlotationCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "BennCalvingModel")
    {
      ptr = new BennCalvingModel(pp);
    }
  else if (type == "VanDerVeenCalvingModel")
    {
      ptr = new VdVCalvingModel(pp);
    }
  else if (type == "ThicknessCalvingModel")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      bool factorMuCoef = false;
      pp.query("factor_mu_coef",factorMuCoef); 
      ptr = new ThicknessCalvingModel
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime, factorMuCoef); 
    }
  else if (type == "CliffCollapseCalvingModel")
    {  
      Real maxCliffHeight = 100.0;
      pp.get("max_cliff_height", maxCliffHeight);
      Real recessionRate = 0.0;
      pp.get("recession_rate", recessionRate );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new CliffCollapseCalvingModel
	(maxCliffHeight, recessionRate, startTime, endTime); 
    }
  else if (type == "MaxiumumExtentCalvingModel")  
    {
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);

      Vector<Real> vect(SpaceDim,0.0);

      pp.getarr("lowLoc",vect,0,SpaceDim);
      RealVect lowLoc(D_DECL(vect[0], vect[1],vect[2]));      

      pp.getarr("highLoc",vect,0,SpaceDim);
      RealVect highLoc(D_DECL(vect[0], vect[1],vect[2]));      

      MaximumExtentCalvingModel* Ptr = new MaximumExtentCalvingModel(highLoc,
                                                                     lowLoc,
                                                                     startTime,
                                                                     endTime);
      ptr = static_cast<CalvingModel*>(Ptr);

    }
  else if (type == "MaskedCalvingModel")
    {
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );

      // masked calving model uses a surfaceFlux as a mask
      std::string mask_prefix(a_prefix);
      mask_prefix += ".mask";
      SurfaceFlux* mask_ptr = SurfaceFlux::parse(mask_prefix.c_str());

      MaskedCalvingModel* Ptr = new MaskedCalvingModel(mask_ptr, minThickness);

      ptr = static_cast<CalvingModel*>(Ptr);

      // MaskedCalvingModel makes a copy of the mask, so clean up here
      if (mask_ptr != NULL)
        {
          delete mask_ptr;
        }
    }
  else if (type == "CompositeCalvingModel")
    {
      int nElements;
      pp.get("nElements",nElements);
     
      std::string elementPrefix(a_prefix);
      elementPrefix += ".element";

      Vector<CalvingModel*> elements(nElements);
      for (int i = 0; i < nElements; i++)
        {
          std::string prefix(elementPrefix);
          char s[32];
          sprintf(s,"%i",i);
          prefix += s;
          ParmParse pe(prefix.c_str());
          elements[i] = parseCalvingModel(prefix.c_str());
          CH_assert(elements[i] != NULL);
        }
      CompositeCalvingModel* compositePtr = new CompositeCalvingModel(elements);
      ptr = static_cast<CalvingModel*>(compositePtr);
    }
  
  else if (type == "VariableRateCalvingModel")
    {
      ptr = new VariableRateCalvingModel(pp);
    }

   else if (type == "RateProportionalToSpeedCalvingModel")
    {
      ptr = new RateProportionalToSpeedCalvingModel(pp);
    }
  

  return ptr;
}


void 
DeglaciationCalvingModelB::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          else if ((mask(iv) == FLOATINGMASKVAL) && (thck(iv) < m_calvingThickness))
            {
	      thck(iv) = m_minThickness;              
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}


void 
ThicknessCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& iceFrac = a_iceFrac[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox effectiveThickness(thck.box(), 1);
      effectiveThickness.copy(thck);

      if (m_factorMuCoef)
	{
	  effectiveThickness *= a_amrIce.muCoef(a_level)[dit];
	}

      
      Box b = thck.box();
      b &= iceFrac.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();          
          // if iceFrac > 0, then rescale effectiveThickness
          // by dividing by iceFrac value, which gives "actual" thickness
          // in the partial cell. Probably eventually want to move this to 
          // fortran
	  Real prevThck = thck(iv);
          if (iceFrac(iv,0) > 0.0)
            {
              effectiveThickness(iv,0) /= iceFrac(iv,0);
            }
            
          if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          // allow ice to spread into open sea regions too, if appropriate
          else if (((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
                   && (effectiveThickness(iv) < m_calvingThickness))
            {
              // note that we're setting thck here, not effectiveThickness, 
              // which is a temporary
              // also set the iceFrac to zero in these cells
	      thck(iv) = m_minThickness; 
              iceFrac(iv,0) = 0.0;
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



  
//alter the thickness field at the end of a time step
void
MaximumExtentCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness,
					  LevelData<FArrayBox>& a_calvedIce,
					  LevelData<FArrayBox>& a_addedIce,
					  LevelData<FArrayBox>& a_removedIce,  
					  LevelData<FArrayBox>& a_iceFrac, 
					  const AmrIce& a_amrIce,
					  int a_level,
					  Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const Real dx = a_amrIce.amrDx()[a_level];
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
          // compute location of cell center
          RealVect loc(iv);          
          loc += 0.5*RealVect::Unit;
          loc *= dx;
          
          // check high and low extents
          if ((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
            {
              if (loc[0] <= m_lowLoc[0])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] <= m_lowLoc[1])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[0] > m_highLoc[0]) 
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] > m_highLoc[1])
                {
                  thck(iv) = 0.0;
                }
            } // end if floating or opensea

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

        } // end loop over cells
  
    }

}




//alter the thickness field at the end of a time step
void 
CompositeCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness, 
				      LevelData<FArrayBox>& a_calvedIce,
				      LevelData<FArrayBox>& a_addedIce,
				      LevelData<FArrayBox>& a_removedIce, 
				      LevelData<FArrayBox>& a_iceFrac, 
				      const AmrIce& a_amrIce,
				      int a_level,
				      Stage a_stage)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
    }
}

  
CompositeCalvingModel::~CompositeCalvingModel()
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      delete m_vectModels[n];
      m_vectModels[n] = NULL;
    }
}

void FlotationCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  m_domainEdgeCalvingModel.applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == FLOATINGMASKVAL)
	    {
	      thck(iv) = 0.0; 
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}

void
CalvingModel::updateCalvedIce(const Real& a_thck, const Real a_prevThck, const int a_mask, Real& a_added, Real& a_calved, Real& a_removed)
{

  if (a_thck > a_prevThck)
    {
      a_added += (a_prevThck-a_thck);
    }
  else 
    {
      if ((a_mask == OPENSEAMASKVAL) || (a_mask == FLOATINGMASKVAL))
	{
	  a_calved += (a_prevThck-a_thck);
	}
      else
	{
	  a_removed += (a_prevThck-a_thck);
	}
    } 

}


VariableRateCalvingModel::VariableRateCalvingModel(ParmParse& a_pp)
{
      Real startTime = -1.2345678e+300;
      a_pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      a_pp.query("end_time",  endTime);
 
      Vector<int> frontLo(2,false); 
      a_pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      a_pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      a_pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      a_pp.query("preserveLand",preserveLand);

      m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);

      std::string prefix (a_pp.prefix());
      m_calvingRate = SurfaceFlux::parse( (prefix + ".CalvingRate").c_str());

}

void VariableRateCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  (*m_domainEdgeCalvingModel).applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox& frac = a_iceFrac[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      Real frac_eps = 1.0e-6;
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  
	  if (frac(iv) < frac_eps)
	    {
	      thck(iv)=0.0;
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



CalvingModel* VariableRateCalvingModel::new_CalvingModel()
  {
    VariableRateCalvingModel* ptr = new VariableRateCalvingModel(*this);
    ptr->m_startTime = m_startTime;
    ptr->m_endTime = m_endTime;
    ptr->m_calvingRate = m_calvingRate->new_surfaceFlux();
    ptr->m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*m_domainEdgeCalvingModel);
    return ptr; 
  }

VariableRateCalvingModel::~VariableRateCalvingModel()
{

  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }

  if (m_calvingRate != NULL)
    {
      delete m_calvingRate;
      m_calvingRate = NULL;
    }

}


void 
CliffCollapseCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness,
					  LevelData<FArrayBox>& a_calvedIce,
					  LevelData<FArrayBox>& a_addedIce,
					  LevelData<FArrayBox>& a_removedIce,  
					  LevelData<FArrayBox>& a_iceFrac, 
					  const AmrIce& a_amrIce,
					  int a_level,
					  Stage a_stage)
{

  // only do this at the end of a timestep
  // (since that's the only time a time-integrated recession rate makes any sense)
  if (a_stage == PostThicknessAdvection)
    {
      Real dt = a_amrIce.dt();
      Real dx = a_amrIce.amrDx()[a_level];
      
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& surfaceHeight = levelCoords.getSurfaceHeight();
      
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& iceFrac = a_iceFrac[dit];
	  FArrayBox& thck = a_thickness[dit];
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  const FArrayBox& surface = surfaceHeight[dit];
	  FArrayBox effectiveSurface(surface.box(),1);
	  effectiveSurface.copy(surface);
	  FArrayBox effectiveThickness(thck.box(), 1);
	  effectiveThickness.copy(thck);
	  Box b = thck.box();
	  b &= iceFrac.box();
	  
	  // keep track of which cells we've already done in order to avoid double-counting
	  BaseFab<int> alreadyDone(b,1);
	  alreadyDone.setVal(0);
	  
	  Real phiNew;
	  
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();          
	      // if iceFrac > 0, then rescale effectiveThickness
	      // by dividing by iceFrac value, which gives "actual" thickness
	      // in the partial cell. Probably eventually want to move this to 
	      // fortran
	      // also compute "effective surface", which is the upper surface height based
	      // on the effective thickness rather than the cell-averaged thickness
	      // Probably eventually want to move this to fortran
	      Real prevThck = thck(iv);
	      if (iceFrac(iv,0) > 0.0)
		{
		  effectiveThickness(iv,0) /= iceFrac(iv,0);
		  effectiveSurface(iv,0) += effectiveThickness(iv,0) - thck(iv,0);	      
		}
	      
	      // if ice is grounded, look at neighbors to see if there are any empty neighbors, then look at
	      // surface differences.
	      if (mask(iv) == GROUNDEDMASKVAL) 
		{
		  // loop over directions
		  for (int dir=0; dir<SpaceDim; dir++)
		    {
		      IntVect shiftVect = BASISV(dir);
		      IntVect ivp = iv + shiftVect;
		      IntVect ivm = iv - shiftVect;
		      
		      // look in both high and low directions at once
		      if (((mask(ivp,0) != GROUNDEDMASKVAL) && (mask(ivp,0) != OPENLANDMASKVAL) &&  ((effectiveSurface(iv,0) - effectiveSurface(ivp,0)) > m_maxCliffThickness)) ||
			  ((mask(ivm,0) != GROUNDEDMASKVAL) && (mask(ivm,0) != OPENLANDMASKVAL) && ((effectiveSurface(iv,0) - effectiveSurface(ivm,0)) > m_maxCliffThickness)))
			{
			  // we have a cliff!  only adjust this cell if we haven't already
			  if (alreadyDone(iv,0) == 0)
			    {
			      alreadyDone(iv,0) = 1;
			      phiNew = iceFrac(iv,0) - m_recessionRate*dt/dx;
			      // don't go below zero
			      phiNew = Max(phiNew, 0.0);
			      
			      // note that we're setting thck here, not effectiveThickness, 
			      // which is a temporary
			      // also modify the iceMask to zero in these cells
			      
			      thck(iv,0) = thck(iv,0)*phiNew/iceFrac(iv,0);
			      iceFrac(iv,0) = phiNew;
			    } // end we haven't already done this one
			} // end if we have a cliff
		    } // end loop over directions
		} // end if this cell is grounded (no floating cliffs)
	      
	      
	      // Record gain/loss of ice
	      if (calved.box().contains(iv))
		{
		  updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		}
	      
	    } // end loop over cells in this box
	} // end loop over grids on this level	  
    } // end if we're at the post-advection stage
}


RateProportionalToSpeedCalvingModel::RateProportionalToSpeedCalvingModel(ParmParse& a_pp)
{
      Real startTime = -1.2345678e+300;
      a_pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      a_pp.query("end_time",  endTime);
 
      Vector<int> frontLo(2,false); 
      a_pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      a_pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      a_pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      a_pp.query("preserveLand",preserveLand);

      m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);

      std::string prefix (a_pp.prefix());
      m_proportion = SurfaceFlux::parse( (prefix + ".proportion").c_str());

}

void RateProportionalToSpeedCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  (*m_domainEdgeCalvingModel).applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox& frac = a_iceFrac[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      Real frac_eps = 1.0e-6;
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  
	  if (frac(iv) < frac_eps)
	    {
	      thck(iv)=0.0;
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



CalvingModel* RateProportionalToSpeedCalvingModel::new_CalvingModel()
  {
    RateProportionalToSpeedCalvingModel* ptr = new RateProportionalToSpeedCalvingModel(*this);
    ptr->m_startTime = m_startTime;
    ptr->m_endTime = m_endTime;
    ptr->m_proportion = m_proportion->new_surfaceFlux();
    ptr->m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*m_domainEdgeCalvingModel);
    return ptr; 
  }

RateProportionalToSpeedCalvingModel::~RateProportionalToSpeedCalvingModel()
{

  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }

  if (m_proportion != NULL)
    {
      delete m_proportion;
      m_proportion = NULL;
    }

}


void
RateProportionalToSpeedCalvingModel::getCalvingRate
(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  m_proportion->evaluate(a_calvingRate, a_amrIce, a_level, 0.0);
  const LevelData<FArrayBox>& vel = *a_amrIce.velocity(a_level); // flux vel might be better  
  for (DataIterator dit(vel.dataIterator()); dit.ok(); ++dit)
    {
      Box b = a_calvingRate[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real usq = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      usq += vel[dit](iv,dir)*vel[dit](iv,dir);
	    }	
	  a_calvingRate[dit](iv) *= (1.0 + std::sqrt(usq));
	}
    }

  int dbg = 0;dbg++; 
}

#include "NamespaceFooter.H"
