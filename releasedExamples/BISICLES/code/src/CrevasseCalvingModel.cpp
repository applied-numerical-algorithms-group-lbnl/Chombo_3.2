#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "CrevasseCalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "CrevasseF_F.H"
#include "CalvingF_F.H"
#include "SigmaCSF_F.H"
#include "NamespaceHeader.H"

/// Calving models based on crevase depth computed from viscous stress 
/// CrevasseCalvingModel implements the methods common to all of these, 
/// and requires its subclasses to compute crevasse depth given the stress

//compute a scalar invariant (e.g an eigenvalue) of the stress tensor
void CrevasseCalvingModel::computeStressMeasure(LevelData<FArrayBox>& a_stressMeasure,
						const AmrIce& a_amrIce,int a_level)
{
  CH_TIME("CrevasseCalvingModel::computeStressMeasure");
  CH_assert(m_stressMeasure < MAX_STRESS_MEASURE);
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const LevelData<FArrayBox>& levelVT  = *a_amrIce.viscousTensor(a_level);
  const LevelData<FArrayBox>& levelVel = *(a_amrIce.velocity(a_level));
  for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& sm = a_stressMeasure[dit];
      const FArrayBox& thck = levelCoords.getH()[dit];
      const FArrayBox& vt = levelVT[dit];
      const FArrayBox& vel = levelVel[dit];
      Box b = levelCoords.grids()[dit];
      b.grow(1); // need one layer of ghosts

      //all this ought to be pushed to fortran kernels
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real oneOnThck = thck(iv) / (1.0e-6 + thck(iv)*thck(iv));
	  Real sxx = oneOnThck*vt(iv,0);
	  Real syy = oneOnThck*vt(iv,3);
	  Real sxy = 0.5*oneOnThck*(vt(iv,2) + vt(iv,1));

	  if (m_stressMeasure == FirstPrincipalStress)
	    {
	      //vertically averaged first  principal stress
	      Real a = 0.5 * (sxx + syy); 
	      Real b = std::sqrt (  std::pow ( 0.5*(sxx - syy), 2) + std::pow(sxy,2));
	      sm(iv) = a + b;
	    }
	  else if  (m_stressMeasure == AlongFlowNormalStress)
	    {
	      //vertically averaged normal stress along flow
	      const Real& ux = vel(iv,0); 
	      const Real& uy = vel(iv,1);
	      Real usq = ux*ux + uy*uy + 1.0e-6;
	      sm(iv) = (sxx*ux*ux + syy*uy*uy + 2.0*sxy*ux*uy)/(usq);
	    } 
	  else if  (m_stressMeasure == Trace)
	    {
	      //Trace of vertically averaged stress tensor
	      sm(iv)  = (sxx + syy);
	    }
	  else
	    {
	      //should never get here
	      CH_assert(m_stressMeasure < MAX_STRESS_MEASURE);
	    }
	}
    }

}

void CrevasseCalvingModel::getWaterDepth(LevelData<FArrayBox>& a_waterDepth, const AmrIce& a_amrIce,int a_level)
{
  m_waterDepth->evaluate(a_waterDepth, a_amrIce, a_level, 0.0);
}
	      		 
void CrevasseCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce, 
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  //domain edge calving model always applies
  m_domainEdgeCalvingModel->applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
  
  if  (a_stage == PostVelocitySolve)
    {
      // stress model only makes sense when the thickness and veolcity are in sync
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      LevelData<FArrayBox> locals(levelCoords.grids(), 6 , IntVect::Unit);

      //compute a (scalar) stress measure;
      LevelData<FArrayBox> stressMeasure;aliasLevelData(stressMeasure, &locals, Interval(0,0));
      computeStressMeasure(stressMeasure, a_amrIce,a_level);

      //compute effective thickness - equal to thickness except in partially covered cells
      //and resulting surface and thickness above flotation
      LevelData<FArrayBox> thcke;aliasLevelData(thcke, &locals, Interval(1,1));
      LevelData<FArrayBox> habe;aliasLevelData(habe, &locals, Interval(2,2));
      LevelData<FArrayBox> usrfe;aliasLevelData(usrfe, &locals, Interval(3,3));

      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  thcke[dit].copy(a_thickness[dit]);
	  FORT_EFFECTIVETHICKNESS( CHF_FRA1(thcke[dit],0), 
				   CHF_CONST_FRA1(a_thickness[dit],0),
				   CHF_CONST_FRA1(a_iceFrac[dit],0), 
				   CHF_BOX(thcke[dit].box()));
	  
	  Real rhoi = levelCoords.iceDensity();
	  Real rhoo = levelCoords.waterDensity(); 
	  Real sl = levelCoords.seaLevel();
	  const FArrayBox& topg = levelCoords.getTopography()[dit];
	  
	  FORT_THICKNESSOVERFLOTATION(CHF_FRA1(habe[dit],0),
				      CHF_CONST_FRA1(thcke[dit],0),
				      CHF_CONST_FRA1(topg,0),
				      CHF_CONST_REAL(rhoi),
				      CHF_CONST_REAL(rhoo),
				      CHF_CONST_REAL(sl),
				      CHF_BOX(habe[dit].box()));
	  
	  FORT_SURFACEHEIGHT(CHF_FRA1(usrfe[dit],0),
			     CHF_CONST_FRA1(thcke[dit],0),
			     CHF_CONST_FRA1(topg,0),
			     CHF_CONST_REAL(rhoi),
			     CHF_CONST_REAL(rhoo),
			     CHF_CONST_REAL(sl),
			     CHF_BOX(usrfe[dit].box()));
	  
	} // end loop over boxes
     
      //compute remaining ice thickness
      LevelData<FArrayBox> remnant;aliasLevelData(remnant, &locals, Interval(4,4));
      LevelData<FArrayBox> waterDepth;aliasLevelData(waterDepth, &locals, Interval(5,5));
      //m_waterDepth->evaluate(waterDepth, a_amrIce, a_level, 0.0);
      getWaterDepth(waterDepth, a_amrIce, a_level);

      computeRemnant(remnant, stressMeasure, thcke, usrfe, habe, waterDepth, levelCoords);
      
      //used to make sure that only cells within some distance of the open sea calve: 
      LevelData<BaseFab<int> > grownOpenSea (levelCoords.grids(), 1 , IntVect::Unit);
      for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
	{
	   if (m_calvingZoneLength > 0.0)
	     {
	       // a finite distance
	       grownOpenSea[dit].copy(levelCoords.getFloatingMask()[dit]);
	     }
	   else
	     {
	       //any distance at all
	       grownOpenSea[dit].setVal( OPENSEAMASKVAL );
	     }
	}

      if (m_calvingZoneLength > 0.0)
	{
	  //grow the opensea mask by nCalve cells. 
	  int niter = std::max( 1,  int( m_calvingZoneLength / a_amrIce.dx(a_level)[0])) ;
	  LevelData<BaseFab<int> > prev (levelCoords.grids(), 1 , IntVect::Unit);
	  
	  for (int iter = 0; iter < niter; iter++)
	    {
	       for (DataIterator dit (levelCoords.grids()); dit.ok(); ++dit)
		 {
		   prev[dit].copy( grownOpenSea[dit]);
		   Box b = levelCoords.grids()[dit];
		   BaseFab<int>& m = grownOpenSea[dit];
		   const BaseFab<int>& mp = prev[dit];
		   for (BoxIterator bit(b); bit.ok(); ++bit)
		     {
		       const IntVect& iv = bit();
		       if (m(iv) != OPENSEAMASKVAL)
			 {
			   bool t = false;
			   t = t || (mp(iv + BASISV(0)) ==  OPENSEAMASKVAL);
			   t = t || (mp(iv - BASISV(0)) ==  OPENSEAMASKVAL);
			   t = t || (mp(iv + BASISV(1)) ==  OPENSEAMASKVAL);
			   t = t || (mp(iv - BASISV(1)) ==  OPENSEAMASKVAL);
			   if (t) m(iv) = OPENSEAMASKVAL;
			 }
		     }
		 }
	       grownOpenSea.exchange();
	    }
	}

      //update thickness and mask.
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  FArrayBox& thck = a_thickness[dit];
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  FArrayBox& iceFrac = a_iceFrac[dit];
	  
	  Box b = thck.box();
	  b &= iceFrac.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit(); 
	      Real prevThck = thck(iv);
	      if ( (grownOpenSea[dit](iv) == OPENSEAMASKVAL ) && ( remnant[dit](iv) < TINY_THICKNESS ))
		{
		  thck(iv) = 0.0;
		  iceFrac(iv,0) = 0.0;
		}

	      // Record gain/loss of ice.
	      if (calved.box().contains(iv))
		{
		  // grownOpenSea can change OPENLAND or GROUNDED mask to OPENSEA
		  updateCalvedIce(thck(iv),prevThck,grownOpenSea[dit](iv),added(iv),calved(iv),removed(iv));
		}

	    } // end loop over cells
	} // end loop over boxes
    } // end (a_stage == PostVelocitySolve || a_stage == PostRegrid)

}

CrevasseCalvingModel::CrevasseCalvingModel(ParmParse& a_pp)
{
  
  //DomainEdgeCalvingModel parameters
  Vector<int> frontLo(2,false); 
  a_pp.getarr("front_lo",frontLo,0,frontLo.size());
  Vector<int> frontHi(2,false);
  a_pp.getarr("front_hi",frontHi,0,frontHi.size());
  bool preserveSea = false;
  a_pp.query("preserveSea",preserveSea);
  bool preserveLand = false;
  a_pp.query("preserveLand",preserveLand);
  m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);
  
  //Crevasse Model parameters
  // m_waterDepth = 0.0;
  //a_pp.get("waterDepth",m_waterDepth);

  std::string prefix (a_pp.prefix());
  m_waterDepth = SurfaceFlux::parse( (prefix + ".CrevasseWaterDepth").c_str());

  if (m_waterDepth == NULL)
    {
      // try the older style waterDepth constant
      Real wd = -1.0;
      a_pp.get("waterDepth", wd);
      constantFlux* ptr = new constantFlux;
      ptr->setFluxVal(wd);
      m_waterDepth = static_cast<SurfaceFlux*>(ptr);
    }


  m_includeBasalCrevasses = false;
  a_pp.query("includeBasalCrevasses",m_includeBasalCrevasses);
  m_calvingZoneLength = -1.0;
  a_pp.query("calvingZoneLength",m_calvingZoneLength);
  
  m_stressMeasure = FirstPrincipalStress; // was the original default
  {
    std::string s = "FirstPrincipalStress"; 
    a_pp.query("stressMeasure",s);
    if (s == "FirstPrincipalStress")
      {
	m_stressMeasure = FirstPrincipalStress;
      }
    else if (s == "AlongFlowNormalStress")
      {
	m_stressMeasure = AlongFlowNormalStress;
      }
    else if (s == "Trace")
      {
	m_stressMeasure = Trace;
      }
    else
      {
	m_stressMeasure = MAX_STRESS_MEASURE;
	CH_assert(m_stressMeasure <  MAX_STRESS_MEASURE);
      }
  }

}



CrevasseCalvingModel::~CrevasseCalvingModel()
{
  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }
  if (m_waterDepth != NULL)
    {
      delete m_waterDepth;
      m_waterDepth = NULL;
    }

 
}

void BennCalvingModel::computeRemnant(LevelData<FArrayBox>& a_remnant,
				      const LevelData<FArrayBox>& a_stress,
				      const LevelData<FArrayBox>& a_thck,
				      const LevelData<FArrayBox>& a_usrf,
				      const LevelData<FArrayBox>& a_hab, 
				      const LevelData<FArrayBox>& a_waterDepth, 
				      const LevelSigmaCS& a_coords)
{
  
  CH_TIME("BennCalvingModel::computeRemnant");
 
  const Real& rhoi = a_coords.iceDensity();
  const Real& rhoo = a_coords.waterDensity();
  const Real& grav = a_coords.gravity();
  
  for (DataIterator dit(a_coords.grids());dit.ok();++dit)
    {
      FArrayBox& remnant = a_remnant[dit];
      const FArrayBox& hab = a_hab[dit];
      const FArrayBox& usrf = a_usrf[dit];
      const FArrayBox& thck = a_thck[dit];
      const FArrayBox& wd = a_waterDepth[dit];
      const FArrayBox& s = a_stress[dit];
      Box b = a_coords.grids()[dit];
      b.grow(1); //need one layer of ghost cells
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  Real Ds = std::max(s(iv),0.0) / (grav*rhoi) + 1000.0/rhoi * wd(iv);

	  if (m_includeBasalCrevasses)
	    {
	      //explicit basal crevasse depth calculation
	      Real Db = ((rhoi/(rhoo-rhoi)) * ( s(iv) /(grav*rhoi) - hab(iv)));
	      remnant(iv) = std::max( 0.0, thck(iv) - Ds -  Db);
	    }
	  else
	    {
	      //assume full thickness fracture if surface crevasses reach sea-level
	      remnant(iv) = std::max( 0.0, usrf(iv) - Ds);
	    }
	}
    }
  a_remnant.exchange();
}

BennCalvingModel::BennCalvingModel(ParmParse& a_pp)  
  : CrevasseCalvingModel(a_pp)
{
  
}

BennCalvingModel::~BennCalvingModel()
{

}


void VdVCalvingModel::computeRemnant(LevelData<FArrayBox>& a_remnant,
				     const LevelData<FArrayBox>& a_stress,
				     const LevelData<FArrayBox>& a_thck,
				     const LevelData<FArrayBox>& a_usrf,
				     const LevelData<FArrayBox>& a_hab, 
				     const LevelData<FArrayBox>& a_waterDepth, 
				     const LevelSigmaCS& a_coords)
{
  
  CH_TIME("VdVCalvingModel::computeRemnant");

  const Real& rhoi = a_coords.iceDensity();
  const Real& rhoo = a_coords.waterDensity();
  const Real& grav = a_coords.gravity();
  // compute the stress intensity for fracture to the waterline
  for (DataIterator dit(a_coords.grids());dit.ok();++dit)
    {
      FArrayBox& remnant = a_remnant[dit];
      const FArrayBox& usrf = a_usrf[dit];
      const FArrayBox& thck = a_thck[dit];
      const FArrayBox& wd = a_waterDepth[dit];
      const FArrayBox& s = a_stress[dit];
      Box b = a_coords.grids()[dit];
      b.grow(1); //need one layer of ghost cells

      FArrayBox depth(b,1); // crevasse depth -
      FArrayBox K(b,2); //stress intensity;

#define TOUGHNESS 0.1e6
      if (m_includeBasalCrevasses)
	{
	  // stress intensity for full thickness basal crevasses
	  FArrayBox thckp(b,1); thckp.setVal(0.0);
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      thckp(iv) = std::max(thck(iv)-usrf(iv),0.0); // peizometric height
	      depth(iv) = thck(iv);
	    }
	  
	  //surface crevasses
	  FORT_VDVSTRESSS( CHF_FRA1(K,0),
			   CHF_CONST_FRA1(thck,0),
			   CHF_CONST_FRA1(wd,0),
			   CHF_CONST_FRA1(depth,0),
			   CHF_CONST_FRA1(s,0),
			   CHF_CONST_REAL(rhoi),
			   CHF_CONST_REAL(rhoo),
			   CHF_CONST_REAL(grav),
			   CHF_BOX(b));
	  //basal crevasses
	  FORT_VDVSTRESSB( CHF_FRA1(K,1),
			   CHF_CONST_FRA1(thck,0),
			   CHF_CONST_FRA1(thckp,0),
			   CHF_CONST_FRA1(depth,0),
			   CHF_CONST_FRA1(s,0),
			   CHF_CONST_REAL(rhoi),
			   CHF_CONST_REAL(rhoo),
			   CHF_CONST_REAL(grav),
			   CHF_BOX(b));


	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      remnant(iv) = std::max(0.0,TOUGHNESS - std::max(K(iv,1),K(iv,0)));
	    }

	}
      else
	{
	  // compute stress intensity for surface crevasse that full thickness or the water-line
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      depth(iv) = std::min(usrf(iv),.999*thck(iv));
	    }
	  FORT_VDVSTRESSS( CHF_FRA1(K,0),
			   CHF_CONST_FRA1(thck,0),
			   CHF_CONST_FRA1(wd,0),
			   CHF_CONST_FRA1(depth,0),
			   CHF_CONST_FRA1(s,0),
			   CHF_CONST_REAL(rhoi),
			   CHF_CONST_REAL(rhoo),
			   CHF_CONST_REAL(grav),
			   CHF_BOX(b));
	  
	 for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      remnant(iv) = std::max(0.0,TOUGHNESS - K(iv,0));
	    }
	} // end if m_includeBasalCrevasses
    } // end loop over boxes

  a_remnant.exchange();
}

VdVCalvingModel::VdVCalvingModel(ParmParse& a_pp)  
  : CrevasseCalvingModel(a_pp)
{
  
}

VdVCalvingModel::~VdVCalvingModel()
{

}

#include "NamespaceFooter.H"
