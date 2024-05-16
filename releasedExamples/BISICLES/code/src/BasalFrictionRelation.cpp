#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "BasalFrictionRelation.H"
#include "BasalFrictionRelationF_F.H"
#include "IceConstants.H"
#include "LevelSigmaCS.H"
#include "IceThermodynamics.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"


/// computes cell-centered \f$ s = - T_b . u_b \f$
/// based on the cell-centered velocity
/// and a coefficient field beta such that \f T_b  = - \alpha (u, C) u_b \f
/** 
    \param a_C -- \f$C\f$ based on the local velocity field. 
    \param a_basalVel -- Cell-centered basal velocity field.
    \param a_scale --   true velocity is a_basalVel/a_scale
    \param a_box: cell-centered box over which to do this computation 
*/
void 
BasalFrictionRelation::computeDissipation(FArrayBox& a_dissipation,  
					  const FArrayBox& a_basalVel,
					  const FArrayBox& a_C, const Real& a_scale,
					  const LevelSigmaCS& a_coords, 
					  const DataIterator& a_dit,
					  int a_level,
					  const Box& a_box) const
{
    
  computeAlpha(a_dissipation, a_basalVel, a_C,  a_scale, a_coords, a_dit, a_level, a_box);
  
  FORT_BFRICTIONAUU(CHF_FRA1(a_dissipation,0),
		    CHF_CONST_FRA(a_basalVel),
		    CHF_BOX(a_box));

}


void
BasalFrictionPowerLaw::computeAlpha(FArrayBox& a_alpha,
				    const FArrayBox& a_basalVel,
				    const FArrayBox& a_C, const Real& a_scale,
				    const LevelSigmaCS& a_coords, 
				    const DataIterator& a_dit,
				    int a_level,
				    const Box& a_box) const
{
  const Real mm1 = m_m - 1.0;
  const FArrayBox& thckOverFlotation = a_coords.getThicknessOverFlotation()[a_dit];
  const BaseFab<int>& mask = a_coords.getFloatingMask()[a_dit];

  FArrayBox Cscale(a_C.box(), a_C.nComp());
  Cscale.copy(a_C); Cscale /= std::pow(a_scale, m_m);
  
  const Real pexp = (m_includeEffectivePressure)?m_m:0.0;
  FORT_BFRICTIONPOWER(CHF_FRA1(a_alpha,0),
		      CHF_CONST_FRA(a_basalVel),
		      CHF_CONST_FRA1(Cscale,0),
		      CHF_CONST_FRA1(thckOverFlotation,0),
		      CHF_CONST_FIA1(mask,0),
		      CHF_CONST_REAL(mm1),
		      CHF_CONST_REAL(pexp),
		      CHF_BOX(a_box));

  if (m_fastSlidingSpeed > TINY_VEL)
    {
      FORT_BFRICTIONJOUGHIN(CHF_FRA1(a_alpha,0),
			    CHF_CONST_FRA(a_basalVel),
			    CHF_CONST_REAL(m_fastSlidingSpeed),
			    CHF_CONST_FRA1(thckOverFlotation,0),
			    CHF_CONST_REAL(m_highThickness),
			    CHF_CONST_REAL(m_m),
			    CHF_BOX(a_box));
    }

  

}
/// computes cell-centered \f$ s = - T_b . u_b \f$
/// based on the cell-centered velocity
/// and a coefficient field beta such that \f T_b  = - \alpha (u, C) u_b \f
/** 
    \param a_C -- \f$C\f$ based on the local velocity field. 
    \param a_basalVel -- Cell-centered basal velocity field.
    \param a_box: cell-centered box over which to do this computation 
*/
void 
PressureLimitedBasalFrictionRelation::computeDissipation(FArrayBox& a_dissipation,  
							 const FArrayBox& a_basalVel,
							 const FArrayBox& a_C, const Real& a_scale,
							 const LevelSigmaCS& a_coords, 
							 const DataIterator& a_dit,
							 int a_level,
							 const Box& a_box) const
{
    
  computeAlpha(a_dissipation, a_basalVel, a_C, a_scale, a_coords, a_dit, a_level, a_box);
  //m_bfr->computeAlpha(a_dissipation, a_basalVel,  a_C, a_coords,  a_dit, a_level, a_box);
  FORT_BFRICTIONAUU(CHF_FRA1(a_dissipation,0),
		    CHF_CONST_FRA(a_basalVel),
		    CHF_BOX(a_box));

}


void
PressureLimitedBasalFrictionRelation::computeAlpha
(FArrayBox& a_alpha,
 const FArrayBox& a_basalVel,
 const FArrayBox& a_C,const Real& a_scale,
 const LevelSigmaCS& a_coords, 
 const DataIterator& a_dit, int a_level,
 const Box& a_box) const
{
  m_bfr->computeAlpha(a_alpha, a_basalVel,  a_C,a_scale, a_coords,  a_dit, a_level, a_box);

  FArrayBox effectivePressure(a_box,1);
  const FArrayBox& thckOverFlotation = a_coords.getThicknessOverFlotation()[a_dit];
  const FArrayBox& thck = a_coords.getH()[a_dit];
  const Real eps = 1.0e-6;
  if ( Abs( m_p - 1.0) < eps)
    {
      // p == 1;  effective pressure = rho * g * (h-hf)
      effectivePressure.copy( thckOverFlotation );
      effectivePressure *= a_coords.iceDensity() * a_coords.gravity();
    }
  else
    {
      // p != 1 : effective pressure formula from Leguy 2014 (could use this for p = too, but
      // since there are extra calculations and we think p == 1 is a common case, don't
      Real rhog = a_coords.iceDensity() * a_coords.gravity();
      FORT_BFRICTIONLEGUYEFFPRES(CHF_FRA1(effectivePressure,0),
				 CHF_CONST_FRA1(thckOverFlotation,0),
				 CHF_CONST_FRA1(thck,0),
				 CHF_CONST_REAL(m_p),
				 CHF_CONST_REAL(rhog),
				 CHF_BOX(a_box));
    }
  
  FArrayBox& N = effectivePressure;
  N += 1.0e-10;
  // optionally, reduce the effective pressure according to the till water depth
  // formula from Van Pelt and Oerlemans 2012 (Eqn 2).
  // \todo allow the Bueler and Ven Pelt 2015 formula?
  if (m_maxTillWaterDepth > 0.0)
    {
      const Real& wtMax = m_maxTillWaterDepth;
      const FArrayBox& wt = (*(*m_tillWaterDepth)[a_level])[a_dit];
      for (BoxIterator bit(N.box()); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real pw = m_tillPressureFactor * N(iv) * std::min(wt(iv),wtMax) / wtMax;
	  N(iv) =  std::max(0.0, std::min(N(iv), N(iv) - pw));
	}
    }
  
  if (m_model == Tsai)
    {

      FArrayBox f(N.box(),1);
      f.setVal(m_a);
      
      FORT_BFRICTIONPLIMITTSAI(CHF_FRA1(a_alpha,0),
			       CHF_CONST_FRA(a_basalVel),
			       CHF_CONST_FRA1(N,0),
			       CHF_CONST_FRA1(f,0),
			       CHF_BOX(a_box));
    }
  else if (m_model == Leguy)
    {
      Real n = 1.0 / power();
      FORT_BFRICTIONPLIMITLEGUY(CHF_FRA1(a_alpha,0),
				CHF_CONST_FRA(a_basalVel),
				CHF_CONST_FRA1(N,0),
				CHF_CONST_REAL(m_a),
				CHF_CONST_REAL(n),
				CHF_BOX(a_box));
    }

}


BasalFrictionRelation* 
BasalFrictionRelation::parse(const char* a_prefix, int a_recursion)
{
  //can't do more than one recursion without changing the input file syntax
  CH_assert(a_recursion < 2);

  BasalFrictionRelation* basalFrictionRelationPtr = NULL;
  std::string type = "powerLaw";
  ParmParse pp(a_prefix);
  pp.query("basalFrictionRelation",type);

  if (type == "powerLaw")
      {
	ParmParse plPP("BasalFrictionPowerLaw");

	Real m = 1.0;
	plPP.query("m",m);
	bool includeEffectivePressure = false;
	plPP.query("includeEffectivePressure",includeEffectivePressure);
	Real fastSlidingSpeed = -HUGE_VEL; // choosing a negative value is equivalent to  fastSlidingSpeed -> infinity
	plPP.query("fastSlidingSpeed", fastSlidingSpeed);
	Real highThickness = -HUGE_THICKNESS; // choosing a negative value is equivalent to highThickness -> infinity
	plPP.query("highThickness", highThickness);
	
	BasalFrictionPowerLaw*  pl = new BasalFrictionPowerLaw(m,fastSlidingSpeed,highThickness,includeEffectivePressure);
	basalFrictionRelationPtr = static_cast<BasalFrictionRelation*>(pl);
      }
  else if (type == "pressureLimitedLaw")
    {
     
      std::string prefix("BasalFrictionPressureLimitedLaw");
      ParmParse ppPLL(prefix.c_str());
     
     
      std::string models = "Tsai";
      ppPLL.query("model",models);

      PressureLimitedBasalFrictionRelation::Model model = 
	PressureLimitedBasalFrictionRelation::MAX_MODEL;
      
      if (models == "Tsai")
	{
	  model = PressureLimitedBasalFrictionRelation::Tsai;
	}
      else if (models == "Leguy")
	{
	  model = PressureLimitedBasalFrictionRelation::Leguy;
	}
      else
	{
	  MayDay::Error("undefined BasalFrictionPressureLimitedLaw.model in inputs");
	}

      Real a = 1.0;
      ppPLL.query("coefficient",a);
      Real p = 1.0;
      ppPLL.query("power",p);

      Real wtMax = -1.0;
      ppPLL.query("maxTillWaterDepth", wtMax);

      Real pwFactor = 0.97; // Corresponds to (alpha) in Van Pelt and Oerlemans 2012 (Eqn 2)
      ppPLL.query("tillPressureFactor", pwFactor);
      
      BasalFrictionRelation* innerLaw = BasalFrictionRelation::parse(prefix.c_str(), a_recursion + 1);
      PressureLimitedBasalFrictionRelation* ptr = 
	new PressureLimitedBasalFrictionRelation(a, p, wtMax, pwFactor,  model, innerLaw);

      basalFrictionRelationPtr = static_cast<BasalFrictionRelation*>(ptr);
    }
  else
    {
      MayDay::Error("undefined basalFrictionRelation in inputs");
    }
  return basalFrictionRelationPtr;
}


#include "NamespaceFooter.H"
