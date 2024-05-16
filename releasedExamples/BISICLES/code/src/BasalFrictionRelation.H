
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BASALFRICTIONRELATION_H_
#define _BASALFRICTIONRELATION_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "NamespaceHeader.H"

class LevelSigmaCS;
class DataIterator;

/// Virtual base class for basal friction relations.
/** Sublclasses must, at least, compute a positive \f$ \alpha \f$
    given basal velocity \f$ u_b \f$ such that basal traction 
    \f$ \tau _b  = \alpha(u_b, ...) u_b \f$. 
*/
class BasalFrictionRelation
{
public:

  virtual ~BasalFrictionRelation() {;}

  

  virtual void updateThermodynamics( Vector<LevelData<FArrayBox>* >* a_tillWaterDepth )
  {
  }

    
  
  /// Computes cell-centered effective basal friction coeffcient \f$ \alpha \f$ such that  \f$ \tau_b  = \alpha (u_b, C) u_b \f$
  /**
      \param a_alpha effective basal friction coefficient based on the local velocity field. 
      \param a_basalVel Cell-centered basal velocity field.
      \param a_C Cell-centered coefficient field
      \param a_coords whole level geometry data
      \param a_dit index needed to retrieve data from a_coords
      \param a_box cell-centered box over which to do this computation
  */

  virtual void computeAlpha(FArrayBox& a_alpha,
			    const FArrayBox& a_basalVel,
			    const FArrayBox& a_C, const Real& a_scale,
			    const LevelSigmaCS& a_coords, 
			    const DataIterator& a_dit,
			    int a_level,
			    const Box& a_box) const = 0;
  
  /// Compute cell-centered dissipation rate \f$ s = T_b . u_b  / k_ice \f$
  /** 
      \param a_dissipation dissipation rate based on the local velocity field. 
      \param a_basalVel Cell-centered basal velocity field.
      \param a_scale assume basalVel is scaled by a factor a_scale
      \param a_C Cell-centered coefficient field
      \param a_coords whole level geometry data
      \param a_dit index needed to retrieve data from a_coords
      \param a_box cell-centered box over which to do this computation
  */
  virtual void computeDissipation(FArrayBox& a_dissipation,  
				  const FArrayBox& a_basalVel,
				  const FArrayBox& a_C, const Real& a_scale,
				  const LevelSigmaCS& a_coords, 
				  const DataIterator& a_dit,
				  int a_level,
				  const Box& a_box) const;

  
  /// creates a new copy of subclass objects.
  virtual BasalFrictionRelation* getNewBasalFrictionRelation() const = 0;

  virtual Real power() const = 0;

  /// Run-time configuration of BasalFrictionRelation objects
  static BasalFrictionRelation* parse(const char* a_prefix, int a_recursion );

};


/// implements BasalFrictionRelation to provide a power law
/**
   \f$ \tau_b  = - C h_{ab}^n |u_b|^{m-1} u_b \f$
   
   where \f$tau_b\f$ is the basal traction, \f$C\f$ is a coefficient 
   not dependent on velocity, \f$h_{ab}\f$ is the thickness above flotation,  
   \f$u_b\f$ is the basal velocity, and \f$n,m\f$ are exponents. Typically \f$m \leq 1\f$.

    - \f$ m = 1 \f$ gives a linear viscous sliding law 
    - \f$ m = 1/3 \f$ is often used (e.g in the MISMIP tests)

    An optional factor \f$ (|u_b|/|u_r| + 1)^{-m} ) \f$ is available (see Joughin 2019)
    that lead \f$ d|T_b|/d|u| \rightarrow -> 0\f$ as \f$|u| \rightarrow \infty \f$ 
	
    An alternate version of the above, \f$ N (N_r^{1/m} |u_b|/|u_r| + N^{1/m} )^{-m} ) \f$ is available

    Run-time configuation is carried out in BasalFrictionRelation::parse.
    If the string `A.B.basalFrictionRelation = powerLaw` is found, the expected
    form is
    
    - `BasalFrictionPowerLaw.m = <float. e,g 0.3333 for m = 1/3>`
    - `BasalFrictionPowerLaw.includeEffectivePressure = <bool>`
    - `BasalFrictionPowerLaw.fastSlidingSpeed = <float. e.g 300 for Joughin 2019>`
    - `BasalFrictionPowerLaw.highThickness = <float. e.g 1e2 for ~100 kPa>`
    
    if includeEffectivePressure, \f$n = m\f$ so that there is a 
    dependence on the thickness over flotation \f$h_ab\f$. Otherwise,
    \f$n = 0\f$. Do not use this option!

 */
class BasalFrictionPowerLaw : public BasalFrictionRelation
{
public:

  BasalFrictionPowerLaw(Real a_m,
			Real a_fastSlidingSpeed = -1.0,
			Real a_highThickness = -1.0  ,
			bool a_includeEffectivePressure = false) 
    : m_m(a_m),  m_includeEffectivePressure(a_includeEffectivePressure) ,
      m_fastSlidingSpeed(a_fastSlidingSpeed), m_highThickness(a_highThickness)
  {;}
  
  virtual ~BasalFrictionPowerLaw() {;}

  virtual void computeAlpha(FArrayBox& a_alpha,
			    const FArrayBox& a_basalVel,
			    const FArrayBox& a_C, const Real& a_scale,
			    const LevelSigmaCS& a_coords, 
			    const DataIterator& a_dit,
			    int a_level,
			    const Box& a_box) const;

  BasalFrictionRelation* getNewBasalFrictionRelation() const{
    return static_cast<BasalFrictionRelation*> (new BasalFrictionPowerLaw(m_m, m_fastSlidingSpeed, m_highThickness, m_includeEffectivePressure));
  }

  /// \callergraph
  virtual Real power() const {return m_m;}
  
private:
  Real m_m;
  bool m_includeEffectivePressure;
  Real m_fastSlidingSpeed; // a scale u_r for fast flow : u >> u_r -> Tb -> a constant
  Real m_highThickness; // a scale for thickness above flotation
  BasalFrictionPowerLaw();

};

/// Basal friction relation which satisfies Tb = min(|Tb'|, a*p) u/|u|
/**
   Tb' is any BasalFrictionRelation
   a is a constant (for now) coefficient, N is the local (hydrostatic) pressure
   For a ~ 1, Tb = Tb' is applied most places, with the Coulomb law 
   Tb = a*p*u/|u| imposed at  flotation

   Three laws are supported, one due to Tsai 2015

   \f$ \tau_b = \min(|\tau _b^{\prime}|, \alpha N^p) u/|u| \f$

   the second due to Schoof 2005, Gagliardini 2007, and Leguy 2014, which
   is smoother but otherwise similar. We use the form from Leguy 2014

   \f$ \tau _b = \tau _b^{\prime} ( N^p / (\alpha |u| + N^p ) )^{1/p}\f$
   

    Run-time configuation is carried out in BasalFrictionRelation::parse.
    If the string `A.B.basalFrictionRelation = pressureLimitedLaw` is found, the expected
    form is
    
     - `BasalFrictionPressureLimitedLaw.coefficient = <float, optional, defualt 1.0>` (\f$\alpha\f$)
     - `BasalFrictionPressureLimitedLaw.power = <float, optional, defualt 1.0>` (\f$p\f$)
     - `BasalFrictionPressureLimitedLaw.model = Tsai|Leguy (optional , default Tsai)`
     - `BasalFrictionPressureLimitedLaw.basalFrictionRelation = <valid BasalFrictionRelation>`


*/ 
class PressureLimitedBasalFrictionRelation : public BasalFrictionRelation
{

public:
  enum Model {Tsai,Leguy,MAX_MODEL};

private:
  Real m_a, m_p, m_maxTillWaterDepth, m_tillPressureFactor;
  Model m_model;
  Vector<LevelData<FArrayBox>* >* m_tillWaterDepth;
  BasalFrictionRelation* m_bfr; 
  
public:


  virtual void computeDissipation(FArrayBox& a_dissipation,  
				  const FArrayBox& a_basalVel,
				  const FArrayBox& a_C, const Real& a_scale,
				  const LevelSigmaCS& a_coords, 
				  const DataIterator& a_dit,
				  int a_level,
				  const Box& a_box) const;
  
  void updateThermodynamics( Vector<LevelData<FArrayBox>* >* a_tillWaterDepth  )
  {
    m_tillWaterDepth = a_tillWaterDepth;
  }
  virtual void computeAlpha(FArrayBox& a_alpha,
			    const FArrayBox& a_basalVel,
			    const FArrayBox& a_C, const Real& a_scale,
			    const LevelSigmaCS& a_coords, 
			    const DataIterator& a_dit,
			    int a_level,
			    const Box& a_box) const;


  PressureLimitedBasalFrictionRelation(Real a_a, Real a_p,  Real a_maxTillWaterDepth,
				       Real a_tillPressureFactor, 
				       Model a_model,  BasalFrictionRelation* a_bfr)
    : m_a(a_a), m_p(a_p), m_maxTillWaterDepth(a_maxTillWaterDepth),
      m_tillPressureFactor(a_tillPressureFactor),  m_model(a_model),  
      m_bfr(a_bfr->getNewBasalFrictionRelation())
  {

  }


  virtual BasalFrictionRelation* getNewBasalFrictionRelation() const
  {
    PressureLimitedBasalFrictionRelation* ptr = 
      new PressureLimitedBasalFrictionRelation(m_a, m_p, m_maxTillWaterDepth,m_tillPressureFactor ,
					       m_model, m_bfr);

    return static_cast<BasalFrictionRelation*>(ptr);
  }
    

  virtual Real power() const
  {
    return m_bfr->power();
  }

};



#include "NamespaceFooter.H"

#endif
