#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#define TEST_WD 0.0;

#include "DamageConstitutiveRelation.H"
#include "NyeCrevasseF_F.H"
#include "SigmaCSF_F.H"
#include "LevelMappedDerivatives.H"
#include "NamespaceHeader.H"

DamageConstitutiveRelation::~DamageConstitutiveRelation()
{
  if ( m_undamagedConstitutiveRelation != NULL)
    {
      delete m_undamagedConstitutiveRelation;
      m_undamagedConstitutiveRelation = NULL;
    }
}

ConstitutiveRelation* 
DamageConstitutiveRelation::getNewConstitutiveRelation() const
{
  DamageConstitutiveRelation* newPtr = 
    new DamageConstitutiveRelation(m_undamagedConstitutiveRelation, m_damageModel);

  return static_cast<ConstitutiveRelation*>(newPtr);
}

void
DamageConstitutiveRelation::computeMu(LevelData<FArrayBox>& a_mu,
                                    const LevelData<FArrayBox>& a_vel, const Real& a_scale,
                                    const LevelData<FArrayBox>* a_crseVelPtr,
                                    int a_nRefCrse,
                                    const LevelData<FArrayBox>& a_A,
                                    const LevelSigmaCS& a_coordSys,
				    const ProblemDomain& a_domain,
                                    const IntVect& a_ghostVect) const
{

  CH_TIME("DamageConstitutiveRelation::computeMu");
  MayDay::Error("DamageConstitutiveRelation::computeMu");
} 

void
DamageConstitutiveRelation::computeDissipation(LevelData<FArrayBox>& a_dissipation,
					     const LevelData<FArrayBox>& a_vel, 
					     const LevelData<FArrayBox>* a_crseVelPtr,
					     int a_nRefCrse,
					     const LevelData<FArrayBox>& a_A,
					     const LevelSigmaCS& a_coordSys,
					     const ProblemDomain& a_domain,
					     const IntVect& a_ghostVect) const
{
  CH_TIME("DamageConstitutiveRelation::computeDissipation");
  MayDay::Error("DamageConstitutiveRelation::computeDissipation");
}

/// modify the effective viscosity to reflect local and transported damage
/** 
    on input a_mu is the undamaged viscosity vertically integrated h*mu
    on output it is damaged h*(1-D)*mu
    
    where hD = max(hD',hD'')

    D' = depth/thickness is the local damage 
    (computed directly from the stress tensor)

    D'' is a transported damage, typically computed by the damage
        transport model 
  
*/
 void  
DamageConstitutiveRelation::modifyMu
(FArrayBox& a_damage,
 FArrayBox& a_mu,
 const FArrayBox& a_transportedDamage, 
 const FArrayBox& a_thck,
 const Box& a_box) 
{
  Real eps = 1.0e-10; // minimum denominator

  CH_assert(a_transportedDamage.max(a_box,0) < 1.0e+5);

  FArrayBox maxdam(a_box,1);
  maxdam.copy(a_damage);
  FORT_FABMAX(CHF_FRA1(maxdam,0), CHF_CONST_FRA1(a_transportedDamage,0), CHF_BOX(a_box));
  
  //limit max damage. should have happened elsewhere, but, belt and braces...
  FArrayBox maxmaxdam(a_box, 1);
  maxmaxdam.copy(a_thck);
  maxmaxdam *= DAMAGE_MAX_VAL;
  FORT_FABMIN(CHF_FRA1(maxdam, 0), CHF_FRA1(maxmaxdam, 0), CHF_BOX(a_box));

  //h * mu *= (1-d/h)
  FORT_CREVASSEMU(CHF_FRA1(a_mu,0),
		  CHF_CONST_FRA1(a_thck,0),
		  CHF_CONST_FRA1(maxdam,0),
		  CHF_CONST_REAL(eps),
		  CHF_BOX(a_box));
}

/// Compute the cell-centred local  damage 
/** local damage depends on the stress tensor, so that
    
    D*h = h*t / (rho * g)
    
*/
 void  DamageConstitutiveRelation::computeLocalDamageVT
(FArrayBox& a_damage,
 const FArrayBox& a_vt, 
 const FArrayBox& a_thck,
 const FArrayBox& a_topg,
 const FArrayBox& a_water,
 const Real& a_rhoi,
 const Real& a_rhoo,
 const Real& a_gravity,
 const Real& a_sealevel,
 const Box& a_box) 
{
  //All components of the velocity gradient;
  Interval DerivDir(0,SpaceDim-1);
  Interval DerivComps(0,SpaceDim-1);

  int dudxComp = derivComponent(0,0);
  int dudyComp = derivComponent(1,0);
  int dvdxComp = derivComponent(0,1);
  int dvdyComp = derivComponent(1,1);
	  
  //principal stresses
  FArrayBox lambda(a_box, SpaceDim);
  FORT_SYMTEIGEN(CHF_FRA(lambda), 
		 CHF_CONST_FRA(a_vt),
		 CHF_INT(dudxComp),
		 CHF_INT(dudyComp),
		 CHF_INT(dvdxComp),
		 CHF_INT(dvdyComp),
		 CHF_BOX(a_box));
  

	  
  Real eps = 1.0e-10; // minimum denominator
  
  FArrayBox thckab(a_box,1);
  FORT_THICKNESSOVERFLOTATION(CHF_FRA1(thckab,0),
			      CHF_CONST_FRA1(a_thck,0),
			      CHF_CONST_FRA1(a_topg,0),
			      CHF_CONST_REAL(a_rhoi),
			      CHF_CONST_REAL(a_rhoo),
			      CHF_CONST_REAL(a_sealevel),
			      CHF_BOX(a_box));
  
  Real maxDamage = DAMAGE_MAX_VAL;

  FORT_NYECREVASSEDEPTHT(CHF_FRA1(a_damage,0), 
		       CHF_CONST_FRA1(a_water,0),
		       CHF_CONST_FRA1(a_thck,0),
		       CHF_CONST_FRA1(lambda,0),
		       CHF_CONST_FRA1(thckab,0),
		       CHF_CONST_REAL(a_rhoi),
		       CHF_CONST_REAL(a_rhoo),
		       CHF_CONST_REAL(a_gravity),
		       CHF_CONST_REAL(eps),
		       CHF_CONST_REAL(maxDamage),
		       CHF_BOX(a_box));

}

/// Compute the local  damage 
/** local damage depends on the deviatoric stress tensor, so that

    h*tij = h*(1-D)*mu*eij

    D*h = h*t / (rho * g)

    mu is undamaged effective viscosity, grad(u), rate factor, etc 
*/
 void  
DamageConstitutiveRelation::computeLocalDamage
(FArrayBox& a_damage,
 const FArrayBox& a_mu,
 const FArrayBox& a_gradU, 
 const FArrayBox& a_thck,
 const FArrayBox& a_topg,
 const FArrayBox& a_water,
 const Real& a_rhoi,
 const Real& a_rhoo,
 const Real& a_gravity,
 const Real& a_sealevel,
 const Box& a_box) 
{
  //All components of the velocity gradient;
  Interval DerivDir(0,SpaceDim-1);
  Interval DerivComps(0,SpaceDim-1);

  int dudxComp = derivComponent(0,0);
  int dudyComp = derivComponent(1,0);
  int dvdxComp = derivComponent(0,1);
  int dvdyComp = derivComponent(1,1);
	  
  //e + I*tr(e)
  FArrayBox ee(a_box, SpaceDim*SpaceDim);
  FORT_EPLUSTRE(CHF_FRA(ee), 
		CHF_CONST_FRA(a_gradU),
		CHF_INT(dudxComp),
		CHF_INT(dudyComp),
		CHF_INT(dvdxComp),
		CHF_INT(dvdyComp),
		CHF_BOX(a_box));
  //2.0* (e + I*tr(e))
  ee *= 2.0;

  //principal strain rate components
  FArrayBox lambda(a_box, SpaceDim);
  FORT_SYMTEIGEN(CHF_FRA(lambda), 
		 CHF_CONST_FRA(ee),
		 CHF_INT(dudxComp),
		 CHF_INT(dudyComp),
		 CHF_INT(dvdxComp),
		 CHF_INT(dvdyComp),
		 CHF_BOX(a_box));

	 
  //principal stress components
  lambda.mult(a_mu,0,0,1); // only need the first principal stress 
  
  //vertical integration (FIXME what do in future 3D problems?)
  lambda.mult(a_thck,0,0,1);
  
  
 
  Real eps = 1.0e-10; // minimum denominator
  
  FArrayBox thckab(a_box,1);
  FORT_THICKNESSOVERFLOTATION(CHF_FRA1(thckab,0),
			      CHF_CONST_FRA1(a_thck,0),
			      CHF_CONST_FRA1(a_topg,0),
			      CHF_CONST_REAL(a_rhoi),
			      CHF_CONST_REAL(a_rhoo),
			      CHF_CONST_REAL(a_sealevel),
			      CHF_BOX(a_box));
  
  Real maxDamage = DAMAGE_MAX_VAL;

  FORT_NYECREVASSEDEPTHTP(CHF_FRA1(a_damage,0), 
		       CHF_CONST_FRA1(a_water,0),
		       CHF_CONST_FRA1(a_thck,0),
		       CHF_CONST_FRA1(lambda,0),
		       CHF_CONST_FRA1(thckab,0),
		       CHF_CONST_REAL(a_rhoi),
		       CHF_CONST_REAL(a_rhoo),
		       CHF_CONST_REAL(a_gravity),
		       CHF_CONST_REAL(eps),
		       CHF_CONST_REAL(maxDamage),
		       CHF_BOX(a_box));
  // {
  //   //sanity check...
  //   FArrayBox coef(a_box,1); coef.copy(a_damage); coef /= a_thck; coef -= 1.0; coef *= -1.0;
  //   lambda.mult(coef,0,0,1);
  //   FArrayBox diff(a_box,1);
  //   FORT_NYECREVASSEDEPTHT(CHF_FRA1(diff,0), 
  //   			   CHF_CONST_FRA1(depthw,0),
  //   			   CHF_CONST_FRA1(a_thck,0),
  //   			   CHF_CONST_FRA1(lambda,0),
  //   			   CHF_CONST_FRA1(thckab,0),
  //   			   CHF_CONST_REAL(a_rhoi),
  //   			   CHF_CONST_REAL(a_rhoo),
  //   			   CHF_CONST_REAL(a_gravity),
  //   			   CHF_CONST_REAL(eps),
  //   			   CHF_CONST_REAL(maxDamage),
  //   			   CHF_BOX(a_box));

  //   diff -= a_damage; diff*=a_thck;
  //   Real err = diff.norm(0);
  //   if (err > 1.0e-5)
  //     pout() << " err = " << diff.norm(0) << std::endl;
  // }

}

void  
DamageConstitutiveRelation::computeFaceMu
(LevelData<FluxBox>& a_mu,
 LevelData<FArrayBox>& a_velocity, const Real& a_scale,
 const LevelData<FArrayBox>* a_crseVelPtr,
 int a_nRefCrse,
 const LevelData<FluxBox>& a_A,
 const LevelSigmaCS& a_coordSys,
 const ProblemDomain& a_domain,
 const IntVect& a_ghostVect) const
{
  CH_TIME("DamageConstitutiveRelation::computeFaceMu");
  const DisjointBoxLayout& grids = a_mu.disjointBoxLayout();

  int lev = level(a_domain);
  LevelData<FluxBox>& a_damage = (*m_damageModel->m_faceMinDamage[lev]);
  computeFaceMuDamage(a_mu,a_damage,a_velocity, a_scale, a_crseVelPtr,
		      a_nRefCrse,a_A,a_coordSys,a_domain,a_ghostVect);
}

void  
DamageConstitutiveRelation::computeFaceMuDamage
(LevelData<FluxBox>& a_mu,
 LevelData<FluxBox>& a_damage,
 LevelData<FArrayBox>& a_velocity, const Real& a_scale,
 const LevelData<FArrayBox>* a_crseVelPtr,
 int a_nRefCrse,
 const LevelData<FluxBox>& a_A,
 const LevelSigmaCS& a_coordSys,
 const ProblemDomain& a_domain,
 const IntVect& a_ghostVect) const
{
  CH_TIME("DamageConstitutiveRelation::computeFaceMuDamage");
  

  const DisjointBoxLayout& grids = a_mu.disjointBoxLayout();
  LevelData<FluxBox> epsSqr(grids, 1, a_ghostVect);
  LevelData<FluxBox> gradU(grids, SpaceDim*SpaceDim, a_ghostVect);
  LevelData<FluxBox> faceTopography(grids, 1, a_ghostVect);

  //we need face-centered topography to compute face-centered thickness above flotation
  CellToEdge(a_coordSys.getTopography(),  faceTopography);

  computeStrainRateInvariantFace
    (epsSqr,gradU,a_velocity, a_crseVelPtr, a_nRefCrse, a_coordSys,
     a_ghostVect);
  
  m_undamagedConstitutiveRelation->computeFaceMu
    ( a_mu, a_velocity, a_scale,  a_crseVelPtr, a_nRefCrse, a_A, a_coordSys, 
      a_domain, a_ghostVect);

  //which level are we on?
  int lev = level(a_domain);
 
  LevelData<FluxBox> faceDamage (grids, 1, a_ghostVect);
  //FIXME also a bit suspect : face-averaging the cell-centred damage when it might 
  //be better to use the face-centered values that get computed for advection
  CellToEdge(*m_damageModel->damage(lev), faceDamage); 

   LevelData<FluxBox> faceWater (grids, 1, a_ghostVect);
  //FIXME also a bit suspect : face-averaging the cell-centred damage when it might 
  //be better to use the face-centered values that get computed for advection
   CellToEdge(*m_damageModel->water(lev), faceWater); 

  
  for (DataIterator dit(grids);dit.ok();++dit)
    {
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
	{
	  Box box = grids[dit];
	  box.surroundingNodes(faceDir);

	  computeLocalDamage(a_damage[dit][faceDir], a_mu[dit][faceDir], gradU[dit][faceDir],
			     a_coordSys.getFaceH()[dit][faceDir], faceTopography[dit][faceDir],
			     faceWater[dit][faceDir],
			     a_coordSys.iceDensity(), a_coordSys.waterDensity(), 
			     a_coordSys.gravity(), a_coordSys.seaLevel(), box);


	  modifyMu(a_damage[dit][faceDir], a_mu[dit][faceDir], faceDamage[dit][faceDir], 
			  a_coordSys.getFaceH()[dit][faceDir], box);

	} //end loop over face directions
    } //end loop over boxes
}



#include "NamespaceFooter.H"
