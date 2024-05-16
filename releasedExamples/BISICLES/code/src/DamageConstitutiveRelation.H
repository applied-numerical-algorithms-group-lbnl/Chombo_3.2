#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _DAMAGECONSTITUTIVERELATION_H_
#define _DAMAGECONSTITUTIVERELATION_H_

#include "ConstitutiveRelation.H"
#include "AMRDamage.H"
#include "NamespaceHeader.H"
/// Computes the effective viscosity \f$\mu\f$ including a factor due to crevassing
/**
   The effective viscosity \f$\mu\f$ is computed by multiplying
   an effective viscosity \f$\mu_0\f$ as calculated by another subclass
   of ConstitutiveRelation by a factor \f$(1-d)/H\f$, where \f$d\f$
   is the total depth of crevasses (surface plus base). \f$d\f$ is
   set to \f$d = max ( d_{\rm NYE} + b d_{\rm OLD}) \f$. \f$d_{\rm NYE}\f$
   is the Nye crevasse depth, and \f$d_{\rm OLD}\f$ is an 'old' crevasse
   depth, computed by a damage transport model

 */

class DamageConstitutiveRelation : public ConstitutiveRelation
{
  friend class AMRDamage;

  ConstitutiveRelation* m_undamagedConstitutiveRelation;
  AMRDamage* m_damageModel;

public:
  
  DamageConstitutiveRelation(ConstitutiveRelation* a_ptr, AMRDamage* a_damage)
  {
    m_undamagedConstitutiveRelation = a_ptr->getNewConstitutiveRelation();
    m_damageModel = a_damage;
  }

  virtual ~DamageConstitutiveRelation();
  
  /// computes cell-centered \f$\mu_{AS}\f$ based on the cell-centered velocity
  /**
     a_mu -- \f$\mu_{AS}\f$ based on the local velocity field.
     a_vel -- Cell-centered velocity field.
     a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
     (NULL if no coarser level)
     a_nRefCrse -- refinement ratio to next coarser level
     a_A: Cell-centered flow law coefficient (glenn's A) field
     a_coordSys:  SigmaCS object containing the geometry of this patch.
     a_ghostVect -- how the boxes over which we want to compute mu relate 
                    to those in the DisjointBoxLayout (can be negative)
  */
  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel, const Real& a_scale,
                         const LevelData<FArrayBox>* a_crseVel,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const;
  
  /// Compute a cell centred bulk dissipation 
  /// \f$\Phi/(\rho _i c _i) = \sigma_ij \epsilon _ji /(\rho _i c _i) \f$ 
  /// (heat source) at the cell centres. This ought to have the same dependencies 
  /// as cell centred mu
  /**
     a_dissipation -- \f$\Phi\f$ based on the local velocity field.
     a_vel -- Cell-centered velocity field.
     a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
     a_nRefCrse -- refinement ratio to next coarser level
     a_A: Cell-centered flow law coefficient (glenn's A) field
     a_coordSys:  SigmaCS object containing the geometry of this patch.
     a_box: cell-centered box over which to do this computation
  **/
  virtual void computeDissipation(LevelData<FArrayBox>& a_dissipation,
				  const LevelData<FArrayBox>& a_vel, 
				  const LevelData<FArrayBox>* a_crseVel,
				  int a_nRefCrse,
				  const LevelData<FArrayBox>& a_A,
				  const LevelSigmaCS& a_coordSys,
				  const ProblemDomain& a_domain,
				  const IntVect& a_ghostVect = IntVect::Zero) const;

  // computes face-centered \f$\mu_{AS}\f$ based on cell-centered velocity
  /** a_mu: face-centered \f$\mu_{AS}\f$ based on the local velocity field.
      a_vel: Cell-centered velocity field.
      a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
      a_nRefCrse -- refinement ratio to next coarser level
      a_A: Cell-centered flow law coefficient (glenn's A) field
      a_coordSys: SigmaCS object containing the geometry of this patch.
      a_ghostVect -- how the boxes over which we want to compute mu relate 
                     to those in the DisjointBoxLayout (can be negative) 
  */  
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel, const Real& a_scale,
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse,
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const;


 


  /// creates a new copy of this ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;

private:

  DamageConstitutiveRelation() : m_undamagedConstitutiveRelation(NULL)
				      
  {;}



  int level(const ProblemDomain&  a_domain) const
  {
    CH_assert(m_damageModel != NULL); // for now, we will say there *has* to be a damage model
    //FIXME this is a bit klugey... we need to work out access various levels from
    // m_damageModel but we typically don't know which level that we are (no level argument), so for
    //now (maybe we will change the ConstitutiveRelation interface later) 
    //work it out by comparing domain boxes
    int lev = 0;
    while ( a_domain.domainBox() != 
	    m_damageModel->damage(lev)->disjointBoxLayout().physDomain().domainBox())
      {
	lev++;
	CH_assert( lev <= m_damageModel->finestLevel());
      }
    return lev;
  }

  /// modify the effective viscosity and the damage to reflect local and transported damage
  static void  modifyMu
  (FArrayBox& a_damage,
   FArrayBox& a_mu,
   const FArrayBox& a_transportedDamage, 
   const FArrayBox& a_thck,
   const Box& a_box) ;

  /// Compute the local  damage given strain rates and undamaged effective viscosity
  static void  computeLocalDamage 
  (FArrayBox& a_damage,
   const FArrayBox& a_mu,
   const FArrayBox& a_gradU, 
   const FArrayBox& a_thck,
   const FArrayBox& a_topg,
   const FArrayBox& a_water,
   const Real& a_rhoi,
   const Real& a_rhow,
   const Real& a_gravity,
   const Real& a_sealevel,
   const Box& a_box);

   /// Compute the local  damage given the vertically integrated deviatoric stress (viscous tensor)
  static void  computeLocalDamageVT
  (FArrayBox& a_damage,
   const FArrayBox& a_vt,
   const FArrayBox& a_thck,
   const FArrayBox& a_topg,
   const FArrayBox& a_water,
   const Real& a_rhoi,
   const Real& a_rhow,
   const Real& a_gravity,
   const Real& a_sealevel,
   const Box& a_box);

  
  // computes face-centered effective viscoity  and a damage coeffcient based on cell-centered velocity
  /** 
      a_mu includes the damage, but certain users (AMRDamage) need to know the damage too

      a_mu: face-centered \f$\mu_{AS}\f$ based on the local velocity field.
      a_vel: Cell-centered velocity field.
      a_crseVel -- coarse-level velocity field (for coarse-fine bc's).
                         (NULL if no coarser level)
      a_nRefCrse -- refinement ratio to next coarser level
      a_A: Cell-centered flow law coefficient (glenn's A) field
      a_coordSys: SigmaCS object containing the geometry of this patch.
      a_ghostVect -- how the boxes over which we want to compute mu relate 
                     to those in the DisjointBoxLayout (can be negative) 
  */  
  virtual void computeFaceMuDamage(LevelData<FluxBox>& a_mu,
				   LevelData<FluxBox>& a_damage,
				   LevelData<FArrayBox>& a_vel, const Real& a_scale,
				   const LevelData<FArrayBox>* a_crseVel,
				   int a_nRefCrse,
				   const LevelData<FluxBox>& a_A, 
				   const LevelSigmaCS& a_coordSys,
				   const ProblemDomain& a_domain,
				   const IntVect& a_ghostVect = IntVect::Zero) const;


};

#include "NamespaceFooter.H"
#endif
