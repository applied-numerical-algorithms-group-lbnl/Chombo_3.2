#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CONSTITUTIVERELATION_H_
#define _CONSTITUTIVERELATION_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelSigmaCS.H"
#include "CellToEdge.H"
#include "NamespaceHeader.H"



/// Abstract class around the englacial constitutive relations for ice
/** 
    Specifies an interface for the calculation of the effecive viscosity
    \f$\mu\f$ such that (stress = \f$\mu\f$*strainrate ) assuming that \f$\mu = \mu(\nabla(\vec{u}),A(T)) \f$. 
    
    Also specifies (and provides a default implementation) for the default rate of strain heating (mu * stress * strainrate)

    Run-time configuration is carried out in ConstitutiveRelation::parse.
  

*/
class 
ConstitutiveRelation
{

public:

  static ConstitutiveRelation* parse(const char* a_prefix);

  ConstitutiveRelation() {;}

  virtual ~ConstitutiveRelation() {;}
 
  /// Compute cell-centered \f$\mu\f$ based on the cell-centered velocity
  /** 
      
     \param a_mu \f$\mu\f$ based on the local velocity field. 
     \param a_vel Cell-centered velocity field.
     \param a_crseVel Next coarser level velocity field (for coarse-fine bc's). (NULL if no coarser level)
     \param a_nRefCrse refinement ratio to next coarser level
     \param a_A Cell- and layer-centered flow law coefficient (Glenn's \f$ A \f$ )  field
     \param a_coordSys: LevelSigmaCS object containing the ice-sheet geometry
     \param a_ghostVect: how the boxes on which we want to compute mu compare
                         with the DisjointBoxLayout. (can be negative)
  */
  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel, const Real& a_scale,
                         const LevelData<FArrayBox>* a_crseVel,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const = 0;

  

  /**

     \brief Compute a cell centred bulk dissipation \f$\Phi/(\rho _i c _i) = \tau_{ij} \epsilon _{ji} /(\rho _i c _i) \f$ 
     (heat source) at the cell centres. This ought to have the same dependencies  as cell centred mu

     \param a_dissipation  \f$\Phi\f$ based on the local velocity field.
     \param a_vel Cell-centered velocity field.
     \param a_crseVel Next coarser level velocity field for coarse-fine bc's (NULL if no coarser level)
     \param a_nRefCrse refinement ratio to next coarser level
     \param a_A Cell- and layer-centered flow law coefficient (Glenn's \f$ A \f$ )  field
     \param a_coordSys  LevelSigmaCS object containing the geometry of this level.
     \param a_box cell-centered box over which to do this computation
  **/
  virtual void computeDissipation(LevelData<FArrayBox>& a_dissipation,
				  const LevelData<FArrayBox>& a_cellVel,
				  const LevelData<FArrayBox>* a_crseVel,
				  int nRefCrse,
				  const LevelData<FArrayBox>& a_A,
				  const LevelSigmaCS& a_coordSys,
				  const ProblemDomain& a_domain,
				  const IntVect& a_ghostVect = IntVect::Zero) const = 0; 

  /// compute face-centered effective viscosity \f$\mu\f$ based on cell-centered velocity
  /** 
      \param a_mu face-centered \f$\mu\f$ based on the local velocity field.
      \param a_vel Cell-centered velocity field.
      \param a_scale compute mu given input a_vel  = m_scale * true_vel
      \param a_A Cell- and layer-centered flow law coefficient (Glenn's \f$ A \f$ )  field
      \param a_crseVel Next coarser level velocity field for coarse-fine bc's (NULL if no coarser level)
                   
      \param a_nRefCrse refinement ratio to next coarser level
      \param a_coordSys LevelSigmaCS object containing the geometry of this level
      \param a_ghostVect how the boxes on which we want to compute mu compare
                         with the DisjointBoxLayout. (can be negative)
  */  
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel, const Real& a_scale, 
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse,
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const = 0;


  /// creates a new copy of this  ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const = 0;

  /// compute cell-centered strainrate invariant epsilon^2
  void computeStrainRateInvariant(LevelData<FArrayBox>& a_epsilonSquared,
                                  const LevelData<FArrayBox>& a_velocity,
                                  const LevelData<FArrayBox>* a_crseVel,
                                  int nRefCrse,
                                  const LevelSigmaCS& a_coordSys,
                                  const IntVect& a_ghostVect = IntVect::Zero) const;

  /// compute face-centered strainrate invariant epsilon^2 based on cell-centered velocity
  void computeStrainRateInvariantFace(LevelData<FluxBox>& a_epsilonSquared,
                                      LevelData<FArrayBox>& a_velocity,
                                      const LevelData<FArrayBox>* a_crseVel,
                                      int a_nRefCrse,
                                      const LevelSigmaCS& a_coordSys,
                                      const IntVect& a_ghostVect = IntVect::Zero) const;  

  /// compute cell-centered strainrate invariant epsilon^2 and velocity gradients. No  more costly than just computing epsilon^2
  void computeStrainRateInvariant(LevelData<FArrayBox>& a_epsilonSquared,
				  LevelData<FArrayBox>& a_gradVelocity,
                                  const LevelData<FArrayBox>& a_velocity,
                                  const LevelData<FArrayBox>* a_crseVel,
                                  int a_nRefCrse,
                                  const LevelSigmaCS& a_coordSys,
                                  const IntVect& a_ghostVect = IntVect::Zero) const;

  /// compute face-centered strainrate invariant epsilon^2 and velocity gradients based on cell-centered velocity. No more costly than just computing epsilon^2
  void computeStrainRateInvariantFace(LevelData<FluxBox>& a_epsilonSquared,
				      LevelData<FluxBox>& a_gradVelocity,
                                      LevelData<FArrayBox>& a_velocity,
                                      const LevelData<FArrayBox>* a_crseVel,
                                      int a_nRefCrse,
                                      const LevelSigmaCS& a_coordSys,
                                      const IntVect& a_ghostVect = IntVect::Zero) const;  

protected:

};
 
/// rate factor A(T) in (e.g) Glen's law  
class RateFactor 
{

public:
  
  virtual ~RateFactor(){;}

  virtual void computeA(FArrayBox& a_A, 
			const FArrayBox& a_thetaStar,
			const FArrayBox& a_pressure,
			const Box& a_box) const = 0;

  virtual RateFactor* getNewRateFactor() const = 0;


} ;

///Constant Rate Factor
/**
   Allows A to be specified as a constant, useful for isothermal
   ice sheets.
*/
class ConstantRateFactor : public RateFactor  
{
  Real m_A; // the constant rate factor
public:
  
  // ConstantRateFactor(Real a_seconds_per_unit_time)
  //   : m_A(9.2e-18 / a_seconds_per_unit_time ) // A measured in Pa^{n} a^{-1}
  // {
  // }
 
  ConstantRateFactor(Real a_A) : m_A(a_A)
  {
  }


  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_theta, 
		const FArrayBox& a_pressure,
		const Box& a_box) const
  {
    a_A.setVal(m_A, a_box, 0);
  }

  RateFactor* getNewRateFactor() const 
  {
    ConstantRateFactor* newPtr = new ConstantRateFactor(m_A);
    return static_cast<RateFactor*>(newPtr);
  }


};

///Arrhenius Rate Factor
/**
   The Arrhenius rate factor is used widely in ice sheet modelling, e.g 
   (Pattyn, 2003)
 */
class ArrheniusRateFactor : public RateFactor  
{
  /// power law exponent
  Real m_n;
  /// Pattyn's `"enhancement factor"
  Real m_enhance;
  /// flow rate factor
  Real m_B0;
  /// limit temperature in flow-rate factor
  Real m_theta_r;
  /// flow rate exponent
  Real m_K;
  /// flow rate factor
  Real m_C;
  /// universal gas constant
  Real m_R;
  /// activation energy for creep 
  // (insert funny comment about lazy creeps here...)
  Real m_Q;

public:
  ArrheniusRateFactor(Real a_seconds_per_unit_time);
 
  void setDefaultParameters(Real a_seconds_per_unit_time);
 
  
  void setParameters(Real a_n,
		     Real a_enhance,
		     Real a_B0,
		     Real a_theta_r,
		     Real a_K,
		     Real a_C,
		     Real a_R,
		     Real a_Q);

  //compute the (temperature-dependent) rate factor A
  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_thetaStar,
		const FArrayBox& a_pressure,
		const Box& a_box) const;

  RateFactor* getNewRateFactor() const; 
  

};

/// Paterson Rate Factor
/**
   Another  Arrhenius law, this one from
   Cuffey and Paterson,  2010, 'The physics of glaciers', p72-74
   (a revision of earlier books by Paterson) 

   A = E*A0*exp(-Qc / R * [1/T* - 1/T0*]

   Th is the usual pressure-corrected temperature
   Th = T + 7e−8 * P
   T0* = T0  + 7e−8 * P (T0 = 263, P = pressure)
   Qc = Qm = if Th < T∗; Qc = Qp if Th > T∗.
  
   Defaults: 

   A0 = 3.5e-25 Pa−3 s−1 * secondsperyear
   Qp = 115 kJ mol−1
   Qm = 6e4 kJ mol−1
   R=8.314 J

 */
class PatersonRateFactor : public RateFactor  
{
  
  /// enhancement factor
  Real m_E;
  /// flow rate factor
  Real m_A0;
  /// limit temperature in flow-rate factor
  Real m_T0;
  /// universal gas constant
  Real m_R;
  /// activation energies for creep 
  // (insert funny comment about two lazy creeps here...)
  Real m_Qm, m_Qp;

public:

  PatersonRateFactor(Real a_seconds_per_unit_time);
 
  void setDefaultParameters(Real a_seconds_per_unit_time);
 
  void setA0(Real a_A0)
  {
    m_A0 = a_A0;
  }

  void setParameters(Real a_E, Real a_A0, Real a_T0, Real a_R, 
		     Real a_Qm, Real a_Qp);

  //compute the (pressure-corrected temperature- 
  //and pressure- dependent) rate factor A
  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_thetaPC, 
		const FArrayBox& a_pressure,
		const Box& a_box) const;

  RateFactor* getNewRateFactor() const; 
  

};


/// Zwinger Rate Factor
/**
   Another  Arrhenius law, this one from
   Zwinger (2007) Annals of Glaciology vol 45 p 29
   and used in (older versions of?) Elmer/Ice. 

   A = E*A0(T)*exp(-Qc / (R * T) ]
   where Qc = Qm = if theta < theta0; Qc = Qp otherwise.
   and   A0 = Am = if theta < theta0; A0 = Ap otherwise.
   Am is computed so that A is continuous with T 

   Th is the usual pressure-corrected temperature
   Th = T + 7e−8 * P
   T0* = T0  + 7e−8 * P (T0 = 263, P = pressure)
   Qc = Qm = if Th < T∗; Qc = Qp if Th > T∗.
  
   Defaults: 

   A0 = 1.916e3 Pa−3 s−1 * secondsperyear
   Qp = 139 kJ mol−1
   Qm = 6e4 kJ mol−1
   R=8.314 J

 */
class ZwingerRateFactor : public RateFactor  
{
  
  /// enhancement factor
  Real m_E;
  /// flow rate factor
  Real m_A0;
  /// limit temperature in flow-rate factor
  Real m_T0;
  /// universal gas constant
  Real m_R;
  /// activation energies for creep 
  // (insert funny comment about two lazy creeps here...)
  Real m_Qm, m_Qp;

public:

  ZwingerRateFactor(Real a_seconds_per_unit_time);
 
  void setDefaultParameters(Real a_seconds_per_unit_time);
 

  void setParameters(Real a_E, Real a_A0, Real a_T0, Real a_R, 
		     Real a_Qm, Real a_Qp);

  //compute the (pressure-corrected temperature- 
  //and pressure- dependent) rate factor A
  void computeA(FArrayBox& a_A, 
		const FArrayBox& a_thetaPC, 
		const FArrayBox& a_pressure,
		const Box& a_box) const;

  RateFactor* getNewRateFactor() const; 
  

};

// $\f \tau^2 \tau _ij = 2 A \dot{\epsilon}_{ij} \f$ 


/// Implements  ConstitutiveRelation to give Glen's flow law.
/** 
    Glen's flow law. Computes an effective viscosity from strain-rate \f$ \dot{\epsilon}_{ij} \f$ 
    such that the stress \f$ \tau _ij \f$ satifies
    
    $\f \tau^2 \tau _ij = 2 A \cdot {\epsilon}_{ij} \f$ 
    
    Run-time configuration is carried out in ConstitutiveRelation::parse
    A  GlensFlowRelation object will be created if the string
    `A.constitutiveRelation =  GlensLaw` is found. Options are:
    
     - GlensLaw.n = <float> - the power in Glen's flow law, default (3.0)
     - GlensLaw.epsSqr0 = <float> - regularization at low strain rates. Default 10^-12
     - GlensLaw.delta = <float> - regularization at high strain rates. Default 0.0

*/
class GlensFlowRelation: public ConstitutiveRelation
{
public:

  GlensFlowRelation();

  virtual ~GlensFlowRelation();



  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                         const LevelData<FArrayBox>* a_crseVel,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const;
  
  virtual void computeDissipation(LevelData<FArrayBox>& a_dissipation,
				  const LevelData<FArrayBox>& a_vel, 
				  const LevelData<FArrayBox>* a_crseVel,
				  int a_nRefCrse,
				  const LevelData<FArrayBox>& a_A,
				  const LevelSigmaCS& a_coordSys,
				  const ProblemDomain& a_domain,
				  const IntVect& a_ghostVect = IntVect::Zero) const;

  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel, const Real& a_scale, 
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse,
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const;


  virtual void setDefaultParameters();

  virtual void setParameters(Real a_n,
                             Real a_epsSqr0, 
			     Real a_delta);

  /// creates a new copy of this ConstitutiveRelation object.
  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;
  
  /// Power law exponent
  Real m_n;

  virtual Real power() const {return m_n;}

  /// Rate factor object
  //RateFactor* m_rateFactor;

  /// small numbers to ensure that viscosity remains finite as
  /// strain rates approach 0 or inifinity
  Real m_epsSqr0, m_delta;

protected:
  /// utility function to compute flow law coefficient dependent part of Glens's
  /// law (ie A^{1/n})
  void computeMu0(FArrayBox& a_mu0,
                  const FArrayBox& a_A,const Real& a_scale,
                  const Box& a_box) const;
  

};


/// implementation of constant-mu constitutive relation
/** The constMuRelation class is publicly derived from the 
    ConstitutiveRelation Class and sets mu to be a constant value
*/
class constMuRelation: public ConstitutiveRelation
{
public:
  ///
  constMuRelation();

  virtual ~constMuRelation();


  virtual void computeMu(LevelData<FArrayBox>& a_mu,
                         const LevelData<FArrayBox>& a_vel, const Real& a_scale,
                         const LevelData<FArrayBox>* a_crseVel,
                         int a_nRefCrse,
                         const LevelData<FArrayBox>& a_A,
                         const LevelSigmaCS& a_coordSys,
			 const ProblemDomain& a_domain,
                         const IntVect& a_ghostVect = IntVect::Zero) const;


  virtual void computeDissipation(LevelData<FArrayBox>& a_dissipation,
				  const LevelData<FArrayBox>& a_vel, 
				  const LevelData<FArrayBox>* a_crseVel,
				  int a_nRefCrse,
				  const LevelData<FArrayBox>& a_A,
				  const LevelSigmaCS& a_coordSys,
				  const ProblemDomain& a_domain,
				  const IntVect& a_ghostVect = IntVect::Zero) const;

 
  virtual void computeFaceMu(LevelData<FluxBox>& a_mu,
                             LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                             const LevelData<FArrayBox>* a_crseVel,
                             int a_nRefCrse, 
                             const LevelData<FluxBox>& a_A, 
                             const LevelSigmaCS& a_coordSys,
			     const ProblemDomain& a_domain,
                             const IntVect& a_ghostVect = IntVect::Zero) const;



  virtual void setDefaultParameters();

  virtual void setConstVal(Real a_mu) {m_mu = a_mu;}

  virtual ConstitutiveRelation* getNewConstitutiveRelation() const;

  // constant value of mu
  Real m_mu;

  virtual Real power() const{return 1.0;}

protected:

};



#include "NamespaceFooter.H"

#endif
