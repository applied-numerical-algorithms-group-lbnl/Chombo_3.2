#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _ICEVELOCITYSOLVER_H_
#define _ICEVELOCITYSOLVER_H_

#include "ConstitutiveRelation.H"
#include "BasalFrictionRelation.H"
#include "RealVect.H"
#include "AMRMultiGrid.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "IceThicknessIBC.H"


#include "NamespaceHeader.H"

/// Abstract class to manage the nonlinear solve for the ice-sheet momentum 
/** This is a virtual base class to facilitate the implementation of
    different solution algorithms.
*/
class IceVelocitySolver
{
public:
  
  IceVelocitySolver() {m_isDefined = false;}
  
  virtual ~IceVelocitySolver() {;}

  /// (Re-)define IceVelocitySolver object for a given mesh and set of physics rules
  /**
     @param a_coarseDomain domain boundaries on the coarsest level
     @param a_constRel physics rule (e.g Glen's flow law) relating englacial stress and strian-rate
     @param a_FrictionRel physics rule (e.g Weertman law) relating basal stress and strain-rate
     @param a_vectGrids mesh layout for each level
     @param a_vectGrids refinement ratio for each level
     @param a_dxCrse mesh spacing on coarsest level
     @param a_bc ice thickness and velocity boundary conditions
     @param a_numLevels number of mesh levels.
   */
  virtual void define(const ProblemDomain& a_coarseDomain,
                      ConstitutiveRelation* a_constRel,
		      BasalFrictionRelation* a_FrictionRel,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_vectRefRatio,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels) = 0;
  
  
  /// Compute the horizontal, cell centered ice velocity ( u(x,y), v(x,y) ) [VERTICALLY INTEGRATED MODELS!]
  /** 
      @param   a_horizontalVel on input an initial guess ( u_0(x,y), v_0(x,y) ), on output ( u(x,y), v(x,y) )
      @param   a_initialResidualNorm  residual norm given the initial guess
      @param   a_finalResidualNorm  residual norm at the end of the calculation
      @param   a_convergenceMetric desired a_finalResidualNorm/a_initialResidualNorm 
      @param   a_rhs right hand side of the stress balance equation: the driving stress (x,y) components  
      @param   a_C main basal friction coefficient : \f$ \tau _b = f_b(C,u) + C_0 u \f$ 
      @param   a_C0 auxilliary basal friction coefficient : \f$ \tau _b =f_b(C,u) + C_0 u \f$  
      @param   a_A flow rate factor: one component per layer
      @param   a_muCoef stiffness factor: one component
      @param   a_coordSys ice sheet geoemtry
      @param   a_time current time
      @param   a_lbase lowest mesh level (probably doesn't work!)
      @param   a_maxLevel highest mesh level
      @return  0 if solver converges successfully
   */
  virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		    Vector<LevelData<FArrayBox>* >& a_calvedIce,
		    Vector<LevelData<FArrayBox>* >& a_addedIce,
		    Vector<LevelData<FArrayBox>* >& a_removedIce,
		    Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		    const Real a_convergenceMetric,
                    const Vector<LevelData<FArrayBox>* >& a_rhs,
                    const Vector<LevelData<FArrayBox>* >& a_C,
		    const Vector<LevelData<FArrayBox>* >& a_C0,
                    const Vector<LevelData<FArrayBox>* >& a_A,
		    const Vector<LevelData<FArrayBox>* >& a_muCoef,
                    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                    Real a_time,
                    int a_lbase, int a_maxLevel) = 0;

  virtual void setVerbosity(int a_verbosity) {m_verbosity = a_verbosity;}

  /// if relevant, set max number of solver iterations
  /** default implementation does nothing -- it's assumed that 
      derived classes will implement this in a meaningful way 
  */
  virtual void setMaxIterations(int a_max_iter) {;}

  /// if relevant, set solver tolerance
  /** default implementation does nothing -- it's assumed that 
      derived classes will implement this in a meaningful way 
  */
  virtual void setTolerance(Real a_tolerance) {;}

  virtual bool isDefined() {return m_isDefined;}


protected:

  bool m_isDefined;
  int m_verbosity;

};

#include "NamespaceFooter.H"

#endif
