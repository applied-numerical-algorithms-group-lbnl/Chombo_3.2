#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _KNOWNVELSOLVER_H_
#define _KNOWNVELSOLVER_H_

#include "IceVelocitySolver.H"
#include "RealFunction.H"
#include "NamespaceHeader.H"
/// Class which implements  IceVelocitySolver for a putative known
/// velocity field 
/// Intended to allow simple tests of (say) temperature advection
class KnownVelocitySolver : public IceVelocitySolver
{

  // Real m_wavenumber, m_phase, m_amplitude;
  //RealVect m_u0,m_k, m_length;
  Vector<RealVect> m_dx;
  RealFunction<RealVect>* m_velFunctionPtr;
  
public:
  
  ~KnownVelocitySolver()
  {
    if (m_velFunctionPtr != NULL)
      delete m_velFunctionPtr;
  }

   virtual void define(const ProblemDomain& a_coarseDomain,
		       ConstitutiveRelation* a_constRel,
		       BasalFrictionRelation* a_basalFrictionRel,
		       const Vector<DisjointBoxLayout>& a_vectGrids,
		       const Vector<int>& a_vectRefRatio,
		       const RealVect& a_dxCrse,
		       IceThicknessIBC* a_bc,
		       int a_numLevels);

  //full solve
  virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
		    Vector<LevelData<FArrayBox>* >& a_calvedIce,
		    Vector<LevelData<FArrayBox>* >& a_addedIce,
		    Vector<LevelData<FArrayBox>* >& a_removedIce,
		    Real& a_initialResidualNorm, Real& a_finalResidualNorm,
		    const Real a_convergenceMetric,
                    const Vector<LevelData<FArrayBox>* >& a_rhs,
                    const Vector<LevelData<FArrayBox>* >& a_beta,
		    const Vector<LevelData<FArrayBox>* >& a_beta0,
                    const Vector<LevelData<FArrayBox>* >& a_A,
		    const Vector<LevelData<FArrayBox>* >& a_muCoef,
                    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                    Real a_time,
                    int a_lbase, int a_maxLevel) ;

};
#include "NamespaceFooter.H"
#endif
