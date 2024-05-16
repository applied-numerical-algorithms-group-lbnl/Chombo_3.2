#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "KnownVelocitySolver.H"
#include "BoxIterator.H"
#include "ParmParse.H"


void KnownVelocitySolver::define(const ProblemDomain& a_coarseDomain,
		       ConstitutiveRelation* a_constRel,
		       BasalFrictionRelation* a_basalFrictionRel,
		       const Vector<DisjointBoxLayout>& a_vectGrids,
		       const Vector<int>& a_vectRefRatio,
		       const RealVect& a_dxCrse,
		       IceThicknessIBC* a_bc,
		       int a_numLevels)
{
 
  m_dx.resize(a_numLevels);
  m_dx[0] = a_dxCrse;
  for (int lev = 1; lev < a_numLevels; lev++){
    m_dx[lev] = m_dx[lev-1] / Real(a_vectRefRatio[lev-1]);
  }
  
  //for (int dir = 0; dir < SpaceDim; ++dir)
  //  m_length[dir] = a_coarseDomain.size(dir)*a_dxCrse[dir];

  ParmParse pp("KnownVelocitySolver");

  std::string functionType = "zero";
  pp.query("function_type",functionType);
  if (functionType == "zero")
    {
      m_velFunctionPtr = new ConstantRealFunction<RealVect>(0.0);
    }
  else if (functionType == "flowline")
    {
      Real dx; 
      std::string file, set;
      pp.get("flowline_dx", dx);
      pp.get("flowline_file", file);
      pp.get("flowline_set", set);
      m_velFunctionPtr = new ExtrudedPieceWiseLinearFlowline(file,set,dx);
    } 
  else
    {
      MayDay::Error("KnownVelocitySolver::define() unknown function_type");
    }

}

int KnownVelocitySolver::solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
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
			       int a_lbase, int a_maxLevel)
{

  
  a_initialResidualNorm = a_finalResidualNorm = 0.0;
  for (int lev = 0; lev <  a_maxLevel + 1; lev++)
    {
      const RealVect& dx = m_dx[lev];
      LevelData<FArrayBox>& levelVel = *a_horizontalVel[lev];
      const DisjointBoxLayout& DBL = levelVel.disjointBoxLayout();
      for (DataIterator dit(DBL); dit.ok(); ++dit)
	{
	  FArrayBox& vel = levelVel[dit];
	  vel.setVal(0.0);
	  const Box& box = vel.box(); 
	  for (BoxIterator bit(box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      RealVect x = iv * dx + 0.5 * dx;
	      vel(iv,0) = (*m_velFunctionPtr)(x);
		
	    }
	}
    }
  return 0;
}
