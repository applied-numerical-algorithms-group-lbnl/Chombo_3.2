#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _REGRESSION_H_
#define _REGRESSION_H_

#include "ParmParse.H"
#include "AMRIO.H"
#include "computeNorm.H" 
//#include "Vector.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5

/// Compute the maximum difference between (a subset of) the AMR data stored in two HDF5 fiels
/**
   @param a_result result hdf5 file
   @param a_ref  reference hdf5 file
   @param a_name component names to compare: if empty compare all components
 */
Real HDF5NormTest(std::string a_result, std::string a_ref)
{

  Vector<LevelData<FArrayBox>* > resultSoln;
  Vector<string> resultVars; // result solution variable names
  Vector<DisjointBoxLayout> resultGrids;
  Box resultDomain;
  Real resultDx, resultDt, resultTime;
  Vector<int> resultRefRatio;
  int resultNumLevels;
  IntVect ghostVect = IntVect::Unit;

  ReadAMRHierarchyHDF5(a_result,
		       resultGrids,
		       resultSoln,
		       resultVars,
		       resultDomain,
		       resultDx,
		       resultDt,
		       resultTime,
		       resultRefRatio,
		       resultNumLevels);

  
  Vector<LevelData<FArrayBox>* > refSoln;
  Vector<string> refVars; // ref soln variable names
  Vector<DisjointBoxLayout> refGrids;
  Box refDomain;
  Real refDx, refDt, refTime;
  Vector<int> refRefRatio;
  int refNumLevels;

  ReadAMRHierarchyHDF5(a_ref,
		       refGrids,
		       refSoln,
		       refVars,
		       refDomain,
		       refDx,
		       refDt,
		       refTime,
		       refRefRatio,
		       refNumLevels);

  CH_assert (refDomain == resultDomain);
  CH_assert (refNumLevels == resultNumLevels);

  Real maxnorm = 0.0;
  for (int j = 0; j < resultSoln[0]->nComp() ; j++)
    {
      int k = 0;
      while (k < refVars.size() && resultVars[j] !=  refVars[k])
	{
	  k++;
	}
      
      if ( k >= refVars.size() )
	{
	  std::string msg = " var " + resultVars[j] + " missing from " + a_ref;
	  MayDay::Warning(msg.c_str());
	}
      else
	{
	  for (int lev = 0; lev < resultNumLevels; lev++)
	    {
	      LevelData<FArrayBox> ref(resultGrids[lev], 1, IntVect::Zero);
	      refSoln[lev]->copyTo(Interval(k,k), ref, Interval(0,0));
	      for (DataIterator dit(resultGrids[lev]); dit.ok(); ++dit)
		{
		  (*resultSoln[lev])[dit].minus(ref[dit], 0, j, 1);
		}
	    }
       
	  Real maxnormj = computeNorm(resultSoln, resultRefRatio, resultDx, Interval(j,j), 0);
	  pout() << "HDF5NormTest: maxnorm difference(" << resultVars[j] << ") = " <<  maxnormj << endl;
	  maxnorm = std::max(maxnorm, maxnormj);
	}
    }

  return maxnorm;

}
#endif


#include "NamespaceFooter.H"
#endif
