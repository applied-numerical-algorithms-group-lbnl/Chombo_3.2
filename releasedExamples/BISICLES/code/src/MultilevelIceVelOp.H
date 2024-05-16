#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


// We can assume that template class T has null construction.

#ifndef _MULTILEVELICEVELOP_H_
#define _MULTILEVELICEVELOP_H_

#include <cmath>

#include "RefCountedPtr.H"
#include "Vector.H"
#include "RealVect.H"
#include "LevelData.H"
#include "MultilevelLinearOp.H"

#include "NamespaceHeader.H"

enum bottomSolverTypes {bicg = 0,
                        PETSC,
                        relax,
                        NUM_BS_TYPES};

///
/**
   Derived class from MultilevelLinearOp to enable specific tweaks for
   ice sheets. Primary reason for this is to overload define to allow for 
   petsc solvers.
 */
class MultilevelIceVelOp : public MultilevelLinearOp<FArrayBox>
{

public:
  /// default constructor (sets defaults for preconditioner)
  MultilevelIceVelOp();

  /// define function -- opFactory should be able to define all required levels
  /**
  */
  virtual void define(const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_refRatios,
                      const Vector<ProblemDomain>& a_domains,
                      const Vector<RealVect>& a_vectDx,
                      RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >& a_opFactory,
                      int a_lBase);


  ///
  virtual ~MultilevelIceVelOp();

  /// 
  int m_bottom_solver_type;
  
};

#include "NamespaceFooter.H"

#endif /*_MULTILEVELICEVELOP_H_*/
