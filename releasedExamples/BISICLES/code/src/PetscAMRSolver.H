#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCAMRSOLVERSOLVER_H_
#define _PETSCAMRSOLVERSOLVER_H_

#include "PetscCompGridVTO.H"
#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#endif 

#include "NamespaceHeader.H"

class LinearizedVTOp;
// solver class for VTO operators
class PetscAMRSolver
{
   static PetscErrorCode ksp_monitor_pout(KSP ksp, PetscInt it, PetscReal rnorm  ,void *ctx)
  {
    pout() << "      KSP:: iteration = " << it << " residual norm = " <<  rnorm << std::endl;
    return 0;
  }

public:
  PetscAMRSolver(int a_verb=0);

  virtual ~PetscAMRSolver() {;}
  ///
  /**
     Solve a_op*a_phi = a_rhs using the PETSC matrix free functions
     The preconditioner (for which a matrix is formed) need not be 
     the same as the actual operator.
  */
  virtual PetscErrorCode solve_mfree( Vector<LevelData<FArrayBox>*>& a_phi, 
				      const Vector<LevelData<FArrayBox>*>& a_rhs, 
				      LinearizedVTOp *a_op );

  virtual PetscErrorCode solve( Vector<LevelData<FArrayBox>*>& a_phi, 
				const Vector<LevelData<FArrayBox>*>& a_rhs );
  
  static PetscErrorCode apply_mfree(Mat A, Vec x, Vec f);  

  void plot(const string, const Vector<LevelData<FArrayBox>*>& a_phi);
  PetscErrorCode setCoords(KSP, const Vector<LevelData<FArrayBox>*>& a_phi);

  LinearizedVTOp *m_op_mfree;
  Vector<LevelData<FArrayBox>*> m_phi_mfree;
  Vector<LevelData<FArrayBox>*> m_Lphi_mfree;
  bool m_mfree_homogeneous;
  PetscCompGridVTO m_petscCompMat; // this is specific for VTO for now, should make PetscCompGrid*
  int m_verbose;
};

#include "NamespaceFooter.H"
#endif
