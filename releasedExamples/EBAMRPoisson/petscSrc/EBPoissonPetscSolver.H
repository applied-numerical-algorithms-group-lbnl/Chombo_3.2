#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPoissonPetscSolver_H_
#define _EBPoissonPetscSolver_H_

#include "LinearSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "LevelData.H"
#include "PetscSolver.H"
#include "EBAMRPoissonOp.H"

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"

#include "NamespaceHeader.H"

template <class T>
class PetscSolverEBCF : public PetscSolver<T>
{
public:
  PetscSolverEBCF() : PetscSolver<LevelData<EBCellFAB> >()
  {
  }

  BaseFab<Real>& getRegFab(LevelData<EBCellFAB>& a_fab, const DataIterator& dit)
  {
    BaseFab<Real>& fabref = a_fab[dit()].getSingleValuedFAB();
    return fabref;
  }

  const BaseFab<Real>& getRegFab(const LevelData<EBCellFAB>& a_fab, const DataIterator& dit) const
  {
    const BaseFab<Real>& fabref = a_fab[dit()].getSingleValuedFAB();
    return fabref;
  }


  BaseFab<Real>& getRegFab(LevelData<EBCellFAB>& a_fab, const DataIndex& dit)
  {
    BaseFab<Real>& fabref = a_fab[dit].getSingleValuedFAB();
    return fabref;
  }

  const BaseFab<Real>& getRegFab(const LevelData<EBCellFAB>& a_fab, const DataIndex& dit) const
  {
    const BaseFab<Real>& fabref = a_fab[dit].getSingleValuedFAB();
    return fabref;
  }

  const BaseFab<Real>& getRegFab(const LevelData<EBCellFAB>& a_fab, DataIterator& a_dit, Box& a_box) const
  {
    const BaseFab<Real>& fabref = a_fab[a_dit()].getSingleValuedFAB();
    a_box = fabref.box();
    return fabref;
  }


  const BaseFab<Real>& getRegFab(const LevelData<EBCellFAB>& a_fab, const DataIndex& a_dit, Box& a_box) const
  {
    const BaseFab<Real>& fabref = a_fab[a_dit].getSingleValuedFAB();
    a_box = fabref.box();
    return fabref;
  }
};

class EBPoissonPetscSolver : public PetscSolverEBCF<LevelData<EBCellFAB> >
{
public:
  EBPoissonPetscSolver(): PetscSolverEBCF<LevelData<EBCellFAB> >()
  {
    m_isDefined = false;
  }
  EBPoissonPetscSolver(RefCountedPtr<EBAMRPoissonOp>& a_op, bool a_homogeneous = true)
  : PetscSolverEBCF<LevelData<EBCellFAB> >()
  {
    define(&(*a_op), a_homogeneous);
  }

  EBPoissonPetscSolver(LinearOp<LevelData<EBCellFAB> >* a_op, bool a_homogeneous = true)
  : PetscSolverEBCF<LevelData<EBCellFAB> >()
  {
    define(a_op, a_homogeneous);
  }

#ifdef CH_USE_PETSC
  virtual void solve( LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs );
#endif
  virtual void define(LinearOp<LevelData<EBCellFAB> > *  a_op, bool a_homogeneous = true )
  {
    m_op = dynamic_cast<EBAMRPoissonOp*>(a_op);
    if(m_op == NULL)
      {
        MayDay::Error("cast failed");
      }
    m_isDefined = true;
  }

  virtual void define(RefCountedPtr<EBAMRPoissonOp>& a_op, bool a_homogeneous = true)
  {
    define(&(*a_op), a_homogeneous);
  }

  virtual void define(EBAMRPoissonOp* a_op, bool a_homogeneous = true)
  {
    CH_assert(a_op != NULL);
    m_op = a_op;
    m_isDefined = true;
  }

  /**
     get an estimate of the number of nnz/row for matrix allocation
  */
  virtual int getNNZPerRow() const
  {
    return 1000;
  }

#ifdef CH_USE_PETSC
  virtual int formMatrix(Mat, const LevelData<EBCellFAB>* =0, PetscInt my0=0, PetscInt nloc=0, PetscInt *d=0, PetscInt *o=0 );
#endif
protected:
  EBAMRPoissonOp*  m_op;
  bool m_isDefined;
};

#include "NamespaceFooter.H"

#ifdef CH_USE_PETSC
#ifndef CH_EXPLICIT_TEMPLATES
#include "PetscSolverI.H"
#endif // CH_EXPLICIT_TEMPLATES
#endif

#endif

#endif /*_PETSCSOLVER_H_*/
