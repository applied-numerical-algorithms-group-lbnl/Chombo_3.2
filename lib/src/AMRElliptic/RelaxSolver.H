#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// BVS, Feb 7, 2005

#ifndef _RELAXSOLVER_H_
#define _RELAXSOLVER_H_

#include "LinearSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

///
/**
   Iterative solver which only uses the operator's preconditioner.
   Probably only useful as a bottom solver.
 */
template <class T>
class RelaxSolver : public LinearSolver<T>
{
public:
  RelaxSolver();
  virtual ~RelaxSolver();

  ///
  /**
     Set whether the solver uses only homogeneous boundary conditions
   */
  virtual void setHomogeneous(bool a_homogeneous)
  {
    m_homogeneous = a_homogeneous;
  }

  virtual void setUseFas(bool a_use_FAS)
  {
    m_use_FAS = a_use_FAS;
  }

  ///
  /**
     Define the solver.   a_op is the linear operator.
   */
  virtual void define(LinearOp<T>* a_op, bool a_homogeneous);

  ///
  /**
     Solve L(phi) = rho
   */
  virtual void solve(T& a_phi, const T& a_rhs);

  ///
  /**
     public member data: whether or not to use inhomogeneous boundary conditions.
   */
  bool m_homogeneous;

  bool m_use_FAS;

  ///
  /**
     public member data:  linear operator.
   */
  LinearOp<T>* m_op;

  ///
  /**
     public member data: maximum number of iterations
   */
  int m_imax;

  ///
  /**
     public member data: which type of norm to use (default is 2-norm)
  */
  int m_normType;

  ///
  /**
     public member data: how much screen output user wants
   */
  int m_verbosity;

  ///
  /**
     public member data:  solver tolerance
   */
  Real m_eps;

  ///
  /**
     public member data:  hang when min(previous residual norms) < m_hang * current residual norm  (default is zero -- don't check for hanging)
  */
  Real m_hang;
  private:
  Real m_minNorm;
  Real m_oldNorm;
};

// *******************************************************
// RelaxSolver Implementation
// *******************************************************

template <class T>
RelaxSolver<T>::RelaxSolver()
  :m_homogeneous(false),
   m_use_FAS(false),
   m_op(NULL),
   m_imax(40),
   m_normType(2),
   m_verbosity(0),
   m_eps(1.0E-6),
   m_hang(0.0)
{
}

template <class T>
RelaxSolver<T>::~RelaxSolver()
{
}

template <class T>
void RelaxSolver<T>::define(LinearOp<T>* a_operator, bool a_homogeneous)
{
  m_homogeneous = a_homogeneous;
  m_op = a_operator;
}

template <class T>
void RelaxSolver<T>::solve(T& a_phi, const T& a_rhs)
{
  CH_assert(m_op != NULL);
  CH_TIME("RelaxSolver::solve");

  T r, e;

  m_op->create(r, a_rhs);
  m_op->create(e, a_phi);

  if (!m_use_FAS) {
      m_op->setToZero(r); // added by petermc, 27 Nov 2013, to zero out ghosts
  }
  m_op->residual(r, a_phi, a_rhs, m_homogeneous);

  Real norm = m_op->norm(r, m_normType);
  m_minNorm = norm;
  m_oldNorm = norm;

  if (m_verbosity >= 3)
    {
      pout() << "      RelaxSolver:: initial Residual norm = "
             << norm << "\n";
      pout() << "      RelaxSolver:: m_eps = " << m_eps << " m_hang = " << m_hang << endl;
    }

  if (norm > 0.) // Do not iterate if norm == 0.
    {
      Real initialNorm = norm;
      int iter =0;
      while (iter <  m_imax)
        {
          if (!m_use_FAS) {
              m_op->setToZero(e);
              m_op->preCond(e, r, r);
              m_op->incr(a_phi, e, 1.0);
          } else {
              m_op->preCond(a_phi, r, a_rhs);
          }
          m_op->residual(r, a_phi, a_rhs, m_homogeneous);
          norm = m_op->norm(r, m_normType);
          iter++;

          if (m_verbosity > 3)
            {
              pout() << "      RelaxSolver:: iteration = " << iter
                     << " residual norm = " << norm
                     << ", rate = " << m_oldNorm/norm << endl;
            }

          if (norm <  m_eps*initialNorm) break;

          if (norm < m_minNorm)
            {
              m_minNorm = norm;
            }

          if ((1- m_hang)*norm >= m_oldNorm)
            {
              if (m_verbosity > 3)
                {
                  pout() << "      RelaxSolver:: exit due to solver hang " << endl;
                }
              break;
            }

            m_oldNorm = norm;
        }
    }

  m_op->clear(r);
  m_op->clear(e);
}

#include "NamespaceFooter.H"
#endif /*_RELAXSOLVER_H_*/
