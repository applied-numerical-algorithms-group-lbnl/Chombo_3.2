#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// SLC, Oct 12, 2011

#ifndef _CGSOLVER_H_
#define _CGSOLVER_H_

#include "LinearSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

///
/**
   Elliptic solver using the preconditioned conjugate gradient algorithm.
 */
template <class T>
class CGSolver : public LinearSolver<T>
{
public:

  CGSolver();

  virtual ~CGSolver();

  virtual void setHomogeneous(bool a_homogeneous)
     {m_homogeneous = a_homogeneous;}

  ///
  /**
     define the solver.   a_op is the linear operator.
     a_homogeneous is whether the solver uses homogeneous boundary
     conditions.
   */
  virtual void define(LinearOp<T>* a_op, bool a_homogeneous);

  ///solve the equation.
  virtual void solve(T& a_phi, const T& a_rhs);

  ///
  virtual void setConvergenceMetrics(Real a_metric,
                                     Real a_tolerance);

  ///
  /**
     public member data: whether the solver is restricted to
     homogeneous boundary conditions
   */
  bool m_homogeneous;

  ///
  /**
     public member data: operator to solve.
   */
  LinearOp<T>* m_op;

  ///
  /**
     public member data:  maximum number of iterations
   */
  int m_imax;

  ///
  /**
     public member data:  how much screen out put the user wants.
     set = 0 for no output.
   */
  int m_verbosity;

  ///
  /**
     public member data:  solver tolerance
   */
  Real m_eps;

  ///
  /**
     public member data:  relative solver tolerance
   */
  Real m_reps;

  ///
  /**
     public member data: solver convergence metric -- if negative, use
     initial residual; if positive, then use m_convergenceMetric
  */
  Real m_convergenceMetric;

  ///
  /**
     public member data:  minium norm of solution should change per iterations
   */
  Real m_hang;

  ///
  /**
     public member data:
     set = -1 if solver exited for an unknown reason
     set =  1 if solver converged to tolerance
     set =  2 if rho = 0
     set =  3 if max number of restarts was reached
   */
  int m_exitStatus;

  ///
  /**
     public member data:  what the algorithm should consider "close to zero"
   */
  Real m_small;

  ///
  /**
     public member data:  number of times the algorithm can restart
   */
  int m_numRestarts;

  ///
  /**
     public member data:  norm to be used when evaluating convergence.
     0 is max norm, 1 is L(1), 2 is L(2) and so on.
   */
  int m_normType;

};

// *******************************************************
// CGSolver Implementation
// *******************************************************

template <class T>
CGSolver<T>::CGSolver()
  :m_homogeneous(false), m_op(NULL), m_imax(80), m_verbosity(3), m_eps(1.0E-6),
   m_reps(1.0E-12), m_convergenceMetric(-1.0), m_hang(1E-8), m_exitStatus(-1),
   m_small(1.0E-30), m_numRestarts(5), m_normType(2){;}

template <class T>
CGSolver<T>::~CGSolver()
{
  m_op = NULL;
}

template <class T>
void CGSolver<T>::define(LinearOp<T>* a_operator, bool a_homogeneous)
{
  m_homogeneous = a_homogeneous;
  m_op = a_operator;
}

template <class T>
void CGSolver<T>::solve(T& a_phi, const T& a_rhs)
{
  CH_TIMERS("CGSolver::solve");

  CH_TIMER("CGSolver::solve::Initialize",timeInitialize);
  CH_TIMER("CGSolver::solve::MainLoop",timeMainLoop);
  CH_TIMER("CGSolver::solve::Cleanup",timeCleanup);

  CH_START(timeInitialize);

  T r,z,p,Lp;

  m_op->create(r,a_rhs);
  m_op->create(z,a_phi);
  m_op->create(p,a_phi);
  m_op->create(Lp,a_rhs);

  CH_assert(m_op != NULL);

  m_op->residual(r, a_phi, a_rhs, m_homogeneous);
 

  Real initial_norm = m_op->norm(r, m_normType);
  Real norm = initial_norm;
  if (m_verbosity >= 3)
    {
      pout() << "      CG:: initial Residual norm = "
             << initial_norm << "\n";
    }
  m_op->setToZero(z);
  m_op->preCond(z, r);
  m_op->assign(p,z);

  // if a convergence metric has been supplied, replace initial residual
  // with the supplied convergence metric...
  if (m_convergenceMetric > 0)
    {
      initial_norm = m_convergenceMetric;
    }

  CH_STOP(timeInitialize);

  CH_START(timeMainLoop);
  int iter = 0;
  while (iter < m_imax )
    {

      m_op->applyOp(Lp, p, true);
      Real rdotz = m_op->dotProduct(r, z);
      Real alpha =   rdotz /  m_op->dotProduct(p , Lp);
      m_op->incr(a_phi, p , alpha);
      m_op->incr(r, Lp , -alpha);
      Real norm_old = norm;
      norm = m_op->norm(r, m_normType);
      if (m_verbosity >= 3)
        {
          pout() << "      CG::     iteration = "  << iter        
		 << ", error norm = "  << norm 
		 << ", rate = " << norm_old/norm << "\n";
        }
      
     
      if(norm <  m_eps*initial_norm)
	{
	  break;
	}
      if (norm_old < m_hang * norm)
	{ 
	 break;
	}
     
      m_op->preCond(z, r);
      Real beta =  m_op->dotProduct(r, z) / rdotz;
      m_op->scale(p,beta);
      m_op->incr(p, z , 1.0);
      iter++;
    }
  CH_STOP(timeMainLoop);

  CH_START(timeCleanup);

  if (m_verbosity >= 3)
    {
      pout() << "      CG:: " << iter 
	     << " iterations, final Residual norm = "
             << norm << "\n";
    }



  m_op->clear(r);
  m_op->clear(z);
  m_op->clear(p);
  m_op->clear(Lp);
  CH_STOP(timeCleanup);
}

template <class T>
void CGSolver<T>::setConvergenceMetrics(Real a_metric,
                                              Real a_tolerance)
{
  m_convergenceMetric = a_metric;
  m_eps = a_tolerance;
}

#include "NamespaceFooter.H"
#endif /*_CGSOLVER_H_*/
