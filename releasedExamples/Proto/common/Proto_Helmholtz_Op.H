
#include "AMRMultiGrid_Proto_Interface.H"

#ifndef __Proto_Helmholtz_Op_H
#define __Proto_Helmholtz_Op_H

/// The namespace everyone wants in on
namespace PrCh_AMR_Elliptic
{
  typedef Chombo::Box                       ch_box;
  typedef Chombo::ProblemDomain             ch_dom;
  typedef Proto::LevelBoxData<double, 1>    pr_lbd;
  typedef Chombo::MayDay                    ch_mayday;
  typedef Chombo::AMRLevelOpFactory<pr_lbd> ch_base_factory;
  
  /// 
  /**
     (2*DIM+1)-point Helmholtz operator with constant coefficients.
  **/
  class Proto_Helmholtz_Op: public   Base_Proto_AMRLevelOp<1>
  {
  private:

    /**
       Because there a factory class, we can dispense with weak construction
    **/
    Proto_Helmholtz_Op();

  public:
    
    virtual ~Proto_Helmholtz_Op()
    {
    }
    
    double    m_alpha;
    double    m_beta;
    double    m_dx;
    
    ///L(phi) = alpha*phi + beta*div(grad(phi));
    Proto_Helmholtz_Op(const double& a_alpha,
                       const double& a_beta,
                       const double& a_dx)
    {
      m_alpha = a_alpha;
      m_beta  = a_beta;
      m_dx    = a_dx;
    }
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(pr_lbd         & a_lph,
                  const pr_lbd   & a_phi)
    {
      ch_mayday::Error("Doh! not implemented");
    }
    
    ///
    virtual void
    relax(pr_lbd       & a_correction,
          const pr_lbd & a_residual, int a_iterations)
    {
      ch_mayday::Error("Doh! not implemented");
    }

    ///
    virtual void
    refluxFineCorrection(pr_lbd       &  a_lphCoar,
                         const pr_lbd &  a_phiFine, 
                         const pr_lbd &  a_phiCoar,
                         Base_Proto_AMRLevelOp<pr_lbd>*  a_finerOp)
    {
      ch_mayday::Error("Doh! not implemented");
    }

    ///
    /**
       correct the fine solution based on coarse correction
       a_phiThisLevel += I[2h->h](a_correctCoarse)
    */
    virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse)
    {
      ch_mayday::Error("Doh! not implemented");
    }
    
    ///a_resCoarse = average(a_resfine) (cannot call this restrict)
    virtual void restriction(T& a_resCoarse, T& a_resFine) = 0;
    {
      ch_mayday::Error("Doh! not implemented");
    }
    
    /// boundary conditions.   We shall restrict ourselves to the big two: Neumann and Dirichlet.
    virtual void
    fillGhostDataOutsideDomain(pr_lbd          & a_phi,
                               const ch_dom    & a_domain)
    {
      ch_mayday::Error("Doh! not implemented");
      return 4586;
    }

    ///
    /**
       return the refinement ratio to next coarser level.
       return 1 when there are no coarser Base_Proto_AMRLevelOp objects.
    */
    virtual int refToCoarser()
    {
      ch_mayday::Error("Doh! not implemented");
      return 4586;
    }
  };

  class Proto_Helmholtz_Op_Factory: public ch_base_factory
  {
    Proto_Helmholtz_Op_Factory()
    {
      ch_mayday::Error("not implemented");
    }
    virtual ~Proto_Helmholtz_Op_Factory()
    {
    }

    ///
    virtual MGLevelOp<T>* MGnewOp(const ch_dom         & a_FineindexSpace,
                                  int                    a_depth,
                                  bool                   a_homoOnly = true)
    {
      ch_mayday::Error("Doh! not implemented");
      return NULL;
    }

    ///
    virtual AMRLevelOp<T>* AMRnewOp(const ch_dom & a_indexSpace)
    {
      ch_mayday::Error("Doh! not implemented");
      return NULL;
    }

    ///
    virtual int refToFiner(const ch_dom & a_indexSpace) const
    {
      ch_mayday::Error("Doh! not implemented");
      return 4586;
    }
  };
}
#endif
