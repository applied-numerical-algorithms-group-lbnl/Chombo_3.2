
#include "AMRMultiGrid_Proto_Interface.H"

#ifndef __Proto_Helmholtz_Op_H
#define __Proto_Helmholtz_Op_H

/****************/
PROTO_KERNEL_START 
unsigned int  gsrbResidF(int              a_pt[DIM],
                         Var<Real, 1>     a_phi,
                         Var<Real, 1>     a_res,
                         Real             a_alpha,
                         Real             a_beta,
                         Real             a_dx,
                         int              a_iredBlack,
                         int              a_iprint)
{
  int sumpt = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    sumpt += a_pt[idir];
  }
  if(sumpt%2 == a_iredBlack)
  {
    static const Real safety = 1.0;
    Real reg_diag_val = -2*DIM/a_dx/a_dx;
    Real reg_realdiag = a_alpha  + a_beta*reg_diag_val;
    Real reg_lambda   = safety/reg_realdiag;
    
    Real phival = a_phi(0);
    Real resval = a_res(0);

    a_phi(0) = phival + reg_lambda*resval;
  }
  return 0;

}
PROTO_KERNEL_END(setRHSF, setRHS) 
/// The namespace everyone wants in on
namespace PrCh_AMR_Elliptic
{
  typedef Proto::LevelBoxData<double, 1>    pr_lbd;
  typedef Proto::Stencil<double>            pr_sten;
  typedef Chombo::Box                       ch_box;
  typedef Chombo::ProblemDomain             ch_dom;
  typedef Chombo::MayDay                    ch_mayday;
  typedef Chombo::AMRLevelOpFactory<pr_lbd> ch_base_factory;
  /// 
  /**
     (2*DIM+1)-point Helmholtz operator with constant coefficients.
  **/
  class Proto_Helmholtz_Op: public   Base_Proto_AMRLevelOp<1>
  {
  private:

    ///Because there a factory class, we can dispense with weak construction
    Proto_Helmholtz_Op();

  public:
    
    virtual ~Proto_Helmholtz_Op()
    {
    }
    
    // identity weight
    double    m_alpha;
    // laplacian weight
    double    m_beta;
    // grid spacing
    double    m_dx;

    ///amr thing
    int m_refToFine;
    int m_refToCoar;
    //we can hide whatever complicated bc we have behind a unique string.    
    string    m_domBC;  
    Stencil<double>  m_sten;
    
    ///L(phi) = alpha*phi + beta*div(grad(phi));
    Proto_Helmholtz_Op(const double& a_alpha,
                       const double& a_beta,
                       const double& a_dx,
                       const int   & a_refToFine,
                       const int   & a_refToCoar,
                       const string& a_domBC)
    {
      CH_TIME("Proto_Helmholtz_Op::Proto_Helmholtz_Op");
      m_alpha     = a_alpha;
      m_beta      = a_beta;
      m_dx        = a_dx;
      m_refToFine = a_refToFine;
      m_refToCoar = a_refToFine;
      m_domBC     = a_domBC;
      
      
      ///in EB, this bit  would be called defineStencil or some such.
      double hsqinv = 1./(m_dx*m_dx);
      // (1/h^2) factor missing 
      Stencil<double>  sten_lapl = Stencil<double>::Laplacian;
      //scaling by (1/h^2) factor    
      sten_lapl *= hsqinv; 
      Stencil<double>  sten_iden = (1.0)*Shift(Point::Zeros());
      m_sten *= m_alpha*sten_iden + m_beta*sten_lapl;
    }
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(pr_lbd         & a_lph,
                  const pr_lbd   & a_phi)
    {
      CH_TIME("Proto_Helmholtz_Op::applyOperator");
      auto grids = a_lph.layout();
      auto dit   = a_lph.begin();
      for(int ibox = 0; ibox <  dit.size(); ibox++)
      {
        auto        & dst = a_lph[dit[ibox]];
        const auto  & src = a_phi[dit[ibox]];

        ///Note: arguments' order is the opposite of usual (we usually do dst, src).
        m_sten.apply(src, dst, bx, initToZero, scale); 
      }
    }
    
    ///phi+= lambda*resid on red black
    virtual void
    relax(pr_lbd       & a_phi,
          const pr_lbd & a_rhs,
          int a_iterations)
    {
      
      pr_lbd resid;
      create(resid, a_rhs);
      auto grids =  a_rhs.layout();
      auto dit   =  a_rhs.begin();
      
      for(int iiter = 0; iiter < a_iterations;  iiter++)
      {
        for(int iredblack = 0; iredblack < 2; iredblack++)
        {
          residual(resid, a_phi, a_rhs);
          int print = 0; //debugging hook
          for(int ibox = 0; ibox < dit.size(); ibox++)
          {
            const auto & resfab =  resid[dit[ibox]];
            const auto & phifab =  a_phi[dit[ibox]];
            forAllInPlace_i(gsrbResid, patch, phifab, resfab, m_alpha, m_beta, m_dx, iredBlack, iprint);
          }
        }
      }
    }

    ///a_phiThisLevel += I[2h->h](a_correctCoarse)
    virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse)
    {
      ch_mayday::Error("Doh! not implemented");
    }
    /**
       correct the fine solution based on coarse correction
       
    */
    
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

    /// amr thing 
    virtual void
    refluxFineCorrection(pr_lbd       &  a_lphCoar,
                         const pr_lbd &  a_phiFine, 
                         const pr_lbd &  a_phiCoar,
                         Base_Proto_AMRLevelOp<pr_lbd>*  a_finerOp)
    {
      ch_mayday::Error("Doh! not implemented");
    }

    ///amr thing
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
