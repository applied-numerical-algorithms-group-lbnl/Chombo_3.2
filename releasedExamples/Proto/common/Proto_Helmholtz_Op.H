
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
/// The new cool namespace that everyone wants in on.   
namespace PrCh_AMR_Elliptic
{
  typedef Proto::LevelBoxData<double, 1>    pr_lbd;
  typedef Proto::Stencil<double>            pr_sten;
  typedef Proto::InterpStencil<double>      pr_interp_sten;
  
  typedef Chombo::Box                       ch_box;
  typedef Chombo::ProblemDomain             ch_dom;
  typedef Chombo::MayDay                    ch_mayday;
  typedef Chombo::AMRLevelOpFactory<pr_lbd> ch_amrlevelop_factory;
  typedef Chombo::AMRLevelOp<pr_lbd> ch_amrlevelop;
  typedef Chombo::MGLevelOp<pr_lbd>  ch_mglevelop;
  
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
    //we can hide whatever complicated bc we have behind a unique string.    
    string    m_domBC;
    ///amr hooks that do not need valid values when they do not apply
    int       m_refToFine;   
    int       m_refToCoar;

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
      m_refToCoar = a_refToCoar;
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

    
    /// domain boundary conditions.
    /** 
        For now we shall restrict ourselves to the big two: Neumann and Dirichlet.
        Of course, anything can be hidden in a string so the interface is general enough.
     **/
    virtual void
    fillGhostDataOutsideDomain(pr_lbd          & a_phi,
                               const ch_dom    & a_domain)
    {
      pr_dom domain = PrCh::getProtoDomain(a_domain);
      for(int idir = 0; idir < DIM; idir++)
      {
        for(int hilo = 0; hilo < 2; hilo++)
        {
          if(m_domBC == string("Dirichlet"))
          {
            Proto_Transition::domainGhostReflectOdd( a_phi, a_domain, idir, hilo);
          }
          else if(m_domBC == string("Neumann"))
          {
            Proto_Transition::domainGhostReflectEven(a_phi, a_domain, idir, hilo);
          }
          else
          {
            ch_mayday::Error("Proto_Helmholtz_Op::fillGhostDataOutsideDomain: unrecognized bc string");
          }
        }
      }
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
      return m_refToCoar;
    }
  };

  class Proto_Helmholtz_Op_Factory: public ch_base_factory
  {
    
  private:
    /// Strong construction, emacs and gdb are  a winning formula for life.
    Proto_Helmholtz_Op_Factory();
    
  public:

    /// no dynamic memory used here
    virtual ~Proto_Helmholtz_Op_Factory() { }

    pr_dom          m_coarsest_dom;
    vector<ch_dbl>  m_grids;
    vector<int>     m_refRatios;
    double          m_coarsest_dx;
    string          m_domBC;
    double          m_alpha;
    double          m_beta;
    
    ///the only constructor
    Proto_Helmholtz_Op_Factory(
      const pr_dom        &        a_coarsest_dom,
      const vector<ch_dbl>&        a_grids,
      const vector<int>   &        a_refRatios,
      const double        &        a_coarsest_dx,
      const string        &        a_domBC,
      const double        &        a_alpha,
      const double        &        a_beta)           
      : Chombo::AMRLevelOp<pr_lbd>()
    {
      m_coarsest_dom  =  a_coarsest_dom;
      m_grids         =  a_grids;
      m_refRatios     =  a_refRatios   ;   
      m_coarsest_dx   =  a_coarsest_dx ;   
      m_domBC         =  a_domBC       ;         
      m_alpha         =  a_alpha       ;   
      m_beta          =  a_beta        ;      
    }

      HERE---remember to read the comment below
    /// remember that this has the interface which returns NULL if not coarsenable
    virtual ch_mglevelop*
    MGnewOp(const ch_dom         & a_FineindexSpace,
            int                    a_depth,
            bool                   a_homoOnly = true)
    {

      Proto_Helmholtz_Op* derivedOp = 
        new Proto_Helmholtz_Op(m_alpha,
                               m_beta,
                               loc_dx,
                               loc_refToFine,
                               loc_refToCoar,
                               m_domBC);

      return static_cast< ch_mglevelop* >(derivedOp);
    }

    ///
    virtual ch_amrlevelop* AMRnewOp(const ch_dom & a_indexSpace)
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
