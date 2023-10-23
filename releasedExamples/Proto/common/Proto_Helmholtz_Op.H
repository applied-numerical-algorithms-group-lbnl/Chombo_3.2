
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

  class Proto_Helmholtz_Op_Factory: public ch_amrlevelop_factory
  {
    
  private:
    /// Strong construction, emacs and gdb are  a winning formula for life.
    Proto_Helmholtz_Op_Factory();
    
  public:

    /// no dynamic memory used here.   
    virtual ~Proto_Helmholtz_Op_Factory() { }

    ch_dom          m_coarsest_dom;
    vector<int>     m_ref_ratios;
    double          m_coarsest_dx;
    string          m_domBC;
    double          m_alpha;
    double          m_beta;
    double          m_max_coarse;

    ///
    /**
       I have found no good way to do this stuff.  
       This is the way amrpoissonopfactory does it.  
    **/
    vector< pair<ch_dom, double> >  m_domain_dx;
    
    ///the only constructor
    Proto_Helmholtz_Op_Factory(
      const ch_dom         &        a_coarsest_dom,
      const vector<ch_dbl> &        a_grids,
      const vector<int>    &        a_ref_ratios,
      const double         &        a_coarsest_dx,
      const string         &        a_domBC,
      const double         &        a_alpha,
      const double         &        a_beta)           
      : Chombo::AMRLevelOp<pr_lbd>()
    {
      m_max_coarse = 2;
      m_coarsest_dom  =  a_coarsest_dom;
      m_ref_ratios    =  a_ref_ratios  ;   
      m_coarsest_dx   =  a_coarsest_dx ;   
      m_domBC         =  a_domBC       ;         
      m_alpha         =  a_alpha       ;   
      m_beta          =  a_beta        ;
      
      //m_domains controls when mg v cycle ends 
      m_domain_dx.resize(a_grids.size());
      m_domain_dx[0] = a_coarsest_dom;
      double fine_dx = a_coarsest_dx;
      for(int ilev = 1; ilev < a_grids.size(); ilev++)
      {
        ch_dom fine_dom =   m_domain[ilev-1];
        fine_dom.refine(a_ref_ratios[ilev-1]);
        fine_dx /= a_ref_ratios[ilev-1];

        m_domain_dx[ilev].first  = fine_dom;
        m_domain_dx[ilev].second = fine_dx;
      }
    }

    ///
    inline int
    get_ref_index(ch_dom a_dom_fine) const
    {
      int iref  = 4586; bool found = false;
      for(int ilev = 0; ilev < m_domains.size(); ilev++)
      {
        if(m_domains[ilev] == a_dom_fine)
        {
          found = true;
          iref  = ilev;
        }
      }
      if(!found)
      {
        ch_mayday::Error("HelmholtzFactory::mgnewop: bogus domain");
      }
      return iref;
    }

    inline void
    get_local_ref_ratios(int    & a_ref_fine,
                         int    & a_ref_coar,
                         const int & a_iref) const
    {
      a_ref_fine = 4586;
      a_ref_coar = 4586;
      if(a_iref > 0)
      {
        a_ref_coar = m_ref_ratios[a_iref];
      }
      //somes folks are picky that the last one does not have to be there.
      if(a_iref < (m_ref_ratios.size()-1))
      {
        a_ref_fine = m_ref_ratios[a_iref+1];
      }
    }
      
    virtual inline ch_mglevelop*
    MGnewOp(const ch_dom         & a_dom_fine,
            int                    a_depth,
            bool                   a_homo_only = true)
    {
      int iref = get_ref_index(a_dom_fine);

      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      
      int factor = 1;
      for(int idep = 0; idep < a_depth; idep++)
      {
        factor*= 2;
      }
      ch_mglevelop* retval  = NULL;
      if((a_depth > 0) && (!a_dom_fine..coarsenable(2*factor)))
      {
        retval = NULL;
      }
      else
      {

        double loc_dx = m_dx[ref]*factor;
        Proto_Helmholtz_Op* derivedOp = 
          new Proto_Helmholtz_Op(m_alpha,
                                 m_beta,
                                 loc_dx,
                                 loc_ref_fine,
                                 loc_ref_coar,
                                 m_domBC);
        retval = static_cast< ch_mglevelop* >(derivedOp);
      }
      return retval;
    }

    ///
    virtual inline ch_amrlevelop* AMRnewOp(const ch_dom & a_indexSpace)
    {
      int iref = get_ref_index(a_indexSpace);
      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      double loc_dx = m_domains_dx[iref].second;

      Proto_Helmholtz_Op* derivedOp = 
        new Proto_Helmholtz_Op(m_alpha,
                               m_beta,
                               loc_dx,
                               loc_ref_fine,
                               loc_ref_coar,
                               m_domBC);
      return  static_cast< ch_amrlevelop* >(derivedOp);
    }

    ///
    virtual inline int refToFiner(const ch_dom & a_indexSpace) const
    {
      int iref = get_ref_index(a_indexSpace);
      int ref_rat = m_ref_ratios[iref];
      return iref;
    }
  };
}
#endif
