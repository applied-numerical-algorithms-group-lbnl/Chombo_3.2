
#include "AMRMultiGrid_Proto_Interface.H"
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Vector.H"

#ifndef __Proto_Helmholtz_Op_H
#define __Proto_Helmholtz_Op_H
/****************/
PROTO_KERNEL_START 

unsigned int axbyHelmF(int              a_pt[DIM],
                            Proto::Var<Real, 1>     a_lph,
                            Proto::Var<Real, 1>     a_x,
                            Proto::Var<Real, 1>     a_y,
                            Real             a_a,
                            Real             a_b,
                            int              a_iprint)
{
  double xval = a_x(0);
  double yval = a_y(0);
  double lph  = a_a*xval + a_b*yval;
  a_lph(0) = lph;
  return 0;

}
PROTO_KERNEL_END(axbyHelmF, axbyHelm) 

/****************/
PROTO_KERNEL_START 
unsigned int  gsrbResidF(int              a_pt[DIM],
                         Proto::Var<Real, 1>     a_phi,
                         Proto::Var<Real, 1>     a_res,
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
PROTO_KERNEL_END(gsrbResidF, gsrbResid) 
/// The new cool namespace that everyone wants in on.   
namespace PrCh_AMR_Elliptic
{
  typedef Proto::LevelBoxData<double, 1>    pr_lbd;
  typedef Proto::Stencil<double>            pr_sten;
  typedef Proto::InterpStencil<double>      pr_interp_sten;
  
  typedef Chombo::Box                       ch_box;
  typedef  Proto::Box                       pr_box;
  typedef Chombo::ProblemDomain             ch_dom;
  typedef Chombo::DisjointBoxLayout         ch_dbl;
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
    typedef Proto::Stencil<double> pr_sten;
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
    pr_sten   m_sten;
    
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
      
      
      // (1/h^2) factor missing  in proto stencil for lapl.
      double hsqinv = 1./(m_dx*m_dx);
      pr_sten  sten_lapl = pr_sten::Laplacian();
      // scaling by (1/h^2) factor    
      sten_lapl *= hsqinv;
      
      pr_sten  sten_iden = (1.0)*Proto::Shift(pr_pt::Zeros());
      m_sten = m_alpha*sten_iden + m_beta*sten_lapl;
      static bool printedOnce = false;
      if(!printedOnce)
      {
        Chombo::pout() << "operator stencil:" << std::endl;
        m_sten.print();
        printedOnce = true;
      }
    }
    
    /// output = aco*xval + bco*yval.   components defeat just doing this in the base class
    /**
       Standard function that I use for all arithmetic.  Because this is not fortran,
       it does not matter if any of the arguments point to overlapping things.
     **/
    virtual void 
    axby(pr_lbd         & a_lphi,
         const pr_lbd   & a_xval,
         const pr_lbd   & a_yval,
         double    a_aco,
         double    a_bco) 
    {
      auto grids =  a_lphi.layout();
      auto dit   =  a_lphi.begin();
      int iprint = 0; //debugging hook
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto       & lphfab =  a_lphi[dit[ibox]];
        const auto &   xfab =  a_xval[dit[ibox]];
        const auto &   yfab =  a_yval[dit[ibox]];
        auto valid          =   grids[dit[ibox]];
        Proto::forallInPlace_i(axbyHelm, valid, lphfab, xfab, yfab,
                               a_aco, a_bco, iprint);
      }
    }
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(pr_lbd         & a_lph,
                  const pr_lbd   & a_phi)
    {
      CH_TIME("Proto_Helmholtz_Op::applyOperator");
      auto grids = a_lph.layout();
      auto dit   = a_lph.begin();
      for(int ibox = 0; ibox <  dit.localSize(); ibox++)
      {
        auto        & dst = a_lph[dit[ibox]];
        const auto  & src = a_phi[dit[ibox]];
        auto valid        = grids[dit[ibox]];
        bool initToZero = true; double scale = 1;
        ///Note: arguments' order is the opposite of usual (we usually do dst, src).
        m_sten.apply(src, dst, valid, initToZero, scale); 
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
      
      int idebug = 0; // gdb debugging hook
      for(int iiter = 0; iiter < a_iterations;  iiter++)
      {
        for(int iredblack = 0; iredblack < 2; iredblack++)
        {
          residual(resid, a_phi, a_rhs);
          for(int ibox = 0; ibox < dit.localSize(); ibox++)
          {
            const auto & resfab =  resid[dit[ibox]];
            const auto & phifab =  a_phi[dit[ibox]];
            auto valid          =  grids[dit[ibox]];
            int iprint = 0; // I left in a debugging hook.
            Proto::forallInPlace_i(gsrbResid, valid, phifab, resfab,
                                   m_alpha, m_beta, m_dx, iredblack, iprint);
          }
          idebug = 0; //stop here for end of dataiter  iteration
        }
        idebug = 0;   //stop here for end of red-black iteration
      }
    }

    
    /// domain boundary conditions.
    /** 
        Anything can be hidden in a string, so the interface is general enough.
        For now, we shall restrict ourselves to the big two: Neumann and Dirichlet.
    **/
    virtual void
    fillGhostDataOutsideDomain(pr_lbd          & a_phi,
                               const ch_dom    & a_domain)
    {
      typedef Proto_Transition::DataUtilities<1> pr_trans_data;
      pr_dom domain = ProtoCh::getProtoDomain(a_domain);
      pr_box dombox = domain.box();
      for(int idir = 0; idir < DIM; idir++)
      {
        for(int hilo = 0; hilo < 2; hilo++)
        {
          if(m_domBC == string("Dirichlet"))
          {
            pr_trans_data::domainGhostReflectOdd( a_phi, dombox, idir, hilo);
          }
          else if(m_domBC == string("Neumann"))
          {
            pr_trans_data::domainGhostReflectEven(a_phi, dombox, idir, hilo);
          }
          else
          {
            ch_mayday::Error("Proto_Helmholtz_Op::fillGhostDataOutsideDomain: unrecognized bc string");
          }
        }
      }
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

    Chombo::Vector<int>                    m_ref_ratios;
    Chombo::Vector<ch_dbl>                 m_grids;
    Chombo::Vector<pair<ch_dom, double> >  m_domain_dx;
    double                                 m_coarsest_dx;
    string                                 m_domBC;
    double                                 m_alpha;
    double                                 m_beta;
    int                                    m_max_coarse;

    
    ///the only constructor
    Proto_Helmholtz_Op_Factory(
      const ch_dom                 &        a_coarsest_dom,
      const Chombo::Vector<ch_dbl> &        a_grids, //just used for coarsenability
      const Chombo::Vector<int>    &        a_ref_ratios,
      const double                 &        a_coarsest_dx,
      const string                 &        a_domBC,
      const double                 &        a_alpha,
      const double                 &        a_beta)           
      : ch_amrlevelop_factory()
    {
      m_max_coarse = 2;
      m_ref_ratios    =  a_ref_ratios  ;
      m_grids         =  a_grids;
      m_coarsest_dx   =  a_coarsest_dx ;   
      m_domBC         =  a_domBC       ;         
      m_alpha         =  a_alpha       ;   
      m_beta          =  a_beta        ;

      m_domain_dx.resize(a_grids.size());
      m_domain_dx[0].first  = a_coarsest_dom;
      m_domain_dx[0].second = a_coarsest_dx;
      double fine_dx = a_coarsest_dx;
      for(int ilev = 1; ilev < a_ref_ratios.size(); ilev++)
      {
        ch_dom fine_dom =   m_domain_dx[ilev-1].first;
        fine_dom.refine(   a_ref_ratios[ilev-1]);
        fine_dx      /=    a_ref_ratios[ilev-1];

        m_domain_dx[ilev].first  = fine_dom;
        m_domain_dx[ilev].second = fine_dx;
      }
    }

    ///
    inline int
    get_ref_index(ch_dom a_dom_fine) const
    {
      int iref  = 4586; bool found = false;
      for(int ilev = 0; ilev < m_domain_dx.size(); ilev++)
      {
        if(m_domain_dx[ilev].first == a_dom_fine)
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
      bool grids_coarsenable = m_grids[iref].coarsenable(2*factor);
      if((a_depth > 0) && (!grids_coarsenable))
      {
        retval = NULL;
      }
      else
      {

        double loc_dx = m_domain_dx[iref].second; 
        loc_dx *= factor;
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
      double loc_dx = m_domain_dx[iref].second;

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
      return ref_rat;
    }
  };
}
#endif
