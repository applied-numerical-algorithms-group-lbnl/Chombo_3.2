
#include "AMRMultiGrid_Proto_Interface.H"
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Vector.H"

#ifndef __Proto_Conductivity_Op_H
#define __Proto_Conductivity_Op_H

///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
PROTO_KERNEL_START 
unsigned int  operaConduF(int              a_pt[DIM],
                          Proto::Var<Real, 1>     a_lph,
                          Proto::Var<Real, 1>     a_phi,
                          Proto::Var<Real, 1>     a_acoef,
                          Proto::Var<Real, 1>     a_xincr,
                          Proto::Var<Real, 1>     a_yincr,
                          Proto::Var<Real, 1>     a_zincr,
                          Real                    a_alpha,
                          Real                    a_beta)
{
  //z is set to zero in 2d
  Real divergence = a_xincr(0) + a_yincr(0) + a_zincr(0);
  Real lphval = a_alpha*a_acoef(0)*a_phi(0) + a_beta*divergence;
  a_lph(0) = lphval;
  return 0;
}
PROTO_KERNEL_END(operaConduF, operaCondu)

///lamdbda = safety/approximate_diagonal
/**
   
   The true diagonal here is 
   diag_true = alpha*acoef 
   - (beta/h^2)*(bcoef_x i+1/2 + bcoef_x(i-1/2 )
   - (beta/h^2)*(bcoef_y j+1/2 + bcoef_y(j-1/2 )
   - (beta/h^2)*(bcoef_z k+1/2 + bcoef_z(k-1/2 )
 
   That is a lot of shifting just for a relaxation coefficient.
   For now, we will use 
   diag_approx = alpha*acoef 
   - (2*beta/h^2)*(bcoef_x(i-1/2 )
   - (2*beta/h^2)*(bcoef_y(j-1/2 )
   - (2*beta/h^2)*(bcoef_z(k-1/2 )
**/
PROTO_KERNEL_START 
unsigned int  lambdaConduF(int              a_pt[DIM],
                           Proto::Var<Real, 1>     a_lambda,
                           Proto::Var<Real, 1>     a_acoef,
                           Proto::Var<Real, 1>     a_bco_x,
                           Proto::Var<Real, 1>     a_bco_y,
                           Proto::Var<Real, 1>     a_bco_z,
                           Real                    a_alpha,
                           Real                    a_beta,
                           Real                    a_dx)
{
  Real safety = 0.9;
  Real factor = -2.*a_beta/(a_dx*a_dx);
  Real diagon = a_alpha*a_acoef(0);
  diagon += factor*a_bco_x(0);
  diagon += factor*a_bco_y(0);
#if DIM==3
  diagon += factor*a_bco_z(0);
#endif
  a_lambda(0) = safety/diag;
  return 0;
}
PROTO_KERNEL_END(lambdaConduF, lambdaCondu)


///
/**
  phi = phi + lambda*(lphi - rho) 
  but only where 
  sum(i+j+k)%2== iredblack
**/
PROTO_KERNEL_START 
unsigned int  gsrbConduF(int              a_pt[DIM],
                         Proto::Var<Real, 1>     a_phi,
                         Proto::Var<Real, 1>     a_res,
                         Proto::Var<Real, 1>     a_lam,
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
    Real lambda = a_lam(0);
    Real phival = a_phi(0);
    Real resval = a_res(0);

    a_phi(0) = phival + lambda*resval;
  }
  return 0;

}
PROTO_KERNEL_END(gsrbConduF, gsrbCondu) 
/// The new cool namespace that everyone wants in on.   
namespace PrCh_AMR_Elliptic
{
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
  typedef Proto::Stencil<double>           pr_sten;
  
  typedef Proto::LevelBoxData<double, 1, Proto::CELL   >   pr_cell_data;
  typedef Proto::LevelBoxData<double, 1, Proto::FACE_0 >   pr_xfac_data;
  typedef Proto::LevelBoxData<double, 1, Proto::FACE_1 >   pr_yfac_data;
  typedef Proto::LevelBoxData<double, 1, Proto::FACE_2 >   pr_zfac_data;
  
  /// Helmholtz operator with spatially varying coefficients
  /**
     L(phi) = alpha*(a(x) phi) + beta*(Div b(x) Grad(phi))A 
     alpha, beta are constants.
     a(x) and b(x) are functions of space only.
  **/
  class Proto_Conductivity_Op: public   Base_Proto_AMRLevelOp<1>
  {
  private:

    ///Because there a factory class, we can dispense with weak construction
    Proto_Conductivity_Op();

  public:
    virtual ~Proto_Conductivity_Op()
    {
    }

    shared_ptr<pr_cell_data>  m_acoef;
    shared_ptr<pr_xfac_data>  m_bco_x;
    shared_ptr<pr_yfac_data>  m_bco_y;
    shared_ptr<pr_zfac_data>  m_bco_z;
    shared_ptr<pr_cell_data>  m_lambda;
    double    m_alpha;
    double    m_beta;
    double    m_dx;
    string    m_domBC;
    ///amr hooks that do not need valid values when they do not apply
    int       m_refToFine;   
    int       m_refToCoar;
    pr_dbl    m_grids;
    
    /// L(phi) = alpha*(a(x) phi) + beta*(Div b(x) Grad(phi))A 
    Proto_Conductivity_Op(
      shared_ptr<pr_cell_data>  a_acoef,
      shared_ptr<pr_xfac_data>  a_bco_x,
      shared_ptr<pr_yfac_data>  a_bco_y,
      shared_ptr<pr_zfac_data>  a_bco_z,
      const pr_dbl&             a_grids,
      const double&             a_alpha,
      const double&             a_beta,
      const double&             a_dx,
      const int   &             a_refToFine,
      const int   &             a_refToCoar,
      const string&             a_domBC)
    {
      CH_TIME("Proto_Conductivity_Op::Proto_Conductivity_Op");
      m_acoef     = a_acoef;
      m_bco_x     = a_bco_x;
      m_bco_y     = a_bco_y;
      m_bco_z     = a_bco_z;
      m_grids     = a_grids;
      m_alpha     = a_alpha;
      m_beta      = a_beta;
      m_dx        = a_dx;
      m_refToFine = a_refToFine;
      m_refToCoar = a_refToCoar;
      m_domBC     = a_domBC;
      defineAndSetRelaxationCoefficient();
    }

    ///
    void defineAndSetRelaxationCoefficient()
    {
      m_lambda = shared_ptr<pr_cell_data>(new pr_cell_data(m_grids, pr_pt::Zero));
      auto dit   = m_lambda.begin();
      ///fill a_lph with div( b(x) grad (phi))
      for(int ibox = 0; ibox <  dit.size(); ibox++)
      {
        const auto  & lamfab = (*m_lambda)[dit[ibox]];
        const auto  & acofab =  (*m_acoef)[dit[ibox]];
        auto valid           =     m_grids[dit[ibox]];
        const auto& bco_x    =  (*m_bco_x)[dit[ibox]];
        const auto& bco_y    =  (*m_bco_x)[dit[ibox]];
#if DIM==3                      
        const auto& bco_z    =  (*m_bco_x)[dit[ibox]];
#else
        //something fake to send
        const auto& bco_z    = (*m_bco_y)[dit[ibox]];
#endif
        ///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
        Proto::forallInPlace_i(lambdaCondu, valid, lamfab,
                               acofab, bco_x, bco_y, bco_z,
                               m_alpha, m_beta, m_dx);
      }

    }

    ///gets (b_i+1/2*grad_phi_i+1/2 - b_i+1/2*grad_phi_i+1/2)/m_dx
    template<Centering cent>
    void
    getDivergenceIncrement()
    {
    }
    
    ///this makes a_lph = div(b(grad phi))
    virtual void
    computeDivBGradPhi(pr_lbd       & a_lph,
                       const pr_lbd & a_phi)
    {
      auto dit   = a_lph.begin();

      ///fill a_lph with div( b(x) grad (phi))
      for(int ibox = 0; ibox <  dit.size(); ibox++)
      {
        auto        & lphfab =      a_lph[dit[ibox]];
        const auto  & phifab =      a_phi[dit[ibox]];
        const auto  & acofab = (*m_acoef)[dit[ibox]];
        auto valid           =    m_grids[dit[ibox]];
        const auto& bco_x    = (*m_bco_x)[dit[ibox]];
        const auto& bco_y    = (*m_bco_x)[dit[ibox]];
#if DIM==3             
        const auto& bco_z    = (*m_bco_x)[dit[ibox]];
#endif        
        BoxData<double, 1> xincr(valid); xincr.setVal(0.);
        BoxData<double, 1> yincr(valid); yincr.setVal(0.);
        BoxData<double, 1> zincr(valid); zincr.setVal(0.);
        
        getDivergenceIncrement<PR_FACE_0>(xincr, phifab, bco_x, valid);
        getDivergenceIncrement<PR_FACE_1>(yincr, phifab, bco_y, valid);
#if DIM==3                                                  
        getDivergenceIncrement<PR_FACE_2>(zincr, phifab, bco_z, valid);
#endif
        ///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
        Proto::forallInPlace_i(operaCondu, valid,
                               lphfab, phifab, acofab,
                               xincr, yincr, zincr,
                               m_alpha, m_beta);
      }
    }
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(pr_lbd         & a_lph,
                  const pr_lbd   & a_phi)
    {
      CH_TIME("Proto_Conductivity_Op::applyOperator");
      auto dit   = a_lph.begin();
      for(int ibox = 0; ibox <  dit.size(); ibox++)
      {
        auto        & dst =   a_lph[dit[ibox]];
        const auto  & src =   a_phi[dit[ibox]];
        auto valid        = m_grids[dit[ibox]];
        //this makes a_lph = div(b(grad phi))
        computeDivBGradPhi(a_lph, a_phi);
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
      auto dit   =  a_rhs.begin();
      
      int idebug = 0; // gdb debugging hook
      for(int iiter = 0; iiter < a_iterations;  iiter++)
      {
        for(int iredblack = 0; iredblack < 2; iredblack++)
        {
          residual(resid, a_phi, a_rhs);
          for(int ibox = 0; ibox < dit.size(); ibox++)
          {
            const auto & resfab =    resid[dit[ibox]];
            const auto & phifab =    a_phi[dit[ibox]];
            auto valid          =  m_grids[dit[ibox]];
            int iprint = 0; // I left in a debugging hook.
            Proto::forallInPlace_i(gsrbCondu, valid, phifab, resfab,
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
            ch_mayday::Error("Proto_Conductivity_Op::fillGhostDataOutsideDomain: unrecognized bc string");
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

  class Proto_Conductivity_Op_Factory: public ch_amrlevelop_factory
  {
    
  private:
    /// Strong construction, emacs and gdb are  a winning formula for life.
    Proto_Conductivity_Op_Factory();
    
  public:

    /// no dynamic memory used here.   
    virtual ~Proto_Conductivity_Op_Factory() { }

    Chombo::Vector<int>                    m_ref_ratios;
    Chombo::Vector<ch_dbl>                 m_grids;
    Chombo::Vector<pair<ch_dom, double> >  m_domain_dx;
    double                                 m_coarsest_dx;
    string                                 m_domBC;
    double                                 m_alpha;
    double                                 m_beta;
    int                                    m_max_coarse;

    
    ///the only constructor
    Proto_Conductivity_Op_Factory(
      const ch_dom                 &        a_coarsest_dom,
      const Chombo::Vector<ch_dbl> &        a_grids,
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
        ch_mayday::Error("ConductivityFactory::mgnewop: bogus domain");
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
        Proto_Conductivity_Op* derivedOp = 
          new Proto_Conductivity_Op(m_alpha,
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

      Proto_Conductivity_Op* derivedOp = 
        new Proto_Conductivity_Op(m_alpha,
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
