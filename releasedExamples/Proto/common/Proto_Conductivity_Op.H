
#include "AMRMultiGrid_Proto_Interface.H"
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Vector.H"

#ifndef __Proto_Conductivity_Op_H
#define __Proto_Conductivity_Op_H

///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
PROTO_KERNEL_START 
unsigned int  operaConduF(int                     a_pt[DIM],
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
  double  divergence = a_xincr(0) + a_yincr(0) + a_zincr(0);
  double  phival     = a_phi(0);  
  double  acoval     = a_acoef(0);  
  double  lphval     = a_alpha*acoval*phival + a_beta*divergence;
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
unsigned int  lambdaConduF(int                     a_pt[DIM],
                           Proto::Var<Real, 1>     a_lambda,
                           Proto::Var<Real, 1>     a_acoef,
                           Proto::Var<Real, 1>     a_bco_x,
                           Proto::Var<Real, 1>     a_bco_y,
                           Proto::Var<Real, 1>     a_bco_z,
                           Real                    a_alpha,
                           Real                    a_beta,
                           Real                    a_dx)
{
  //we are approximating the diagonal.   This may have to go lower.
  // set == 1 for direct comparisons with Proto_Helmholtz_Op 
  Real safety = 1.0;

  Real factor = -2.*a_beta/(a_dx*a_dx);
  Real diagon = a_alpha*a_acoef(0);
  diagon += factor*a_bco_x(0);
  diagon += factor*a_bco_y(0);
#if DIM==3
  diagon += factor*a_bco_z(0);
#endif
  a_lambda(0) = safety/diagon;
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
  ///trying to stick with mostly Proto data holders here
  typedef Proto::Stencil<double>                                                       pr_sten;
  typedef Proto::InterpStencil<double>                                                 pr_interp_sten;
  typedef Proto::Point                                                                 pr_pt;
  typedef Proto::Box                                                                   pr_box;
  typedef Proto::BoxData<double, 1>                                                    pr_box_data;
  typedef Proto::DisjointBoxLayout                                                     pr_dbl;
  typedef Proto::Shift                                                                 pr_shift;
  typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data;
  typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >   pr_xfac_data;
  typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >   pr_yfac_data;
  typedef Proto::LevelBoxData<double, 1, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >   pr_zfac_data;
                                                              
  ///AMRMultigrid does force some of these.
  typedef Chombo::ProblemDomain                                                        ch_dom;
  typedef Chombo::MayDay                                                               ch_mayday;
  typedef Chombo::AMRLevelOpFactory<pr_cell_data>                                      ch_amrlevelop_factory;
  typedef Chombo::AMRLevelOp<       pr_cell_data>                                      ch_amrlevelop;
  typedef Chombo::MGLevelOp<        pr_cell_data>                                      ch_mglevelop;

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
      CH_TIME("Pr_Conductivity::defineAndSetRelaxationCoefficient");
      m_lambda = shared_ptr<pr_cell_data>(new pr_cell_data(m_grids, pr_pt::Zero() ));
      auto dit   = m_lambda->begin();
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

    ///sets a_incr =  (b_i+1/2*grad_phi_i+1/2 - b_i+1/2*grad_phi_i+1/2)/m_dx
    int
    getDivergenceIncrement(pr_box_data           & a_incr,
                           const pr_box_data     & a_phi,
                           const pr_box_data     & a_bcoef,
                           const pr_box          & a_valid,
                           const int             & a_face_dir)
    {
      CH_TIME("Pr_Conductivity::getDivergenceIcrement");
      

      //Proto::Stencil<double> home(Proto::Shift(Proto::Point::Zeros()), 1.0);
      
      pr_pt left_point = (-1)*pr_pt::Basis(a_face_dir);
      pr_shift left(left_point);
      pr_shift home(Proto::Point::Zeros());

      pr_sten hi(home,  1./m_dx);
      pr_sten lo(left, -1./m_dx);

      pr_sten grad_sten = hi - lo;
      pr_box face_box = a_valid;
      face_box.m_high[a_face_dir] += 1;
      face_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      pr_box_data bco_grad_phi(face_box);
      bool initToZero = true; double applyScale = 1;
      pr_box_data& phicast  = (pr_box_data &) a_phi;
      ///remember that this goes src, dst instead of chombo standard dst,src
      grad_sten.apply(phicast, bco_grad_phi, face_box, initToZero, applyScale);
      bco_grad_phi *= a_bcoef;
      ///remember that this goes src, dst instead of chombo standard dst,src
      grad_sten.apply(bco_grad_phi, a_incr , a_valid , initToZero, applyScale);
      
      return 0;
    }
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(pr_lbd         & a_lph,
                  const pr_lbd   & a_phi)
    {
      CH_TIME("Proto_Conductivity_Op::applyOperator");
      auto dit   = a_lph.begin();

      ///fill a_lph with div( b(x) grad (phi))
      for(int ibox = 0; ibox <  dit.size(); ibox++)
      {
        auto        & lphfab =      a_lph[dit[ibox]];
        const auto  & phifab =      a_phi[dit[ibox]];
        const auto  & acofab = (*m_acoef)[dit[ibox]];
        auto valid           =    m_grids[dit[ibox]];
        const auto& bco_x    = (*m_bco_x)[dit[ibox]];
        const auto& bco_y    = (*m_bco_y)[dit[ibox]];
#if DIM==3             
        const auto& bco_z    = (*m_bco_z)[dit[ibox]];
#endif        
        pr_box_data xincr(valid); xincr.setVal(0.);
        pr_box_data yincr(valid); yincr.setVal(0.);
        pr_box_data zincr(valid); zincr.setVal(0.);
        
        getDivergenceIncrement(xincr, phifab, bco_x, valid, 0);
        getDivergenceIncrement(yincr, phifab, bco_y, valid, 1);
#if DIM==3                                       
        getDivergenceIncrement(zincr, phifab, bco_z, valid, 2);
#endif
        ///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
        Proto::forallInPlace_i(operaCondu, valid,
                               lphfab, phifab, acofab,
                               xincr, yincr, zincr,
                               m_alpha, m_beta);
      }
    }
    
    ///phi+= lambda*resid on red black
    virtual void
    relax(pr_lbd       & a_phi,
          const pr_lbd & a_rhs,
          int a_iterations)
    {
      CH_TIME("Pr_Conductivity::relax");
      
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
            const auto & resfab =       resid[dit[ibox]];
            const auto & phifab =       a_phi[dit[ibox]];
            const auto & lamfab = (*m_lambda)[dit[ibox]];
            auto valid          =     m_grids[dit[ibox]];
            int iprint = 0; // I left in a debugging hook.
            Proto::forallInPlace_i(gsrbCondu, valid, phifab, resfab, lamfab,
                                   iredblack, iprint);
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
      CH_TIME("Pr_Conductivity::fillGhostdataOutsideDomain");
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
  }; ///end class proto_conductivity_op

  class Proto_Conductivity_Op_Factory: public ch_amrlevelop_factory
  {
    
  private:
    /// Strong construction, emacs, and gdb are a winning formula for life.
    Proto_Conductivity_Op_Factory();
    
  public:

    virtual ~Proto_Conductivity_Op_Factory() { }

    vector<int>                            m_ref_ratios;
    vector<pr_dbl>                         m_grids;
    vector<pair<ch_dom, double> >          m_domain_dx;
    double                                 m_coarsest_dx;
    string                                 m_domBC;
    double                                 m_alpha;
    double                                 m_beta;
    int                                    m_max_coarse;
    vector< shared_ptr<pr_cell_data> >     m_acoef;
    vector< shared_ptr<pr_xfac_data> >     m_bco_x;
    vector< shared_ptr<pr_yfac_data> >     m_bco_y;
    vector< shared_ptr<pr_zfac_data> >     m_bco_z;

    
    ///the only constructor
    Proto_Conductivity_Op_Factory(
      const ch_dom                             &  a_coarsest_dom,
      const vector<pr_dbl>                     &  a_grids,
      const vector<int>                        &  a_ref_ratios,
      const double                             &  a_coarsest_dx,
      const string                             &  a_domBC,
      const double                             &  a_alpha,
      const double                             &  a_beta,
      const vector< shared_ptr<pr_cell_data> > &  a_acoef,
      const vector< shared_ptr<pr_xfac_data> > &  a_bco_x,
      const vector< shared_ptr<pr_yfac_data> > &  a_bco_y,
      const vector< shared_ptr<pr_zfac_data> > &  a_bco_z)           
    : ch_amrlevelop_factory()
    {
      CH_TIME("Pr_Conductivity_Factory_Constructor");
      m_max_coarse = 2;
      m_ref_ratios    =  a_ref_ratios  ;
      m_grids         =  a_grids;
      m_coarsest_dx   =  a_coarsest_dx ;   
      m_domBC         =  a_domBC       ;         
      m_alpha         =  a_alpha       ;   
      m_beta          =  a_beta        ;

      m_acoef = a_acoef;
      m_bco_x = a_bco_x;
      m_bco_y = a_bco_y;
      m_bco_z = a_bco_z;

      
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
      } //end loop over AMR levels
    }   //end factory constructor

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
    } //end function get_ref_index

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
    }//end function get_local_ref_ratios

    pr_sten getStencilFromRefBox(pr_box a_ref_box,
                                 int    a_ref_rat)
    {
      int numPts = a_ref_box.size();
      double weight = 1.0/double(numPts);
      pr_pt src_ref_ratio = pr_pt::Ones(a_ref_rat);
      pr_pt dst_ref_ratio = pr_pt::Ones(1);
      pr_pt dst_shift     = pr_pt::Zeros();
#if 0      
      /**
         This breaks due to a proto quirk that will probably be fixed soon.
         Right now, adding a stencil with non unity refinement ratio 
         to an empty stencil causes a runtime error.
         dtg 11-27-2023
      **/
      pr_sten whole_sten;
      for(auto boxit = a_ref_box.begin(); boxit != a_ref_box.end(); ++boxit)
      {
        pr_pt ref_pt = *boxit;
        pr_sten pt_sten(pr_shift(ref_pt), weight,
                        dst_ref_ratio, dst_shift, src_ref_ratio);
        whole_sten += pt_sten;
      }
    
#else
      ///This is the ugly workaround.
      vector<pr_sten> pt_stens(a_ref_box.size());
      int iref = 0;
      for(auto boxit = a_ref_box.begin(); boxit != a_ref_box.end(); ++boxit)
      {
        pr_pt ref_pt = *boxit;
        pr_sten pt_sten(pr_shift(ref_pt), weight,
                        dst_ref_ratio, dst_shift, src_ref_ratio);
        pt_stens[iref] = pt_sten;
        iref++;
      }
      pr_sten whole_sten(pt_stens[0]);
      for(int ivec = 1 /* important that this is = 1*/;
          ivec < pt_stens.size(); ivec++) 
      {
        whole_sten += pt_stens[ivec];
      }
#endif      
      return whole_sten;
    }

    /// for averaging stuff based in cells
    /**
       set coarse value to arithmetic average of all fine cell values.
     **/
    static void
    averageDownCell(pr_box_data       & a_coar,
                    const pr_box_data & a_fine,
                    const pr_box      & a_coar_valid,
                    int                 a_ref_rat)
    {
      pr_box unit_box(pr_pt::Zeros(), pr_pt::Ones());
      pr_pt pt_ref_rat = pr_pt::Ones( a_ref_rat);
      pr_box ref_box = unit_box.refine(pt_ref_rat);
      pr_sten whole_sten = getStencilFromRefBox(ref_box);

      
      pr_box_data& fine_cast = (pr_box_data&) a_fine;
      bool initToZero = true; double applyScale = 1.;
      whole_sten.apply(fine_cast, a_coar, a_coar_valid, initToZero, applyScale);
    }
    /// for averaging stuff based in faces 
    /**
       set coarse value to arithmetic average of all fine face values.   Only faces in the same facedir plane are used.
       The valid box is cell centered.
     **/
    static void
    averageDownFace(pr_box_data       & a_coar,
                    const pr_box_data & a_fine,
                    const pr_box      & a_coar_valid,
                    int                 a_ref_rat,
                    int                 a_face_dir)
    {
      pr_box unit_box(pr_pt::Zeros(), pr_pt::Ones());
      pr_pt pt_ref_rat = pr_pt::Ones( a_ref_rat);
      pr_box ref_box = unit_box.refine(pt_ref_rat);
      // only include faces on the same facedir plane
      ref_box.m_high[a_face_dir] = 0;
      ref_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      
      int numPts = ref_box.size();
      double weight = 1.0/double(numPts);
      pr_sten whole_sten = getStencilFromRefBox(ref_box);

      pr_box_data& fine_cast = (pr_box_data&) a_fine;
      bool initToZero = true; double applyScale = 1.;
      pr_box coar_face_valid = a_coar_valid;
      coar_face_valid.m_high[a_face_dir] += 1;
      coar_face_valid.recomputeSize();  // required to keep box data consistent after above shenannigans.
      
      whole_sten.apply(fine_cast, a_coar, coar_face_valid, initToZero, applyScale);

    }

    virtual void
    getCoefficients(
      shared_ptr<pr_cell_data>  & a_acoef,
      shared_ptr<pr_xfac_data>  & a_bco_x,
      shared_ptr<pr_yfac_data>  & a_bco_y,
      shared_ptr<pr_zfac_data>  & a_bco_z,
      int                         a_iref_index,
      int                         a_mg_ref_rat)
    {
      CH_TIME("Pr_Conductivity_Factory_getCoefficients");
      if(a_mg_ref_rat == 1)
      {
        a_acoef = m_acoef[a_iref_index];
        a_bco_x = m_bco_x[a_iref_index];
        a_bco_y = m_bco_y[a_iref_index];
        a_bco_z = m_bco_z[a_iref_index];
      }
      else
      {
        auto fine_dbl = m_grids[a_iref_index];
        //auto coar_dbl = fine_dbl.coarsen(pr_pt::Ones(a_mg_ref_rat));
        pr_dbl coar_dbl = fine_dbl.coarsen(pr_pt::Ones(a_mg_ref_rat));
        a_acoef = shared_ptr<pr_cell_data>(new pr_cell_data(coar_dbl,  pr_pt::Zeros()));
        a_bco_x = shared_ptr<pr_xfac_data>(new pr_xfac_data(coar_dbl,  pr_pt::Zeros()));
        a_bco_y = shared_ptr<pr_yfac_data>(new pr_yfac_data(coar_dbl,  pr_pt::Zeros()));
#if DIM==3        
        a_bco_z = shared_ptr<pr_zfac_data>(new pr_zfac_data(coar_dbl,  pr_pt::Zeros()));
#endif        
        auto dit   =  a_acoef->begin();
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          averageDownCell((*a_acoef)[dit[ibox]], (*m_acoef[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat);
          averageDownFace((*a_bco_x)[dit[ibox]], (*m_bco_x[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 0);
          averageDownFace((*a_bco_y)[dit[ibox]], (*m_bco_y[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 1);
#if DIM==3
          averageDownFace((*a_bco_z)[dit[ibox]], (*m_bco_z[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 2);
#endif        
        }//end loop over boxes
      }//  end averaging down case
    }//end function getCoefficients
      
    virtual inline ch_mglevelop*
    MGnewOp(const ch_dom         & a_dom_fine,
            int                    a_depth,
            bool                   a_homo_only = true)
    {
      CH_TIME("Pr_Conductivity_Factory_MGnewOp");
      int iref = get_ref_index(a_dom_fine);

      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      
      int mgRefRat = 1;
      for(int idep = 0; idep < a_depth; idep++)
      {
        mgRefRat*= 2;
      }
      ch_mglevelop* retval  = NULL;
      bool grids_coarsenable = m_grids[iref].coarsenable(pr_pt::Ones(2*mgRefRat));
      if((a_depth > 0) && (!grids_coarsenable))
      {
        retval = NULL;
      }
      else
      {
        shared_ptr<pr_cell_data>  acoef;
        shared_ptr<pr_xfac_data>  bco_x;
        shared_ptr<pr_yfac_data>  bco_y;
        shared_ptr<pr_zfac_data>  bco_z;
        
        getCoefficients(acoef,bco_x,bco_y,bco_z, iref, mgRefRat);
        auto dblMG = acoef->layout();
        double loc_dx = m_domain_dx[iref].second; 
        loc_dx *= mgRefRat;
        Proto_Conductivity_Op* derivedOp = 
          new Proto_Conductivity_Op(acoef,
                                    bco_x,
                                    bco_y,
                                    bco_z,
                                    dblMG,
                                    m_alpha,
                                    m_beta,
                                    loc_dx,
                                    loc_ref_fine,
                                    loc_ref_coar,
                                    m_domBC);
        retval = static_cast< ch_mglevelop* >(derivedOp);
      } //end case (we can return an operator)
      return retval;
    } //end function mgnewop

    ///
    virtual inline ch_amrlevelop* AMRnewOp(const ch_dom & a_indexSpace)
    {
      CH_TIME("Pr_Conductivity_Factory_AMRnewOp");
      int iref = get_ref_index(a_indexSpace);
      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      double loc_dx = m_domain_dx[iref].second;

      Proto_Conductivity_Op* derivedOp = 
        new Proto_Conductivity_Op(m_acoef[iref],
                                  m_bco_x[iref],
                                  m_bco_y[iref],
                                  m_bco_z[iref],
                                  m_grids[iref],
                                  m_alpha,
                                  m_beta,
                                  loc_dx,
                                  loc_ref_fine,
                                  loc_ref_coar,
                                  m_domBC);
      return  static_cast< ch_amrlevelop* >(derivedOp);
    } //end function AMRnewOp

    ///
    virtual inline int refToFiner(const ch_dom & a_indexSpace) const
    {
      int iref = get_ref_index(a_indexSpace);
      int ref_rat = m_ref_ratios[iref];
      return ref_rat;
    }
  };
} //end namespace Pr_Ch_AMR_AMRElliptic
#endif
