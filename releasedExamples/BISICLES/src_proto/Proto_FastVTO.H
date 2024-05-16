
#include "AMRMultiGrid_Proto_Interface.H"
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Proto_LoHiCenterStencil.H"
#include "Vector.H"

#ifndef __Proto_FastVTO_H__
#define __Proto_FastVTO_H__
/****************/
PROTO_KERNEL_START 

unsigned int axbyFastVTOF(int                       a_pt[DIM],
                        Proto::Var<Real, DIM>     a_lph,
                        Proto::Var<Real, DIM>     a_x,
                        Proto::Var<Real, DIM>     a_y,
                        Real                      a_a,
                        Real                      a_b,
                        int                       a_iprint)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    double xval = a_x(idir);
    double yval = a_y(idir);
    double lph  = a_a*xval + a_b*yval;
    a_lph(idir) = lph;
  }
  return 0;

}
PROTO_KERNEL_END(axbyFastVTOF, axbyFastVTO) 
///a_flux = eta( grad u + grad u T) + lambda I div U;
PROTO_KERNEL_START 
unsigned int  fluxFuncFastVTOF(int                     a_pt[DIM],
                             Proto::Var<Real, 1>     a_flux,
                             Proto::Var<Real, 1>     a_grad_u,
                             Proto::Var<Real, 1>     a_grad_u_t,
                             Proto::Var<Real, 1>     a_div_u,
                             Proto::Var<Real, 1>     a_eta,
                             Proto::Var<Real, 1>     a_lambda,
                             int                     a_face_dir,
                             int                     a_vel_comp)
{
  a_flux(0) = a_eta(0)*(a_grad_u(0) + a_grad_u_t(0));
  if(a_face_dir == a_vel_comp)
  {
    a_flux(0) += a_lambda(0)*a_div_u(0);
  }
  return 0;
}
PROTO_KERNEL_END(fluxFuncFastVTOF, fluxFuncFastVTO)


///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
PROTO_KERNEL_START 
unsigned int  operaFastVTOF(int                       a_pt[DIM],
                          Proto::Var<Real, DIM>     a_lph,
                          Proto::Var<Real, DIM>     a_phi,
                          Proto::Var<Real,  1 >     a_acoef,
                          Proto::Var<Real, DIM>     a_xincr,
                          Proto::Var<Real, DIM>     a_yincr,
                          Proto::Var<Real, DIM>     a_zincr,
                          Real                      a_alpha,
                          Real                      a_beta)
{
  //zincr is set to zero in 2d
  for(int idir = 0; idir < DIM; idir++)
  {
    double  divergence = a_xincr(idir) + a_yincr(idir) + a_zincr(idir);
    double  phival     =   a_phi(idir);  
    double  acoval     = a_acoef(0);   //same acoef for all variables
    double  lphval     = a_alpha*acoval*phival + a_beta*divergence;
    a_lph(idir) = lphval;
  }
  return 0;
}
PROTO_KERNEL_END(operaFastVTOF, operaFastVTO)

/// 
/**
   I am modelling this on lib/src/ViscousTensorOp but I will
   not account for differences in coefficients from one side of the cell
   to the other.   
   The VTO fortran amounts to 
   invlam(icomp) = acoef*a_alpha;
   for(int face_dir = 0; face_dir < SpaceDim; face_dir++)
   {
   double incr = 0;
   if(icomp == face_dir)
   {
   incr = beta*(2*eta(i+1/2) + lambda(i+1/2)
   +2*eta(i-1/2) + lambda(i-1/2))/(dx*dx);
   }
   else
   {
   incr = beta*(eta(i+1/2) + 2*eta(i-1/2));
   }
      
   inv_lam(idir) -= incr;
   }
   lam(idir) = -1/invlam(idir)

   Forall likes things to be at the same point so we will
   We will approximate this by;
   invlam(icomp) = acoef*a_alpha;
   for(int face_dir = 0; face_dir < SpaceDim; face_dir++)
   {
   double incr = 0;
   if(icomp == face_dir)
   {
   incr = beta*(4*eta(i-1/2) + 2*lambda(i-1/2))/(dx*dx);
                     
   }
   else
   {
   incr = 2*beta*eta(i-1/2);
   }
      
   inv_lam(idir) -= incr;
   }
   lam(idir) = 1/invlam(idir)
**/
PROTO_KERNEL_START 
unsigned int  relaxCoefFastVTOF(int                        a_pt[DIM],
                              Proto::Var<Real,  DIM>     a_relCoef,
                              Proto::Var<Real,  1  >     a_acoef,
                              Proto::Var<Real,  1  >     a_eta_x,
                              Proto::Var<Real,  1  >     a_eta_y,
                              Proto::Var<Real,  1  >     a_eta_z,
                              Proto::Var<Real,  1  >     a_lam_x,
                              Proto::Var<Real,  1  >     a_lam_y,
                              Proto::Var<Real,  1  >     a_lam_z,
                              double                     a_alpha,
                              double                     a_beta,
                              double                     a_dx)
{
  //this this is approximate, might have to have a real safety here
  double safety = 0.9; //same as ViscousTensorOpF.ChF
  double inv_lam[DIM];
  double eta_val[DIM];
  double lam_val[DIM];
  eta_val[0] = a_eta_x(0) ;
  lam_val[0] = a_lam_x(0) ;
  
  eta_val[1] = a_eta_y(0) ;
  lam_val[1] = a_lam_y(0) ;
#if DIM==3
  eta_val[2] = a_eta_z(0) ;
  lam_val[2] = a_lam_z(0) ;
#endif
  for(int icomp = 0; icomp < DIM; icomp++)
  {
    inv_lam[icomp] = a_alpha*a_acoef(0);
    for(int face_dir = 0; face_dir < DIM; face_dir++)
    {
      double incr = 0;
      if(icomp == face_dir)
      {
        incr = a_beta*(4.*eta_val[face_dir] + 2.*lam_val[face_dir])/(a_dx*a_dx);
      }
      else
      {
        incr = a_beta*(2.*eta_val[face_dir])/(a_dx*a_dx);
      }
      inv_lam[icomp] -= incr;
    }//end loop over face_dir
    a_relCoef(icomp) = safety/inv_lam[icomp];
  }

  return 0;
}
PROTO_KERNEL_END(relaxCoefFastVTOF, relaxCoefFastVTO)


///
/**
   phi = phi + lambda*(lphi - rho) 
   but only where 
   sum(i+j+k)%2== iredblack
**/
PROTO_KERNEL_START 
unsigned int  gsrbFastVTOF(int                       a_pt[DIM],
                         Proto::Var<Real, DIM>     a_phi,
                         Proto::Var<Real, DIM>     a_res,
                         Proto::Var<Real, DIM>     a_lam,
                         int                       a_iredBlack)
{
  int sumpt = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    sumpt += a_pt[idir];
  }
  if(sumpt%2 == a_iredBlack)
  {
    for(int icomp = 0; icomp < DIM; icomp++)
    {
      Real lambda = a_lam(icomp);
      Real phival = a_phi(icomp);
      Real resval = a_res(icomp);

      a_phi(icomp) = phival + lambda*resval;
    }
  }
  return 0;

}
PROTO_KERNEL_END(gsrbFastVTOF, gsrbFastVTO) 
/// The new namespace that all the cool operators want in on.   
namespace PrCh_AMR_Elliptic
{
  ///trying to stick with mostly Proto data holders here
  typedef Proto::Stencil<double>                                                              pr_sten;
  typedef Proto::InterpStencil<double>                                                        pr_interp_sten;
  typedef Proto::Point                                                                        pr_pt;
  typedef Proto::Box                                                                          pr_box;
  typedef Proto::BoxData<double, 1  >                                                         pr_box_data_sca;
  typedef Proto::BoxData<double, DIM>                                                         pr_box_data_vec;
  typedef Proto::BoxData<double, DIM*DIM>                                                     pr_box_data_ten;
  typedef Proto::DisjointBoxLayout                                                            pr_dbl;
  typedef Proto::Shift                                                                        pr_shift;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data_sca;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >   pr_xfac_data_sca;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >   pr_yfac_data_sca;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >   pr_zfac_data_sca;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data_vec;
  typedef Proto::LevelBoxData<double,  DIM*DIM, Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data_ten;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >   pr_xfac_data_vec;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >   pr_yfac_data_vec;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >   pr_zfac_data_vec;
                                                                 
  ///AMRMultigrid forces me to use some of these.
  typedef Chombo::ProblemDomain                                                               ch_dom; // 
  typedef Chombo::MayDay                                                                      ch_mayday;
  typedef Chombo::AMRLevelOpFactory<pr_cell_data_vec>                                         ch_amrlevelop_fact_vec;
  typedef Chombo::AMRLevelOp<       pr_cell_data_vec>                                         ch_amrlevelop_vec;
  typedef Chombo::MGLevelOp<        pr_cell_data_vec>                                         ch_mglevelop_vec;

  typedef PrCh_Lo_Hi_Center::LHC_MetaData     lo_hi_center_t;

  /// Viscous Tensor operator with variable coefficients.
  /**
     vel is a vector.
     alpha, beta are constants.
     lam = lam(x), eta = eta(x)
     I is the DIMxDIM identity matrix.
     L(vel) = alpha*a(x) vel + beta*(Div F(vel))
     F(vel) = eta*(grad(vel) + grad(vel)^T)  + lam* I div(vel)
  **/
  class Proto_FastVTO: public   Base_Proto_AMRLevelOp<DIM>
  {
  private:

    ///Because there a factory class, we can dispense with weak construction
    Proto_FastVTO();

  public:
    virtual ~Proto_FastVTO()
    {
    }


    
    shared_ptr<pr_cell_data_sca>  m_acoef;
    shared_ptr<pr_xfac_data_sca>  m_eta_x;
    shared_ptr<pr_yfac_data_sca>  m_eta_y;
    shared_ptr<pr_zfac_data_sca>  m_eta_z;
    shared_ptr<pr_xfac_data_sca>  m_lam_x;
    shared_ptr<pr_yfac_data_sca>  m_lam_y;
    shared_ptr<pr_zfac_data_sca>  m_lam_z;
    shared_ptr<pr_cell_data_vec>  m_relcoef;
    shared_ptr<pr_cell_data_vec>  m_resid;
    shared_ptr<pr_cell_data_vec>  m_lphi;

    double    m_alpha;
    double    m_beta;
    double    m_dx;
    string    m_domBC;
    ///amr hooks that do not need valid values when they do not apply
    int       m_refToFine;   
    int       m_refToCoar;
    pr_dbl    m_grids;
    
    ///
    Proto_FastVTO(
      shared_ptr<pr_cell_data_sca>  a_acoef,
      shared_ptr<pr_xfac_data_sca>  a_eta_x,
      shared_ptr<pr_yfac_data_sca>  a_eta_y,
      shared_ptr<pr_zfac_data_sca>  a_eta_z,
      shared_ptr<pr_xfac_data_sca>  a_lam_x,
      shared_ptr<pr_yfac_data_sca>  a_lam_y,
      shared_ptr<pr_zfac_data_sca>  a_lam_z,
      const pr_dbl&                 a_grids,
      const double&                 a_alpha,
      const double&                 a_beta,
      const double&                 a_dx,
      const int   &                 a_refToFine,
      const int   &                 a_refToCoar,
      const string&                 a_domBC)
    {
      CH_TIME("Proto_FastVTO::Proto_FastVTO");
      m_acoef     = a_acoef;

      m_eta_x     = a_eta_x;
      m_eta_y     = a_eta_y;
      m_eta_z     = a_eta_z;

      m_lam_x     = a_lam_x;
      m_lam_y     = a_lam_y;
      m_lam_z     = a_lam_z;

      m_grids     = a_grids;
      m_alpha     = a_alpha;
      m_beta      = a_beta;
      m_dx        = a_dx;
      m_refToFine = a_refToFine;
      m_refToCoar = a_refToCoar;
      m_domBC     = a_domBC;
      defineAuxilliaries();
      setRelaxationCoefficient();
    }

    ///
    virtual void
    homogeneousCFInterp(pr_lbd & a_phi)
    {
      CH_TIME("Proto_FastVTO::homogeneousCFInterp");
      Proto_Transition::DataUtilities<DIM>::homogeneousCFInterp<2>(a_phi, pr_pt::Ones(m_refToCoar));
    }
    
    /// amr thing
    virtual void
    fillGhostDataOverCoarseFineInterface(pr_lbd       &  a_data_fine, 
                                         const pr_lbd &  a_data_coar)
    {
      CH_TIME("Proto_FastVTO::fillDataOverCoarseFineInterface");
      int   reflect_sign_lo[DIM];
      int   reflect_sign_hi[DIM];  
      for(int idir = 0; idir < DIM; idir++)
      {
        if(m_domBC == string("Dirichlet"))
        {
          reflect_sign_lo[idir] = -1;
          reflect_sign_hi[idir] = -1;
        }
        else if(m_domBC == string("Neumann"))
        {
          reflect_sign_lo[idir] = 1;
          reflect_sign_hi[idir] = 1;
        }
        else
        {
          ch_mayday::Error("Proto_FastVTO::fillCFGhost: unrecognized bc string");
        }
      }
      typedef Proto_Transition::DataUtilities<DIM> pr_trans_data;

      pr_box domain_coar = a_data_coar.layout().domain().box();
      pr_pt ref_pt(pr_pt::Ones(m_refToCoar));
      pr_trans_data::fillGhostOverCoarseFineBoundaries<2>(a_data_fine,
                                                          a_data_coar,
                                                          domain_coar,
                                                          ref_pt,
                                                          reflect_sign_lo,
                                                          reflect_sign_hi);
      return;
      
    }
    /// output = aco*xval + bco*yval.   components defeat just doing this in the base class
    /**
       Standard function that I use for all arithmetic.  Because this is not fortran,
       it does not matter if any of the arguments point to overlapping things.
    **/
    virtual void 
    axby(pr_cell_data_vec         & a_lphi,
         const pr_cell_data_vec   & a_xval,
         const pr_cell_data_vec   & a_yval,
         double    a_aco,
         double    a_bco) 
    {
      CH_TIME("Proto_Fast_VTO::axby");
      auto grids =  a_lphi.layout();
      auto dit   =  a_lphi.begin();
      int iprint = 0; //debugging hook
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto       & lphfab =  a_lphi[dit[ibox]];
        const auto &   xfab =  a_xval[dit[ibox]];
        const auto &   yfab =  a_yval[dit[ibox]];
        auto valid          =   grids[dit[ibox]];
        Proto::forallInPlace_i(axbyFastVTO, valid, lphfab, xfab, yfab,
                               a_aco, a_bco, iprint);
      }
    }
    
    ///
    void defineAuxilliaries()
    {
      CH_TIME("Proto_FastVTO::defineAuxilliaries");
      m_relcoef = shared_ptr<pr_cell_data_vec>(new pr_cell_data_vec(m_grids, pr_pt::Zero()  ));
      m_resid   = shared_ptr<pr_cell_data_vec>(new pr_cell_data_vec(m_grids, pr_pt::Zero()  ));
      m_lphi    = shared_ptr<pr_cell_data_vec>(new pr_cell_data_vec(m_grids, pr_pt::Zero()  ));

    }
    
    ///
    void setRelaxationCoefficient()
    {
      CH_TIME("Pr_FastVTO::defineAndSetRelaxationCoefficient");
      auto dit   = m_relcoef->begin();

      auto& relco_lbd = *m_relcoef;
      ///fill a_lph with div( b(x) grad (phi))
      for(int ibox = 0; ibox <  dit.localSize(); ibox++)
      {
        const auto  & relcofab =   relco_lbd[dit[ibox]];
        const auto  & acoeffab =  (*m_acoef)[dit[ibox]];
        auto valid             =     m_grids[dit[ibox]];
        const auto& eta_x      =  (*m_eta_x)[dit[ibox]];
        const auto& eta_y      =  (*m_eta_y)[dit[ibox]];
        const auto& lam_x      =  (*m_lam_x)[dit[ibox]];
        const auto& lam_y      =  (*m_lam_y)[dit[ibox]];
#if DIM==3                        
        const auto& eta_z      =  (*m_eta_z)[dit[ibox]];
        const auto& lam_z      =  (*m_lam_z)[dit[ibox]];
#else
        //something fake to send in 2d
        const auto& eta_z      =  (*m_eta_x)[dit[ibox]];
        const auto& lam_z      =  (*m_lam_x)[dit[ibox]];
#endif
        ///compute relaxation coefficient
        Proto::forallInPlace_i(relaxCoefFastVTO, valid, 
                               relcofab, acoeffab,
                               eta_x, eta_y, eta_z,
                               lam_x, lam_y, lam_z,
                               m_alpha, m_beta, m_dx);
      }
      
      return;
    }

    
    /// find which var you care about in grad
    static int gradIndex(int a_velocity_component,
                         int a_derivative_direction)
    {
      int varnum = a_velocity_component*DIM + a_derivative_direction;
      return varnum;
    }//end function gradIndex
    
  
      
    ///
    /**
       Stencil for partial(u_facedir)/partial(x_facedir)
       This is a centered difference of velocity  only and relies upon 
       the velocity  having correct ghost cell values.
    **/
    int getNormalFaceGradStencil(lo_hi_center_t &  a_diff_velo_lhc, 
                                 const pr_box   &  a_lo_box,
                                 const pr_box   &  a_hi_box,
                                 const pr_box   &  a_ce_box,
                                 const int      &  a_has_lo,
                                 const int      &  a_has_hi,
                                 const int      &  a_face_dir,
                                 const int      &  a_diff_dir,
                                 const pr_box   &  a_flux_box)
    {
      ///this is just to make sure I am keeping this stuff straight.
      if(a_face_dir != a_diff_dir)
      {
        ch_mayday::Error("getNormStencil only for facedir==diffdir");
      }
      a_diff_velo_lhc.setBoxes(a_has_lo,
                              a_has_hi,
                              a_lo_box,
                              a_hi_box,
                              a_ce_box);

      pr_pt home     =  pr_pt::Zeros();
      //pr_pt righ_one = (+1)*pr_pt::Basis(a_diff_dir); //arg must be positive
      pr_pt   left_one = (-1)*pr_pt::Basis(a_diff_dir); //arg must be positive
      //center stencil
      pr_shift home_s(home);
      //pr_shift righ_one_s(righ_one);
      pr_shift   left_one_s(left_one);
      pr_sten home_sten(    home_s    , 1.);
      //pr_sten righ_one_sten(righ_one_s, 1.);
      pr_sten   left_one_sten(left_one_s, 1.);
      a_diff_velo_lhc.m_ce_sten = (1.0/m_dx)*(home_sten - left_one_sten);
      /**
         This stencil gets applied to velocity
         and relies upon velo having correct ghost cells 
         so only centered diff.   
      **/
      a_diff_velo_lhc.m_ce_box = a_flux_box;
      a_diff_velo_lhc.m_has_lo = 0;
      a_diff_velo_lhc.m_has_hi = 0;
      return 0;
    } ///end function getNormalFaceGradStencil

   ///
    /**
       Stencil for partial(u_veldir)/partial(x_diffdir) | facedir
       where facedir != diffdir.
       We average the nearby centered differences in the interior.
       We follow previous operators (VTO and Proto_VTO) which 
       extrapolate the cell centered gradient to ghost cells outside
       the domain.    Since we are not using m_grad anymore, I had
       to make the stencils just from phi.  Hence the wackiness at domain boundaries here.
       Normal gradients just use the ghost cell value of phi.
    **/
    int getTangentialFaceGradStencil(lo_hi_center_t &  a_diff_velo_lhc, 
                                     const pr_box   &  a_lo_box,
                                     const pr_box   &  a_hi_box,
                                     const pr_box   &  a_ce_box,
                                     const int      &  a_has_lo,
                                     const int      &  a_has_hi,
                                     const int      &  a_face_dir,
                                     const int      &  a_diff_dir,
                                     const pr_box   &  a_flux_box)
    {
      ///this is just to make sure I am keeping this stuff straight.
      if(a_face_dir == a_diff_dir)
      {
        ch_mayday::Error("getTanStencil one only for facedir!=diffdir.");
      }

      a_diff_velo_lhc.setBoxes(a_has_lo,
                               a_has_hi,
                               a_lo_box,
                               a_hi_box,
                               a_ce_box);

      pr_sten  minu_two_ef_plus_one_ed_sten(pr_shift(-2*pr_pt::Basis(a_face_dir) +  pr_pt::Basis(a_diff_dir)), 1.);
      pr_sten  minu_two_ef_minu_one_ed_sten(pr_shift(-2*pr_pt::Basis(a_face_dir) -  pr_pt::Basis(a_diff_dir)), 1.);
      pr_sten  minu_one_ef_plus_one_ed_sten(pr_shift(-1*pr_pt::Basis(a_face_dir) +  pr_pt::Basis(a_diff_dir)), 1.);
      pr_sten  minu_one_ef_minu_one_ed_sten(pr_shift(-1*pr_pt::Basis(a_face_dir) -  pr_pt::Basis(a_diff_dir)), 1.);
      pr_sten  plus_zer_ef_plus_one_ed_sten(pr_shift( 0*pr_pt::Basis(a_face_dir) +  pr_pt::Basis(a_diff_dir)), 1.);
      pr_sten  plus_zer_ef_minu_one_ed_sten(pr_shift( 0*pr_pt::Basis(a_face_dir) -  pr_pt::Basis(a_diff_dir)), 1.);
      pr_sten  plus_one_ef_plus_one_ed_sten(pr_shift( 1*pr_pt::Basis(a_face_dir) +  pr_pt::Basis(a_diff_dir)), 1.);
      pr_sten  plus_one_ef_minu_one_ed_sten(pr_shift( 1*pr_pt::Basis(a_face_dir) -  pr_pt::Basis(a_diff_dir)), 1.);


      /**
      ///For gradient based stuff:
      pr_sten     home_sten(pr_shift(        pr_pt::Zeros()         ), 1.0);
      pr_sten righ_one_sten(pr_shift( = ( 1)*pr_pt::Basis(a_face_dir), 1.0);
      pr_sten left_one_sten(pr_shift( = (-1)*pr_pt::Basis(a_face_dir), 1.0);
      pr_sten left_two_sten(pr_shift( = (-2)*pr_pt::Basis(a_face_dir), 1.0);

      ///Differences here are always in graddir direction
     
      //    ce_sten = (1./2.)*(home_sten + left_one_sten);
      //
      //  home_sten     = ((home+ Basisv(graddir) )                
      //                  -(home- Basisv(graddir) )                    )/2dx
      //  left_one_sten = ((home+ Basisv(graddir)  -  Basisv(face_dir))  
      //                  -(home+ Basisv(graddir)  -  Basisv(face_dir)))/2dx
      //
      //  ce_sten       = ((home+ Basisv(graddir)                     )                
      //                  -(home- Basisv(graddir)                     )
      //                  +(home+ Basisv(graddir)  -  Basisv(face_dir))  
      //                  -(home+ Basisv(graddir)  -  Basisv(face_dir)))/4dx
      **/
      a_diff_velo_lhc.m_ce_sten = (1.0/(4.*m_dx))*( plus_zer_ef_plus_one_ed_sten
                                                   -plus_zer_ef_minu_one_ed_sten
                                                  + minu_one_ef_plus_one_ed_sten 
                                                  - minu_one_ef_minu_one_ed_sten);
                                                    

      if(a_has_lo == 1)
      {
        /**
        /// in terms of gradient:
        /// lo_sten = (1./2.)*(3.*home_sten - righ_one_sten);
        //  home_sten     = ((home+ Basisv(graddir) )                
        //                  -(home- Basisv(graddir) )                   )/2dx
        //  righ_onesten = ((home+ Basisv(graddir)  + Basisv(face_dir))          
        //                 -(home- Basisv(graddir)  + Basisv(face_dir)) )/2dx
        /// lo_sten = (3/2)*home_sten - (1/2)righ_one_sten
        ///
        /// lo_sten = (3/(4dx))*(home+ Basisv(graddir)                   )
        ///          -(3/(4dx)) (home- Basisv(graddir)                   )
        ///          -(1/(4dx))*(home+ Basisv(graddir) + Basisv(face_dir))  
        ///          +(1/(4dx)) (home- Basisv(graddir) + Basisv(face_dir))
        ///
        ///
        ///
        **/
        
        a_diff_velo_lhc.m_lo_sten =
           (3.0/(4.*m_dx))*( plus_zer_ef_plus_one_ed_sten )
          -(3.0/(4.*m_dx))*( plus_zer_ef_minu_one_ed_sten )
          -(1.0/(4.*m_dx))*( plus_one_ef_plus_one_ed_sten )
          +(1.0/(4.*m_dx))*( plus_one_ef_minu_one_ed_sten );

      }
      if(a_has_hi == 1)
      {
        /**
        /// in terms of gradient:
           left_one_sten = ((home+ Basisv(graddir)  -   Basisv(face_dir))          
                           -(home- Basisv(graddir)  -   Basisv(face_dir)) )/2dx
           left_two_sten = ((home+ Basisv(graddir)  - 2*Basisv(face_dir))          
                           -(home- Basisv(graddir)  - 2*Basisv(face_dir)) )/2dx
           hi_sten = (1./2.)*(3.*left_one_sten - left_two_sten);
           hi_sten = (3/(4dx))*((home+ Basisv(graddir)  -   Basisv(face_dir)))          
                    -(3/(4dx))*((home- Basisv(graddir)  -   Basisv(face_dir)))          
                    -(1/(4dx))*((home+ Basisv(graddir)  - 2*Basisv(face_dir)))          
                    +(1/(4dx))*((home+ Basisv(graddir)  - 2*Basisv(face_dir)))          
        **/
        
        a_diff_velo_lhc.m_hi_sten =
           (3.0/(4.*m_dx))*( minu_one_ef_plus_one_ed_sten )
          -(3.0/(4.*m_dx))*( minu_one_ef_minu_one_ed_sten )
          -(1.0/(4.*m_dx))*( minu_two_ef_plus_one_ed_sten )
          +(1.0/(4.*m_dx))*( minu_two_ef_minu_one_ed_sten );
      }

      return 0;
      
      
    } ///end function getTangentialFaceGradStencil

    
    ///
    /**
       Average cell-centered gradients to the face centers
       Overwrite this with phi difference if face_dir == diff_dir.
    **/
    int
    getFaceDivAndGrad(pr_box_data_sca          & a_face_div,
                      pr_box_data_ten          & a_face_grad,
                      const pr_box_data_vec    & a_phi,
                      const pr_box             & a_flux_box,
                      const pr_box             & a_valid,
                      const int                & a_face_dir)
    {
      CH_TIME("Proto_FastVTO::getFaceDivAndGrad");

      pr_dom domain = m_grids.domain();
      ///first divergence 
      a_face_div.setVal(0.);
      pr_box_data_sca div_u_incr(a_flux_box);
      pr_box lo_box, hi_box, ce_box;
      int has_lo=4586, has_hi=4586;
      PrCh_Lo_Hi_Center::loHiCenterFace(lo_box, hi_box, ce_box,
                                        has_lo, has_hi, a_flux_box,
                                        a_valid, a_face_dir, domain);

      {
        CH_TIME("div_u_calculation");
        for(int diff_dir = 0; diff_dir < DIM; diff_dir++)
        {
          div_u_incr.setVal(0.);  //this is important.
          lo_hi_center_t diff_velo_lhc;
          if(diff_dir == a_face_dir)
          {
            /**
               Stencil for partial(u_facedir)/partial(x_facedir)
               This is a centered difference of velocity  only and relies upon 
               the velocity  having correct ghost cell values.
            **/
            getNormalFaceGradStencil(diff_velo_lhc, lo_box, hi_box, ce_box, has_lo, has_hi,
                                     a_face_dir, diff_dir, a_flux_box);
          
          }
          else
          {
            /**
               Stencil for partial(u_veldir)/partial(x_diffdir) | facedir
               where facedir != diffdir.
               We average the nearby centered differences.   This means 4
               nearby phis in the interior, 2 at domain boundaries.
            **/
            getTangentialFaceGradStencil(diff_velo_lhc, lo_box, hi_box, ce_box, has_lo, has_hi,
                                         a_face_dir, diff_dir, a_flux_box);
          }
          /**
             Component matches diff_dir because this is divergence
          **/
          pr_box_data_sca velo_comp = Proto::slice(a_phi, diff_dir);
          /**
             This is standard chombo dst, src so div_u_incr is  the output.
          **/
          PrCh_Lo_Hi_Center::loHiCenterApply(div_u_incr, velo_comp, &velo_comp, &velo_comp,
                                             diff_velo_lhc);
          a_face_div += div_u_incr;
        }
      }
      {
        CH_TIME("gradient calculation");
        //now for the gradient
        int ibreak = 4586;
        for(int diff_dir = 0; diff_dir < DIM; diff_dir++)
        {
          for(int velo_dir = 0; velo_dir < DIM; velo_dir++)
          {
            int grad_ind = gradIndex(velo_dir, diff_dir); // indicies match because this is divergence
            pr_box_data_sca velo_comp      = Proto::slice(a_phi      , velo_dir);
            pr_box_data_sca face_grad_comp = Proto::slice(a_face_grad, grad_ind);
            lo_hi_center_t diff_velo_lhc;
            if( (a_face_dir == diff_dir) )
            {
              getNormalFaceGradStencil(diff_velo_lhc, lo_box, hi_box, ce_box, has_lo, has_hi, // 
                                       a_face_dir, diff_dir, a_flux_box);
            }
            else
            {
              getTangentialFaceGradStencil(diff_velo_lhc, lo_box, hi_box, ce_box, has_lo, has_hi, // 
                                           a_face_dir, diff_dir, a_flux_box);
            }

            //this is standard chombo dest, src so face_grad_comp is the output
            PrCh_Lo_Hi_Center::loHiCenterApply(face_grad_comp, velo_comp, &velo_comp, &velo_comp,
                                               diff_velo_lhc);
            ibreak = 0;

          }   // end velo_dir loop
        }     // end diff_dir loop         
      }     // timer boundary
      return 0;
    }///end function getFaceDivAndGrad
       
    ///gets flux vector for a given face direction
    int
    getFluxDotNormal(pr_box_data_vec           & a_flux,
                     const pr_box_data_vec     & a_phi,
                     const pr_box_data_sca     & a_eta,
                     const pr_box_data_sca     & a_lam,
                     const pr_box              & a_valid,
                     const int                 & a_face_dir)
    {
      CH_TIME("Proto_Fast_VTO::getFluxDotNormal");

      //get face centered velocity divergence
      pr_box flux_box = a_valid;
      flux_box.m_high[a_face_dir] += 1;
      flux_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      pr_box_data_sca face_divu(flux_box);
      pr_box_data_ten face_grad(flux_box);
      getFaceDivAndGrad(face_divu, face_grad, a_phi, flux_box, a_valid,a_face_dir);
      
      ///F(u) = eta*(grad(u) + grad(u)^T) + lambda*I div(u)
      int ibreak = 4586;
      for(int icomp = 0; icomp < DIM; icomp++)
      {
        pr_box_data_sca flux_comp = Proto::slice(a_flux, icomp);

        int ind_gradu  = gradIndex(a_face_dir, icomp     );
        int ind_graduT = gradIndex(icomp     , a_face_dir);
        pr_box_data_sca grad_u   = Proto::slice(face_grad, ind_gradu);
        pr_box_data_sca grad_u_t = Proto::slice(face_grad, ind_graduT);
        
        ///a_flux = eta( grad u + grad u T) + lambda I div U;
        Proto::forallInPlace_i(fluxFuncFastVTO, flux_box, flux_comp,
                               grad_u, grad_u_t, face_divu,
                               a_eta, a_lam, a_face_dir, icomp);
        ibreak = 0;
      }
      return 0;
    }///end function getFluxDotNormal
    ///sets a_incr =  (flux_i+1/2 - flux_i-1/2)/m_dx
    int
    getDivergenceIncrement(pr_box_data_vec           & a_incr,
                           const pr_box_data_vec     & a_phi,
                           const pr_box_data_sca     & a_eta,
                           const pr_box_data_sca     & a_lam,
                           const pr_box              & a_valid,
                           const int                 & a_face_dir)
    {
      CH_TIME("Pr_FastVTO::getDivergenceIcrement");
      pr_box face_box = a_valid;
      face_box.m_high[a_face_dir] += 1;
      face_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      
      pr_box_data_vec  flux_dot_normal(face_box);
      getFluxDotNormal(flux_dot_normal, a_phi, a_eta, a_lam, a_valid, a_face_dir);
      

      //Proto::Stencil<double> home(Proto::Shift(Proto::Point::Zeros()), 1.0);
      
      pr_pt left_point = (-1)*pr_pt::Basis(a_face_dir);
      pr_pt righ_point = (+1)*pr_pt::Basis(a_face_dir);
      pr_shift righ(righ_point);
      pr_shift home(Proto::Point::Zeros());

      pr_sten home_sten(home, 1./m_dx);
      pr_sten righ_sten(righ, 1./m_dx);

      pr_sten righ_minus_home = righ_sten - home_sten;
      
      bool initToZero = true; double applyScale = 1;
      a_incr.setVal(-4586.);
      int ibreak = 0;
      for(int icomp = 0; icomp < DIM; icomp++)
      {
        pr_box_data_sca flux_comp = Proto::slice(flux_dot_normal, icomp);
        pr_box_data_sca incr_comp = Proto::slice(a_incr         , icomp);
        ///this arg list is src, dst so  incr_comp is the output
        righ_minus_home.apply(flux_comp, incr_comp,
                              a_valid , initToZero, applyScale);
        ibreak = 4586;
      }
      
      return 0;
    }///end function getDivergenceIncrement
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(pr_cell_data_vec         & a_lph,
                  const pr_cell_data_vec   & a_phi)
    {
      CH_TIME("Proto_FastVTO::applyOperator");
      auto dit   = a_lph.begin();

      //m_grad is no longer a thing
      //all stencils are directly from velo
      ///fill a_lph with div( b(x) grad (phi))
      for(int ibox = 0; ibox <  dit.localSize(); ibox++)
      {
        auto        & lphfab =      a_lph[dit[ibox]];
        const auto  & phifab =      a_phi[dit[ibox]];
        const auto  & acofab = (*m_acoef)[dit[ibox]];
        auto valid           =    m_grids[dit[ibox]];
        const auto& eta_x    = (*m_eta_x)[dit[ibox]];
        const auto& eta_y    = (*m_eta_y)[dit[ibox]];
        const auto& lam_x    = (*m_lam_x)[dit[ibox]];
        const auto& lam_y    = (*m_lam_y)[dit[ibox]];
#if DIM==3                         
        const auto& eta_z    = (*m_eta_z)[dit[ibox]];
        const auto& lam_z    = (*m_lam_z)[dit[ibox]];
#endif        

        //the setvals are meaningful here. 
        pr_box_data_vec xincr(valid); xincr.setVal(0.);
        pr_box_data_vec yincr(valid); yincr.setVal(0.);
        pr_box_data_vec zincr(valid); zincr.setVal(0.); 
        
        getDivergenceIncrement(xincr, phifab, eta_x, lam_x,  valid, 0);
        getDivergenceIncrement(yincr, phifab, eta_y, lam_y,  valid, 1);
#if DIM==3                                                        
        getDivergenceIncrement(zincr, phifab, eta_z, lam_z,  valid, 2);
#endif
        ///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
        Proto::forallInPlace_i(operaFastVTO, valid,
                               lphfab, phifab, acofab,
                               xincr, yincr, zincr,
                               m_alpha, m_beta);
      }
    }///end function applyOperator
    
    ///phi+= lambda*resid on red black
    virtual void
    relax(pr_cell_data_vec       & a_phi,
          const pr_cell_data_vec & a_rhs,
          int a_iterations)
    {
      CH_TIME("Pr_FastVTO::relax");
      
      pr_cell_data_vec& resid = *m_resid;

      auto  dit       =  a_rhs.begin();
      auto& relco_lbd = (*m_relcoef);
      // gdb debugging hook.  semantically does nothing.
      int idebug = 0; 
      for(int iiter = 0; iiter < a_iterations;  iiter++)
      {
        for(int iredblack = 0; iredblack < 2; iredblack++)
        {
          //begin debug
          //Chombo::pout() << "iiter = " << iiter << ", iredblack = " << iredblack << endl;
          //end debug
          residual(resid, a_phi, a_rhs);
          for(int ibox = 0; ibox < dit.localSize(); ibox++)
          {
            const auto & resfab =       resid[dit[ibox]];
            const auto & phifab =       a_phi[dit[ibox]];
            const auto & lamfab =   relco_lbd[dit[ibox]];
            auto valid          =     m_grids[dit[ibox]];

            Proto::forallInPlace_i(gsrbFastVTO, valid,
                                   phifab, resfab, lamfab,
                                   iredblack);
          }
          idebug = 0; //stop here for end of dataiter  iteration
        }
        idebug = 0;   //stop here for end of red-black iteration
      }
      idebug = 0;     //stop here for end of relax 
    } //end functtion relax

    
    /// precond is just relax here.   
    virtual void
    preCond(pr_cell_data_vec        & a_cor,
            const pr_cell_data_vec  & a_res)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::preCond");
      /**
         ViscousTensorOp only relaxes once.
      **/
      int num_relax = 1;

      this->relax(a_cor, a_res, num_relax);
    }
    /// domain boundary conditions.
    /** 
        Anything can be hidden in a string, so the interface is general enough.
        For now, we shall restrict ourselves to the big two: Neumann and Dirichlet.
    **/
    virtual void
    fillGhostDataOutsideDomain(pr_cell_data_vec          & a_phi,
                               const ch_dom              & a_domain)
    {
      CH_TIME("Pr_FastVTO::fillGhostdataOutsideDomain");
      typedef Proto_Transition::DataUtilities<DIM> pr_trans_data_vec;
      pr_dom domain = ProtoCh::getProtoDomain(a_domain);
      pr_box dombox = domain.box();
      int force_num_ghost = 1;
      for(int idir = 0; idir < DIM; idir++)
      {
        for(int hilo = 0; hilo < 2; hilo++)
        {
          if(m_domBC == string("Dirichlet"))
          {
            pr_trans_data_vec::domainGhostReflectOdd( a_phi, dombox, idir, hilo, force_num_ghost);
          }
          else if(m_domBC == string("Neumann"))
          {
            pr_trans_data_vec::domainGhostReflectEven(a_phi, dombox, idir, hilo, force_num_ghost);
          }
          else
          {
            ch_mayday::Error("Proto_FastVTO::fillGhostDataOutsideDomain: unrecognized bc string");
          }
        }
      }
    } //end functtion fillGhostDataOutsideDomain


    ///amr thing
    virtual int refToCoarser()
    {
      return m_refToCoar;
    }

    ///used in refluxing
    /**
       This is the interface that refluxFineCorrection uses. 
       Fill in the values. memory is already allocated.
       ibox is in relation to the operator's grids
       (so you can find your variable coefficients if you have them).
     **/
    virtual void getFlux(pr_box_data       & a_flux,
                         const pr_box_data & a_phi,
                         const pr_box      & a_face_box,
                         int a_face_dir,  int a_ibox)
    {
      CH_TIME("Proto_FastVTO::getFlux");
      auto dit   = m_acoef->begin();
      pr_box            valid =   m_grids[dit[a_ibox]];
      pr_box_data_sca * eta_ptr = NULL;
      pr_box_data_sca * lam_ptr = NULL;
      if(a_face_dir == 0)
      {
        eta_ptr = &( (*m_eta_x)[dit[a_ibox]] );
        lam_ptr = &( (*m_lam_x)[dit[a_ibox]] );
      }
      else if(a_face_dir == 1)
      {
        eta_ptr = &( (*m_eta_y)[dit[a_ibox]] );
        lam_ptr = &( (*m_lam_y)[dit[a_ibox]] );
      }
#if DIM==3
      else if(a_face_dir == 2)
      {
        eta_ptr = &( (*m_eta_z)[dit[a_ibox]] );
        lam_ptr = &( (*m_lam_z)[dit[a_ibox]] );
      }
#endif
      else
      {
        ch_mayday::Error("Proto_VTO::getFlux: bogus face_dir");
      }

      getFluxDotNormal(a_flux,a_phi, *eta_ptr, *lam_ptr, valid, a_face_dir);
    }

    ///
    double getDx() const
    {
      return m_dx;
    }
  }; ///end class Proto_FastVTO

  class Proto_FastVTO_Factory: public ch_amrlevelop_fact_vec
  {
    
  private:
    /// Strong construction, emacs, gdb, and Oxford commas are a winning formula for life.
    Proto_FastVTO_Factory();
    
  public:

    virtual ~Proto_FastVTO_Factory() { }

    vector<int>                                m_ref_ratios;
    vector<pr_dbl>                             m_grids;
    vector<pair<ch_dom, double> >              m_domain_dx;
    double                                     m_coarsest_dx;
    string                                     m_domBC;
    double                                     m_alpha;
    double                                     m_beta;
    int                                        m_max_coarse;
    vector< shared_ptr<pr_cell_data_sca> >     m_acoef;
    vector< shared_ptr<pr_xfac_data_sca> >     m_eta_x;
    vector< shared_ptr<pr_yfac_data_sca> >     m_eta_y;
    vector< shared_ptr<pr_zfac_data_sca> >     m_eta_z;
    vector< shared_ptr<pr_xfac_data_sca> >     m_lam_x;
    vector< shared_ptr<pr_yfac_data_sca> >     m_lam_y;
    vector< shared_ptr<pr_zfac_data_sca> >     m_lam_z;

    
    ///the only constructor
    Proto_FastVTO_Factory(
      const ch_dom                                 &  a_coarsest_dom,
      const vector<pr_dbl>                         &  a_grids,
      const vector<int>                            &  a_ref_ratios,
      const double                                 &  a_coarsest_dx,
      const string                                 &  a_domBC,
      const double                                 &  a_alpha,
      const double                                 &  a_beta,
      const vector< shared_ptr<pr_cell_data_sca> > &  a_acoef,
      const vector< shared_ptr<pr_xfac_data_sca> > &  a_eta_x,
      const vector< shared_ptr<pr_yfac_data_sca> > &  a_eta_y,
      const vector< shared_ptr<pr_zfac_data_sca> > &  a_eta_z,
      const vector< shared_ptr<pr_xfac_data_sca> > &  a_lam_x,
      const vector< shared_ptr<pr_yfac_data_sca> > &  a_lam_y,
      const vector< shared_ptr<pr_zfac_data_sca> > &  a_lam_z)
    : ch_amrlevelop_fact_vec()
    {
      CH_TIME("Pr_FastVTO_Factory_Constructor");
      m_max_coarse = 2;
      m_ref_ratios    =  a_ref_ratios  ;
      m_grids         =  a_grids;
      m_coarsest_dx   =  a_coarsest_dx ;   
      m_domBC         =  a_domBC       ;         
      m_alpha         =  a_alpha       ;   
      m_beta          =  a_beta        ;

      m_acoef = a_acoef;
      m_eta_x = a_eta_x;
      m_eta_y = a_eta_y;
      m_eta_z = a_eta_z;
      m_lam_x = a_lam_x;
      m_lam_y = a_lam_y;
      m_lam_z = a_lam_z;

      
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
        ch_mayday::Error("Viscous_TensorFactory::mgnewop: bogus domain");
      }
      return iref;
    } //end function get_ref_index

    ///
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


    /// for averaging stuff based in cells
    /**
       set coarse value to arithmetic average of all fine cell values.
    **/
    static void
    averageDownCell(pr_box_data_sca       & a_coar,
                    const pr_box_data_sca & a_fine,
                    const pr_box          & a_coar_valid,
                    int                     a_ref_rat)
    {
      ///got more generally useful
      PrChTools::AveragingUtilities::averageDownCell(a_coar, a_fine, a_coar_valid, a_ref_rat);
    }
    /// for averaging stuff based in faces 
    /**
       set coarse value to arithmetic average of all fine face values.   Only faces in the same facedir plane are used.
       The valid box is cell centered.
    **/
    static void
    averageDownFace(pr_box_data_sca       & a_coar,
                    const pr_box_data_sca & a_fine,
                    const pr_box      & a_coar_valid,
                    int                 a_ref_rat,
                    int                 a_face_dir)
    {
      PrChTools::AveragingUtilities::averageDownFace(a_coar, a_fine, a_coar_valid, a_ref_rat, a_face_dir);
    }

    virtual void
    getCoefficients(
      shared_ptr<pr_cell_data_sca>  & a_acoef,
      shared_ptr<pr_xfac_data_sca>  & a_eta_x,
      shared_ptr<pr_yfac_data_sca>  & a_eta_y,
      shared_ptr<pr_zfac_data_sca>  & a_eta_z,
      shared_ptr<pr_xfac_data_sca>  & a_lam_x,
      shared_ptr<pr_yfac_data_sca>  & a_lam_y,
      shared_ptr<pr_zfac_data_sca>  & a_lam_z,
      int                         a_iref_index,
      int                         a_mg_ref_rat)
    {
      CH_TIME("Pr_FastVTO_Factory_getCoefficients");
      if(a_mg_ref_rat == 1)
      {
        a_acoef = m_acoef[a_iref_index];
        a_eta_x = m_eta_x[a_iref_index];
        a_eta_y = m_eta_y[a_iref_index];
        a_eta_z = m_eta_z[a_iref_index];

        a_lam_x = m_lam_x[a_iref_index];
        a_lam_y = m_lam_y[a_iref_index];
        a_lam_z = m_lam_z[a_iref_index];

      }
      else
      {
        auto fine_dbl = m_grids[a_iref_index];
        //auto coar_dbl = fine_dbl.coarsen(pr_pt::Ones(a_mg_ref_rat));
        pr_dbl coar_dbl = fine_dbl.coarsen(pr_pt::Ones(a_mg_ref_rat));
        a_acoef = shared_ptr<pr_cell_data_sca>(new pr_cell_data_sca(coar_dbl,  pr_pt::Zeros()));

        a_eta_x = shared_ptr<pr_xfac_data_sca>(new pr_xfac_data_sca(coar_dbl,  pr_pt::Zeros()));
        a_eta_y = shared_ptr<pr_yfac_data_sca>(new pr_yfac_data_sca(coar_dbl,  pr_pt::Zeros()));
#if DIM==3
        a_eta_z = shared_ptr<pr_zfac_data_sca>(new pr_zfac_data_sca(coar_dbl,  pr_pt::Zeros()));
#endif        
        a_lam_x = shared_ptr<pr_xfac_data_sca>(new pr_xfac_data_sca(coar_dbl,  pr_pt::Zeros()));
        a_lam_y = shared_ptr<pr_yfac_data_sca>(new pr_yfac_data_sca(coar_dbl,  pr_pt::Zeros()));
#if DIM==3
        a_lam_z = shared_ptr<pr_zfac_data_sca>(new pr_zfac_data_sca(coar_dbl,  pr_pt::Zeros()));
#endif        
        auto dit   =  a_acoef->begin();
        for(int ibox = 0; ibox < dit.localSize(); ibox++)
        {
          averageDownCell((*a_acoef)[dit[ibox]], (*m_acoef[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat);
          averageDownFace((*a_eta_x)[dit[ibox]], (*m_eta_x[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 0);
          averageDownFace((*a_eta_y)[dit[ibox]], (*m_eta_y[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 1);
#if DIM==3           
          averageDownFace((*a_eta_z)[dit[ibox]], (*m_eta_z[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 2);
#endif        
          averageDownFace((*a_lam_x)[dit[ibox]], (*m_lam_x[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 0);
          averageDownFace((*a_lam_y)[dit[ibox]], (*m_lam_y[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 1);
#if DIM==3     
          averageDownFace((*a_lam_z)[dit[ibox]], (*m_lam_z[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 2);
#endif        
        }//end loop over boxes
      }//  end averaging down case
    }//end function getCoefficients
      
    virtual inline ch_mglevelop_vec*
    MGnewOp(const ch_dom         & a_dom_fine,
            int                    a_depth,
            bool                   a_homo_only = true)
    {
      CH_TIME("Pr_FastVTO_Factory_MGnewOp");
      int iref = get_ref_index(a_dom_fine);

      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      
      int mgRefRat = 1;
      for(int idep = 0; idep < a_depth; idep++)
      {
        mgRefRat*= 2;
      }
      ch_mglevelop_vec* retval  = NULL;
      bool grids_coarsenable = m_grids[iref].coarsenable(pr_pt::Ones(2*mgRefRat));
      if((a_depth > 0) && (!grids_coarsenable))
      {
        retval = NULL;
      }
      else
      {
        shared_ptr<pr_cell_data_sca>  acoef;
        shared_ptr<pr_xfac_data_sca>  eta_x;
        shared_ptr<pr_yfac_data_sca>  eta_y;
        shared_ptr<pr_zfac_data_sca>  eta_z;
        shared_ptr<pr_xfac_data_sca>  lam_x;
        shared_ptr<pr_yfac_data_sca>  lam_y;
        shared_ptr<pr_zfac_data_sca>  lam_z;

        
        getCoefficients(acoef,
                        eta_x,eta_y,eta_z,
                        lam_x,lam_y,lam_z,
                        iref, mgRefRat);
        
        auto dblMG = acoef->layout();
        double loc_dx = m_domain_dx[iref].second; 
        loc_dx *= mgRefRat;
        Proto_FastVTO* derivedOp = 
          new Proto_FastVTO(acoef,
                                      eta_x,
                                      eta_y,
                                      eta_z,
                                      lam_x,
                                      lam_y,
                                      lam_z,
                                      dblMG,
                                      m_alpha,
                                      m_beta,
                                      loc_dx,
                                      loc_ref_fine,
                                      loc_ref_coar,
                                      m_domBC);
        retval = static_cast< ch_mglevelop_vec* >(derivedOp);
      } //end case (we can return an operator)
      return retval;
    } //end function mgnewop

    ///
    virtual inline ch_amrlevelop_vec* AMRnewOp(const ch_dom & a_indexSpace)
    {
      CH_TIME("Pr_FastVTO_Factory_AMRnewOp");
      int iref = get_ref_index(a_indexSpace);
      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      double loc_dx = m_domain_dx[iref].second;

      Proto_FastVTO* derivedOp = 
        new Proto_FastVTO(m_acoef[iref],
                                    m_eta_x[iref],
                                    m_eta_y[iref],
                                    m_eta_z[iref],
                                    m_lam_x[iref],
                                    m_lam_y[iref],
                                    m_lam_z[iref],
                                    m_grids[iref],
                                    m_alpha,
                                    m_beta,
                                    loc_dx,
                                    loc_ref_fine,
                                    loc_ref_coar,
                                    m_domBC);
      return  static_cast< ch_amrlevelop_vec* >(derivedOp);
    } //end function AMRnewOp

    ///
    virtual inline int refToFiner(const ch_dom & a_indexSpace) const
    {
      int iref = get_ref_index(a_indexSpace);
      int ref_rat = m_ref_ratios[iref];
      return ref_rat;
    }//end function refToFiner
  }; //end class Proto_FastVTO_Factory
}    //end namespace Pr_Ch_AMR_AMRElliptic
#endif
