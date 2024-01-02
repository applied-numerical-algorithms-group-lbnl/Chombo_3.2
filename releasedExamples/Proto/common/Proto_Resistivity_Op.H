
#include "AMRMultiGrid_Proto_Interface.H"
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Proto_LoHiCenterStencil.H"
#include "Vector.H"

#ifndef __Proto_Resistivity_Op_H__
#define __Proto_Resistivity_Op_H__
/****************/
PROTO_KERNEL_START 
///always 3  variables here
unsigned int axbyResistF(int                       a_pt[DIM],
                         Proto::Var<Real,   3>     a_lph,
                         Proto::Var<Real,   3>     a_x,
                         Proto::Var<Real,   3>     a_y,
                         Real                      a_a,
                         Real                      a_b,
                         int                       a_iprint)
{
  for(int idir = 0; idir < 3; idir++)
  {
    double xval = a_x(idir);
    double yval = a_y(idir);
    double lph  = a_a*xval + a_b*yval;
    a_lph(idir) = lph;
  }
  return 0;

}
PROTO_KERNEL_END(axbyResistF, axbyResist) 


///a_lph = alpha*phi + beta(xincr + yincr + zincr)
PROTO_KERNEL_START 
unsigned int  operaResistF(int                     a_pt[DIM],
                           Proto::Var<Real, 3>     a_lph,
                           Proto::Var<Real, 3>     a_phi,
                           Proto::Var<Real, 3>     a_xincr,
                           Proto::Var<Real, 3>     a_yincr,
                           Proto::Var<Real, 3>     a_zincr,
                           Real                    a_alpha,
                           Real                    a_beta)
{
  //zincr is set to zero in 2d
  for(int idir = 0; idir < 3; idir++)
  {
    double  divergence = a_xincr(idir) + a_yincr(idir) + a_zincr(idir);
    double  phival     =   a_phi(idir);  
    double  lphval     = a_alpha*phival + a_beta*divergence;
    a_lph(idir) = lphval;
  }
  return 0;
}
PROTO_KERNEL_END(operaResistF, operaResist)

/// 
/**
   Ignoring differences within the cell.  ResistivityOp has lambda =
   alpha - beta*sum_dir(ehi + elo)/dx/dx We will set ehi = elo =
   e_i+/2 because forall does not do stencils and I do not want to do
   a bunch of shifting nonsense.  Note: You could send in ehi and elo
   as shifted and unshifted versions.  It would not be that
   complicated since there is only eta to send.  I think this will
   work, however, and if coefficients are varying that strongly, I
   doubt that the other version would save convergence.
   dtg 12/30/2023.
 **/
PROTO_KERNEL_START 
unsigned int  relaxCoefResistF(int                        a_pt[DIM],
                               Proto::Var<Real,  3  >     a_relCoef,
                               Proto::Var<Real,  1  >     a_eta_x,
                               Proto::Var<Real,  1  >     a_eta_y,
                               Proto::Var<Real,  1  >     a_eta_z,
                               double                     a_alpha,
                               double                     a_beta,
                               double                     a_dx)
{
  //this this is approximate, might have to have a real safety here
  double safety = 1.0;// matches resistivityOp
  double inv_lam[3];
  double eta_val[DIM];

  eta_val[0] = a_eta_x(0) ;
  eta_val[1] = a_eta_y(0) ;
#if DIM==3
  eta_val[2] = a_eta_z(0) ;
#endif
  for(int icomp = 0; icomp < 3; icomp++)
  {
    inv_lam[icomp] = a_alpha; //acoef always one here
    for(int face_dir = 0; face_dir < DIM; face_dir++)
    {
      double incr = a_beta*(2.*eta_val[face_dir])/(a_dx*a_dx);
      inv_lam[icomp] -= incr;
    }//end loop over face_dir
    a_relCoef(icomp) = safety/inv_lam[icomp];
  }

  return 0;
}
PROTO_KERNEL_END(relaxCoefResistF, relaxCoefResist)


///
/**
   phi = phi + lambda*(lphi - rho) 
   but only where 
   sum(i+j+k)%2== iredblack
**/
PROTO_KERNEL_START 
unsigned int  gsrbResistF(int                     a_pt[DIM],
                          Proto::Var<Real, 3>     a_phi,
                          Proto::Var<Real, 3>     a_res,
                          Proto::Var<Real, 3>     a_lam,
                          int                     a_iredBlack)
{
  int sumpt = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    sumpt += a_pt[idir];
  }
  if(sumpt%2 == a_iredBlack)
  {
    for(int icomp = 0; icomp < 3; icomp++)
    {
      Real lambda = a_lam(icomp);
      Real phival = a_phi(icomp);
      Real resval = a_res(icomp);

      a_phi(icomp) = phival + lambda*resval;
    }
  }
  return 0;

}
PROTO_KERNEL_END(gsrbResistF, gsrbResist) 
/// The hip new namespace that all the groovy operators want in on.   
namespace PrCh_AMR_Elliptic
{
  ///trying to stick with mostly Proto data holders here
  typedef Proto::Stencil<double>                                                        pr_sten;
  typedef Proto::InterpStencil<double>                                                  pr_interp_sten;
  typedef Proto::Point                                                                  pr_pt;
  typedef Proto::Box                                                                    pr_box;
  typedef Proto::BoxData<double, 1  >                                                   pr_box_data_scalar;
  //always 3 variables. can be two derivatives
  typedef Proto::BoxData<double, 3>                                                     pr_box_data_three;
  typedef Proto::BoxData<double, 3*DIM>                                                 pr_box_data_thrdim;
  typedef Proto::DisjointBoxLayout                                                      pr_dbl;
  typedef Proto::Shift                                                                  pr_shift;
  typedef Proto::LevelBoxData<double, 1    , Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   > pr_cell_data_scalar;
  typedef Proto::LevelBoxData<double, 1    , Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_0 > pr_xfac_data_scalar;
  typedef Proto::LevelBoxData<double, 1    , Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_1 > pr_yfac_data_scalar;
  typedef Proto::LevelBoxData<double, 1    , Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_2 > pr_zfac_data_scalar;
  //always 3 variables. can be two derivatives
  typedef Proto::LevelBoxData<double, 3    , Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   > pr_cell_data_three;
  typedef Proto::LevelBoxData<double, 3*DIM, Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   > pr_cell_data_thrdim;
  typedef Proto::LevelBoxData<double, 3    , Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_0 > pr_xfac_data_three;
  typedef Proto::LevelBoxData<double, 3    , Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_1 > pr_yfac_data_three;
  typedef Proto::LevelBoxData<double, 3    , Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_2 > pr_zfac_data_three;
                                                                 
  ///AMRMultigrid forces me to use some of these.
  typedef Chombo::ProblemDomain                                                         ch_dom; // 
  typedef Chombo::MayDay                                                                ch_mayday;
  typedef Chombo::AMRLevelOpFactory<pr_cell_data_three>                                 ch_amrlevelop_fact_three;
  typedef Chombo::AMRLevelOp<       pr_cell_data_three>                                 ch_amrlevelop_three;
  typedef Chombo::MGLevelOp<        pr_cell_data_three>                                 ch_mglevelop_three;

  typedef PrCh_Lo_Hi_Center::LHC_MetaData     lo_hi_center_t;

  /// Resistivity operator with variable coefficients.
  /**
     Define alpha, beta to be constant.
     Define eta = eta(x) to be the resisitivity coefficient.
     Given a magnetic field (always 3 variables) B, our operator L is defined as
     L(B)  = alpha B + beta*(divF),  where 
     F = eta((grad B - grad B^T) + I  (div B)).
 
     B always has 3 variables (B = (Bx, By, Bz)).  
     DIM  only affects what  derivatives are taken.
     Tensors (fluxes, gradients) are  of dimension 3*DIM.
     I is the 3xDIM identity matrix.
  **/
  class Proto_Resistivity_Op: public   Base_Proto_AMRLevelOp<3>
  {
  private:

    ///Because there a factory class, we can dispense with weak construction
    Proto_Resistivity_Op();

  public:
    virtual ~Proto_Resistivity_Op()
    {
    }

    shared_ptr<pr_xfac_data_scalar>  m_eta_x;
    shared_ptr<pr_yfac_data_scalar>  m_eta_y;
    shared_ptr<pr_zfac_data_scalar>  m_eta_z;
    shared_ptr<pr_cell_data_three>  m_relcoef;
    shared_ptr<pr_cell_data_thrdim>  m_grad;
    double    m_alpha;
    double    m_beta;
    double    m_dx;
    string    m_domBC;
    int       m_refToFine;   
    int       m_refToCoar;
    pr_dbl    m_grids;
    
    ///
    Proto_Resistivity_Op(
      shared_ptr<pr_xfac_data_scalar>  a_eta_x,
      shared_ptr<pr_yfac_data_scalar>  a_eta_y,
      shared_ptr<pr_zfac_data_scalar>  a_eta_z,
      const pr_dbl   &                 a_grids,
      const double   &                 a_alpha,
      const double   &                 a_beta,
      const double   &                 a_dx,
      const int      &                 a_refToFine,
      const int      &                 a_refToCoar,
      const string   &                 a_domBC)
    {
      CH_TIME("Proto_Resistivity_Op::Proto_Resistivity_Op");
      m_eta_x     = a_eta_x;
      m_eta_y     = a_eta_y;
      m_eta_z     = a_eta_z;

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

    /// output = aco*xval + bco*yval.   components defeat just doing this in the base class
    /**
       Standard function that I use for all arithmetic.  Because this is not fortran,
       it does not matter if any of the arguments point to overlapping things.
    **/
    virtual void 
    axby(pr_cell_data_three         & a_lphi,
         const pr_cell_data_three   & a_xval,
         const pr_cell_data_three   & a_yval,
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
        Proto::forallInPlace_i(axbyResist, valid, lphfab, xfab, yfab,
                               a_aco, a_bco, iprint);
      }
    }
    
    ///
    void defineAuxilliaries()
    {
      m_relcoef = shared_ptr<pr_cell_data_three>(new pr_cell_data_three(m_grids, pr_pt::Zero()  ));
      m_grad    = shared_ptr<pr_cell_data_thrdim>(new pr_cell_data_thrdim(m_grids, pr_pt::Ones(1) ));
    }
    
    ///
    void setRelaxationCoefficient()
    {
      CH_TIME("Pr_Resistivity::defineAndSetRelaxationCoefficient");
      auto dit   = m_relcoef->begin();

      auto& relco_lbd = *m_relcoef;
      ///fill a_lph with div( b(x) grad (phi))
      for(int ibox = 0; ibox <  dit.localSize(); ibox++)
      {
        const auto  & relcofab =   relco_lbd[dit[ibox]];
        auto valid             =     m_grids[dit[ibox]];
        const auto& eta_x      =  (*m_eta_x)[dit[ibox]];
        const auto& eta_y      =  (*m_eta_y)[dit[ibox]];
#if DIM==3                        
        const auto& eta_z      =  (*m_eta_z)[dit[ibox]];
#else
        //something fake to send in 2d
        const auto& eta_z      =  (*m_eta_x)[dit[ibox]];
#endif
        ///compute relaxation coefficient
        Proto::forallInPlace_i(relaxCoefResist, valid, 
                               relcofab, 
                               eta_x, eta_y, eta_z,
                               m_alpha, m_beta, m_dx);
      }
      
      return;
    }

    
    /// find which var you care about in m_grad
    static int gradIndex(int a_magnetic_component,
                         int a_derivative_direction)
    {
      //I think DIM (not 3) is correct here. 
      int varnum = a_magnetic_component*DIM + a_derivative_direction;
      return varnum;
    }//end function gradIndex
    
    ///
    virtual void
    cellGrad(pr_box_data_thrdim        &  a_grad,
             const pr_box_data_three   &  a_phi,
             pr_box                       a_valid)
    {
      
      //always.
      bool initToZero = true; double applyScale = 1.;
      
      pr_sten gradient_stencils[DIM];
      for(int deriv_dir = 0;  deriv_dir < DIM; deriv_dir++)
      {        
        pr_pt left_point = (-1)*pr_pt::Basis(deriv_dir);
        pr_pt righ_point = (+1)*pr_pt::Basis(deriv_dir);
        pr_shift righ(righ_point);
        pr_shift left(left_point);
        double weight = 1.0/(2.*m_dx);
        pr_sten  left_sten(left, weight);
        pr_sten  righ_sten(righ, weight);
        pr_sten  cent_grad = (righ_sten - left_sten);
        gradient_stencils[deriv_dir] = cent_grad;
      }

      int ibreak = 4586;  //gdb hook
      ///always 3 mag comps
      for(int comp_dir = 0; comp_dir < 3; comp_dir++)
      {
        pr_box_data_scalar mag_comp = Proto::slice(a_phi, comp_dir);
        //derivs only get dim
        for(int deriv_dir = 0;  deriv_dir < DIM; deriv_dir++)
        {
          int grad_index = gradIndex(comp_dir, deriv_dir);
          pr_box_data_scalar grad_comp = Proto::slice(a_grad, grad_index);

          //rememeber that this goes src, dst so grad_comp is the output
          gradient_stencils[deriv_dir].apply(mag_comp, grad_comp, a_valid, initToZero, applyScale);
          ibreak++;
        }
      }
    }
    ///Assumes that the ghost cells of a_phi are already set correctly.
    virtual void
    fillGrad(const pr_cell_data_three   & a_phi)
    {
      auto dit   = a_phi.begin();
      auto dbl   = a_phi.layout();
      int ibreak = 4586; //gdb hook
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto valid = dbl[dit[ibox]];
        cellGrad((*m_grad)[dit[ibox]], a_phi[dit[ibox]], valid);
        ibreak++;
      }
      m_grad->exchange(); //this is why grad is over the whole level
    } //end function co mputeCCGradient
  
      
    ///
    /**
       Stencil for partial(u_facedir)/partial(x_facedir)
       This is a centered difference of magnetism  only and relies upon 
       the magnetism having correct ghost cell values.
    **/
    int getFaceGradFromMagnStencil(lo_hi_center_t &  a_div_magn_lhc, 
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
        ch_mayday::Error("This one only should be when facedir==diffdir. Calling wrong stencil, maybe?");
      }
      a_div_magn_lhc.setBoxes(a_has_lo,
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
      a_div_magn_lhc.m_ce_sten = (1.0/m_dx)*(home_sten - left_one_sten);
      /**
         This stencil gets applied to magnetism
         and relies upon magn having correct ghost cells 
         so only centered diff.   
      **/
      a_div_magn_lhc.m_ce_box = a_flux_box;
      a_div_magn_lhc.m_has_lo = 0;
      a_div_magn_lhc.m_has_hi = 0;
      return 0;
    } ///end function getFaceGradFromMagnStencil

    ///
    /**
       This gets a face centered gradient where face_dir != div_dir.
       Where there a gradient on both sides of the face, it is a simple average 
       from cell centers.
       At boundaries it is extrapolated from neighbors.
    **/
    int getFaceGradFromGradStencil(lo_hi_center_t &  a_div_grad_lhc, 
                                   const pr_box   &  a_lo_box,
                                   const pr_box   &  a_hi_box,
                                   const pr_box   &  a_ce_box,
                                   const int      &  a_has_lo,
                                   const int      &  a_has_hi,
                                   const int      &  a_face_dir,
                                   const int      &  a_diff_dir,
                                   const pr_box   &  a_flux_box)
    {
      ///trying to keep stuff straight
      if(a_face_dir == a_diff_dir)
      {
        ch_mayday::Error("This one only should be when facedir!=diffdir. Calling wrong stencil, maybe?");
      }
      a_div_grad_lhc.setBoxes(a_has_lo,
                              a_has_hi,
                              a_lo_box,
                              a_hi_box,
                              a_ce_box);
      pr_pt home     =  pr_pt::Zeros();
      ///this extrapolates in face dir
      pr_pt righ_one = ( 1)*pr_pt::Basis(a_face_dir);
      pr_pt left_one = (-1)*pr_pt::Basis(a_face_dir);
      pr_pt left_two = (-2)*pr_pt::Basis(a_face_dir);
      //center stencil
      pr_shift righ_one_s(righ_one);
      pr_shift left_one_s(left_one);
      pr_shift left_two_s(left_two);
      pr_sten  righ_one_sten(righ_one_s, 1.);
      pr_sten  left_one_sten(left_one_s, 1.);
      pr_sten  left_two_sten(left_two_s,1.);
      
      pr_shift   home_s(home);
      pr_sten home_sten(home_s, 1.);
      //these get applied to gradients.
      a_div_grad_lhc.m_ce_sten = (1./2.)*(home_sten + left_one_sten);
       
      if(a_has_lo == 1)
      {
        a_div_grad_lhc.m_lo_sten = (1./2.)*(3.*home_sten - righ_one_sten);
      }
      if(a_has_hi == 1)
      {
        a_div_grad_lhc.m_hi_sten = (1./2.)*(3.*left_one_sten - left_two_sten);
      }
      return 0;
    }/// end function getFaceGradFromGradStencil
    

    ///
    int getFaceGradFromCellGradStencil(lo_hi_center_t &  a_grad_grad_lhc, 
                                       const pr_box   &  a_lo_box,
                                       const pr_box   &  a_hi_box,
                                       const pr_box   &  a_ce_box,
                                       const int      &  a_has_lo,
                                       const int      &  a_has_hi,
                                       const int      &  a_face_dir,
                                       const int      &  a_diff_dir,
                                       const pr_box   &  a_flux_box)
    {
      if(a_face_dir == a_diff_dir)
      {
        ch_mayday::Error("This one only should be when facedir!=diffdir. Calling wrong stencil, maybe?");
      }

      a_grad_grad_lhc.setBoxes(a_has_lo,
                               a_has_hi,
                               a_lo_box,
                               a_hi_box,
                               a_ce_box);
      pr_pt home     =  pr_pt::Zeros();
      ///extrapolations in the facedir direction.  
      pr_pt left_one = (-1)*pr_pt::Basis(a_face_dir);
      pr_pt left_two = (-2)*pr_pt::Basis(a_face_dir);
      pr_pt righ_one = ( 1)*pr_pt::Basis(a_face_dir);
      pr_shift home_s(home);
      pr_shift left_one_s(  left_one);
      pr_shift left_two_s(  left_two);
      pr_shift righ_one_s(  righ_one);
      
      pr_sten left_one_sten(left_one_s, 1.);
      pr_sten left_two_sten(left_two_s, 1.);
      pr_sten righ_one_sten(righ_one_s, 1.);
      pr_sten home_sten(    home_s    , 1.);
      
      //these get applied to gradients.
      //center stencil
      a_grad_grad_lhc.m_ce_sten = (1./2.)*(home_sten + left_one_sten);
       
      if(a_has_lo == 1)
      {
        a_grad_grad_lhc.m_lo_sten = (1.0/2.)*(3.*home_sten - righ_one_sten);
      }
      if(a_has_hi == 1)
      {
        a_grad_grad_lhc.m_hi_sten = (1.0/2.)*(3.*left_one_sten - left_two_sten);
      }
      return 0;
    }///end function getFaceGradFromCellGradStencil
    
    ///
    /**
       Average cell-centered gradients to the face centers
       Overwrite this with phi difference if face_dir == diff_dir.
    **/
    int
    getFaceDivAndGrad(pr_box_data_scalar          & a_face_div,
                      pr_box_data_thrdim          & a_face_grad,
                      const pr_box_data_three     & a_phi,
                      const pr_box_data_thrdim    & a_grad,
                      const pr_box                & a_flux_box,
                      const pr_box                & a_valid,
                      const int                   & a_face_dir)
    {

      pr_dom domain = m_grids.domain();
      ///first divergence 
      a_face_div.setVal(0.);
      pr_box_data_scalar div_b_incr(a_flux_box);
      pr_box lo_box, hi_box, ce_box;
      int has_lo=4586, has_hi=4586;
      PrCh_Lo_Hi_Center::loHiCenterFace(lo_box, hi_box, ce_box,
                                        has_lo, has_hi, a_flux_box,
                                        a_valid, a_face_dir, domain);
      
      for(int diff_dir = 0; diff_dir < DIM; diff_dir++)
      {
        div_b_incr.setVal(0.);  //this is important.
        if(diff_dir == a_face_dir)
        {
          /**
             Stencil for partial(u_facedir)/partial(x_facedir)
             This is a centered difference of magnetism  only and relies upon 
             the magn having correct ghost cell values.
             srcLo and srcHi are set to NULL because they will not be used.
          **/
          lo_hi_center_t div_magn_lhc;
          getFaceGradFromMagnStencil(div_magn_lhc, lo_box, hi_box, ce_box, has_lo, has_hi,
                                     a_face_dir, diff_dir, a_flux_box);
          
          pr_box_data_scalar magn_comp = Proto::slice(a_phi, diff_dir);
          /**
             This is standard chombo dst, src so div_b_incr is  the output.
          **/
          PrCh_Lo_Hi_Center::loHiCenterApply(div_b_incr, magn_comp, NULL, NULL,
                                             div_magn_lhc);
        }
        else
        {
          /**
             This gets a face centered gradient where face_dir != div_dir.
             Where there a gradient on both sides of the face, it is a simple average 
             from cell centers.
             At boundaries it is extrapolated from neighbors.
          **/
          lo_hi_center_t grad_grad_lhc;
          getFaceGradFromGradStencil(grad_grad_lhc, lo_box, hi_box, ce_box, has_lo, has_hi,
                                     a_face_dir, diff_dir, a_flux_box);
          int grad_ind = gradIndex(diff_dir, diff_dir); // indicies match because this is divergence
          pr_box_data_scalar grad_comp = Proto::slice(a_grad, grad_ind);
          /**
             This is standard chombo dst, src so div_b_incr is  the output.
          **/
          PrCh_Lo_Hi_Center::loHiCenterApply(div_b_incr, grad_comp, &grad_comp, &grad_comp,
                                             grad_grad_lhc);
        }
        a_face_div += div_b_incr;
      }

      //now for the gradient
      int ibreak = 4586;
      for(int diff_dir = 0; diff_dir < DIM; diff_dir++)
      {
        //magnetism has 3 vars
        for(int magn_dir = 0; magn_dir < 3; magn_dir++)
        {
          int grad_ind = gradIndex(magn_dir, diff_dir); // indicies match because this is divergence
          pr_box_data_scalar magn_comp      = Proto::slice(a_phi      , magn_dir);
          pr_box_data_scalar cell_grad_comp = Proto::slice(a_grad     , grad_ind);
          pr_box_data_scalar face_grad_comp = Proto::slice(a_face_grad, grad_ind);
          if( (a_face_dir == diff_dir) )
          {
            lo_hi_center_t div_magn_lhc;
            getFaceGradFromMagnStencil(div_magn_lhc, lo_box, hi_box, ce_box, has_lo, has_hi,
                                       a_face_dir, diff_dir, a_flux_box);
            //this is standard chombo dest, src so face_grad_comp is the output
            PrCh_Lo_Hi_Center::loHiCenterApply(face_grad_comp, magn_comp, NULL, NULL, 
                                               div_magn_lhc);
          }
          else
          {
            lo_hi_center_t grad_grad_lhc;
            getFaceGradFromGradStencil(grad_grad_lhc, lo_box, hi_box, ce_box, has_lo, has_hi,
                                       a_face_dir, diff_dir, a_flux_box);
            ///this is standard chombo dst, src so face_grad_comp is  the output
            PrCh_Lo_Hi_Center::loHiCenterApply(face_grad_comp, cell_grad_comp, &cell_grad_comp, &cell_grad_comp,
                                               grad_grad_lhc);
            ibreak = 0; //sometimes scoping is confusing in gdb
          }
          ibreak = 0; //stop here for //end of if() branch
        }   // end magn_dir loop
        ibreak = 0;   //stop here for //end of magn_dir loop
      }     // end diff_dir loop         
      return 0;
    }///end function getFaceDivAndGrad
       
    ///gets flux vector for a given face direction
    int
    getFlux(pr_box_data_three           & a_flux,
            const pr_box_data_three     & a_phi,
            const pr_box_data_scalar    & a_eta,
            const pr_box_data_thrdim    & a_grad,
            const pr_box                & a_valid,
            const int                   & a_face_dir)
    {

      //get face centered velocity divergence
      pr_box flux_box = a_valid;
      flux_box.m_high[a_face_dir] += 1;
      flux_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      pr_box_data_scalar face_divu(flux_box);
      pr_box_data_thrdim face_grad(flux_box);
      
      getFaceDivAndGrad(face_divu, face_grad, a_phi, a_grad, flux_box, a_valid,a_face_dir);
      
      //this is needed
      a_flux.setVal(0.);
      //copy the divergence into the diagonal component of the flux
      //diagonal component of the flux is the a_face_dir comp
      pr_box_data_scalar flux_diag = Proto::slice(a_flux, a_face_dir);
      flux_diag += face_divu;
      //now add in grad v - (grad v)^T
      //There are always 3 components.  There are SpaceDim directions.
      for (int component = 0; component < 3; component++)
      {
        pr_box_data_scalar flux_comp = Proto::slice(a_flux, component);
        //otherwise, grad(B)-grad^T(B)   = 0
        if(component != a_face_dir)
        {
          int ind_gradb  = gradIndex(component, a_face_dir);
          pr_box_data_scalar grad_b    = Proto::slice(face_grad, ind_gradb);
          flux_comp +=  grad_b;
          //do not do 3d derivs in 2d
          if(component < DIM)
          {
            int ind_gradb_T  = gradIndex(a_face_dir, component);
            pr_box_data_scalar grad_b_T    = Proto::slice(face_grad, ind_gradb_T);
            flux_comp -=  grad_b_T;
          }
        }
        flux_comp *= a_eta;
        flux_comp *= m_beta;
      }
      return 0;
    }///end function getFlux
    ///sets a_incr =  (flux_i+1/2 - flux_i-1/2)/m_dx
    int
    getDivergenceIncrement(pr_box_data_three           & a_incr,
                           const pr_box_data_three     & a_phi,
                           const pr_box_data_scalar    & a_eta,
                           const pr_box_data_thrdim    & a_grad,
                           const pr_box                & a_valid,
                           const int                   & a_face_dir)
    {
      CH_TIME("Pr_Resistivity::getDivergenceIcrement");
      pr_box face_box = a_valid;
      face_box.m_high[a_face_dir] += 1;
      face_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      
      pr_box_data_three  flux_dot_normal(face_box);
      getFlux(flux_dot_normal, a_phi, a_eta, a_grad, a_valid, a_face_dir);
      
      pr_pt left_point = (-1)*pr_pt::Basis(a_face_dir);
      pr_pt righ_point = (+1)*pr_pt::Basis(a_face_dir);
      pr_shift righ(righ_point);
      pr_shift home(Proto::Point::Zeros());

      pr_sten home_sten(home, 1./m_dx);
      pr_sten righ_sten(righ, 1./m_dx);

      pr_sten righ_minus_home = righ_sten - home_sten;
      
      bool initToZero = true; double applyScale = 1;
      a_incr.setVal(-4586.);
      ///always 3 for magnetism
      for(int icomp = 0; icomp < 3; icomp++)
      {
        pr_box_data_scalar flux_comp = Proto::slice(flux_dot_normal, icomp);
        pr_box_data_scalar incr_comp = Proto::slice(a_incr         , icomp);
        ///this arg list is src, dst so  incr_comp is the output
        righ_minus_home.apply(flux_comp, incr_comp,
                              a_valid , initToZero, applyScale);
      }
      
      return 0;
    }///end function getDivergenceIncrement
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(pr_cell_data_three         & a_lph,
                  const pr_cell_data_three   & a_phi)
    {
      CH_TIME("Proto_Resistivity_Op::applyOperator");
      auto dit   = a_lph.begin();

      //this fills m_grad with Grad(a_phi)
      fillGrad(a_phi);
      ///fill a_lph with div( b(x) grad (phi))
      for(int ibox = 0; ibox <  dit.localSize(); ibox++)
      {
        auto        & lphfab =      a_lph[dit[ibox]];
        const auto  & phifab =      a_phi[dit[ibox]];
        auto valid           =    m_grids[dit[ibox]];
        const auto& eta_x    = (*m_eta_x)[dit[ibox]];
        const auto& eta_y    = (*m_eta_y)[dit[ibox]];
        const auto& gradfab  =  (*m_grad)[dit[ibox]];
#if DIM==3                         
        const auto& eta_z    = (*m_eta_z)[dit[ibox]];
#endif        

        //the setvals are meaningful here. 
        pr_box_data_three xincr(valid); xincr.setVal(0.);
        pr_box_data_three yincr(valid); yincr.setVal(0.);
        pr_box_data_three zincr(valid); zincr.setVal(0.); 
        
        getDivergenceIncrement(xincr, phifab, eta_x, gradfab, valid, 0);
        getDivergenceIncrement(yincr, phifab, eta_y, gradfab, valid, 1);
#if DIM==3                                                         
        getDivergenceIncrement(zincr, phifab, eta_z, gradfab, valid, 2);
#endif
        ///a_lph = alpha**phi + beta(xincr + yincr + zincr)
        Proto::forallInPlace_i(operaResist, valid,
                               lphfab, phifab,
                               xincr, yincr, zincr,
                               m_alpha, m_beta);
      }
    }///end function applyOperator
    
    ///phi+= lambda*resid on red black
    virtual void
    relax(pr_cell_data_three       & a_phi,
          const pr_cell_data_three & a_rhs,
          int a_iterations)
    {
      CH_TIME("Pr_Resistivity::relax");
      
      pr_cell_data_three resid;
      create(resid, a_rhs);
      auto  dit       =  a_rhs.begin();
      auto& relco_lbd = (*m_relcoef);
      // gdb debugging hook.  semantically does nothing.
      int idebug = 0; 
      for(int iiter = 0; iiter < a_iterations;  iiter++)
      {
        for(int iredblack = 0; iredblack < 2; iredblack++)
        {
          residual(resid, a_phi, a_rhs);
          for(int ibox = 0; ibox < dit.localSize(); ibox++)
          {
            const auto & resfab =       resid[dit[ibox]];
            const auto & phifab =       a_phi[dit[ibox]];
            const auto & lamfab =   relco_lbd[dit[ibox]];
            auto valid          =     m_grids[dit[ibox]];

            Proto::forallInPlace_i(gsrbResist, valid,
                                   phifab, resfab, lamfab,
                                   iredblack);
          }
          idebug = 0; //stop here for end of dataiter  iteration
        }
        idebug = 0;   //stop here for end of red-black iteration
      }
      idebug = 0;     //stop here for end of relax 
    } //end functtion relax

    
    /// domain boundary conditions.
    /** 
        Anything can be hidden in a string, so the interface is general enough.
        For now, we shall restrict ourselves to the big two: Neumann and Dirichlet.
    **/
    virtual void
    fillGhostDataOutsideDomain(pr_cell_data_three          & a_phi,
                               const ch_dom                & a_domain)
    {
      CH_TIME("Pr_Resistivity::fillGhostdataOutsideDomain");
      typedef Proto_Transition::DataUtilities<3> pr_trans_data_three;
      pr_dom domain = ProtoCh::getProtoDomain(a_domain);
      pr_box dombox = domain.box();
      for(int idir = 0; idir < DIM; idir++)
      {
        for(int hilo = 0; hilo < 2; hilo++)
        {
          if(m_domBC == string("Dirichlet"))
          {
            pr_trans_data_three::domainGhostReflectOdd( a_phi, dombox, idir, hilo);
          }
          else if(m_domBC == string("Neumann"))
          {
            pr_trans_data_three::domainGhostReflectEven(a_phi, dombox, idir, hilo);
          }
          else
          {
            ch_mayday::Error("Proto_Resistivity_Op::fillGhostDataOutsideDomain: unrecognized bc string");
          }
        }
      }
    } //end functtion fillGhostDataOutsideDomain


    ///amr thing
    virtual int refToCoarser()
    {
      return m_refToCoar;
    }
  }; ///end class Proto_Resistivity_Op

  class Proto_Resistivity_Op_Factory: public ch_amrlevelop_fact_three
  {
    
  private:
    /// Strong construction, emacs, gdb, and Oxford commas are a winning formula for life.
    Proto_Resistivity_Op_Factory();
    
  public:

    virtual ~Proto_Resistivity_Op_Factory() { }

    vector<int>                                   m_ref_ratios;
    vector<pr_dbl>                                m_grids;
    vector<pair<ch_dom, double> >                 m_domain_dx;
    double                                        m_coarsest_dx;
    string                                        m_domBC;
    double                                        m_alpha;
    double                                        m_beta;
    int                                           m_max_coarse;
    vector< shared_ptr<pr_xfac_data_scalar> >     m_eta_x;
    vector< shared_ptr<pr_yfac_data_scalar> >     m_eta_y;
    vector< shared_ptr<pr_zfac_data_scalar> >     m_eta_z;

    
    ///the only constructor
    Proto_Resistivity_Op_Factory(
      const ch_dom                                    &  a_coarsest_dom,
      const vector<pr_dbl>                            &  a_grids,
      const vector<int>                               &  a_ref_ratios,
      const double                                    &  a_coarsest_dx,
      const string                                    &  a_domBC,
      const double                                    &  a_alpha,
      const double                                    &  a_beta,
      const vector< shared_ptr<pr_xfac_data_scalar> > &  a_eta_x,
      const vector< shared_ptr<pr_yfac_data_scalar> > &  a_eta_y,
      const vector< shared_ptr<pr_zfac_data_scalar> > &  a_eta_z)
    : ch_amrlevelop_fact_three()
    {
      CH_TIME("Pr_Resistivity_Factory_Constructor");
      m_max_coarse = 2;
      m_ref_ratios    =  a_ref_ratios;
      m_grids         =  a_grids;
      m_coarsest_dx   =  a_coarsest_dx;   
      m_domBC         =  a_domBC;         
      m_alpha         =  a_alpha;   
      m_beta          =  a_beta;
      m_eta_x         =  a_eta_x;
      m_eta_y         =  a_eta_y;
      m_eta_z         =  a_eta_z;

      
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
        ch_mayday::Error("ResistivityFactory::mgnewop: bogus domain");
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

    static pr_sten
    getStencilFromRefBox(pr_box a_ref_box,
                         int    a_ref_rat)
    {
      int numPts = a_ref_box.size();
      double weight = 1.0/double(numPts);
      pr_pt src_ref_ratio = pr_pt::Ones(a_ref_rat);
      pr_pt dst_ref_ratio = pr_pt::Ones(1);
      pr_pt dst_shift     = pr_pt::Zeros();
      /**
         Yes, I know I am setting things on an  empty stencil.
         Right now, adding a stencil with non unity refinement ratio 
         to an empty stencil without doing this causes a runtime error.
         dtg 11-28-2023
      **/
      pr_sten whole_sten;
      whole_sten.srcRatio()  = src_ref_ratio;
      whole_sten.destRatio() = dst_ref_ratio;
      for(auto boxit = a_ref_box.begin(); boxit != a_ref_box.end(); ++boxit)
      {
        pr_pt ref_pt = *boxit;
        pr_sten pt_sten(pr_shift(ref_pt), weight, dst_ref_ratio, dst_shift, src_ref_ratio);
        whole_sten += pt_sten;
      }
    
      return whole_sten;
    }

    /// for averaging stuff based in cells
    /**
       set coarse value to arithmetic average of all fine cell values.
    **/
    static void
    averageDownCell(pr_box_data_scalar       & a_coar,
                    const pr_box_data_scalar & a_fine,
                    const pr_box      & a_coar_valid,
                    int                 a_ref_rat)
    {
      pr_box unit_box(pr_pt::Zeros(), pr_pt::Ones());
      pr_pt pt_ref_rat = pr_pt::Ones( a_ref_rat);
      pr_box ref_box = unit_box.refine(pt_ref_rat);
      pr_sten whole_sten = getStencilFromRefBox(ref_box, a_ref_rat);

      
      pr_box_data_scalar& fine_cast = (pr_box_data_scalar&) a_fine;
      bool initToZero = true; double applyScale = 1.;
      //rememeber that this goes src, dst so a_coar is the output
      whole_sten.apply(fine_cast, a_coar, a_coar_valid, initToZero, applyScale);
    }
    /// for averaging stuff based in faces 
    /**
       set coarse value to arithmetic average of all fine face values.   Only faces in the same facedir plane are used.
       The valid box is cell centered.
    **/
    static void
    averageDownFace(pr_box_data_scalar       & a_coar,
                    const pr_box_data_scalar & a_fine,
                    const pr_box             & a_coar_valid,
                    int                        a_ref_rat,
                    int                        a_face_dir)
    {
      pr_box unit_box(pr_pt::Zeros(), pr_pt::Ones());
      pr_pt pt_ref_rat = pr_pt::Ones( a_ref_rat);
      pr_box ref_box = unit_box.refine(pt_ref_rat);
      // only include faces on the same facedir plane
      ref_box.m_high[a_face_dir] = 0;
      ref_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      
      pr_sten whole_sten = getStencilFromRefBox(ref_box, a_ref_rat);

      pr_box_data_scalar& fine_cast = (pr_box_data_scalar&) a_fine;
      bool initToZero = true; double applyScale = 1.;
      pr_box coar_face_valid = a_coar_valid;
      coar_face_valid.m_high[a_face_dir] += 1;
      coar_face_valid.recomputeSize();  // required to keep box data consistent after above shenannigans.
      
      whole_sten.apply(fine_cast, a_coar, coar_face_valid, initToZero, applyScale);

    }

    virtual void
    getCoefficients(
      shared_ptr<pr_xfac_data_scalar>  & a_eta_x,
      shared_ptr<pr_yfac_data_scalar>  & a_eta_y,
      shared_ptr<pr_zfac_data_scalar>  & a_eta_z,
      int                                a_iref_index,
      int                                a_mg_ref_rat)
    {
      CH_TIME("Pr_Resistivity_Factory_getCoefficients");
      if(a_mg_ref_rat == 1)
      {
        a_eta_x = m_eta_x[a_iref_index];
        a_eta_y = m_eta_y[a_iref_index];
        a_eta_z = m_eta_z[a_iref_index];
      }
      else
      {
        auto fine_dbl = m_grids[a_iref_index];
        //auto coar_dbl = fine_dbl.coarsen(pr_pt::Ones(a_mg_ref_rat));
        pr_dbl coar_dbl = fine_dbl.coarsen(pr_pt::Ones(a_mg_ref_rat));

        a_eta_x = shared_ptr<pr_xfac_data_scalar>(new pr_xfac_data_scalar(coar_dbl,  pr_pt::Zeros()));
        a_eta_y = shared_ptr<pr_yfac_data_scalar>(new pr_yfac_data_scalar(coar_dbl,  pr_pt::Zeros()));
#if DIM==3
        a_eta_z = shared_ptr<pr_zfac_data_scalar>(new pr_zfac_data_scalar(coar_dbl,  pr_pt::Zeros()));
#endif        
        auto dit   =  a_eta_x->begin();
        for(int ibox = 0; ibox < dit.localSize(); ibox++)
        {
          averageDownFace((*a_eta_x)[dit[ibox]], (*m_eta_x[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 0);
          averageDownFace((*a_eta_y)[dit[ibox]], (*m_eta_y[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 1);
#if DIM==3           
          averageDownFace((*a_eta_z)[dit[ibox]], (*m_eta_z[a_iref_index])[dit[ibox]], coar_dbl[dit[ibox]], a_mg_ref_rat, 2);
#endif        
        }//end loop over boxes
      }//  end averaging down case
    }//end function getCoefficients
      
    virtual inline ch_mglevelop_three*
    MGnewOp(const ch_dom         & a_dom_fine,
            int                    a_depth,
            bool                   a_homo_only = true)
    {
      CH_TIME("Pr_Resistivity_Factory_MGnewOp");
      int iref = get_ref_index(a_dom_fine);

      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      
      int mgRefRat = 1;
      for(int idep = 0; idep < a_depth; idep++)
      {
        mgRefRat*= 2;
      }
      ch_mglevelop_three* retval  = NULL;
      bool grids_coarsenable = m_grids[iref].coarsenable(pr_pt::Ones(2*mgRefRat));
      if((a_depth > 0) && (!grids_coarsenable))
      {
        retval = NULL;
      }
      else
      {
        shared_ptr<pr_xfac_data_scalar>  eta_x;
        shared_ptr<pr_yfac_data_scalar>  eta_y;
        shared_ptr<pr_zfac_data_scalar>  eta_z;

        
        getCoefficients(eta_x,eta_y,eta_z,
                        iref, mgRefRat);
        
        auto dblMG = eta_x->layout();
        double loc_dx = m_domain_dx[iref].second; 
        loc_dx *= mgRefRat;
        Proto_Resistivity_Op* derivedOp = 
          new Proto_Resistivity_Op(eta_x,
                                   eta_y,
                                   eta_z,
                                   dblMG,
                                   m_alpha,
                                   m_beta,
                                   loc_dx,
                                   loc_ref_fine,
                                   loc_ref_coar,
                                   m_domBC);
        retval = static_cast< ch_mglevelop_three* >(derivedOp);
      } //end case (we can return an operator)
      return retval;
    } //end function mgnewop

    ///
    virtual inline ch_amrlevelop_three* AMRnewOp(const ch_dom & a_indexSpace)
    {
      CH_TIME("Pr_Resistivity_Factory_AMRnewOp");
      int iref = get_ref_index(a_indexSpace);
      int loc_ref_fine, loc_ref_coar;
      get_local_ref_ratios(loc_ref_fine, loc_ref_coar, iref);
      double loc_dx = m_domain_dx[iref].second;

      Proto_Resistivity_Op* derivedOp = 
        new Proto_Resistivity_Op(m_eta_x[iref],
                                 m_eta_y[iref],
                                 m_eta_z[iref],
                                 m_grids[iref],
                                 m_alpha,
                                 m_beta,
                                 loc_dx,
                                 loc_ref_fine,
                                 loc_ref_coar,
                                 m_domBC);
      return  static_cast< ch_amrlevelop_three* >(derivedOp);
    } //end function AMRnewOp

    ///
    virtual inline int refToFiner(const ch_dom & a_indexSpace) const
    {
      int iref = get_ref_index(a_indexSpace);
      int ref_rat = m_ref_ratios[iref];
      return ref_rat;
    }//end function refToFiner
  }; //end class Proto_Resistivity_Op_Factory
}    //end namespace Pr_Ch_AMR_AMRElliptic
#endif
