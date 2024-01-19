
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Vector.H"

#ifndef _PrChFluxRegister__
#define _PrChFluxRegister__
/// 
/**
   // his is a very simple flux register.   
   // t is more memory heavy than it has to be to achieve that clarity.
   // his is rigged so that if the coarse flux values and the fine flux values
   // re equal, the flux register does nothing.   This means that fine fluxes 
   // re applied nref^(DIM-1) times to account for the fact that 
   // here  nref^(DIM-1) fine fluxes  occupy same area as a coarse flux. 
   //  A 3D face looks like this:   
   //                 |---------------|---------------|           
   //                 |               |               |           
   //                 |               |               |           
   //                 |        F_FINE |        F_FINE |           
   //                 |      /\       |      /\       |           
   //                 |      \/       |      \/       |           
   //                 |               |               |           
   //                 |               |               |           
   //                 |              /|\ F_COARSE     |           
   //                 |---------------|---------------|           
   //                 |              \|/              |          
   //                 |               |               |           
   //                 |        F_FINE |        F_FINE |           
   //                 |      /\       |      /\       |           
   //                 |      \/       |      \/       |           
   //                 |               |               |           
   //                 |               |               |           
   //                 |               |               |           
   //                 |---------------|---------------|           
   // 
   // If F_FINE=F_COARSE, this class does the right thing.
   // Because it is easier for me (and I am the only user), 
   // I am abandoning the box by box interface traditionally used 
   // by flux registers and sending in the fluxes a level at a time.
   // Proto supports this very well.
   // I am also abandoning the incrementFine/Coarse interface. 
   // That semantic requires temporaries to program easily.
   // The usage pattern here is:
   // setFineFlux({fine_fluxes});  //set m_cofi_flux = Ave(fine_fluxes) 
   // setCoarFlux({coar_fluxes});  //set m_coar_flux =     coar_fluxes
   // reflux(solution, dx);        // solution -= (1/dx)(Div(cofi_flux - coar_flux))
   // Where there are no fine fluxes cofi_flux==coar_flux
**/
namespace PrCh_Tools
{
  ///
  template <int ncomp>
  class FluxRegister
  {
  private:
    
    FluxRegister();

  public:
    typedef Chombo::MayDay ch_mayday;
    typedef Proto::Stencil<double>                                                              pr_sten;
    typedef Proto::InterpStencil<double>                                                        pr_interp_sten;
    typedef Proto::Point                                                                        pr_pt;
    typedef Proto::Box                                                                          pr_box;
    typedef Proto::BoxData<double, 1  >                                                         pr_box_data_sca;
    typedef Proto::BoxData<double,ncomp>                                                        pr_box_data_ncomp;
    typedef Proto::DisjointBoxLayout                                                            pr_dbl;
    typedef Proto::Shift                                                                        pr_shift;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >     pr_cell_data;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >     pr_xfac_data;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >     pr_yfac_data;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >     pr_zfac_data;
    //meta data
    pr_dbl m_grids_coar;
    pr_dbl m_grids_cofi;
    pr_dbl m_grids_fine;
    int    m_ref_ratio;
    ///The data buffers.
    shared_ptr<pr_xfac_data> m_xflux_coar;
    shared_ptr<pr_xfac_data> m_xflux_cofi;
    shared_ptr<pr_yfac_data> m_yflux_coar;
    shared_ptr<pr_yfac_data> m_yflux_cofi;
    shared_ptr<pr_zfac_data> m_zflux_coar;
    shared_ptr<pr_zfac_data> m_zflux_cofi;

    ///scratch space on coarse layout
    shared_ptr<pr_xfac_data> m_xbuff_coar;
    shared_ptr<pr_yfac_data> m_ybuff_coar;
    shared_ptr<pr_zfac_data> m_zbuff_coar;

    
    /// 
    ~FluxRegister()
    {
      CH_TIME("PrChFluxRegister destructor");
    }
      
    /// 
    FluxRegister(const pr_dbl& a_grids_coar,
                 const pr_dbl& a_grids_fine,
                 int           a_ref_ratio)
    {
      CH_TIME("PrChFluxRegister constructor");
      m_grids_coar = a_grids_coar;
      m_grids_fine = a_grids_fine;
      m_ref_ratio  = a_ref_ratio;
      m_grids_cofi = a_grids_fine.coarsen(pr_pr::Ones(a_ref_ratio));
      m_xflux_coar = shared_ptr<pr_xfac_data>(new pr_xfac_data(m_grids_coar, pr_pt::Zeros()));
      m_xbuff_coar = shared_ptr<pr_xfac_data>(new pr_xfac_data(m_grids_coar, pr_pt::Zeros()));
      m_xflux_cofi = shared_ptr<pr_xfac_data>(new pr_xfac_data(m_grids_cofi, pr_pt::Zeros()));
      m_yflux_coar = shared_ptr<pr_yfac_data>(new pr_yfac_data(m_grids_coar, pr_pt::Zeros()));
      m_ybuff_coar = shared_ptr<pr_yfac_data>(new pr_yfac_data(m_grids_coar, pr_pt::Zeros()));
      m_yflux_cofi = shared_ptr<pr_yfac_data>(new pr_yfac_data(m_grids_cofi, pr_pt::Zeros()));
#if DIM==3      
      m_zflux_coar = shared_ptr<pr_zfac_data>(new pr_zfac_data(m_grids_coar, pr_pt::Zeros()));
      m_zbuff_coar = shared_ptr<pr_zfac_data>(new pr_zfac_data(m_grids_coar, pr_pt::Zeros()));
      m_zflux_cofi = shared_ptr<pr_zfac_data>(new pr_zfac_data(m_grids_cofi, pr_pt::Zeros()));
#endif      

    }


    ///
    void setCoarFlux(shared_ptr<pr_xfac_data> a_xflux_coar,
                     shared_ptr<pr_yfac_data> a_yflux_coar,
                     shared_ptr<pr_zfac_data> a_zflux_coar)
    {
      
      CH_TIME("PrChFluxRegister::setCoarFlux");
      a_xflux_coar->copyTo(*m_xflux_coar);
      a_yflux_coar->copyTo(*m_yflux_coar);
#if DIM==3      
      a_zflux_coar->copyTo(*m_zflux_coar);
#endif      
    }
    
    /// sets cofi buffers = average of fine fluxes
    /***
        cofi_flux = sum(fine_fluxes)/nfine_per_coar
    **/
    void setFineFlux(shared_ptr<pr_xfac_data> a_xflux_fine_ptr,
                     shared_ptr<pr_yfac_data> a_yflux_fine_ptr,
                     shared_ptr<pr_zfac_data> a_zflux_fine_ptr)
    {
      CH_TIME("PrChFluxRegister::setFineFlux");
      auto dit   =  a_xflux_fine_ptr->begin();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto & fine_xfab = (*a_xflux_fine_ptr)[dit[ibox]];
        auto & cofi_xfab = (*m_xflux_cofi_ptr)[dit[ibox]];
        auto & fine_yfab = (*a_yflux_fine_ptr)[dit[ibox]];
        auto & cofi_yfab = (*m_yflux_cofi_ptr)[dit[ibox]];
#if DIM==3        
        auto & fine_zfab = (*a_zflux_fine_ptr)[dit[ibox]];
        auto & cofi_zfab = (*a_zflux_cofi_ptr)[dit[ibox]];
#endif
        pr_box cofi_valid = m_grids_cofi[dit[ibox]];
        for(int icomp = 0; icomp < ncomp++)
        {
          pr_box_data_sca xfine_comp = Proto::slice(fine_xfab, icomp);
          pr_box_data_sca xcofi_comp = Proto::slice(cofi_xfab, icomp);
          pr_box_data_sca yfine_comp = Proto::slice(fine_yfab, icomp);
          pr_box_data_sca ycofi_comp = Proto::slice(cofi_yfab, icomp);
#if DIM==3        
          pr_box_data_sca zfine_comp = Proto::slice(fine_zfab, icomp);
          pr_box_data_sca zcofi_comp = Proto::slice(cofi_zfab, icomp);
#endif
          PrChUtilities<1>::averageDownFace(xcofi_comp, xfine_comp, cofi_valid, m_ref_ratio, 0);
          PrChUtilities<1>::averageDownFace(ycofi_comp, yfine_comp, cofi_valid, m_ref_ratio, 1);
#if DIM==3        
          PrChUtilities<1>::averageDownFace(zcofi_comp, zfine_comp, cofi_valid, m_ref_ratio, 2);
#endif
        }
      }
    }

    /// 
    /**
       // U -= Div(Ave(F_fine) - Fcoar)
     **/
    void reflux(pr_cell_data & a_soln_coar, const double& a_dx)
    {
      CH_TIME("PrChFluxRegister::reflux");
      
      //always.
      bool initToZero = true; double applyScale = 1.;
      
      pr_sten div_stencils[DIM];
      for(int deriv_dir = 0;  deriv_dir < DIM; deriv_dir++)
      {        
        double weight = 1.0/(a_dx);
        pr_pt righ_point = (+1)*pr_pt::Basis(deriv_dir);
        pr_pt home_point =      pr_pt::Zeros();
        pr_shift righ(righ_point);
        pr_shift home(home_point);
        pr_sten  home_sten(home, weight);
        pr_sten  righ_sten(righ, weight);
        pr_sten  iph_minu_imh = (righ_sten - home_sten);
        div_stencils[deriv_dir] = iph_minu_imh;
      }
      /**
         // Set buff = flux_coar.   
      **/
      m_xflux_coar->copyTo(*m_xbuff_coar);
      m_yflux_coar->copyTo(*m_ybuff_coar);
#if DIM==3                          
      m_zflux_coar->copyTo(*m_zbuff_coar);
#endif      
      /**
         //  Copy over from the flux_cofi
         //  so now buff = f_composite(flux_coar, flux_cofi)
         //    f_composite_i = flux_cofi_i  if i is covered by finer grid
         //                    flux_coar_i  otherwise
         //  These are the only all-to-all copies.  
      **/
      m_xflux_cofi->copyTo(*m_xbuff_coar);
      m_yflux_cofi->copyTo(*m_ybuff_coar);
#if DIM==3                          
      m_zflux_cofi->copyTo(*m_zbuff_coar);
#endif      

      /**
         //  Refluxing dance.   Should be able to do all this within the dataiterator loop
         //  to avoid yet more temporary data.
         //  1. Subtract off coarse flux everywhere  to make buff = Fdiff
                Fdiff = Fcofi - Fcoar
         //  2. Compute reflux divergence D_R= Div(Fdiff) 
         //     Div(Fdiff) = (1/dxcoar)*sum(i<DIM)(Fdiff_i+1/2 - Fdiff_i-1/2)
         //  3. a_lph -= Div(Fdiff)
       **/
    }

    auto dit   =  a_soln_coar.begin();
    for(int ibox = 0; ibox < dit.localSize(); ibox++)
    {
      auto & buff_x = (*m_xbuff_coar)[dit[ibox]];
      auto & flux_x = (*m_xflux_coar)[dit[ibox]];
      auto & buff_y = (*m_ybuff_coar)[dit[ibox]];
      auto & flux_y = (*m_yflux_coar)[dit[ibox]];
#if DIM==3                          
      auto & buff_z = (*m_zbuff_coar)[dit[ibox]];
      auto & flux_z = (*m_zbuff_coar)[dit[ibox]];
#endif        

      //  1. Subtract off coarse flux everywhere  to make buff = Fdiff
      buff_x -= flux_x;
      buff_y -= flux_z;
#if DIM==3                          
      buff_z -= flux_z;
#endif      
      //  2. Compute reflux divergence D_R= Div(Fdiff)
      pr_box valid_coar = m_grids_coar[dit[ibox]];
      pr_box_data_ncomp div_reflux(valid_coar);
      for(int icomp = 0; icomp < ncomp; icomp++)
      {
        auto output_comp = Proto::slice(div_reflux, icomp);
        
        pr_box_data_sca div_incr(valid_coar);
        auto xflux_comp = Proto::slice(buff_x, icomp);
        auto yflux_comp = Proto::slice(buff_y, icomp);
        auto zflux_comp = Proto::slice(buff_z, icomp);

        //rememeber that this goes src, dst so div_incr is the output
        div_stencils[0].apply(xflux_comp, div_incr, a_valid, initToZero, applyScale);
        output_comp+= div_incr;
        //rememeber that this goes src, dst so div_incr is the output
        div_stencils[1].apply(yflux_comp, div_incr, a_valid, initToZero, applyScale);
        output_comp+= div_incr;
#if DIM==3        
        //rememeber that this goes src, dst so div_incr is the output
        div_stencils[2].apply(zflux_comp, div_incr, a_valid, initToZero, applyScale);
        output_comp+= div_incr;
#endif              
      } // end loop over components

      //  3. a_lph -= Div(Fdiff)
      a_soln_coar[dit[ibox]] -= div_reflux;
    } //end loop over boxes
  }//end function reflux
        
  };
} //end namespace Pr_Ch_AMR_AMRElliptic
#endif
