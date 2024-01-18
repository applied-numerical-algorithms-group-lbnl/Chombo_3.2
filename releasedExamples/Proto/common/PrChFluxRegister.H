
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
    typedef Proto::BoxData<double, DIM>                                                         pr_box_data_vec;
    typedef Proto::BoxData<double, DIM*DIM>                                                     pr_box_data_ten;
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
    ///The data buffers.
    shared_ptr<pr_xfac_data> m_xflux_coar;
    shared_ptr<pr_xfac_data> m_xflux_cofi;
    shared_ptr<pr_yfac_data> m_yflux_coar;
    shared_ptr<pr_yfac_data> m_yflux_cofi;
    shared_ptr<pr_zfac_data> m_zflux_coar;
    shared_ptr<pr_zfac_data> m_zflux_cofi;
    /// 
    ~FluxRegister() {}
      
    ///
    FluxRegister(const pr_dbl& a_dbl_coar,
                 const pr_dbl& a_dbl_fine,
                 int           a_ref_ratio)
    {
      ch_mayday::Error("not implemented");
    }


    ///
    /**
       copies coarse fluxes into internal buffer.
     **/
    void incrementCoar(shared_ptr<pr_xfac_data> a_xflux_coar_ptr,
                       shared_ptr<pr_yfac_data> a_yflux_coar_ptr,
                       shared_ptr<pr_zfac_data> a_zflux_coar_ptr)
    {
      ch_mayday::Error("not implemented");
    }
    
    /// increments cofi buffers by average of fine fluxes
    /***
        cofi_flux = sum(fine_fluxes)/nfine_per_coar
    **/
    void incrementFine(shared_ptr<pr_xfac_data> a_xflux_fine_ptr,
                       shared_ptr<pr_yfac_data> a_yflux_fine_ptr,
                       shared_ptr<pr_zfac_data> a_zflux_fine_ptr)
    {
      ch_mayday::Error("not implemented");
    }

    /// 
    /**
       flux_diff = flux_cofi - flux_coar
       U -= weight*flux_diff
       typically, weight = 1./dx
     **/
    void reflux(pr_cell_data & a_coarseData, const double& a_weight)
    {
      ch_mayday::Error("not implemented");
    }
        
    /// sets all buffers to zero
    void setToZero()
    {
      ch_mayday::Error("not implemented");
    }
  };
} //end namespace Pr_Ch_AMR_AMRElliptic
#endif
