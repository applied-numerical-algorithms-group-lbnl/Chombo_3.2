
#include "PrChUtilities.H"
#include "PrChAveraging.H"
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
   // reflux(solution, dx);        // solution += (1/dx)(Div(cofi_flux - coar_flux))
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
    typedef Proto::LevelIterator                                                                pr_level_it;
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
    shared_ptr<pr_xfac_data>   m_xbuff_coar;
    shared_ptr<pr_yfac_data>   m_ybuff_coar;
    shared_ptr<pr_zfac_data>   m_zbuff_coar;
    vector< vector<pr_box > >  m_boxes_CFI[2*DIM];

    /// this is way more complicated with periodic
    void setCoarseLocations()
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vector< vector<pr_box> > & lo_cfi = m_boxes_CFI[idir];
        vector< vector<pr_box> > & hi_cfi = m_boxes_CFI[idir+ DIM];
        auto dit_coar   =  a_xflux_coar->begin();
        lo_cfi.resize(dit.localSize());
        hi_cfi.resize(dit.localSize());
        for(int ibox = 0; ibox < dit.localSize(); ibox++)
        {
          pr_box grid_coar = m_grids_coar[dit[ibox]];
          for(auto level_it = m_grids_cofi.begin(); level_it != m_grids_cofi.end(); ++level_it)
          {
            pr_box grid_cofi = m_grids_cofi[level_it];
            pr_box lo_box = grid_cofi;  lo_box.m_high[idir] = lo_box.m_low [idir]; lo_box.recomputeSize();
            pr_box hi_box = grid_cofi;  hi_box.m_low [idir] = hi_box.m_high[idir]; hi_box.recomputeSize();
            pr_box buff_box_lo = lo_box.shift(pr_pt::Basis(idir, -1));
            pr_box buff_box_hi = hi_box.shift(pr_pt::Basis(idir,  1));
            if(grid_coar.contains(buff_box_lo))
            {
              lo_cfi[ibox].push_back(lo_box);
            }
            if(grid_coar.contains(buff_box_hi))
            {
              hi_cfi[ibox].push_back(hi_box);
            }
        }
    }
    
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
      m_grids_cofi = a_grids_fine.coarsen(pr_pt::Ones(a_ref_ratio));
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
      setCoarseLocations();
    }


    ///
    void setCoarFlux(shared_ptr<pr_xfac_data> a_xflux_coar,
                     shared_ptr<pr_yfac_data> a_yflux_coar,
                     shared_ptr<pr_zfac_data> a_zflux_coar,
                     bool a_verbose)
    {
      
      CH_TIME("PrChFluxRegister::setCoarFlux");
      a_xflux_coar->copyTo(*m_xflux_coar);
      a_yflux_coar->copyTo(*m_yflux_coar);
#if DIM==3      
      a_zflux_coar->copyTo(*m_zflux_coar);
#endif
      if(a_verbose)
      {
        Chombo::pout() << " after setCoarFlux: " << endl;
        Chombo::pout() << " m_xflux_coar->max() = " << m_xflux_coar->max() << " m_xflux_coar->min() = " << m_xflux_coar->min()  << endl;
        Chombo::pout() << " m_yflux_coar->max() = " << m_yflux_coar->max() << " m_yflux_coar->min() = " << m_yflux_coar->min()  << endl;
#if DIM==3                                                                                                                     
        Chombo::pout() << " m_zflux_coar->max() = " << m_zflux_coar->max() << " m_zflux_coar->min() = " << m_zflux_coar->min()  << endl;
#endif
      } // end verbosity check
    }   // end function setCoarFlux
    
    /// sets internal cofi buffers = average of fine fluxes
    /***
        cofi_flux = sum(fine_fluxes)/nfine_per_coar
    **/
    void setFineFlux(shared_ptr<pr_xfac_data> a_xflux_fine,
                     shared_ptr<pr_yfac_data> a_yflux_fine,
                     shared_ptr<pr_zfac_data> a_zflux_fine,
                     bool a_verbose)
    {
      CH_TIME("PrChFluxRegister::setFineFlux");
      auto dit   =  a_xflux_fine->begin();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto & fine_xfab = (*a_xflux_fine)[dit[ibox]];
        auto & cofi_xfab = (*m_xflux_cofi)[dit[ibox]];
        auto & fine_yfab = (*a_yflux_fine)[dit[ibox]];
        auto & cofi_yfab = (*m_yflux_cofi)[dit[ibox]];
#if DIM==3        
        auto & fine_zfab = (*a_zflux_fine)[dit[ibox]];
        auto & cofi_zfab = (*m_zflux_cofi)[dit[ibox]];
#endif
      if(a_verbose)
      {
        Chombo::pout() << " going into setFineFlux: " << endl;
        Chombo::pout() << " a_xflux_fine->max() = " << a_xflux_fine->max() << " a_xflux_fine->min() = " << a_xflux_fine->min()  << endl;
        Chombo::pout() << " a_yflux_fine->max() = " << a_yflux_fine->max() << " a_yflux_fine->min() = " << a_yflux_fine->min()  << endl;
#if DIM==3        
        Chombo::pout() << " a_zflux_fine->max() = " << a_zflux_fine->max() << " a_zflux_fine->min() = " << a_zflux_fine->min()  << endl;
#endif
      }   // end verbosity check
        pr_box cofi_valid = m_grids_cofi[dit[ibox]];
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          pr_box_data_sca xfine_comp = Proto::slice(fine_xfab, icomp);
          pr_box_data_sca xcofi_comp = Proto::slice(cofi_xfab, icomp);
          pr_box_data_sca yfine_comp = Proto::slice(fine_yfab, icomp);
          pr_box_data_sca ycofi_comp = Proto::slice(cofi_yfab, icomp);
#if DIM==3        
          pr_box_data_sca zfine_comp = Proto::slice(fine_zfab, icomp);
          pr_box_data_sca zcofi_comp = Proto::slice(cofi_zfab, icomp);
#endif
          PrChTools::AveragingUtilities::averageDownFace(xcofi_comp, xfine_comp, cofi_valid, m_ref_ratio, 0);
          PrChTools::AveragingUtilities::averageDownFace(ycofi_comp, yfine_comp, cofi_valid, m_ref_ratio, 1);
#if DIM==3    
          PrChTools::AveragingUtilities::averageDownFace(zcofi_comp, zfine_comp, cofi_valid, m_ref_ratio, 2);
#endif
        } // end loop over components
      }   // end loop over boxes
      if(a_verbose)
      {
        Chombo::pout() << " after setFineFlux: " << endl;
        Chombo::pout() << " m_xflux_cofi->max() = " << m_xflux_cofi->max() << " m_xflux_cofi->min() = " << m_xflux_cofi->min()  << endl;
        Chombo::pout() << " m_yflux_cofi->max() = " << m_yflux_cofi->max() << " m_yflux_cofi->min() = " << m_yflux_cofi->min()  << endl;
#if DIM==3                                                                                                                     
        Chombo::pout() << " m_zflux_cofi->max() = " << m_zflux_cofi->max() << " m_zflux_cofi->min() = " << m_zflux_cofi->min()  << endl;
#endif
      }   // end verbosity check
    }     // end function setFineFlux

    /// 
    /**
       // U += Div(Ave(F_fine) - Fcoar)
     **/
    void reflux(pr_cell_data & a_soln_coar, const double& a_dx, bool a_verbose)

    {
      CH_TIME("PrChFluxRegister::reflux");
      
      /**
       // Stage 1: 
       //  Set buffers = averaged fine flux (where it exists) and zero elsewhere.
       //  This makes all the data live on the coarse layout
      **/
      m_xbuff_coar->setVal(0.);
      m_ybuff_coar->setVal(0.);
      m_xflux_cofi->copyTo(*m_xbuff_coar);
      m_yflux_cofi->copyTo(*m_ybuff_coar);
#if DIM==3
      m_zbuff_coar.setVal(0.);
      m_zflux_cofi->copyTo(*m_zbuff_coar);
#endif
      if(a_verbose)
      {
        double max_x = m_xbuff_coar->max(); double min_x = m_xbuff_coar->min();
        double max_y = m_ybuff_coar->max(); double min_y = m_xbuff_coar->min();
#if DIM==3                                                                  
        double max_z = m_zbuff_coar->max(); double min_z = m_zbuff_coar->min();
#endif
        Chombo::pout() << "reflux: After stage 1 (coarsend fine flux only): " << endl;
        Chombo::pout() <<  " max_x = " << max_x << ", min_x = " << min_x << endl;
        Chombo::pout() <<  " max_y = " << max_y << ", min_y = " << min_y << endl;
#if DIM==3                                                               
        Chombo::pout() <<  " max_z = " << max_z << ", min_z = " << min_z << endl;
#endif
      }
      /**
         //  Stage 2:
         //  Loop through coarse/fine boxes and update solution 
      **/
      if(a_verbose)
      {
        Chombo::pout() << "reflux: After stage 2 (composite coar/cofi): " << endl;

        double max_x = m_xbuff_coar->max(); double min_x = m_xbuff_coar->min();
        double max_y = m_ybuff_coar->max(); double min_y = m_xbuff_coar->min();
#if DIM==3                                                                  
        double max_z = m_zbuff_coar->max(); double min_z = m_zbuff_coar->min();
#endif
        Chombo::pout() <<  " max_x = " << max_x << ", min_x = " << min_x << endl;
        Chombo::pout() <<  " max_y = " << max_y << ", min_y = " << min_y << endl;
#if DIM==3                                                               
        Chombo::pout() <<  " max_z = " << max_z << ", min_z = " << min_z << endl;
#endif
      }

      for(int idir = 0; idir < DIM; idir++)
      {
        vector<pr_box>& lo =
      auto dit   =  a_soln_coar.begin();
      ch_dit DataIterator = m_grids
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        
      }
      /** Stage 3:
       //  Refluxing dance.   Should be able to do all this within the dataiterator loop
       //  to avoid yet more temporary data.
       //  1. Subtract off coarse flux everywhere  to make buff = Fdiff
       Fdiff = Fcofi - Fcoar
       //  2. Compute reflux divergence D_R= Div(Fdiff) 
       //     Div(Fdiff) = (1/dxcoar)*sum(i<DIM)(Fdiff_i+1/2 - Fdiff_i-1/2)
       //  3. a_lph += Div(Fdiff)
       **/
      if(a_verbose)
      {
        Chombo::pout() << "reflux stage 3: reflux divergence and increment: " << endl;
      }
      int ibox_print = 0; ///only print tons of stuff on box this one box
        ///whether to print stuff for this particular box
        bool verb_box = a_verbose && (ibox == ibox_print);
        
        auto & buff_x = (*m_xbuff_coar)[dit[ibox]];
        auto & flux_x = (*m_xflux_coar)[dit[ibox]];
        auto & buff_y = (*m_ybuff_coar)[dit[ibox]];
        auto & flux_y = (*m_yflux_coar)[dit[ibox]];
#if DIM==3                          
        auto & buff_z = (*m_zbuff_coar)[dit[ibox]];
        auto & flux_z = (*m_zbuff_coar)[dit[ibox]];
#endif        

        if(verb_box)
        {
          Chombo::pout() << " reflux: Buffers should hold composite coar/cofi flux: " << endl;
          
          Chombo::pout() << " buff_x.max() = " << buff_x.max() << ", buff_x.min() = " << buff_x.min() << endl;
          Chombo::pout() << " buff_y.max() = " << buff_y.max() << ", buff_y.min() = " << buff_y.min() << endl;
#if DIM==3                                    
          Chombo::pout() << " buff_z.max() = " << buff_z.max() << ", buff_z.min() = " << buff_z.min() << endl;
#endif          
          Chombo::pout() << " flux_x.max() = " << flux_x.max() << ", flux_x.min() = " << flux_x.min() << endl;
          Chombo::pout() << " flux_y.max() = " << flux_y.max() << ", flux_y.min() = " << flux_y.min() << endl;
#if DIM==3                                                                                           
          Chombo::pout() << " flux_z.max() = " << flux_z.max() << ", flux_z.min() = " << flux_z.min() << endl;
#endif          
        }
        //  1. Subtract off coarse flux everywhere  to make buff = Fdiff
        buff_x -= flux_x;
        buff_y -= flux_y;
#if DIM==3                          
        buff_z -= flux_z;
#endif
        if(verb_box)
        {
          Chombo::pout() << " reflux: Buffers should hold Fdiff = Fcomposite - Fcoar (zero where they match): " << endl;
          Chombo::pout() << " buff_x.max() = " << buff_x.max() << ", buff_x.min() = " << buff_x.min() << endl;
          Chombo::pout() << " buff_y.max() = " << buff_y.max() << ", buff_y.min() = " << buff_y.min() << endl;
#if DIM==3                                    
          Chombo::pout() << " buff_z.max() = " << buff_z.max() << ", buff_z.min() = " << buff_z.min() << endl;
#endif          
        }
        
        //  2. Compute reflux divergence D_R= Div(Fdiff)
        pr_box valid = m_grids_coar[dit[ibox]];
        pr_box_data_ncomp div_reflux(valid);
        if(verb_box)
        {
          Chombo::pout() << "reflux going into component loop: " << endl;
          Chombo::pout() << " valid =   " << valid << endl;
        }
        ///only printing first component
        int icomp_print = 0;
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          bool verb_comp = verb_box && (icomp == icomp_print);
          auto output_comp = Proto::slice(div_reflux, icomp);

          /// the setvals are not strictly necessary
          pr_box_data_sca div_incr_x(valid);  div_incr_x.setVal(0.);
          pr_box_data_sca div_incr_y(valid);  div_incr_y.setVal(0.);
          pr_box_data_sca div_incr_z(valid);  div_incr_z.setVal(0.);

          auto xflux_comp = Proto::slice(buff_x, icomp);
          auto yflux_comp = Proto::slice(buff_y, icomp);
#if DIM==3          
          auto zflux_comp = Proto::slice(buff_z, icomp);
#endif          
          //rememeber that this goes src, dst so div_incr is the output
          div_stencils[0].apply(xflux_comp, div_incr_x, valid, initToZero, applyScale);
          div_stencils[1].apply(yflux_comp, div_incr_y, valid, initToZero, applyScale);
#if DIM==3        
          div_stencils[2].apply(zflux_comp, div_incr_z, valid, initToZero, applyScale);
#endif
          ///initialize divergence to zero then add increment from each direction
          output_comp.setVal(0.);
          output_comp += div_incr_x;
          output_comp += div_incr_y;
#if DIM==3        
          output_comp += div_incr_z;
#endif          
          if(verb_comp)
          {
            Chombo::pout() << " reflux: big incr bit: " << endl;

            Chombo::pout() << " div_incr_x.max() = " << div_incr_x.max()<< " div_incr_x.min() = " << div_incr_x.min()<< endl;
            Chombo::pout() << " div_incr_y.max() = " << div_incr_y.max()<< " div_incr_y.min() = " << div_incr_y.min()<< endl;
#if DIM==3                                                                                                           
            Chombo::pout() << " div_incr_z.max() = " << div_incr_z.max()<< " div_incr_z.min() = " << div_incr_z.min()<< endl;
#endif                                                                                                               
            Chombo::pout() << " divergence.max() = " <<output_comp.max()<< " divergence.min() = " <<output_comp.min() << endl;
          }
        } // end loop over components

        //  3. a_lph += Div(Fdiff)
        a_soln_coar[dit[ibox]] += div_reflux;
      } //end loop over boxes
    }//end function reflux
  };
} //end namespace Pr_Ch_AMR_AMRElliptic
#endif
