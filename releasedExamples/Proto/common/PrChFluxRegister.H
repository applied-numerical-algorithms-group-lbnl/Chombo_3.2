
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
    typedef Chombo::Side   ch_side;
    //meta data
    pr_dbl m_grids_coar;
    pr_dbl m_grids_cofi;
    pr_dbl m_grids_fine;
    int    m_ref_ratio;
    /**
     // Flux buffers.   Have to do this flux-based because 
     // proto does not have reverse copier.
     //
     // All fluxes live on the coarse layout.  
     // cofi/fine indicates from where the flux originated from
     **/
    shared_ptr<pr_xfac_data> m_xflux_coar;
    shared_ptr<pr_xfac_data> m_xflux_fine;
    shared_ptr<pr_yfac_data> m_yflux_coar;
    shared_ptr<pr_yfac_data> m_yflux_fine;
    shared_ptr<pr_zfac_data> m_zflux_coar;
    shared_ptr<pr_zfac_data> m_zflux_fine;

    vector< vector<pr_box > >  m_boxes_CFI[2*DIM];

    ///
    /**
       This is way more complicated with periodic.
    **/
    void setCoarseLocations()
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vector< vector<pr_box> > & lo_cfi = m_boxes_CFI[idir];
        vector< vector<pr_box> > & hi_cfi = m_boxes_CFI[idir+ DIM];
        auto dit_coar   =  m_xflux_coar->begin();
        lo_cfi.resize(dit_coar.localSize());
        hi_cfi.resize(dit_coar.localSize());
        for(int ibox = 0; ibox < dit_coar.localSize(); ibox++)
        {
          pr_box valid_coar = m_grids_coar[dit_coar[ibox]];
          for(auto level_it = m_grids_cofi.begin(); level_it != m_grids_cofi.end(); ++level_it)
          {
            pr_box grid_cofi = m_grids_cofi[*level_it];
            pr_box lo_box = grid_cofi;  lo_box.m_high[idir] = lo_box.m_low [idir]; lo_box.recomputeSize();
            pr_box hi_box = grid_cofi;  hi_box.m_low [idir] = hi_box.m_high[idir]; hi_box.recomputeSize();
            pr_box buff_box_lo = lo_box.shift(pr_pt::Basis(idir, -1));
            pr_box buff_box_hi = hi_box.shift(pr_pt::Basis(idir,  1));
            if(valid_coar.contains(buff_box_lo))
            {
              lo_cfi[ibox].push_back(lo_box);
            }
            if(valid_coar.contains(buff_box_hi))
            {
              hi_cfi[ibox].push_back(hi_box);
            }
          } //end loop over coarsened fine grids
        }   //end loop over coarse grids
      }     //end loop over directions
    }//end function setCoarseLocations
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
      /**
       // Flux buffers.   Have to do this flux-based because 
       // proto does not have reverse copier.
       //
       // All fluxes live on the coarse layout.  
       // fine/coar indicates from where the flux originated from
       **/
      m_grids_cofi = a_grids_fine.coarsen(pr_pt::Ones(a_ref_ratio));
      m_xflux_coar = shared_ptr<pr_xfac_data>(new pr_xfac_data(m_grids_coar, pr_pt::Zeros()));
      m_xflux_fine = shared_ptr<pr_xfac_data>(new pr_xfac_data(m_grids_coar, pr_pt::Zeros()));
      m_yflux_coar = shared_ptr<pr_yfac_data>(new pr_yfac_data(m_grids_coar, pr_pt::Zeros()));
      m_yflux_fine = shared_ptr<pr_yfac_data>(new pr_yfac_data(m_grids_coar, pr_pt::Zeros()));
#if DIM==3      
      m_zflux_coar = shared_ptr<pr_zfac_data>(new pr_zfac_data(m_grids_coar, pr_pt::Zeros()));
      m_zflux_fine = shared_ptr<pr_zfac_data>(new pr_zfac_data(m_grids_coar, pr_pt::Zeros()));
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
      if(a_verbose)
      {
        Chombo::pout() << " Entering setCoarFlux: " << endl;
      } // end verbosity check
      a_xflux_coar->copyTo(*m_xflux_coar);
      a_yflux_coar->copyTo(*m_yflux_coar);
#if DIM==3      
      a_zflux_coar->copyTo(*m_zflux_coar);
#endif
      if(a_verbose)
      {
        Chombo::pout() << " Leaving setCoarFlux: " << endl;
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
      if(a_verbose)
      {
        Chombo::pout() << " Entering setFineFlux: " << endl;
      } // end verbosity check
      //create temporary cofi flux buffers
      shared_ptr<pr_xfac_data> t_xflux_cofi(new pr_xfac_data(m_grids_cofi, pr_pt::Zeros()));
      shared_ptr<pr_yfac_data> t_yflux_cofi(new pr_yfac_data(m_grids_cofi, pr_pt::Zeros()));
#if DIM==3                               
      shared_ptr<pr_zfac_data> t_zflux_cofi(new pr_zfac_data(m_grids_cofi, pr_pt::Zeros()));
#endif
      auto dit   =  a_xflux_fine->begin();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto & fine_xfab = (*a_xflux_fine)[dit[ibox]];
        auto & cofi_xfab = (*t_xflux_cofi)[dit[ibox]];
        auto & fine_yfab = (*a_yflux_fine)[dit[ibox]];
        auto & cofi_yfab = (*t_yflux_cofi)[dit[ibox]];
#if DIM==3        
        auto & fine_zfab = (*a_zflux_fine)[dit[ibox]];
        auto & cofi_zfab = (*t_zflux_cofi)[dit[ibox]];
#endif
        
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

      /**
       // Copy over to coarse distribution.  There will be many faces that are on  the coarse
       // level that are not on cofi.   These fluxes get initialized to zero so max/min and so on will make sense.
       // the _fine is just about where the flux came from.   It is on the coarse layout.
       **/
      m_xflux_fine->setVal(0.);
      m_yflux_fine->setVal(0.);
#if DIM==3                                                                                                         
      m_zflux_fine->setVal(0.);
#endif    

      if(a_verbose)
      {
        Chombo::pout() << " setFineFlux: before copyTo" << endl;
      } // end verbosity check
      
      t_xflux_cofi->copyTo(*m_xflux_fine);
      t_yflux_cofi->copyTo(*m_yflux_fine);
#if DIM==3                                                                                                         
      t_zflux_cofi->copyTo(*m_zflux_fine);
#endif
      if(a_verbose)
      {
        Chombo::pout() << " Leaving setFineFlux: " << endl;
      } // end verbosity check
    }     // end function setFineFlux

    ///
    /**
     **/
    pr_box_data_ncomp* getFluxPtrFine(int a_idir, int a_ibox)
    {
      pr_box_data_ncomp* retval = NULL;
      auto dit = m_xflux_fine->begin();
      if(a_idir == 0)
      {
        retval = &((*m_xflux_fine)[dit[a_ibox]]);
      }
      else if(a_idir == 1)
      {
        retval = &((*m_yflux_fine)[dit[a_ibox]]);
      }
#if DIM==3       
      else if(a_idir == 2)
      {
        retval = &((*m_zflux_fine)[dit[a_ibox]]);
      }
#endif
          
      return retval;
    }
    ///
    /**
     **/
    pr_box_data_ncomp* getFluxPtrCoar(int a_idir, int a_ibox)
    {
      pr_box_data_ncomp* retval = NULL;
      auto dit = m_xflux_coar->begin();
      if(a_idir == 0)
      {
        retval = &((*m_xflux_coar)[dit[a_ibox]]);
      }
      else if(a_idir == 1)
      {
        retval = &((*m_yflux_coar)[dit[a_ibox]]);
      }
#if DIM==3       
      else if(a_idir == 2)
      {
        retval = &((*m_zflux_coar)[dit[a_ibox]]);
      }
#endif        
      return retval;
    }
    ///
    /**
       All fluxes input on coarse resolution and coarse layout
       Side refers to which side of the coarse cell it is on.
    **/
    void refluxSide(pr_box_data_ncomp            & a_solu_coar,
                    const pr_box_data_ncomp      & a_flux_fine,
                    const pr_box_data_ncomp      & a_flux_coar,
                    const double                 & a_dx,
                    const vector<pr_box>         & a_cfi_boxes,
                    int a_idir,  ch_side::LoHiSide    a_coarCellSide)
    {
      /**
       // Stencils weights for fine and cofi will have opposite signs.
       // The shifts  will be the same for both.
       **/
      double wgt_fine = 4586,wgt_coar = 4586;
      pr_shift shift_val;
      if(a_coarCellSide == ch_side::Lo)
      {
        wgt_fine = -1/a_dx;
        wgt_coar =  1/a_dx;
        shift_val = pr_shift(pr_pt::Zeros());
      }
      else
      {
        wgt_fine =  1/a_dx;
        wgt_coar = -1/a_dx;
        shift_val = pr_shift(pr_pt::Basis(a_idir, 1));
      }
      pr_sten sten_fine(shift_val, wgt_fine);
      pr_sten sten_coar(shift_val, wgt_coar);
      for(int icomp = 0; icomp < ncomp; icomp++)
      {
        pr_box_data_sca solu_coar = Proto::slice(a_solu_coar, icomp);
        pr_box_data_sca flux_coar = Proto::slice(a_flux_coar, icomp);
        pr_box_data_sca flux_fine = Proto::slice(a_flux_fine, icomp);

        ///these are all done incrementally
        bool initToZero = false;   double applyScale = 1.;
        for(int ibox = 0; ibox < a_cfi_boxes.size(); ibox++)
        {
          ///this arg list is src, dst  so solu_coar is the output
          sten_fine.apply(flux_fine, solu_coar, a_cfi_boxes[ibox], 
                          initToZero, applyScale);
          sten_coar.apply(flux_coar, solu_coar, a_cfi_boxes[ibox], 
                          initToZero, applyScale);
        } //end loop over cfi boxes for this coarse grid
      }   //end loop over components
      return;
    }    //end function refluxSide
    /// 
    /**
     // U += Div(Ave(F_fine) - Fcoar)
     **/
    void reflux(pr_cell_data & a_solut_coar, const double& a_dx, bool a_verbose)

    {
      CH_TIME("PrChFluxRegister::reflux");
      
      if(a_verbose)
      {
        Chombo::pout() << " entering reFlux: " << endl;
      } // end verbosity check
      /**
       // Stage 1: 
       //  Set buffers = averaged fine flux (where it exists) and zero elsewhere.
       //  This makes all the data live on the coarse layout
       **/

      for(int idir = 0; idir < DIM; idir++)
      {
        auto dit_coar   =  a_solut_coar.begin();
        for(int ibox = 0; ibox < dit_coar.localSize(); ibox++)
        {
          pr_box_data_ncomp& coar_solut_fab = a_solut_coar[dit_coar[ibox]];
          pr_box_data_ncomp* flux_ptr_fine = getFluxPtrFine(idir, ibox);
          pr_box_data_ncomp* flux_ptr_coar = getFluxPtrCoar(idir, ibox);
          {
            /** 
                Boxes on the low side of the fine grid 
                (which makes the update on the high side of the coarse cell)
            **/
            const vector<pr_box> & cfi_boxes = m_boxes_CFI[idir][ibox];
            refluxSide(coar_solut_fab, *flux_ptr_fine, *flux_ptr_coar, a_dx, cfi_boxes, idir, ch_side::Hi);
          }  // end low side (of the fine grid) reflux
          {
            /** 
                Boxes on the high side of the fine grid 
                (which makes the update on the low side of the coarse cell)
            **/
            const vector<pr_box> & cfi_boxes = m_boxes_CFI[idir+ DIM][ibox];
            refluxSide(coar_solut_fab, *flux_ptr_fine, *flux_ptr_coar, a_dx, cfi_boxes, idir, ch_side::Lo);
          } //end high side (of the fine grid) reflux 
        }   //end loop over coarse local boxes
      }     //end loop over face directions

      if(a_verbose)
      {
        Chombo::pout() << " leaving reFlux " << endl;
      } // end verbosity check
      return;
    } //end function reflux
  };  //end class FluxRegister
}     //end namespace PrChTools
#endif
