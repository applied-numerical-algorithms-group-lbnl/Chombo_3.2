#include "Base_Simplified_AMRLevelOp.H"
#include "ProtoInterface.H"
#include "PrChFluxRegister.H"
#include "PrChAveraging.H"
#include "Proto.H"
#include "SPMD.H"
#ifndef ___AMRMultiGrid_Proto_Interface___
#define ___AMRMultiGrid_Proto_Interface___

///
/**
   Namespace PrCh_AMR_Elliptic

   Infrastructurue for using AMRMultiGrid as a solver front end and a
   proto back end.  I intend to rewrite in this fashion all of the
   some of the non-EB elliptic operators in Chombo3.2, including
   Poisson (constant coefficient and variable), ViscousTensor and
   Resistivity (variable coefficient for both).  I ignore all the
   use_homogeneous flags in the base classes because all domain bcs
   here should be homogeneous.  If you have inhomogeneous boundary
   conditions, slip them into the RHS.  Your data holder needs to be
   able to give its domain.  

   Brian van Straalen's wildly successful original templated
   elliptic solver to do the Martin-Cartwright algorithm has become
   a very large interface.  The base class hierarchy
   LevelOp/MGLevelOp/AMRLevelOp has (as of October 2023) 83 virtual
   functions.  Production codes can get large as more and more
   aggressive steps are taken to improve performance.  The
   AMRMultiGrid_Proto_Interface classes are meant to

   (1) Simplify using AMRMultiGrid.  
   (2) Provide the Proto backend for people who want to run on the device.

   The class structure reflects those priorities. 
   (1) Base_Simplified_AMRLevelOp provides reasonable defaults for many of these functions (independent of the proto interface).
   (2) Base_Proto_AMRLevelOp uses the proto interface to further winnow away at 83.
     
   In the end, all operators have to implement 5 functions.
**/
///
/**
   Base_Proto_AMRLevelOp
   All the proto-dependent default implementations go here.
   The distinction is probably academic but it does make for  smaller classes.

   This layer deals with proto<---->chombo conversions and will deal
   only with the purely virtual functions of
   Base_Simplified_AMRLevelOp.  This layer leans heavily upon the
   proto interface.  Only operator-dependent stuff is left out.

   The operator class has to implement the following:
   (1) applyOperator
   (2) relax 
   (3) fillGhostCellsOutsideDomain
   (4) getFlux
   (5) refToCoarser.

   Some of these could be reasonably hardwired (force gsrb, first
   order prolongation, etc.), If this gets reused enough, we can
   do that. But for now, five shall be the number (and 5 << 83).  
   There are many default function
   implementations here.  All will have timers.  If, for a
   particular operator, one of these functions become odious, that
   operator can supercede it (all these functions are virtual).

   An important design point here is that both base classes have no
   data.  I prefer base classes to just be base classes with all
   data living in the derived class.  
**/
namespace PrCh_AMR_Elliptic
{
  template<int ncomp>
  class Base_Proto_AMRLevelOp: public Chombo::Base_Simplified_AMRLevelOp<Proto::LevelBoxData<double, ncomp> >
  {
  public:
    
    typedef Chombo::Base_Simplified_AMRLevelOp< Proto::LevelBoxData<double, ncomp> >    ch_simple_op;
    typedef Chombo::AMRLevelOp<                 Proto::LevelBoxData<double, ncomp> >    ch_amrlevelop;
    typedef                                     Proto::LevelBoxData<double, ncomp  >    pr_lbd;
    typedef                                    Chombo::Box                              ch_box;
    typedef                                    Chombo::MayDay                           ch_mayday;
    typedef                                    Chombo::ProblemDomain                    ch_dom;
    typedef                                     Proto::ProblemDomain                    pr_dom;
    typedef                                     Proto::Point                            pr_pt;
    typedef                                     Proto::Shift                            pr_shift;
    typedef                                     Proto::Stencil<double >                 pr_sten;
    typedef                                     Proto::BoxData<double, ncomp>           pr_box_data;
    typedef                                     Proto::BoxData<double, 1    >           pr_box_data_sca;
    typedef                                     Proto::Box                              pr_box;
    typedef                                     Proto::DisjointBoxLayout                pr_dbl;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >     pr_cell_data;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >     pr_xfac_data;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >     pr_yfac_data;
    typedef Proto::LevelBoxData<double,  ncomp, Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >     pr_zfac_data;

    ///
    Base_Proto_AMRLevelOp():ch_simple_op()
    {
    }

    ///This class has no data
    virtual ~Base_Proto_AMRLevelOp(){}

    /// These can be on different layouts
    virtual void assign(    pr_lbd& a_lhs, const pr_lbd& a_rhs)
    {
      CH_TIME("Base_Proto_AMRLevelOp::assign");
      a_rhs.copyTo(a_lhs);
    }
    ///
    virtual ch_dom
    getDomainFromData(const pr_lbd & a_data)
    {
      CH_TIME("Base_Proto_AMRLevelOp::getDomainFromData");
      auto   dbl_pr = a_data.layout();
      pr_dom dom_pr = dbl_pr.domain();
      ch_dom retval = ProtoCh::getChomboDomain(dom_pr);
      return retval;
    }      
    
    ///exchange
    virtual void
    fillGhostDataOverFineFineInterface(pr_lbd       &  a_phiFine)
    {
      CH_TIME("Base_Proto_AMRLevelOp::fillGhostDataOverFineFineInterface");
      a_phiFine.exchange();
    }

    ///these two can be on different layouts
    virtual void
    copyTo(pr_lbd       &  a_dst,
           const pr_lbd &  a_src)
    {
      CH_TIME("Base_Proto_AMRLevelOp::copyTo");
      a_src.copyTo(a_dst);
    }

    ///
    virtual void
    homogeneousCFInterp(pr_lbd & a_phi)  = 0;
    
    ///
    virtual
    void createCoarser(pr_lbd                   & a_coar,
                       const pr_lbd             & a_fine,
                       bool ghosted)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarser");
      auto   dbl_fine = a_fine.layout();
      //in proto, dbl_fine.coarsen(pt) is a const function that returns the coarse layout.
      //in chombo, the equivalent changes dbl_fine and would be a bug here.
      auto   dbl_coar = dbl_fine.coarsen(2*pr_pt::Ones());
      pr_pt ghost = a_fine.ghost();
      a_coar.define(dbl_coar, ghost);
    }
    
    /// Cheers.
    virtual double  norm(const pr_lbd & a_rhs, int a_ord)
    {
      CH_TIME("Base_Proto_AMRLevelOp::norm");
      return a_rhs.absMax();
    }
    
    ///From EBLevelDataOps (credit to NASA's Mike Barad)
    double parallelSum(const double& a_value)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarser");
      // Find the sum of all a_value's
      Real sum = a_value;
#ifdef CH_MPI
      Real sendBuf = a_value;
      int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_CH_REAL,MPI_SUM, MPI_COMM_WORLD);

      if (result != MPI_SUCCESS)
      {
        Chombo::MayDay::Error("Base_Proto_AMRLevelOp::parallelSum communication error");
      }
#endif
      return sum;
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
         double    a_bco) = 0;
    
    ///
    /**
       These two have to on the same layout for this to work.
       Fortunately, I think that is forced by the calling function.
    **/
    virtual double
    dotProduct(const pr_lbd & a_fir,
               const pr_lbd & a_sec)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarser");
      auto dbl = a_fir.layout();
      auto dit = a_fir.begin();
      double localSum = 0.;
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        const auto& firfab = a_fir[dit[ibox]];
        const auto& secfab = a_sec[dit[ibox]];
        auto multfab = firfab*secfab;
        localSum += multfab.sum();
      }
      double mpisum = parallelSum(localSum);
      return mpisum;
    }
    
    ///
    virtual void
    create(pr_lbd       & a_lhs,
           const pr_lbd & a_rhs)
    {
      CH_TIME("Base_Proto_AMRLevelOp::create");
      auto grids = a_rhs.layout();
      auto ghost = a_rhs.ghost();
      a_lhs.define(grids, ghost);
    }
    
    ///
    virtual void createCoarsened(pr_lbd       & a_phiCoar,
                                 const pr_lbd & a_phiFine,
                                 const int    & a_refRat)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarsened (the AMR one)");
      const auto& gridsFine = a_phiFine.layout();
      auto        ghostFine = a_phiFine.ghost();
      //in proto, gridsFine.coarsen(pt) is a const function that returns the coarse layout.
      //in chombo, the equivalent changes gridsFine and would be a bug here.
      auto gridsCoar = gridsFine.coarsen( pr_pt::Ones(a_refRat));
      a_phiCoar.define(gridsCoar, ghostFine);
    }
    
    ///a_phiThisLevel += I[2h->h](a_correctCoarse).  factor of 2.  On proc.
    virtual void
    prolongIncrement(pr_lbd              & a_phi_fine,
                     const pr_lbd        & a_cor_coar)
    {
      int ref_rat = 2;
      vector<pr_sten> zeroth_interp =
        PrChTools::InterpolationUtilities::getZerothInterpStencils(ref_rat);
      auto grids =  a_phi_fine.layout();
      auto dit   =  a_phi_fine.begin();
      pr_lbd & cor_coar_cast = (pr_lbd &) a_cor_coar;
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        pr_box_data & cor_coar_fab =  cor_coar_cast[dit[ibox]];
        pr_box_data & phi_fine_fab =     a_phi_fine[dit[ibox]];
        pr_box          valid_fine =          grids[dit[ibox]];
        pr_box valid_coar = valid_fine.coarsen(2);
        for(int isten  = 0; isten < zeroth_interp.size(); isten++)
        {
          pr_sten color_sten = zeroth_interp[isten];
          /**
           // Doing this as an incr was too much trouble.   
           // Might have to change the sign of applyScale
           **/
          bool initToZeroFalse = false;  double applyScalePos = 1.;  
          color_sten.apply(cor_coar_fab   , phi_fine_fab, ///this arg list is src, dst  so phi_fine_fab is the output
                           valid_coar, initToZeroFalse, applyScalePos);
        }// end loop over colors
      }//   end loop over boxes
    } //    end function prolongIncrement


    ///  a_resCoar = I[h-2h]( a_rhsFine - L(phiFine, phiCoar)  )
    virtual void
    AMRRestrict(pr_lbd      & a_res_coar,
                const pr_lbd& a_rhs_fine,
                const pr_lbd& a_phi_fine,
                const pr_lbd& a_phi_coar,
                bool a_skip_res)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMRRestrict");
      int ref_rat = this->refToCoarser();
      pr_lbd res_fine, res_cofi;
      this->create(         res_fine, a_rhs_fine);
      //makes resfine = rhs-lphi
      this->AMRResidualNF(res_fine, a_phi_fine, a_phi_coar, a_rhs_fine, false);
     
      this->createCoarsened(res_cofi, res_fine, ref_rat);
      this->restriction(    res_cofi, res_fine, ref_rat);
      copyTo(a_res_coar, res_cofi);
    }

    ///a_phiFine += I[2h->h](a_phiCoar) 
    virtual void
    AMRProlong(pr_lbd      & a_phi_fine,
               const pr_lbd& a_cor_coar)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMRProlong");
      pr_lbd cor_cofi;
      int ref_rat = this->refToCoarser();
      this->createCoarsened(cor_cofi, a_phi_fine, ref_rat);
      copyTo(cor_cofi, a_cor_coar);
    
      vector<pr_sten> zeroth_interp =
        PrChTools::InterpolationUtilities::getZerothInterpStencils(ref_rat);
      auto grids =  a_phi_fine.layout();
      auto dit   =  a_phi_fine.begin();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        pr_box_data & cor_cofi_fab =    cor_cofi[dit[ibox]];
        pr_box_data & phi_fine_fab =  a_phi_fine[dit[ibox]];
        pr_box        valid_fine =       grids[dit[ibox]];
        pr_box valid_coar = valid_fine.coarsen(ref_rat);
        for(int isten  = 0; isten < zeroth_interp.size(); isten++)
        {
          pr_sten color_sten = zeroth_interp[isten];
          /**
           // Doing this as an incr was too much trouble.   
           // Might have to change the sign of applyScale
           **/
          bool initToZeroFalse = false;  double applyScalePos = 1.;  
          color_sten.apply(cor_cofi_fab   , phi_fine_fab, ///this arg list is src, dst  so phi_fine_fab is the output
                           valid_coar, initToZeroFalse, applyScalePos);
        } //end loop over colors
      }   //end loop over boxes
    } //end function AMRProlong
    
    ///a_resCoarse = average(a_resfine) on processor. 
    virtual void
    restriction(pr_lbd & a_resCoar,
                pr_lbd & a_resFine,
                int a_ref_rat)
    {
      auto grids =  a_resCoar.layout();
      auto dit   =  a_resCoar.begin();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        pr_box coar_valid =       grids[dit[ibox]];
        auto & finefab    =   a_resFine[dit[ibox]];
        auto & coarfab    =   a_resCoar[dit[ibox]];
        for(int icomp = 0;icomp <  ncomp; icomp++)
        {
          pr_box_data_sca fine_sca = Proto::slice(finefab, icomp);
          pr_box_data_sca coar_sca = Proto::slice(coarfab, icomp);
          PrChTools::AveragingUtilities::averageDownCell(coar_sca, fine_sca, coar_valid, a_ref_rat);
        }
      }
    }

    ///
    virtual void setToZero(pr_lbd & a_lhs)
    {
      a_lhs.setToZero();
    }

    /// 
    virtual void
    refluxFineCorrection(pr_lbd       &  a_lphCoar,
                         const pr_lbd &  a_phiFine, 
                         const pr_lbd &  a_phiCoar,
                         ch_amrlevelop*  a_finerOp)
    {
      Base_Proto_AMRLevelOp* cast_finer_op = dynamic_cast<Base_Proto_AMRLevelOp*>(a_finerOp) ;
      double  dxCoar = getDx();

      const pr_dbl& dbl_fine = a_phiFine.layout();
      const pr_dbl& dbl_coar = a_phiCoar.layout();
      int ref_ratio   = a_finerOp->refToCoarser();
      typedef PrCh_Tools::FluxRegister<ncomp> pr_ch_flux_reg;
      pr_ch_flux_reg flux_register(dbl_coar, dbl_fine, ref_ratio);

      shared_ptr<pr_xfac_data> xflux_coar(new pr_xfac_data(dbl_coar, pr_pt::Zero()));
      shared_ptr<pr_yfac_data> yflux_coar(new pr_yfac_data(dbl_coar, pr_pt::Zero()));
      shared_ptr<pr_xfac_data> xflux_fine(new pr_xfac_data(dbl_fine, pr_pt::Zero()));
      shared_ptr<pr_yfac_data> yflux_fine(new pr_yfac_data(dbl_fine, pr_pt::Zero()));
#if DIM==3                                                                      
      shared_ptr<pr_zfac_data> zflux_coar(new pr_zfac_data(dbl_coar, pr_pt::Zero()));
      shared_ptr<pr_zfac_data> zflux_fine(new pr_zfac_data(dbl_fine, pr_pt::Zero()));
#else
      shared_ptr<pr_zfac_data> zflux_coar, zflux_fine;
#endif
      /**
         // Because fillFluxes does not do anything about ghost cells,
         //  The finer solution its ghost cells over the C/F interface to be correct.
         //  Both resolutions need to be exchanged.
         // AMRPoissonOp.cpp contains the following cryptic comment which applies here:
         // const cast:  OK because we're changing ghost cells only
         // Translation.  We going to a bunch of evil casting of const objects but
         // you should trust us since we are just changing ghost cells.   
      **/
      pr_lbd& phi_fine_cast = (pr_lbd &) a_phiFine;
      pr_lbd& phi_coar_cast = (pr_lbd &) a_phiCoar;
      ch_dom domain_fine = this->getDomainFromData(phi_fine_cast);
      ch_dom domain_coar = this->getDomainFromData(phi_coar_cast);
      this->fillGhostDataOutsideDomain(         phi_coar_cast, domain_coar);
#if 1
      //fine boundary conditions are not done here in AMRPoissonOp.
      cast_finer_op->fillGhostDataOutsideDomain(phi_fine_cast, domain_fine);
      cast_finer_op->fillGhostDataOverCoarseFineInterface(phi_fine_cast, phi_coar_cast);
#endif      
      phi_fine_cast.exchange();
      phi_coar_cast.exchange();
      
      this->fillFluxes(         xflux_coar, yflux_coar, zflux_coar, a_phiCoar);
      cast_finer_op->fillFluxes(xflux_fine, yflux_fine, zflux_fine, a_phiFine);
      bool verbose = false; //turn on for flux register diagnostics
      flux_register.setCoarFlux(xflux_coar, yflux_coar, zflux_coar, verbose);
      flux_register.setFineFlux(xflux_fine, yflux_fine, zflux_fine, verbose);
      flux_register.reflux(a_lphCoar, dxCoar, verbose);

      ///  end reflux bit
    }    //end function refluxFineCorrection
      

    virtual 
    void fillFluxes(shared_ptr<pr_xfac_data> a_xflux,
                    shared_ptr<pr_yfac_data> a_yflux,
                    shared_ptr<pr_zfac_data> a_zflux,
                    const pr_lbd          &  a_phi)
    {
      auto grids =  a_phi.layout();
      auto dit   =  a_phi.begin();
      int ibreak = 4586; //just for getting gdb to stop in the right place
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        pr_box valid = grids[dit[ibox]];
        const pr_box_data & phi_fab = a_phi[dit[ibox]];
        {
          pr_box_data & xflux_fab =   (*a_xflux)[dit[ibox]];
          getFlux(xflux_fab, phi_fab, valid, 0, ibox);
          ibreak = 4586; //just for getting gdb to stop in the right place
        }
        {
          pr_box_data & yflux_fab =   (*a_yflux)[dit[ibox]];
          getFlux(yflux_fab, phi_fab, valid, 1, ibox);
          ibreak = 4586; //just for getting gdb to stop in the right place
        }
#if DIM==3        
        {
          pr_box_data & zflux_fab =   (*a_zflux)[dit[ibox]];
          getFlux(zflux_fab, phi_fab, valid, 2, ibox);
          ibreak = 4586; //just for getting gdb to stop in the right place
        }
#endif        
      }
    }
    ///
    /**
       This is the getFlux interface that refluxFineCorrection uses. 
       Fill in the values over a_face_box (memory is already allocated).
       ibox is in relation to the operator's grids
       (so you can find your variable coefficients if you have them).
    **/
    virtual void getFlux(pr_box_data       & a_flux,
                         const pr_box_data & a_phi,
                         const pr_box      & a_face_box,
                         int a_face_dir, int a_ibox) = 0;

    ///used in refluxing.
    virtual double getDx() const = 0;
  };
} //end namespace PrCh_AMR_Elliptic

#undef MG_NUM_COLORS

#endif
