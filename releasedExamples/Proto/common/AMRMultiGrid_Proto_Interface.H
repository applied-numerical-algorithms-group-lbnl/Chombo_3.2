#include "Base_Simplified_AMRLevelOp.H"
#include "ProtoInterface.H"

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
     
   In the end, all operators have to implement 7 functions.
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
     (4) refluxFineCorrection
     (5) refToCoarser.

     Some of these could be reasonably hardwired (force gsrb, first
     order prolongation, etc.), If this gets reused enough, we can
     do that. But for now, five shall be the number (and 5 << 83).  
     (4) and (5) can be delayed because they are amr things.

     An important design point here is that both base classes have no
     data.  I prefer base classes to just be base classes with all
     data living in the derived class.  There are many default function
     implementations here.  All will have timers.  If, for a
     particular operator, one of these functions become odious, that
     operator can supercede it (all these functions are virtual),
  **/
namespace PrCh_AMR_Elliptic
{
  template<int ncomp>
  class Base_Proto_AMRLevelOp: public Chombo::Base_Simplified_AMRLevelOp<Proto::LevelBoxData<double, ncomp> >
  {
  public:
    typedef                                     Proto::LevelBoxData<double, ncomp  >    pr_lbd;
    typedef Chombo::Base_Simplified_AMRLevelOp< Proto::LevelBoxData<double, ncomp> >    ch_amrlevelop;
    typedef                                    Chombo::Box                              ch_box;
    typedef                                    Chombo::MayDay                           ch_mayday;
    typedef                                    Chombo::ProblemDomain                    ch_dom;
    typedef                                     Proto::ProblemDomain                    pr_dom;
    typedef                                     Proto::Point                            pr_pt;
    
    ///This class has no data.  
    Base_Proto_AMRLevelOp():ch_amrlevelop() {}

    ///This class has no data
    virtual ~Base_Proto_AMRLevelOp(){}

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
      int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_CH_REAL,MPI_SUM, Chombo_MPI::comm);

      if (result != MPI_SUCCESS)
      {
        MayDay::Error("Base_Proto_AMRLevelOp::parallelSum communication error");
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
    
    ///a_phiThisLevel += I[2h->h](a_correctCoarse).  factor of 2.  indifferent to layout.
    virtual void
    prolongIncrement(pr_lbd              & a_phiMedi,
                     const pr_lbd        & a_corCoar)
    {
      //borrowed from Phil's levelmultigrid example in proto/example
      typedef Proto::BoxData<double, ncomp> pr_bd;
      pr_lbd localCor;
      this->createCoarser(   localCor, a_phiMedi, true);
      a_corCoar.copyTo(localCor);
      typedef Proto::InterpStencil<double> pr_interp_sten;
      pr_interp_sten interp_stencil = pr_interp_sten::Constant(pr_pt::Ones(2));
      auto grids =  a_phiMedi.layout();
      auto dit   =  a_phiMedi.begin();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        pr_bd & corfab =   localCor[dit[ibox]];
        pr_bd & phifab =  a_phiMedi[dit[ibox]];

        pr_bd fineCorr = interp_stencil(corfab);
        phifab += fineCorr;
      }
    }

    
    ///a_resCoarse = average(a_resfine) on processor. 
    virtual void
    restriction(pr_lbd & a_resCoar,
                pr_lbd & a_resFine)
    {
      ///taken from Phil's levelmultigrid example
      typedef Proto::Stencil<double> pr_sten;
      auto A = pr_sten::AvgDown(2);
      auto grids =  a_resCoar.layout();
      auto dit   =  a_resCoar.begin();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto & finefab =   a_resFine[dit[ibox]];
        auto & coarfab =   a_resCoar[dit[ibox]];
        coarfab.setVal(0.);
        coarfab |= A(finefab);
      }
    }

    ///
    virtual void setToZero(pr_lbd & a_lhs)
    {
      a_lhs.setToZero();
    }

    /// amr thing
    virtual void
    fillGhostDataOverCoarseFineInterface(pr_lbd       &  a_phiFine, 
                                         const pr_lbd &  a_phiCoar)
    {
      CH_TIME("Base_Proto_AMRLevelOp::fillGhostDataOverCoarseFineInterface");
      ch_mayday::Error("I think we can leave this until we are doing AMR");
    }
    

  };
} //end namespace PrCh_AMR_Elliptic

#endif
